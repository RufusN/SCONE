module adjointLSRRPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential, expG, expG2 
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use physicsPackage_inter,           only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat            => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init         => init, &
                                             ndReg_getMatNames  => getMatNames, &
                                             ndReg_activate     => activate, &
                                             ndReg_kill         => kill, &
                                             ndReg_getNeutronMG => getNeutronMG
  use materialHandle_inter,           only : materialHandle
  use mgNeutronDatabase_inter,        only : mgNeutronDatabase
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  ! Visualisation
  use visualiser_class,               only : visualiser

  ! Tally map for fission rate
  use tallyMap_inter,                 only : tallyMap
  use tallyMapFactory_func,           only : new_tallyMap

  ! Random ray - or a standard particle
  ! Also particleState for easier output
  use particle_class,                      only : ray => particle, particleState

  ! For locks
  use omp_lib

  implicit none
  private

  ! Parameter for when to skip a tiny volume
  real(defReal), parameter :: volume_tolerance = 1.0E-10_defReal

  ! Parameters for indexing into matrices and spatial moments
  integer(shortInt), parameter :: x = 1, y = 2, z = 3, nDim = 3, &
                                  xx = 1, xy = 2, xz = 3, &
                                  yy = 4, yz = 5, zz = 6, &
                                  matSize = 6

  ! Convenient arithmetic parameters
  real(defFlt), parameter :: one_two = real(HALF,defFlt), &
                             two_three = real(2.0_defFlt/3.0_defFlt,defFlt)

  !!
  !! Physics package to perform The Random Ray Method (TRRM) eigenvalue calculations
  !! Uses linear sources
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes, sources, and keff
  !! are updated.
  !!
  !! Both inactive and active cycles occur, as in Monte Carlo. These can be terminated
  !! after a specified number of iterations or on reaching some chosen convergence
  !! criterion (though the latter hasn't been implemented yet).
  !!
  !! Calculates relative volume of different materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised such that the total domain
  !! volume is 1.0.
  !!
  !! IMPORTANT N.B.: Geometry type must be extended! Won't run if shrunk.
  !! This is because spatial discretisation is determined by the number of unique cells in the
  !! geometry.
  !! Also, this is obviously for multi-group calculations only.
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type adjointLSRRPhysicsPackage;
  !!     dead 10;              // Dead length where rays do not score to scalar fluxes
  !!     termination 100;      // Length a ray travels before it is terminated
  !!     rays 1000;            // Number of rays to sample per iteration
  !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;           // Number of scoring cycles (would use eps otherwise)
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #rho  0.6;#           // Optional stabilisation factor - default is 0, no stabilisation
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #fluxMap {<map>}#     // Optionally output one-group fluxes according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition>}
  !!   }
  !!
  !! Private Members
  !!   geom        -> Pointer to the geometry.
  !!   geomIdx     -> Index of the geometry in geometry Registry.
  !!   top         -> Top co-ordinates of the geometry bounding box.
  !!   bottom      -> Bottom co-ordinates of the geometry bounding box.
  !!   rand        -> Random number generator.
  !!   timerMain   -> Index of the timer defined to measure calculation time.
  !!   mgData      -> MG database. Calculation obviously cannot be run in CE.
  !!   nG          -> Number of energy groups, kept for convenience.
  !!   nCells      -> Number of unique cells in the geometry, kept for convenience.
  !!   nMat        -> Number of unique materials in the geometry - for convenience.
  !!   lengthPerIt -> Distance all rays travel in a single iteration - for convenience.
  !!
  !!   termination -> Distance a ray can travel before it is terminated
  !!   dead        -> Distance a ray must travel before it becomes active
  !!   pop         -> Number of rays to track per cycle
  !!   inactive    -> Number of inactive cycles to perform
  !!   active      -> Number of active cycles to perform
  !!   cache       -> Logical check whether to use distance caching
  !!   rho         -> Diagonal stabilisation factor, following Gunow
  !!   outputFile  -> Output file name
  !!   outputFormat-> Output file format
  !!   plotResults -> Plot results?
  !!   printFluxes -> Print fluxes?
  !!   printVolume -> Print volumes?
  !!   printCells  -> Print cell positions?
  !!   viz         -> Output visualiser
  !!   mapFission  -> Output fission rates across a given map?
  !!   resultsMap  -> The map across which to output fission rate results
  !!   mapFission  -> Output 1G flux across a given map?
  !!   resultsMap  -> The map across which to output 1G flux results
  !!
  !!   sigmaT      -> Local total cross section vector
  !!   nuSigmaF    -> Local nuSigmaF vector
  !!   sigmaS      -> Local flattened scattering matrix
  !!   chi         -> Local chi vector
  !!
  !!   keff        -> Estimated value of keff
  !!   keffScore   -> Vector holding cumulative keff score and keff^2 score
  !!   scalarFlux  -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux    -> Array of previous scalar flux values of length = nG * nCells
  !!   scalarX     -> Array of scalar flux X moments  of length = nG * nCells
  !!   scalarY     -> Array of scalar flux Y moments  of length = nG * nCells
  !!   scalarZ     -> Array of scalar flux Z moments  of length = nG * nCells
  !!   prevX       -> Array of previous scalar flux X moments of length = nG * nCells
  !!   prevY       -> Array of previous scalar flux Y moments of length = nG * nCells
  !!   prevZ       -> Array of previous scalar flux Z moments of length = nG * nCells
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!   sourceX     -> Array of neutron source X gradients of length = nG * nCells
  !!   sourceY     -> Array of neutron source Y gradients of length = nG * nCells
  !!   sourceZ     -> Array of neutron source Z gradients of length = nG * nCells
  !!   volume      -> Array of stochastically estimated cell volumes of length = nCells
  !!   centroid    -> Array of stochastically estimated cell centroids of length = 3 * nCells
  !!   momMat      -> Array of stochastically estimated cell spatial moments of length 6 * nCells
  !!                  (6 due to symmetry: xx, xy, xz, yy, yz, zz)
  !!   volumeTracks   -> Array of sum of length ever traversing a cell. Size = nCells
  !!   centroidTracks -> Array of sum of scores to the centroid of a cell. Size  = 3 * nCells
  !!   momTracks      -> Array of sum of scores to moments of a cell. Size = 6 * nCells
  !!   cellHit     -> Array tracking whether given cells have been hit during tracking
  !!   cellFound   -> Array tracking whether a cell was ever found
  !!   cellPos     -> Array of cell positions, populated once they are found
  !!   intersectionsTotal -> Total number of intersections performed during the calculation.
  !!
  !!   locks       -> OpenMP locks used when writing to scalar flux during transport
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: adjointLSRRPhysicsPackage
    private
    ! Components
    class(geometryStd), pointer           :: geom
    integer(shortInt)                     :: geomIdx     = 0
    real(defReal), dimension(3)           :: top         = ZERO
    real(defReal), dimension(3)           :: bottom      = ZERO
    type(RNG)                             :: rand
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    integer(shortInt)                     :: nMat        = 0
    real(defReal)                         :: lengthPerIt = ZERO

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
    real(defReal)      :: rho         = ZERO
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    logical(defBool)   :: plotResults = .false.
    logical(defBool)   :: printFlux   = .false.
    logical(defBool)   :: printVolume = .false.
    logical(defBool)   :: printCells  = .false.
    type(visualiser)   :: viz
    logical(defBool)   :: mapFission  = .false.
    class(tallyMap), allocatable :: resultsMap
    logical(defBool)   :: mapFlux     = .false.
    class(tallyMap), allocatable :: fluxMap

    !Perturbation
    integer(shortInt), dimension(:), allocatable  :: energyId 
    integer(shortInt)  :: NSegMax     = 100
    integer(shortInt)  :: NSegMaxTemp = 100
    integer(shortInt)  :: XStype      = 0
    integer(shortInt)  :: matPert     = 0
    real(defReal)      :: XSchange    = ZERO
    real(defReal)      :: numSum      = 0
    real(defReal)      :: denSum      = 0

! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: fission
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi
    real(defFlt), dimension(:), allocatable     :: sigmaC
    real(defFlt), dimension(:), allocatable     :: adjNuSigmaF
    real(defFlt), dimension(:), allocatable     :: adjSigmaS
    real(defFlt), dimension(:), allocatable     :: adjChi

    ! Results space
    real(defFlt)                               :: keff
    real(defReal), dimension(2)                :: keffScore
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defFlt), dimension(:), allocatable    :: scalarX
    real(defFlt), dimension(:), allocatable    :: scalarY
    real(defFlt), dimension(:), allocatable    :: scalarZ
    real(defFlt), dimension(:), allocatable    :: prevX
    real(defFlt), dimension(:), allocatable    :: prevY
    real(defFlt), dimension(:), allocatable    :: prevZ
    real(defReal), dimension(:,:), allocatable :: fluxScores
    real(defFlt), dimension(:), allocatable    :: source
    real(defFlt), dimension(:), allocatable    :: sourceX
    real(defFlt), dimension(:), allocatable    :: sourceY
    real(defFlt), dimension(:), allocatable    :: sourceZ
    real(defReal), dimension(:), allocatable   :: volume
    real(defReal), dimension(:), allocatable   :: volumeTracks
    real(defReal), dimension(:), allocatable   :: momMat
    real(defReal), dimension(:), allocatable   :: momTracks
    real(defReal), dimension(:), allocatable   :: centroid
    real(defReal), dimension(:), allocatable   :: centroidTracks

    real(defFlt)                               :: adjkeff
    real(defReal)                              :: sensitivity 
    real(defReal)                              :: deltaKeff 
    real(defReal), dimension(2)                :: adjkeffScore
    real(defReal), dimension(2)                :: sensitivityScore
    real(defReal), dimension(2)                :: deltaKeffScore
    real(defFlt), dimension(:), allocatable    :: adjScalarFlux
    real(defFlt), dimension(:), allocatable    :: adjPrevFlux
    real(defReal), dimension(:), allocatable   :: angularIP
    real(defFlt), dimension(:), allocatable    :: adjSource

    real(defFlt), dimension(:), allocatable    :: adjScalarX
    real(defFlt), dimension(:), allocatable    :: adjScalarY
    real(defFlt), dimension(:), allocatable    :: adjScalarZ
    real(defFlt), dimension(:), allocatable    :: adjPrevX
    real(defFlt), dimension(:), allocatable    :: adjPrevY
    real(defFlt), dimension(:), allocatable    :: adjPrevZ
    real(defFlt), dimension(:), allocatable    :: adjSourceX
    real(defFlt), dimension(:), allocatable    :: adjSourceY
    real(defFlt), dimension(:), allocatable    :: adjSourceZ

    ! Tracking cell properites
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos
    integer(longInt)                             :: intersectionsTotal = 0

    ! OMP locks
    integer(kind=omp_lock_kind), dimension(:), allocatable :: locks

    ! Timer bins
    integer(shortInt) :: timerMain
    integer(shortInt) :: timerTransport
    real (defReal)    :: time_transport = ZERO
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: initPerturbation
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sourceUpdateKernel
    procedure, private :: adjointSourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: adjointNormaliseFluxAndVolume
    procedure, private :: normaliseInnerProduct
    procedure, private :: resetFluxes
    procedure, private :: sensitivityCalculation
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings

  end type adjointLSRRPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    integer(shortInt)                             :: seed_temp, i, g, g1, m
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    class(dictionary),pointer                     :: tempDict, graphDict
    class(mgNeutronDatabase),pointer              :: db
    character(nameLen)                            :: geomName, graphType, nucData
    class(geometry), pointer                      :: geom
    type(outputFile)                              :: test_out
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
    character(100), parameter :: Here = 'init (adjointLSRRPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % active, 'active')
    call dict % get(self % inactive, 'inactive')
    
    ! Perform distance caching?
    call dict % getOrDefault(self % cache, 'cache', .false.)

    ! Stabilisation factor for negative in-group scattering
    call dict % getOrDefault(self % rho, 'rho', ZERO)

    ! Print fluxes?
    call dict % getOrDefault(self % printFlux, 'printFlux', .false.)

    ! Print volumes?
    call dict % getOrDefault(self % printVolume, 'printVolume', .false.)

    ! Print cell positions?
    call dict % getOrDefault(self % printCells, 'printCells', .false.)

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Check settings
    if (self % termination <= ZERO) call fatalError(Here, 'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError(Here, 'Must have 1 or more rays (pop).')

    ! Dead length can be less than zero but will be reset to zero if so
    if (self % dead < ZERO) then
      self % dead = ZERO
      print *,'Warning: Dead length of rays (dead) was negative. This has been set to 0 instead.'
    end if

    ! Ensure termination length is longer than dead length
    if (self % termination <= self % dead) call fatalError(Here,&
            'Ray termination length must be greater than ray dead length')

    ! Check whether there is a map for outputting fission rates
    ! If so, read and initialise the map to be used
    if (dict % isPresent('fissionMap')) then
      self % mapFission = .true.
      tempDict => dict % getDictPtr('fissionMap')
      call new_tallyMap(self % resultsMap, tempDict)
    else
      self % mapFission = .false.
    end if

    tempDict => dict % getDictPtr("Perturbation")
    call self % initPerturbation(tempDict)
    
    ! Check whether there is a map for outputting one-group fluxes
    ! If so, read and initialise the map to be used
    if (dict % isPresent('fluxMap')) then
      self % mapFlux = .true.
      tempDict => dict % getDictPtr('fluxMap')
      call new_tallyMap(self % fluxMap, tempDict)
    else
      self % mapFlux = .false.
    end if

    ! Register timer
    self % timerMain = registerTimer('simulationTime')
    self % timerTransport = registerTimer('transportTime')

    ! Initialise RNG
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')
    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % rand % init(seed)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))
    
    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'randomRayGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    geom    => gr_geomPtr(self % geomIdx)

    ! Ensure geometry is geometryStd
    select type(geom)
      type is (geometryStd)
        self % geom => geom
      class default
        call fatalError(Here,'Unrecognised geometry type')
    end select

    ! Ensure that geometry graph is extended
    graphDict => tempDict % getDictPtr('graph')
    call graphDict % get(graphType,'type')
    if (graphType /= 'extended') call fatalError(Here,&
            'Geometry graph type must be "extended" for random ray calculations.')

    ! Activatee nuclear data
    call ndReg_activate(P_NEUTRON_MG, nucData, self % geom % activeMats())

    ! Ensure that nuclear data is multi-group
    db => ndReg_getNeutronMG()
    if (.not. associated(db)) call fatalError(Here,&
            'No MG nuclear database was constructed')

    ! Ensure nuclear data is baseMgNeutronDatabase
    select type(db)
      type is (baseMgNeutronDatabase)
        self % mgData => db
      class default
        call fatalError(Here,'Unrecognised MG database type')
    end select

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()

    ! Get lower and upper corner of bounding box
    associate (aabb => self % geom % bounds())
      self % bottom = aabb(1:3)
      self % top    = aabb(4:6)
    end associate
    
    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing visualisation"
      call self % viz % makeViz()
      call self % viz % kill()
    endif
    
    ! Check for results plotting and initialise VTK
    call dict % getOrDefault(self % plotResults,'plot',.false.)
    if (self % plotResults) then
      ! Initialise a visualiser to be used when results are available
      print *, "Initialising results visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing geometry visualisation"
      call self % viz % initVTK()
    end if

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()

    ! Allocate results space
    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % scalarX(self % nCells * self % nG))
    allocate(self % scalarY(self % nCells * self % nG))
    allocate(self % scalarZ(self % nCells * self % nG))
    allocate(self % prevX(self % nCells * self % nG))
    allocate(self % prevY(self % nCells * self % nG))
    allocate(self % prevZ(self % nCells * self % nG))
    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % sourceX(self % nCells * self % nG))
    allocate(self % sourceY(self % nCells * self % nG))
    allocate(self % sourceZ(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % momMat(self % nCells * matSize))
    allocate(self % momTracks(self % nCells * matSize))
    allocate(self % centroid(self % nCells * nDim))
    allocate(self % centroidTracks(self % nCells * nDim))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))

    allocate(self % adjScalarFlux(self % nCells * self % nG))
    allocate(self % adjPrevFlux(self % nCells * self % nG))
    allocate(self % adjSource(self % nCells * self % nG)) 
    allocate(self % adjSourceX(self % nCells * self % nG))
    allocate(self % adjSourceY(self % nCells * self % nG))
    allocate(self % adjSourceZ(self % nCells * self % nG))
    allocate(self % adjScalarX(self % nCells * self % nG))
    allocate(self % adjScalarY(self % nCells * self % nG))
    allocate(self % adjScalarZ(self % nCells * self % nG))
    allocate(self % adjPrevX(self % nCells * self % nG))
    allocate(self % adjPrevY(self % nCells * self % nG))
    allocate(self % adjPrevZ(self % nCells * self % nG))

    allocate(self % angularIP(self % nCells * self % nG * self % nG))
    
    ! Set active length traveled per iteration
    self % lengthPerIt = (self % termination - self % dead) * self % pop
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
      call omp_init_lock(self % locks(i))
    end do

    ! Initialise local nuclear data
    ! TODO: clean nuclear database afterwards! It is no longer used
    !       and takes up memory.
    self % nMat = mm_nMat()
    allocate(self % sigmaT(self % nMat * self % nG))
    allocate(self % sigmaC(self % nMat * self % nG))
    allocate(self % nuSigmaF(self % nMat * self % nG))
    allocate(self % fission(self % nMat * self % nG))
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG))

    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, self % rand),defFlt)
        self % sigmaC(self % nG * (m - 1) + g) = real(mat % getCaptureXS(g, self % rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, self % rand),defFlt)
        self % fission(self % nG * (m - 1) + g) = real(mat % getFissionXS(g, self % rand),defFlt)
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, self % rand),defFlt)
        ! Include scattering multiplicity
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                  real(mat % getScatterXS(g1, g, self % rand) * mat % scatter % prod(g1, g) , defFlt)
        end do
      end do
    end do

    ! allocate(self % adjSigmaT(self % nMat * self % nG))
    allocate(self % adjNuSigmaF(self % nMat * self % nG))
    allocate(self % adjChi(self % nMat * self % nG))
    allocate(self % adjSigmaS(self % nMat * self % nG * self % nG))

    ! Adjoint fission source
    self % adjNuSigmaF = self % chi
    self % adjChi      = self % nuSigmaF

    do m = 1, self % nMat
      do g = 1, self % nG
        do g1 = 1, self % nG
          self % adjSigmaS(self % nG * self % nG * (m - 1) + self % nG * (g1 - 1) + g)  = &
                  self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)
        end do
      end do
    end do

  end subroutine init

  subroutine initPerturbation(self,dict)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(in)                   :: dict
    class(dictionary), pointer                      :: tempDict
    integer(shortInt)                               :: i, g
    character(100), parameter :: Here = 'initPerturbation (adjointLSRRPhysicsPackage_class.f90)'

    call dict % getOrDefault(self % XStype, 'XStype', 0)

    call dict % getOrDefault(self % matPert, 'material', 0)

    call dict % getOrDefault(self % XSchange, 'XSchange', ZERO)

    call dict % get(self % energyId, 'energyGroups')

    if (self % XStype == 3 .and. mod(size(self % energyId), 2) /= 0) then
      call fatalError(Here, 'energyGroups for scattering XS Perturbation must be given in pairs.')
    end if 

  
  end subroutine initPerturbation

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles()
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of The Random Ray Method.
  !!
  !! Randomly places the ray starting point and direction uniformly.
  !! Rays are tracked until they reach some specified termination length.
  !! During tracking, fluxes are attenuated (and adjusted according to BCs),
  !! scoring to fluxes and volume estimates when the ray has surpassed its 
  !! specified dead length.
  !!
  !! Inactive and active iterations occur, terminating subject either to 
  !! given criteria or when a fixed number of iterations has been passed.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  subroutine cycles(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defReal)                                 :: numTotal, denTotal
    real(defReal),save                            :: numSum, denSum
    real(defFlt)                                  :: hitRate, ONE_KEFF, adj_keff
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints, numSum, denSum)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % keff       = 1.0_defFlt
    self % scalarFlux = 0.0_defFlt
    self % scalarX  = 0.0_defFlt
    self % scalarY  = 0.0_defFlt
    self % scalarZ  = 0.0_defFlt
    self % prevFlux = 1.0_defFlt
    self % prevX    = 0.0_defFlt
    self % prevY    = 0.0_defFlt
    self % prevZ    = 0.0_defFlt
    self % fluxScores = 0.0_defReal
    self % keffScore  = 0.0_defReal
    self % source     = 0.0_defFlt
    self % sourceX = 0.0_defFlt
    self % sourceY = 0.0_defFlt
    self % sourceZ = 0.0_defFlt

    self % adjkeff    = 1.0_defFlt
    self % adjScalarFlux = 0.0_defFlt
    self % adjPrevFlux   = 1.0_defFlt
    self % adjSource     = 0.0_defFlt
    self % angularIP     = 0.0_defReal
    self % adjScalarX  = 0.0_defFlt
    self % adjScalarY  = 0.0_defFlt
    self % adjScalarZ  = 0.0_defFlt
    self % adjSourceX  = 0.0_defFlt
    self % adjSourceY  = 0.0_defFlt
    self % adjSourceZ  = 0.0_defFlt
    self % adjPrevX    = 0.0_defFlt
    self % adjPrevY    = 0.0_defFlt
    self % adjPrevZ    = 0.0_defFlt
    self % adjKeffScore     = ZERO
    self % deltaKeffScore   = ZERO
    self % sensitivityScore = ZERO

    ! Initialise other results
    self % cellHit        = 0
    self % volume         = 0.0_defReal
    self % volumeTracks   = 0.0_defReal
    self % momMat         = 0.0_defReal
    self % momTracks      = 0.0_defReal
    self % centroid       = 0.0_defReal
    self % centroidTracks = 0.0_defReal
    self % intersectionsTotal  = 0

    ! Initialise cell information
    self % cellFound = .false.
    self % cellPos = -INFINITY

    ! Stopping criterion is initially on flux convergence or number of convergence iterations.
    ! Will be replaced by RMS error in flux or number of scoring iterations afterwards.
    itInac = 0
    itAct  = 0
    isActive = .false.
    stoppingCriterion = .true.
    
    ! Power iteration
    do while( stoppingCriterion )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct

      ONE_KEFF = 1.0_defFlt / self % keff
      adj_KEFF = 1.0_defFlt / self % adjkeff
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        call self % sourceUpdateKernel(i, ONE_KEFF, it)
        call self % adjointSourceUpdateKernel(i, adj_KEFF, it)
      end do
      !$omp end parallel do

      ! Reset and start transport timer
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
      intersections = 0
      
      !$omp parallel do schedule(dynamic) reduction(+: intersections)
      do i = 1, self % pop

        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 

        ! Set ray attributes
        call self % initialiseRay(r)

        ! Transport ray until termination criterion met
        call self % transportSweep(r,ints, it,isActive)
        intersections = intersections + ints

      end do
      !$omp end parallel do

      self % intersectionsTotal = self % intersectionsTotal + intersections
      
      call timerStop(self % timerTransport)

      ! Update RNG on master thread
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)
      call self % adjointNormaliseFluxAndVolume(it)

      ! Calculate new k
      call self % calculateKeff()

      if (isActive) then
        call self % normaliseInnerProduct()
        ! Accumulate flux scores
        numTotal = ZERO
        denTotal = ZERO 

        !$omp parallel do schedule(static) reduction(+: numTotal, denTotal)
        do i = 1, self % nCells
          call self % sensitivityCalculation(i, ONE_KEFF, it, numSum, denSum)
          numTotal = numTotal + numSum
          denTotal = denTotal + denSum
        end do
        !$omp end parallel do

        self % deltaKeff  = self % keff * self % keff * (numTotal / denTotal) 
        self % sensitivity = (self % deltaKeff / ( self % keff * self % XSchange ) ) 
        call self % accumulateFluxAndKeffScores()
      end if

      ! Calculate proportion of cells that were hit
      hitRate = real(sum(self % cellHit),defFlt) / self % nCells
      self % cellHit = 0

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        stoppingCriterion = (itAct < self % active)
      else
        isActive = (itInac >= self % inactive)
      end if

      ! Set previous iteration flux to scalar flux
      ! and zero scalar flux
      call self % resetFluxes()

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)
      transport_T = timerTime(self % timerTransport)
      self % time_transport = self % time_transport + transport_T

      ! Predict time to end
      end_T = real(self % active + self % inactive, defReal) * elapsed_T / it
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(it)
      print *
      print *, 'Iteration: ', numToChar(it), ' of ', numToChar(self % active + self % inactive)
      if(isActive) then
        print *,'Active iterations'
      else
        print *,'Inactive iterations'
      end if
      print *, 'Cell hit rate: ', trim(numToChar(real(hitRate,defReal)))
      print *, 'keff: ', trim(numToChar(real(self % keff,defReal)))
      print *, 'adjoint_keff: ', trim(numToChar(real(self % adjkeff,defReal)))
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      print *, 'Time per integration (ns): ', &
              trim(numToChar(transport_T*10**9/(self % nG * intersections)))


    end do

    ! Finalise flux scores
    call self % finaliseFluxAndKeffScores(itAct)

  end subroutine cycles

  !!
  !! Initialises rays: samples initial position and direction,
  !! and performs the build operation
  !!
  subroutine initialiseRay(self, r)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, coord
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (adjointLSRRPhysicsPackage_class.f90)'

    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      coord = self % bottom + (self % top - self % bottom) * rand3

      ! Exit if point is inside the geometry
      call self % geom % whatIsAt(matIdx, cIdx, coord, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(coord, u, 1, ONE)
    call self % geom % placeCoord(r % coords)

    if (.not. self % cellFound(cIdx)) then
      !$omp critical 
      self % cellFound(cIdx) = .true.
      self % cellPos(cIdx,x) = coord(x)
      self % cellPos(cIdx,y) = coord(y)
      self % cellPos(cIdx,z) = coord(z)
      !$omp end critical
    end if

  end subroutine initialiseRay

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux, volume, flux moments, centroid, and moment matrix
  !! Records the number of integrations/ray movements.
  !!
  !! Note: when scoring moment matrices and centroids, assumes rays cross the
  !! entire cell in a flight: this will not always happen and hence there will
  !! be a bias - this should be small provided rays travel relatively large 
  !! distances compared to cell sizes
  !!
  subroutine transportSweep(self, r, ints, it, isActive)
    class(adjointLSRRPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: it
    logical(defBool), intent(in)                          :: isActive
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, &
                                                             matIdx0, baseIdx, centIdx, momIdx
    real(defReal)                                         :: totalLength, length, len2_12
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt, lenFlt2_2
    real(defFlt), dimension(self % nG)                    :: F1, F2, G1, G2, H, tau, delta, fluxVec, &
                                                             flatQ, gradQ, xInc, yInc, zInc, fluxVec0, &
                                                              flatQ0, Gn
    real(defReal), dimension(matSize)                     :: matScore
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec, &
                                                             xGradVec, yGradVec, zGradVec, &
                                                             xMomVec, yMomVec, zMomVec, &
                                                             xMomVecFW, yMomVecFW, zMomVecFW, &
                                                             xIncFW, yIncFW, zIncFW, scalarVecFW, &
                                                             sourceVecFW, deltaFW, avgFluxVecFW
    real(defReal), pointer, dimension(:)                  :: mid, momVec, centVec
    real(defReal), pointer                                :: volTrack
    real(defReal), dimension(3)                           :: r0, mu0, rC, r0Norm, rNorm
    real(defFlt), dimension(3)                            :: muFlt, rNormFlt, r0NormFlt
    integer(shortInt)                                     :: segCount, i, segCountCrit, segIdx, gIn, pIdx
    logical(defBool)                                      :: tally,  oversized
    real(defFlt), dimension(self % nG)                    :: avgFluxVec
    real(defReal), pointer, dimension(:)                  :: angularProd
    real(defFlt), dimension(self % NSegMax*2)             :: lenBack
    real(shortInt), dimension(self % NSegMax*2)           :: cIdxBack
    logical(defBool), dimension(self % NSegMax*2)         :: vacBack
    real(defReal), allocatable, dimension(:)              :: lenBackOversized, lenBackBuffer
    real(shortInt), allocatable, dimension(:)             :: cIdxBackOversized, cIdxBackBuffer
    logical(defBool), allocatable, dimension(:)           :: vacBackOversized, vacBackBuffer 
    real(defFlt), target, allocatable, dimension(:)       :: deltaRecordOversized
    real(defFlt), allocatable, dimension(:)               :: deltaRecordBuffer
    real(defFlt), target, dimension(self % nG * self % NSegMax*2) :: deltaRecord
    real(defFlt), target, dimension(self % nG * self % NSegMax*2) :: xIncRecord, yIncRecord,zIncRecord
    real(defFlt), target, allocatable, dimension(:)               :: xIncRecordOversized, yIncRecordOversized, zIncRecordOversized
    real(defFlt), allocatable, dimension(:)                       :: xIncRecordBuffer, yIncRecordBuffer, zIncRecordBuffer
    real(defFlt), target, dimension(3, self % NSegMax*2)          :: muFltRecord, rNormFltRecord, r0NormFltRecord
    real(defFlt), target, allocatable, dimension(:,:)             :: muFltOversized, rNormFltOversized, r0NormFltOversized
    real(defFlt), allocatable, dimension(:,:)                     :: muFltBuffer, rNormFltBuffer, r0NormFltBuffer
    real(defFlt), target, allocatable, dimension(:)       :: avgRecordOversized
    real(defFlt), allocatable, dimension(:)               :: avgRecordBuffer
    real(defFlt), target, dimension(self % nG * self % NSegMax*2) :: avgRecord

    
    
    ! Set initial angular flux to angle average of cell source
    ! Perhaps this should be different for linear source?
    ! Could in principle use the position and source gradient
    cIdx = r % coords % uniqueID
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % source(idx)
    end do

    ints = 0
    matIdx0 = 0
    segCount = 0
    totalLength = ZERO
    activeRay = .false.
    tally = .true.
    oversized = .false.

    do while (totalLength < self % termination + self % dead)

      ! Get ray coords for LS calculations
      mu0 = r % dirGlobal()
      r0  = r % rGlobal()

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        ! Cache total cross section
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
      end if

       ! Set maximum flight distance and ensure ray is active
      if (totalLength >= self % dead) then
        length = self % termination - totalLength 
        activeRay = .true.
      else
        length = self % dead - totalLength
      end if

      !second dead length
      if (totalLength >= self % termination) then
        length = self % termination + self % dead - totalLength
        tally = .false.
      end if

      ! Move ray
      ! Use distance caching or standard ray tracing
      ! Distance caching seems a little bit more unstable
      ! due to FP error accumulation, but is faster.
      ! This can be fixed by resetting the cache after X number
      ! of distance calculations.
      if (self % cache) then
        if (mod(ints,100_longInt) == 0)  cache % lvl = 0
        call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
      else
        call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)
      end if
      totalLength = totalLength + length

      ! Calculate the track centre
      rC = r0 + length * HALF * mu0
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = rC(:)
        !$omp end critical
      end if

      ints = ints + 1

      baseIdx = (cIdx - 1) * self % nG
      xGradVec => self % sourceX((baseIdx + 1):(baseIdx + self % nG))
      yGradVec => self % sourceY((baseIdx + 1):(baseIdx + self % nG))
      zGradVec => self % sourceZ((baseIdx + 1):(baseIdx + self % nG))
      sourceVec => self % source((baseIdx + 1):(baseIdx + self % nG))
      mid => self % centroid(((cIdx - 1) * nDim + 1):(cIdx * nDim))

      if (self % volume(cIdx) > volume_tolerance) then
        ! Compute the track centroid in local co-ordinates
        rNorm = rC - mid(1:nDim)
        ! Compute the entry point in local co-ordinates
        r0Norm = r0 - mid(1:nDim)
      else
        rNorm = 0
        r0Norm = - mu0 * HALF * length
      end if 

      ! Convert to floats for speed
      r0NormFlt = real(r0Norm,defFlt)
      rNormFlt = real(rNorm,defFlt)
      muFlt = real(mu0,defFlt)
      lenFlt = real(length,defFlt)

      ! Calculate source terms
      !$omp simd aligned(xGradVec, yGradVec, zGradVec, sourceVec)
      do g = 1, self % nG
        flatQ(g) = rNormFlt(x) * xGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(y) * yGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(z) * zGradVec(g)
        flatQ(g) = flatQ(g) + sourceVec(g)

        gradQ(g) = muFlt(x) * xGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(y) * yGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(z) * zGradVec(g)
      end do
      
      ! Compute exponentials necessary for angular flux update
      !$omp simd
      do g = 1, self % nG
        tau(g) = totVec(g) * lenFlt
        if (tau(g) < 1E-8_defFlt) then 
          tau(g) = 0.0_defFlt
        end if
      end do

      !$omp simd
      do g = 1, self % nG
        Gn(g) = expG(tau(g))
      end do

     !$omp simd
      do g = 1, self % nG
        F1(g)  = 1.0_defFlt - tau(g) * Gn(g) !expTau(tau(g)) * lenFlt
      end do

      !$omp simd
      do g = 1, self % nG
        F2(g) = (2.0_defFlt * Gn(g) - F1(g)) * lenFlt 
      end do

      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - flatQ(g)) * F1(g)  & 
                    - one_two * gradQ(g) * F2(g) 
      end do

      ! Intermediate flux variable creation or update
      !$omp simd
      do g = 1, self % nG
        fluxVec0(g) = fluxVec(g)
      end do

      ! Flux vector update
      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g) * tau(g) 
      end do

      !$omp simd
      do g = 1, self % nG
        flatQ0(g) = flatQ(g)
      end do

      ! Accumulate to scalar flux
      if (activeRay) then
        segCount = segCount + 1
        ! Precompute geometric info to keep it out of the lock
        len2_12 = length * length / 12
        matScore(xx) = length * (rnorm(x) * rnorm(x) + mu0(x) * mu0(x) * len2_12)
        matScore(xy) = length * (rnorm(x) * rnorm(y) + mu0(x) * mu0(y) * len2_12)
        matScore(xz) = length * (rnorm(x) * rnorm(z) + mu0(x) * mu0(z) * len2_12)
        matScore(yy) = length * (rnorm(y) * rnorm(y) + mu0(y) * mu0(y) * len2_12)
        matScore(yz) = length * (rnorm(y) * rnorm(z) + mu0(y) * mu0(z) * len2_12)
        matScore(zz) = length * (rnorm(z) * rnorm(z) + mu0(z) * mu0(z) * len2_12)
        centIdx = nDim * (cIdx - 1)
        momIdx = matSize * (cIdx - 1)
        
        rC = rC * length
        
        ! Compute necessary exponential quantities
        ! Follows those in Gunow
        lenFlt2_2 = lenFlt * lenFlt * one_two
        
       !$omp simd
        do g = 1, self % nG
          H(g) = ( F1(g) - Gn(g) ) !expH(tau(g))
        end do
      
        !$omp simd
        do g = 1, self % nG
            G1(g) = one_two - H(g)
        end do

        !$omp simd 
        do g = 1, self % nG
          G2(g) = expG2(tau(g)) 
        end do

        do g = 1, self % nG
          G1(g) = G1(g) * flatQ(g) * lenFlt
          G2(g) = G2(g) * gradQ(g) * lenFlt2_2
          H(g)  = H(g) * fluxVec0(g) * lenFlt
          H(g) = (G1(g) + G2(g) + H(g)) * lenFlt
          flatQ(g) = flatQ(g) * lenFlt + delta(g) * lenFlt
        end do

        !$omp simd
        do g = 1, self % nG
          xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) 
          yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) 
          zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g)
        end do

        if (tally) then
            segCountCrit = segCount
        end if
 
        if ( (segCount >= (self % NSegMax - 1)) .and. .not. oversized) then

          oversized = .true.
          allocate(lenBackOversized(size(lenBack) * 2))
          lenBackOversized(1:size(lenBack)) = lenBack

          allocate(cIdxBackOversized(size(cIdxBack) * 2))
          cIdxBackOversized(1:size(cIdxBack)) = cIdxBack

          allocate(vacBackOversized(size(vacBack) * 2))
          vacBackOversized(1:size(vacBack)) = vacBack

          allocate(deltaRecordOversized(size(deltaRecord) * 2))
          deltaRecordOversized(1:size(deltaRecord)) = deltaRecord

          allocate(xIncRecordOversized(size(xIncRecord) * 2))
          xIncRecordOversized(1:size(xIncRecord)) = xIncRecord

          allocate(yIncRecordOversized(size(yIncRecord) * 2))
          yIncRecordOversized(1:size(yIncRecord)) = yIncRecord

          allocate(zIncRecordOversized(size(zIncRecord) * 2))
          zIncRecordOversized(1:size(zIncRecord)) = zIncRecord

          allocate(rNormFltOversized(3, size(rNormFltRecord(1,:)) * 2))
          rNormFltOversized(:, 1:size(rNormFltRecord(1,:))) = rNormFltRecord

          allocate(r0NormFltOversized(3, size(r0NormFltRecord(1,:)) * 2))
          r0NormFltOversized(:, 1:size(r0NormFltRecord(1,:))) = r0NormFltRecord

          allocate(muFltOversized(3, size(muFltRecord(1,:)) * 2))
          muFltOversized(:, 1:size(muFltRecord(1,:))) = muFltRecord

          allocate(avgRecordOversized(size(avgRecord) * 2))
          avgRecordOversized(1:size(avgRecord)) = avgRecord

        elseif (oversized) then

          if  (segCount >= (size(vacBackOversized) - 1) ) then
            allocate(lenBackBuffer(size(lenBackOversized) * 2)) 
            lenBackBuffer(1:size(lenBackOversized)) = lenBackOversized
            call move_alloc(lenBackBuffer, lenBackOversized)

            allocate(cIdxBackBuffer(size(cIdxBackOversized) * 2)) 
            cIdxBackBuffer(1:size(cIdxBackOversized)) = cIdxBackOversized
            call move_alloc(cIdxBackBuffer, cIdxBackOversized)

            allocate(vacBackBuffer(size(vacBackOversized) * 2))  
            vacBackBuffer(1:size(vacBackOversized)) = vacBackOversized
            call move_alloc(vacBackBuffer, vacBackOversized)

            allocate(deltaRecordBuffer(size(deltaRecordOversized) * 2)) 
            deltaRecordBuffer(1:size(deltaRecordOversized)) = deltaRecordOversized
            call move_alloc(deltaRecordBuffer, deltaRecordOversized)

            allocate(xIncRecordBuffer(size(xIncRecordOversized) * 2))  
            xIncRecordBuffer(1:size(xIncRecordOversized)) = xIncRecordOversized 
            call move_alloc(xIncRecordBuffer, xIncRecordOversized)

            allocate(yIncRecordBuffer(size(yIncRecordOversized) * 2))  
            yIncRecordBuffer(1:size(yIncRecordOversized)) = yIncRecordOversized
            call move_alloc(yIncRecordBuffer, yIncRecordOversized)

            allocate(zIncRecordBuffer(size(zIncRecordOversized) * 2))   
            zIncRecordBuffer(1:size(zIncRecordOversized)) = zIncRecordOversized
            call move_alloc(zIncRecordBuffer, zIncRecordOversized)

            allocate(rNormFltBuffer(3, size(rNormFltOversized(1,:)) * 2))
            rNormFltBuffer(:, 1:size(rNormFltOversized,2)) = rNormFltOversized
            call move_alloc(rNormFltBuffer, rNormFltOversized)

            allocate(r0NormFltBuffer(3, size(r0NormFltOversized(1,:)) * 2))
            r0NormFltBuffer(:, 1:size(r0NormFltOversized,2)) = r0NormFltOversized
            call move_alloc(r0NormFltBuffer, r0NormFltOversized)

            allocate(muFltBuffer(3, size(muFltOversized(1,:)) * 2))
            muFltBuffer(:, 1:size(muFltOversized,2)) = muFltOversized
            call move_alloc(muFltBuffer, muFltOversized)

            allocate(avgRecordBuffer(size(avgRecordOversized) * 2))
            avgRecordBuffer(1:size(avgRecordOversized)) = avgRecordOversized
            call move_alloc(avgRecordBuffer, avgRecordOversized)

          end if

        end if

        if (oversized) then
          cIdxBackOversized(segCount) = cIdx
          lenBackOversized(segCount) = length
          vacBackOversized(segCount + 1) = hitVacuum
          rNormFltOversized(:, segCount) = rNormFlt
          r0NormFltOversized(:, segCount) = r0NormFlt
          muFltOversized(:, segCount) = muFlt
        else
          cIdxBack(segCount) = cIdx
          lenBack(segCount) = length
          vacBack(segCount + 1) = hitVacuum
          rNormFltRecord(:, segCount) = rNormFlt
          r0NormFltRecord(:, segCount) = r0NormFlt
          muFltRecord(:, segCount) = muFlt
        end if


        if (tally) then
          if (oversized) then
            !$omp simd
            do g = 1, self % nG
              avgRecordOversized((segCount - 1) * self % nG + g) = (flatQ0(g) + delta(g))
              deltaRecordOversized((segCount - 1) * self % nG + g) = delta(g) 
              xIncRecordOversized((segCount - 1) * self % nG + g) = xInc(g)
              yIncRecordOversized((segCount - 1) * self % nG + g) = yInc(g)
              zIncRecordOversized((segCount - 1) * self % nG + g) = zInc(g)
            end do
          else
            !$omp simd
            do g = 1, self % nG
              avgRecord((segCount - 1) * self % nG + g) = (flatQ0(g) + delta(g))
              deltaRecord((segCount - 1) * self % nG + g) = delta(g)  
              xIncRecord((segCount - 1) * self % nG + g) = xInc(g)
              yIncRecord((segCount - 1) * self % nG + g) = yInc(g)
              zIncRecord((segCount - 1) * self % nG + g) = zInc(g)
            end do
          end if

          centVec => self % centroidTracks((centIdx + 1):(centIdx + nDim))
          momVec => self % momTracks((momIdx + 1):(momIdx + matSize))
          ! Update centroid
          do g = 1, nDim
            !$omp atomic
            centVec(g) = centVec(g) + rC(g)
          end do

          ! Update spatial moment scores
          do g = 1, matSize
            !$omp atomic
            momVec(g) = momVec(g) + matScore(g)
          end do

          if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
        
        end if

      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

    end do

    ! Flux guess for new adjoint deadlength
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % adjSource(idx)
    end do

    !iterate over segments
    do i = segCount, 1, -1

      if (oversized) then
        lenFlt = real(lenBackOversized(i),defFlt)
        length = lenBackOversized(i)
        hitVacuum = vacBackOversized(i)
        cIdx = cIdxBackOversized(i)
        rNormFlt = rNormFltOversized(:,i)
        muFlt = muFltOversized(:,i)
        r0NormFlt = r0NormFltOversized(:,i)

        avgFluxVecFW => avgRecordOversized((i - 1) * self % nG + 1 : i * self % nG)
        deltaFW => deltaRecordOversized((i - 1) * self % nG + 1 : i * self % nG)
        xIncFW  => xIncRecordOversized((i - 1) * self % nG + 1 : i * self % nG)
        yIncFW  => yIncRecordOversized((i - 1) * self % nG + 1 : i * self % nG)
        zIncFW  => zIncRecordOversized((i - 1) * self % nG + 1 : i * self % nG)
      else
        lenFlt = real(lenBack(i),defFlt)
        length = lenBack(i)
        hitVacuum = vacBack(i)
        cIdx = cIdxBack(i)
        rNormFlt = rNormFltRecord(:,i)
        muFlt = muFltRecord(:,i)
        r0NormFlt = r0NormFltRecord(:,i)

        avgFluxVecFW => avgRecord((i - 1) * self % nG + 1 : i * self % nG)
        deltaFW => deltaRecord((i - 1) * self % nG + 1 : i * self % nG)
        xIncFW  => xIncRecord((i - 1) * self % nG + 1 : i * self % nG)
        yIncFW  => yIncRecord((i - 1) * self % nG + 1 : i * self % nG)
        zIncFW  => zIncRecord((i - 1) * self % nG + 1 : i * self % nG)
      end if

      lenFlt2_2 = lenFlt * lenFlt * one_two

      matIdx = self % geom % geom % graph % getMatFromUID(cIdx)
      baseIdx = (cIdx - 1) * self % nG
      segIdx =  (i - 1) * self % nG

      sourceVec => self % adjSource((baseIdx + 1):(baseIdx + self % nG))
      xGradVec  => self % adjSourceX((baseIdx + 1):(baseIdx + self % nG))
      yGradVec  => self % adjSourceY((baseIdx + 1):(baseIdx + self % nG))
      zGradVec  => self % adjSourceZ((baseIdx + 1):(baseIdx + self % nG))

      !sourceVecFW => self % source((baseIdx + 1):(baseIdx + self % nG))
      
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        ! Cache total cross section
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
      end if

      !$omp simd aligned(xGradVec, yGradVec, zGradVec)
      do g = 1, self % nG
        flatQ(g) = rNormFlt(x) * xGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(y) * yGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(z) * zGradVec(g)
        flatQ(g) = flatQ(g) + sourceVec(g)

        gradQ(g) = muFlt(x) * xGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(y) * yGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(z) * zGradVec(g)
      end do

      !$omp simd
      do g = 1, self % nG
        tau(g) = totVec(g) * lenFlt
        if (tau(g) < 1E-8_defFlt) then 
          tau(g) = 0.0_defFlt
        end if
      end do

      !$omp simd
      do g = 1, self % nG
        Gn(g) = expG(tau(g))
      end do

     !$omp simd
      do g = 1, self % nG
        F1(g)  = 1.0_defFlt - tau(g) * Gn(g) 
      end do

      !$omp simd
      do g = 1, self % nG
        F2(g) = (2.0_defFlt * Gn(g) - F1(g)) * lenFlt
      end do

      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - flatQ(g)) * F1(g) & 
                    - one_two * gradQ(g) * F2(g) 
      end do

      ! Intermediate flux variable creation or update
      !$omp simd
      do g = 1, self % nG
        fluxVec0(g) = fluxVec(g)
      end do

      ! Flux vector update
      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g) * tau(g) 
      end do

      if (i <= segCountCrit) then

        !$omp simd
        do g = 1, self % nG
          avgFluxVec(g) = (flatQ(g) + delta(g))
        end do

        lenFlt2_2 = lenFlt * lenFlt * one_two
        
        !$omp simd
        do g = 1, self % nG
          H(g) = ( F1(g) - Gn(g) ) 
        end do
      
        !$omp simd
        do g = 1, self % nG
            G1(g) = one_two - H(g)
        end do

        !$omp simd 
        do g = 1, self % nG
          G2(g) = expG2(tau(g)) 
        end do

        do g = 1, self % nG
          G1(g) = G1(g) * flatQ(g) * lenFlt
          G2(g) = G2(g) * gradQ(g) * lenFlt2_2
          H(g)  = H(g) * fluxVec0(g) * lenFlt
          H(g) = (G1(g) + G2(g) + H(g)) * lenFlt
          flatQ(g) = flatQ(g) * lenFlt + delta(g) * lenFlt
        end do

        !$omp simd
        do g = 1, self % nG
          xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) 
          yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) 
          zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g)
        end do

        call OMP_set_lock(self % locks(cIdx))

        scalarVec => self % adjScalarFlux((baseIdx + 1):(baseIdx + self % nG))
        xMomVec => self % adjScalarX((baseIdx + 1):(baseIdx + self % nG))
        yMomVec => self % adjScalarY((baseIdx + 1):(baseIdx + self % nG)) 
        zMomVec => self % adjScalarZ((baseIdx + 1):(baseIdx + self % nG))

        angularProd => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG )

        scalarVecFW => self % scalarFlux((baseIdx + 1):(baseIdx + self % nG))
        xMomVecFW => self % scalarX((baseIdx + 1):(baseIdx + self % nG))
        yMomVecFW => self % scalarY((baseIdx + 1):(baseIdx + self % nG))
        zMomVecFW => self % scalarZ((baseIdx + 1):(baseIdx + self % nG))  

        !$omp simd 
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + delta(g) * lenFlt
            xMomVec(g) = xMomVec(g) + xInc(g)
            yMomVec(g) = yMomVec(g) + yInc(g)
            zMomVec(g) = zMomVec(g) + zInc(g)

            scalarVecFW(g) = scalarVecFW(g) + deltaFW(g) * lenFlt
            xMomVecFW(g) = xMomVecFW(g) + xIncFW(g)
            yMomVecFW(g) = yMomVecFW(g) + yIncFW(g)
            zMomVecFW(g) = zMomVecFW(g) + zIncFW(g)
        end do
        if (isActive) then
          !$omp simd 
          do g = 1, self % nG
              do gIn = 1, self % nG
                pIdx = self % nG * (g - 1) + gIn
                angularProd(pIdx) = angularProd(pIdx) + avgFluxVec(g) * avgFluxVecFW(gIn)
              end do
          end do
        end if
        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
        if (segCount > self % NSegMax) then
          self % NSegMaxTemp = segCount
        end if
        call OMP_unset_lock(self % locks(cIdx))

      end if  

      if (hitVacuum) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

     end do

  end subroutine transportSweep

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defReal), save                           :: invVol
    real(defFlt), save                            :: total, vol, NTV, D, SigGG
    integer(shortInt)                             :: cIdx
    integer(shortInt), save                       :: g, matIdx, idx, dIdx, mIdx
    !$omp threadprivate(total, vol, idx, mIdx, dIdx, g, matIdx, invVol, NTV, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)
    self % NSegMax = self % NSegMaxTemp

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)
      dIdx = (cIdx - 1) * nDim
      mIdx = (cIdx - 1) * matSize
      
      ! Update volume
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx),defFlt)
      
      if (self % volume(cIdx) > volume_tolerance) then
        invVol = ONE / self % volumeTracks(cIdx)
        
        ! Update centroids
        self % centroid(dIdx + x) =  self % centroidTracks(dIdx + x) * invVol
        self % centroid(dIdx + y) =  self % centroidTracks(dIdx + y) * invVol
        self % centroid(dIdx + z) =  self % centroidTracks(dIdx + z) * invVol
      
        ! Update spatial moments
        self % momMat(mIdx + xx) = self % momTracks(mIdx + xx) * invVol
        self % momMat(mIdx + xy) = self % momTracks(mIdx + xy) * invVol
        self % momMat(mIdx + xz) = self % momTracks(mIdx + xz) * invVol
        self % momMat(mIdx + yy) = self % momTracks(mIdx + yy) * invVol
        self % momMat(mIdx + yz) = self % momTracks(mIdx + yz) * invVol
        self % momMat(mIdx + zz) = self % momTracks(mIdx + zz) * invVol

      else
        self % centroid(dIdx + x) =  ZERO
        self % centroid(dIdx + y) =  ZERO
        self % centroid(dIdx + z) =  ZERO

        self % momMat(mIdx + xx) = ZERO
        self % momMat(mIdx + xy) = ZERO
        self % momMat(mIdx + xz) = ZERO
        self % momMat(mIdx + yy) = ZERO
        self % momMat(mIdx + yz) = ZERO
        self % momMat(mIdx + zz) = ZERO

      end if

      do g = 1, self % nG

        idx   = self % nG * (cIdx - 1) + g
        NTV = norm / vol
        total = self % sigmaT((matIdx - 1) * self % nG + g)
        sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)

        if (vol > volume_tolerance) then
          self % scalarFlux(idx) =  self % scalarFlux(idx) * NTV
          self % scalarX(idx) = self % scalarX(idx) * NTV 
          self % scalarY(idx) = self % scalarY(idx) * NTV 
          self % scalarZ(idx) = self % scalarZ(idx) * NTV 
        end if

        ! Presumes non-zero total XS
        if ((sigGG < 0.0_defFlt) .and. (total > 0.0_defFlt)) then
          D = -real(self % rho, defFlt) * sigGG / total
        else
          D = 0.0_defFlt
        end if

        self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)
        self % scalarX(idx) =  (self % scalarX(idx) + D * self % prevX(idx) ) / (1 + D)
        self % scalarY(idx) =  (self % scalarY(idx) + D * self % prevY(idx) ) / (1 + D)
        self % scalarZ(idx) =  (self % scalarZ(idx) + D * self % prevZ(idx) ) / (1 + D)
      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume

  subroutine adjointNormaliseFluxAndVolume(self, it)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal), save                           :: invVol
    real(defFlt), save                            :: total, vol, NTV, D, SigGG
    integer(shortInt)                             :: cIdx
    integer(shortInt), save                       :: g, matIdx, idx, dIdx, mIdx
    !$omp threadprivate(total, vol, idx, mIdx, dIdx, g, matIdx, invVol, NTV, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)
      vol = real(self % volume(cIdx),defFlt)
      
      do g = 1, self % nG

        idx   = self % nG * (cIdx - 1) + g
        NTV = norm / vol
        total = self % sigmaT((matIdx - 1) * self % nG + g)
        sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)

        if (vol > volume_tolerance) then
          self % adjScalarFlux(idx) = self % adjScalarFlux(idx) * NTV
          self % adjScalarX(idx) = self % adjScalarX(idx) * NTV 
          self % adjScalarY(idx) = self % adjScalarY(idx) * NTV 
          self % adjScalarZ(idx) = self % adjScalarZ(idx) * NTV 
        end if

        ! Presumes non-zero total XS
        if ((sigGG < 0.0_defFlt) .and. (total > 0.0_defFlt)) then
          D = -real(self % rho, defFlt) * sigGG / total
        else
          D = 0.0_defFlt
        end if

        self % adjScalarFlux(idx) =  (self % adjScalarFlux(idx) + self % adjSource(idx) + D * self % adjPrevFlux(idx) ) / (1 + D)
        self % adjScalarX(idx) =  (self % adjScalarX(idx) + D * self % adjPrevX(idx) ) / (1 + D)
        self % adjScalarY(idx) =  (self % adjScalarY(idx) + D * self % adjPrevY(idx) ) / (1 + D)
        self % adjScalarZ(idx) =  (self % adjScalarZ(idx) + D * self % adjPrevZ(idx) ) / (1 + D)
      end do

    end do
    !$omp end parallel do

  end subroutine adjointNormaliseFluxAndVolume


  !!
  !! Normalise inner product by volume
  !!
  subroutine normaliseInnerProduct(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defFlt), save                            :: total, vol
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells

      do g = 1, self % nG * self % nG
        idx = (cIdx - 1) * self % nG * self % nG + g
        self % angularIP(idx) = self % angularIP(idx) * (self % volume(cIdx))
      end do

    end do

    !$omp end parallel do
  end subroutine normaliseInnerProduct
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF, it)
    class(adjointLSRRPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    integer(shortInt), intent(in)                         :: it
    real(defFlt)                                          :: scatter, xScatter, yScatter, zScatter, &
                                                             fission, xFission, yFission, zFission, &
                                                             xSource, ySource, zSource
    real(defFlt)                                          :: invMxx, invMxy, invMxz, invMyy, invMyz, invMzz
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec, xFluxVec, yFluxVec, zFluxVec
    real(defReal), pointer, dimension(:)                  :: momVec
    real(defReal)                                         :: det, one_det 
    
    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx)

    ! Guard against void cells
    if (matIdx >= VOID_MAT - 1) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
        self % sourceX(idx) = 0.0_defFlt
        self % sourceY(idx) = 0.0_defFlt
        self % sourceZ(idx) = 0.0_defFlt
      end do
      return
    end if

    momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize))

    det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
    - momVec(yy) * momVec(xz) * momVec(xz) - momVec(zz) * momVec(xy) * momVec(xy) &
    + 2 * momVec(xy) * momVec(xz) * momVec(yz)

    if ((abs(det) > volume_tolerance) .and. self % volume(cIdx) > 1E-10) then ! maybe: vary volume check depending on avg cell size..and. (self % volume(cIdx) > 1E-6)
      one_det = ONE/det
      invMxx = real(one_det * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)),defFlt)
      invMxy = real(one_det * (momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz)),defFlt)
      invMxz = real(one_det * (momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz)),defFlt)
      invMyy = real(one_det * (momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)),defFlt)
      invMyz = real(one_det * (momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz)),defFlt)
      invMzz = real(one_det * (momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)),defFlt)
    else
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt
      det = ONE 
    end if
     
    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx+1):(baseIdx + self % nG))
    xFluxVec => self % prevX((baseIdx + 1):(baseIdx + self % nG))
    yFluxVec => self % prevY((baseIdx + 1):(baseIdx + self % nG))
    zFluxVec => self % prevZ((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !$omp simd reduction(+:fission, xFission, yFission, zFission) aligned(fluxVec, xFluxVec, yFluxVec, zFluxVec, nuFission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
      xFission = xFission + xFluxVec(gIn) * nuFission(gIn)
      yFission = yFission + yFluxVec(gIn) * nuFission(gIn)
      zFission = zFission + zFluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG))

      ! Calculate scattering source
      scatter = 0.0_defFlt
      xScatter = 0.0_defFlt
      yScatter = 0.0_defFlt
      zScatter = 0.0_defFlt

      ! Sum contributions from all energies
      !$omp simd reduction(+:scatter, xScatter, yScatter, zScatter) aligned(fluxVec, xFluxVec, yFluxVec, zFluxVec, scatterVec)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        xScatter = xScatter + xFluxVec(gIn) * scatterVec(gIn)
        yScatter = yScatter + yFluxVec(gIn) * scatterVec(gIn)
        zScatter = zScatter + zFluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      self % source(idx) = self % source(idx) / total(g)
      xSource = chi(g) * xFission + xScatter
      xSource = xSource / total(g)
      ySource = chi(g) * yFission + yScatter
      ySource = ySource / total(g)
      zSource = chi(g) * zFission + zScatter
      zSource = zSource / total(g)
        
      
      if (it > 10) then
        self % sourceX(baseIdx + g) = invMxx * xSource + &
                invMxy * ySource + invMxz * zSource
        self % sourceY(baseIdx + g) = invMxy * xSource + &
                invMyy * ySource + invMyz * zSource
        ! self % sourceZ(baseIdx + g) = invMxz * xSource + &
        !      invMyz * ySource + invMzz * zSource
        self % sourceZ(baseIdx + g) = 0.0_defFlt
      else
        self % sourceX(baseIdx + g) = 0.0_defFlt
        self % sourceY(baseIdx + g) = 0.0_defFlt
        self % sourceZ(baseIdx + g) = 0.0_defFlt
      end if
    end do

  end subroutine sourceUpdateKernel


  subroutine adjointSourceUpdateKernel(self, cIdx, ONE_KEFF, it)
    class(adjointLSRRPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    integer(shortInt), intent(in)                         :: it
    real(defFlt)                                          :: scatter, xScatter, yScatter, zScatter, &
                                                             fission, xFission, yFission, zFission, &
                                                             xSource, ySource, zSource
    real(defFlt)                                          :: invMxx, invMxy, invMxz, invMyy, invMyz, invMzz
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec, xFluxVec, yFluxVec, zFluxVec
    real(defReal), pointer, dimension(:)                  :: momVec
    real(defReal)                                         :: det, one_det 
    
    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx)

    ! Guard against void cells
    if (matIdx >= VOID_MAT - 1) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % adjSource(idx) = 0.0_defFlt
        self % adjSourceX(idx) = 0.0_defFlt
        self % adjSourceY(idx) = 0.0_defFlt
        self % adjSourceZ(idx) = 0.0_defFlt
      end do
      return
    end if

    momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize))

    det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
    - momVec(yy) * momVec(xz) * momVec(xz) - momVec(zz) * momVec(xy) * momVec(xy) &
    + 2 * momVec(xy) * momVec(xz) * momVec(yz)

    if ((abs(det) > volume_tolerance) .and. self % volume(cIdx) > 1E-10) then ! maybe: vary volume check depending on avg cell size..and. (self % volume(cIdx) > 1E-6)
      one_det = ONE/det
      invMxx = real(one_det * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)),defFlt)
      invMxy = real(one_det * (momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz)),defFlt)
      invMxz = real(one_det * (momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz)),defFlt)
      invMyy = real(one_det * (momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)),defFlt)
      invMyz = real(one_det * (momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz)),defFlt)
      invMzz = real(one_det * (momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)),defFlt)
    else
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt
      det = ONE 
    end if
     
    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % adjSigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
    nuFission => self % adjNuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % adjChi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % adjPrevFlux((baseIdx+1):(baseIdx + self % nG))
    xFluxVec => self % adjPrevX((baseIdx + 1):(baseIdx + self % nG))
    yFluxVec => self % adjPrevY((baseIdx + 1):(baseIdx + self % nG))
    zFluxVec => self % adjPrevZ((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !$omp simd reduction(+:fission, xFission, yFission, zFission) aligned(fluxVec, xFluxVec, yFluxVec, zFluxVec, nuFission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
      xFission = xFission + xFluxVec(gIn) * nuFission(gIn)
      yFission = yFission + yFluxVec(gIn) * nuFission(gIn)
      zFission = zFission + zFluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG))

      ! Calculate scattering source
      scatter = 0.0_defFlt
      xScatter = 0.0_defFlt
      yScatter = 0.0_defFlt
      zScatter = 0.0_defFlt

      ! Sum contributions from all energies
      !$omp simd reduction(+:scatter, xScatter, yScatter, zScatter) aligned(fluxVec, xFluxVec, yFluxVec, zFluxVec, scatterVec)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        xScatter = xScatter + xFluxVec(gIn) * scatterVec(gIn)
        yScatter = yScatter + yFluxVec(gIn) * scatterVec(gIn)
        zScatter = zScatter + zFluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % adjSource(idx) = chi(g) * fission + scatter
      self % adjSource(idx) = self % adjSource(idx) / total(g)
      xSource = chi(g) * xFission + xScatter
      xSource = xSource / total(g)
      ySource = chi(g) * yFission + yScatter
      ySource = ySource / total(g)
      zSource = chi(g) * zFission + zScatter
      zSource = zSource / total(g)
        
      if (it > 10) then
        !Calculate source gradients by inverting the moment matrix
        self % adjSourceX(baseIdx + g) = invMxx * xSource + &
                invMxy * ySource + invMxz * zSource
        self % adjSourceY(baseIdx + g) = invMxy * xSource + &
                invMyy * ySource + invMyz * zSource
        ! self % adjSourceZ(baseIdx + g) = invMxz * xSource + &
        !      invMyz * ySource + invMzz * zSource
        self % adjSourceZ(baseIdx + g) = 0.0_defFlt
      else
        self % adjSourceX(baseIdx + g) = 0.0_defFlt
        self % adjSourceY(baseIdx + g) = 0.0_defFlt
        self % adjSourceZ(baseIdx + g) = 0.0_defFlt
      end if
    end do

  end subroutine adjointSourceUpdateKernel

  subroutine sensitivityCalculation(self, cIdx, ONE_KEFF, it, numSum, denSum)
    class(adjointLSRRPhysicsPackage), target, intent(inout) :: self
    real(defFlt), intent(in)                      :: ONE_KEFF
    integer(shortInt), intent(in)                 :: it
    integer(shortInt), intent(in)                 :: cIdx
    real(defReal), intent(out)                    :: numSum, denSum
    real(defReal)                                 :: delta, fission, fission_pert, scatter_pert
    integer(shortInt)                             :: baseIdx, idx, matIdx, g, gIn, mat, g1Pert, g2pert, i
    real(defFlt), dimension(:), pointer           :: nuFission, total, chi, capture, scatterXS, &
                                                      fissVec, scatterVec
    real(defReal), dimension(:), pointer          :: IPVec

    mat  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    baseIdx = (cIdx - 1) * self % nG 
    matIdx = (mat - 1) * self % nG

    IPVec => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG)
    ! total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    fissVec => self % fission((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))
    capture => self % sigmaC((matIdx + 1):(matIdx + self % nG)) 
    scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))

    numSum = ZERO
    denSum = ZERO

    do g = 1, self % nG
      idx = baseIdx + g
      fission = 0.0_defFlt
      !$omp simd
      do gIn = 1, self % nG 
        fission = fission + IPVec((g - 1) * self % nG + gIn) * nuFission(gIn)
      end do
      fission  = fission * chi(g)
      denSum = denSum + fission
    end do

    if (mat == self % matPert .and. self % XStype == 1) then ! capture - complete 
      do i = 1, size(self % energyId)
        g1Pert = self % energyId(1)
        delta = 0.0_defFlt
        delta = self % XSchange * capture(g1Pert) * IPVec((g1Pert - 1) * self % nG + g1Pert)
        numSum = numSum - delta
      end do

    elseif ( mat == self % matPert .and. self % XStype == 2) then ! fission - complete
      do i = 1, size(self % energyId)
        g1Pert = self % energyId(i)
        do g = 1, self % nG 
          delta = 0.0_defFlt
          fission_pert = 0.0_defFlt
          fission_pert = fission_pert + IPVec(self % nG * (g - 1) + g1Pert) * nuFission(g1Pert) * self % XSchange
          if ( g == g1Pert ) then 
            delta = IPVec((g - 1) * self % nG + g) * fissVec(g) * self % XSchange
          end if
          delta = fission_pert * chi(g) * ONE_KEFF - delta
          numSum = numSum + delta
        end do
      end do

    elseif (mat == self % matPert .and. self % XStype == 3) then !scatter - complete
      do i = 1, size(self % energyId), 2
        g1Pert = self % energyId(i)
        g2Pert = self % energyId(i+1)
        delta = IPVec((g1Pert - 1) * self % nG + g1Pert) * scatterXS((g2Pert - 1) * self % nG + g1Pert) * &
                      self % XSchange   
           
        scatter_pert = scatterXS( (g2Pert - 1) * self % nG + g1Pert ) * &
                                    IPVec((g2Pert - 1 ) * self % nG + g1Pert) * self % XSchange
        delta = scatter_pert - delta
        numSum = numSum + delta
      end do
    end if


  end subroutine sensitivityCalculation

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    real(defFlt)                                  :: fissionRate, prevFissionRate,adjfissionRate, adjprevFissionRate
    real(defFlt), save                            :: fissLocal, prevFissLocal, vol, adjfissLocal, adjprevFissLocal
    integer(shortInt), save                       :: matIdx, g, idx, mIdx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, matIdx, g, idx, mIdx, vol)

    fissionRate        = 0.0_defFlt
    prevFissionRate    = 0.0_defFlt
    adjFissionRate     = 0.0_defFlt
    adjPrevFissionRate = 0.0_defFlt
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate, adjFissionRate, adjPrevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      if (matIdx >= VOID_MAT) cycle

      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      if (.not. mat % isFissile()) cycle

      vol = real(self % volume(cIdx), defFlt)

      if (vol <= volume_tolerance) cycle

      fissLocal = 0.0_defFlt
      prevFissLocal = 0.0_defFlt
      adjfissLocal = 0.0_defFlt
      adjprevFissLocal = 0.0_defFlt
      mIdx = (matIdx - 1) * self % nG
      do g = 1, self % nG
        
        ! Source index
        idx = self % nG * (cIdx - 1) + g
        fissLocal     = fissLocal     + self % scalarFlux(idx) * self % nuSigmaF(mIdx + g)
        prevFissLocal = prevFissLocal + self % prevFlux(idx) * self % nuSigmaF(mIdx + g)

        adjfissLocal     = adjfissLocal     + self % adjscalarFlux(idx) * self % adjNuSigmaF(mIdx + g)
        adjprevFissLocal = adjprevFissLocal + self % adjprevFlux(idx) * self % adjNuSigmaF(mIdx + g)

      end do

      fissionRate     = fissionRate     + fissLocal * vol
      prevFissionRate = prevFissionRate + prevFissLocal * vol

      adjfissionRate     = adjfissionRate     + adjfissLocal * vol
      adjprevFissionRate = adjprevFissionRate + adjprevFissLocal * vol

    end do
    !$omp end parallel do

    ! Update k
    self % keff = self % keff * fissionRate / prevFissionRate

    self % adjkeff = self % adjkeff * adjfissionRate / adjprevFissionRate

  end subroutine calculateKeff
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! Same for flux moments
  !!
  subroutine resetFluxes(self)
    class(adjointLSRRPhysicsPackage), intent(inout)    :: self
    integer(shortInt)                                  :: idx, g

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
      self % prevX(idx) = self % scalarX(idx)
      self % scalarX(idx) = 0.0_defFlt
      self % prevY(idx) = self % scalarY(idx)
      self % scalarY(idx) = 0.0_defFlt
      self % prevZ(idx) = self % scalarZ(idx)
      self % scalarZ(idx) = 0.0_defFlt

      self % adjPrevFlux(idx) = self % adjScalarFlux(idx)
      self % adjScalarFlux(idx) = 0.0_defFlt
      self % adjPrevX(idx) = self % adjScalarX(idx)
      self % adjScalarX(idx) = 0.0_defFlt
      self % adjPrevY(idx) = self % adjScalarY(idx)
      self % adjScalarY(idx) = 0.0_defFlt
      self % adjPrevZ(idx) = self % adjScalarZ(idx)
      self % adjScalarZ(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

    !$omp parallel do schedule(static)
    do g = 1, size(self % angularIP)
      self % angularIP(g) = 0.0_defFlt
    end do
    !$omp end parallel do
    
  end subroutine resetFluxes

  subroutine accumulateFluxAndKeffScores(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    real(defReal), save                            :: flux
    integer(shortInt)                              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux * flux
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % adjkeffScore(1) = self % adjkeffScore(1) + self % adjkeff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff
    self % adjkeffScore(2) = self % adjkeffScore(2) + self % adjkeff * self % adjkeff

    self % sensitivityScore(1) = self % sensitivityScore(1) + self % sensitivity
    self % sensitivityScore(2) = self % sensitivityScore(2) + self % sensitivity * self % sensitivity

    self % deltaKeffScore(1) = self % deltaKeffScore(1) + self % deltaKeff
    self % deltaKeffScore(2) = self % deltaKeffScore(2) + self % deltaKeff * self % deltaKeff

  end subroutine accumulateFluxAndKeffScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt)                             :: idx
    real(defReal)                                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) * N1
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) * N1
      self % fluxScores(idx,2) = Nm1 *(self % fluxScores(idx,2) - &
            self % fluxScores(idx,1) * self % fluxScores(idx,1)) 
      if (self % fluxScores(idx,2) <= ZERO) then
        self % fluxScores(idx,2) = ZERO
      else
        self % fluxScores(idx,2) = sqrt(self % fluxScores(idx,2))
      end if
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) * N1
    self % keffScore(2) = self % keffScore(2) * N1
    self % keffScore(2) = (Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 
    self % keffScore(2) = sqrt(abs(self % keffScore(2)))

    self % adjkeffScore(1) = self % adjkeffScore(1) * N1
    self % adjkeffScore(2) = self % adjkeffScore(2) * N1
    self % adjkeffScore(2) = (Nm1*(self % adjkeffScore(2) - &
            self % adjkeffScore(1) * self % adjkeffScore(1))) 
    self % adjkeffScore(2) = sqrt(abs(self % adjkeffScore(2)))

    self % sensitivityScore(1) = self % sensitivityScore(1) * N1
    self % sensitivityScore(2) = self % sensitivityScore(2) * N1
    self % sensitivityScore(2) = (Nm1*(self % sensitivityScore(2) - &
            self % sensitivityScore(1) * self % sensitivityScore(1))) 
    self % sensitivityScore(2) = sqrt(abs(self % sensitivityScore(2)))

    self % deltaKeffScore(1) = self % deltaKeffScore(1) * N1
    self % deltaKeffScore(2) = self % deltaKeffScore(2) * N1
    self % deltaKeffScore(2) = (Nm1*(self % deltaKeffScore(2) - &
            self % deltaKeffScore(1) * self % deltaKeffScore(1))) 
    self % deltaKeffScore(2) = sqrt(abs(self % deltaKeffScore(2)))

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    integer(shortInt)                             :: cIdx, g1
    integer(shortInt), save                       :: idx, matIdx, i, g
    real(defReal), save                           :: vol, SigmaF
    type(particleState), save                     :: s
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: groupFlux, fiss, fissSTD
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g)

    call out % init(self % outputFormat)
    
    name = 'seed'
    call out % printValue(self % rand % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)
    
    name = 'Time_Per_Integration'
    call out % printValue(self % time_transport/(self % intersectionsTotal * self % nG),name)
    
    name = 'Clock_Time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print keff
    name = 'keff'
    call out % startBlock(name)
    call out % printResult(real(self % keffScore(1),defReal), real(self % keffScore(2),defReal), name)
    call out % endBlock()

    ! Print keff
    name = 'Adjoint_keff'
    call out % startBlock(name)
    call out % printResult(real(self % adjkeffScore(1),defReal), real(self % adjkeffScore(2),defReal), name)
    call out % endBlock()

    name = 'delta_keff'
    call out % startBlock(name)
    call out % printResult(real(self % deltaKeffScore(1),defReal), &
          real(self % deltaKeffScore(2),defReal), name)
    call out % endBlock()

    name = 'pert_keff'
    call out % startBlock(name)
    call out % printResult(real(self % keffScore(1) + self % deltaKeffScore(1),defReal), &
          real(self % keffScore(1) - self % deltaKeffScore(1),defReal), name)
    call out % endBlock()

    name = 'sensitivity'
    call out % startBlock(name)
    call out % printResult(real(self % sensitivityScore(1),defReal), real(self % sensitivityScore(2)&
                                        /self % sensitivityScore(1),defReal), name)
    call out % endBlock()

    ! Print cell volumes
    if (self % printVolume) then
      name = 'volume'
      call out % startBlock(name)
      resArrayShape = [size(self % volume)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(self % volume(cIdx), ZERO)
      end do
      call out % endArray()
      call out % endBlock()
    end if

    ! Print cell positions
    if (self % printCells) then
      name = 'position'
      call out % startBlock(name)
      resArrayShape = [size(self % cellPos)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(self % cellPos(cIdx,1), ZERO)
        call out % addResult(self % cellPos(cIdx,2), ZERO)
        call out % addResult(self % cellPos(cIdx,3), ZERO)
      end do
      call out % endArray()
      call out % endBlock()
    end if

    ! Print fluxes
    if (self % printFlux) then
      resArrayShape = [size(self % volume)]
      do g = 1, self % nG
        name = 'flux_g'//numToChar(g)
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          call out % addResult(self % fluxScores(idx,1), self % fluxScores(idx,2))
        end do
        call out % endArray()
        call out % endBlock()
      end do
    end if

    ! Send fission rates to map output
    if (self % mapFission) then
      resArrayShape = self % resultsMap % binArrayShape()
      allocate(fiss(self % resultsMap % bins(0)))
      allocate(fissSTD(self % resultsMap % bins(0)))
      fiss    = ZERO
      fissSTD = ZERO

      ! Find whether cells are in map and sum their contributions
      !$omp parallel do reduction(+: fiss, fissSTD)
      do cIdx = 1, self % nCells
        
        ! Identify material
        matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
        matPtr => self % mgData % getMaterial(matIdx)
        mat    => baseMgNeutronMaterial_CptrCast(matPtr)
        vol    =  self % volume(cIdx)

        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % resultsMap % map(s)

        if (i > 0) then
          do g = 1, self % nG
            SigmaF = mat % getFissionXS(g, self % rand)
            idx = (cIdx - 1)* self % nG + g
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1) * SigmaF
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    vol * vol * self % fluxScores(idx,2)*self % fluxScores(idx,2) * SigmaF * SigmaF
          end do
        end if

      end do
      !$omp end parallel do

      do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
      end do

      name = 'fissionRate'
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      ! Add all map elements to results
      do idx = 1, self % resultsMap % bins(0)
        call out % addResult(fiss(idx), fissSTD(idx))
      end do
      call out % endArray()
      ! Output tally map
      call self % resultsMap % print(out)
      call out % endBlock()
      
      deallocate(fiss)
      deallocate(fissSTD)
    end if

    ! Send fluxes to map output
    ! Will the results still be correct with linear source?
    ! Is some extra logic required?
    if (self % mapFlux) then
      resArrayShape = self % fluxMap % binArrayShape()
      allocate(fiss(self % fluxMap % bins(0)))
      allocate(fissSTD(self % fluxMap % bins(0)))
      fiss    = ZERO
      fissSTD = ZERO

      ! Find whether cells are in map and sum their contributions
      !$omp parallel do reduction(+: fiss, fissSTD)
      do cIdx = 1, self % nCells
        
        vol    =  self % volume(cIdx)
        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % fluxMap % map(s)

        if (i > 0) then
          do g = 1, self % nG
            idx = (cIdx - 1)* self % nG + g
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1)
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    self % fluxScores(idx,2)*self % fluxScores(idx,2) * vol * vol
          end do
        end if

      end do
      !$omp end parallel do

      do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
      end do

      name = 'flux1G'
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      ! Add all map elements to results
      do idx = 1, self % fluxMap % bins(0)
        call out % addResult(fiss(idx), fissSTD(idx))
      end do
      call out % endArray()
      ! Output tally map
      call self % fluxMap % print(out)
      call out % endBlock()
      
      deallocate(fiss)
      deallocate(fissSTD)
    end if
    
    call out % writeToFile(self % outputFile)

    ! Send all fluxes and stds to VTK
    ! TODO: can sources be interpolated more finely using
    !       moment and centroid info?
    if (self % plotResults) then
      allocate(groupFlux(self % nCells))
      do g1 = 1, self % nG
        name = 'flux_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'std_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,2) /self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      call self % viz % finaliseVTK
    end if


  end subroutine printResults

  !!
  !! Print settings of the random ray calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(adjointLSRRPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ LINEAR SOURCE RANDOM RAY EIGENVALUE CALCULATION /\/\"
    print *, "Using "//numToChar(self % inactive)// " iterations for "&
              //"the inactive cycles"
    print *, "Using "//numToChar(self % active)// " iterations for "&
              //"the active cycles"
    print * 
    print *, "Rays per cycle: "// numToChar(self % pop)
    print *, "Ray dead length: "//numToChar(self % dead)
    print *, "Ray termination length: "//numToChar(self % termination)
    print *, "Initial RNG Seed:   "// numToChar(self % rand % getSeed())
    print *
    print *, "Number of cells in the geometry: "// numToChar(self % nCells)
    print *, "Number of energy groups: "// numToChar(self % nG)
    if (self % cache) print *, "Accelerated with distance caching"
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(adjointLSRRPhysicsPackage), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean Nuclear Data, Geometry and visualisation
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0
    self % timerTransport = 0

    self % top       = ZERO
    self % bottom    = ZERO
    self % mgData    => null()
    self % nG        = 0
    
    if(allocated(self % locks)) then
      do i = 1, self % nCells
        call OMP_destroy_lock(self % locks(i))
      end do
      deallocate(self % locks)
    end if
    self % nCells    = 0
    self % nMat      = 0
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nusigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % chi)) deallocate(self % chi)

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    self % rho         = ZERO
    self % cache       = .false.
    self % mapFission  = .false.
    self % mapFlux     = .false.
    self % plotResults = .false.
    self % printFlux   = .false.
    self % printVolume = .false.
    self % printCells  = .false.

    self % keff        = ZERO
    self % keffScore   = ZERO
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % scalarX)) deallocate(self % scalarX)
    if(allocated(self % scalarY)) deallocate(self % scalarY)
    if(allocated(self % scalarZ)) deallocate(self % scalarZ)
    if(allocated(self % prevX)) deallocate(self % prevX)
    if(allocated(self % prevY)) deallocate(self % prevY)
    if(allocated(self % prevZ)) deallocate(self % prevZ)
    if(allocated(self % sourceX)) deallocate(self % sourceX)
    if(allocated(self % sourceY)) deallocate(self % sourceY)
    if(allocated(self % sourceZ)) deallocate(self % sourceZ)
    if(allocated(self % momMat)) deallocate(self % momMat)
    if(allocated(self % momTracks)) deallocate(self % momTracks)
    if(allocated(self % centroid)) deallocate(self % centroid)
    if(allocated(self % centroidTracks)) deallocate(self % centroidTracks)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % resultsMap)) then

    if(allocated(self % scalarX)) deallocate(self % adjScalarX)
    if(allocated(self % scalarY)) deallocate(self % adjScalarY)
    if(allocated(self % scalarZ)) deallocate(self % adjScalarZ)
    if(allocated(self % prevX)) deallocate(self % adjPrevX)
    if(allocated(self % prevY)) deallocate(self % adjPrevY)
    if(allocated(self % prevZ)) deallocate(self % adjPrevZ)
    if(allocated(self % sourceX)) deallocate(self % adjSourceX)
    if(allocated(self % sourceY)) deallocate(self % adjSourceY)
    if(allocated(self % sourceZ)) deallocate(self % adjSourceZ)

    if(allocated(self % adjScalarX)) deallocate(self % adjScalarX)
    if(allocated(self % adjScalarY)) deallocate(self % adjScalarY)
    if(allocated(self % adjScalarZ)) deallocate(self % adjScalarZ)

    if(allocated(self % adjPrevX)) deallocate(self % adjPrevX)
    if(allocated(self % adjPrevY)) deallocate(self % adjPrevY)
    if(allocated(self % adjPrevZ)) deallocate(self % adjPrevZ)

    if(allocated(self % adjsourceX)) deallocate(self % adjsourceX)
    if(allocated(self % adjsourceY)) deallocate(self % adjsourceY)
    if(allocated(self % adjsourceZ)) deallocate(self % adjsourceZ)



      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if
    if(allocated(self % fluxMap)) then
      call self % resultsMap % kill()
      deallocate(self % fluxMap)
    end if

  end subroutine kill

end module adjointLSRRPhysicsPackage_class
