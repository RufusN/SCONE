module linearAnisotropicScattering_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential
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
  real(defReal), parameter :: volume_tolerance = 1.0E-10, &
                                SQRT3 = sqrt(3._defReal), &
                                SQRT5_8 = sqrt(5._defReal/8._defReal), &
                                SQRT3_8 = sqrt(3._defReal/8._defReal), &
                                SQRT15  = sqrt(15._defReal)

  ! Parameters for indexing into matrices and spatial moments
  integer(shortInt), parameter :: x = 1, y = 2, z = 3, nDim = 3, &
                                  xx = 1, xy = 2, xz = 3, &
                                  yy = 4, yz = 5, zz = 6, &
                                  matSize = 6

  ! Parameters for deciding how to invert the moment matrix
  integer(shortInt), parameter :: invertXYZ = 7, invertXY = 6, &
                                  invertXZ = 5, invertYZ = 3, &
                                  invertX = 4, invertY = 2, &
                                  invertZ = 1

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
  !!     type linearAnisotropicScattering;
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
  type, public, extends(physicsPackage) :: linearAnisotropicScattering
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
    integer(shortInt)  :: harmonicOrder = ZERO
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

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi

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

    !real(defReal), dimension(:), allocatable   :: sphericalHarmonic

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
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings
    procedure, private :: sphericalHarmonicCalculator

  end type linearAnisotropicScattering

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(linearAnisotropicScattering), intent(inout) :: self
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
    character(100), parameter :: Here = 'init (linearAnisotropicScattering_class.f90)'

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

    ! Stabilisation factor for negative in-group scattering
    call dict % getOrDefault(self % harmonicOrder, 'harmonicOrder', ZERO)

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

    !allocate(self % sphericalHarmonic((harmonicOrder + 1) ** 2 ))
    
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
    allocate(self % nuSigmaF(self % nMat * self % nG))
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG))

    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, self % rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, self % rand),defFlt)
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, self % rand),defFlt)
        ! Include scattering multiplicity in sigmaS directly
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                  real(mat % getScatterXS(g1, g, self % rand) * mat % scatter % prod(g1, g) , defFlt)
        end do
      end do
    end do

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(linearAnisotropicScattering), intent(inout) :: self

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
    class(linearAnisotropicScattering), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints)

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
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        call self % sourceUpdateKernel(i, ONE_KEFF)
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
        call self % transportSweep(r,ints, it)
        intersections = intersections + ints

      end do
      !$omp end parallel do

      self % intersectionsTotal = self % intersectionsTotal + intersections
      
      call timerStop(self % timerTransport)

      ! Update RNG on master thread
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)

      ! Calculate new k
      call self % calculateKeff()

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxAndKeffScores()

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
    class(linearAnisotropicScattering), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, coord
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (linearAnisotropicScattering_class.f90)'

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
  subroutine transportSweep(self, r, ints, it)
    class(linearAnisotropicScattering), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: it
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, &
                                                             matIdx0, baseIdx, centIdx, momIdx
    real(defReal)                                         :: totalLength, length, len2_12
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt, lenFlt2_2
    real(defFlt), dimension(self % nG)                    :: F1, F2, G1, G2, H, tau, delta, fluxVec, &
                                                             flatQ, gradQ, xInc, yInc, zInc, fluxVec0
    real(defReal), dimension(matSize)                     :: matScore
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec, &
                                                             xGradVec, yGradVec, zGradVec, &
                                                             xMomVec, yMomVec, zMomVec
    real(defReal), pointer, dimension(:)                  :: mid, momVec, centVec
    real(defReal), pointer                                :: volTrack
    real(defReal), dimension(3)                           :: r0, mu0, rC, r0Norm, rNorm
    real(defFlt), dimension(3)                            :: muFlt, rNormFlt, r0NormFlt
    
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
    totalLength = ZERO
    activeRay = .false.
    do while (totalLength < self % termination)

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
        self % cellPos(cIdx,x) = rC(x)
        self % cellPos(cIdx,y) = rC(y)
        self % cellPos(cIdx,z) = rC(z)
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

      ! Calculate source terms
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

      lenFlt = real(length,defFlt)
      
      ! Compute exponentials necessary for angular flux update
      !$omp simd
      do g = 1, self % nG
        tau(g) = totVec(g) * lenFlt
      end do

      !$omp simd
      do g = 1, self % nG
        F1(g)  = exponential(tau(g))
      end do

      !$omp simd
      do g = 1, self % nG
        F2(g)  = 2 * (tau(g) - F1(g)) - tau(g) * F1(g)
      end do
      
      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - flatQ(g)) * F1(g) - &
                one_two * gradQ(g) * F2(g) / totVec(g)
      end do

      ! Create an intermediate flux variable for use in LS scores
      !$omp simd
      do g = 1, self % nG
        fluxVec0(g) = fluxVec(g)
      end do
      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g)
      end do

      ! Accumulate to scalar flux
      if (activeRay) then
      
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
          if (length > 1E-6) then
            G1(g) = 1 + one_two * tau(g) - (1 + 1.0_defFlt/tau(g)) * F1(g)
          else
            G1(g) = tau(g)*tau(g)/3
          end if
        end do

        !$omp simd
        do g = 1, self % nG
          if (length > 1E-6) then
            G2(g) = (two_three * tau(g) - (1 + 2.0_defFlt/tau(g)) * G1(g) )
          else
            G2(g) = -tau(g)*tau(g)/12
          end if
        end do
     
        !$omp simd
        do g = 1, self % nG
          H(g) = (one_two * tau(g) - G1(g))
        end do

        ! Make some more condensed variables to help vectorisation
        !$omp simd 
        do g = 1, self % nG
          G1(g) = G1(g) * flatQ(g) * lenFlt
          G2(g) = G2(g) * gradQ(g) * lenFlt2_2
          H(g)  = H(g) * fluxVec0(g) * lenFlt
          H(g) = G1(g) + G2(g) + H(g)
          flatQ(g) = flatQ(g) * tau(g) + delta(g)
        end do
        
        ! Attempt based on Fitzgerald... Going to be very explicitly written
        !!$omp simd 
        !do g = 1, self % nG
        !  xInc(g) = real(r0Norm(x),defFlt) * (delta(g) + totVec(g) * flatQ(g) * lenFlt) + &
        !        real(mu0(x),defFlt) * (G1(g) * flatQ(g) * lenFlt + gradQ(g) * G2(g) * lenFlt2_2 + &
        !       lenFlt * fluxVec0(g) * H(g)) 
        !  yInc(g) = real(r0Norm(y),defFlt) * (delta(g) + totVec(g) * flatQ(g) * lenFlt) + &
        !        real(mu0(y),defFlt) * (G1(g) * flatQ(g) * lenFlt + gradQ(g) * G2(g) * lenFlt2_2 + &
        !       lenFlt * fluxVec0(g) * H(g)) 
        !  zInc(g) = real(r0Norm(z),defFlt) * (delta(g) + totVec(g) * flatQ(g) * lenFlt) + &
        !       real(mu0(z),defFlt) * (G1(g) * flatQ(g) * lenFlt + gradQ(g) * G2(g) * lenFlt2_2 + &
        !       lenFlt * fluxVec0(g) * H(g))
        !end do
        !$omp simd
        do g = 1, self % nG
          xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) 
          yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) 
          zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g) 
        end do
      

        call OMP_set_lock(self % locks(cIdx))
        
        scalarVec => self % scalarFlux((baseIdx + 1):(baseIdx + self % nG))
        xMomVec => self % scalarX((baseIdx + 1):(baseIdx + self % nG))
        yMomVec => self % scalarY((baseIdx + 1):(baseIdx + self % nG))
        zMomVec => self % scalarZ((baseIdx + 1):(baseIdx + self % nG))

        ! Update flux moments
        !$omp simd aligned(scalarVec, xMomVec, yMomVec, zMomVec)
        do g = 1, self % nG
          scalarVec(g) = scalarVec(g) + delta(g) 
          xMomVec(g) = xMomVec(g) + xInc(g) 
          yMomVec(g) = yMomVec(g) + yInc(g)
          zMomVec(g) = zMomVec(g) + zInc(g) 
        end do
        
        centVec => self % centroidTracks((centIdx + 1):(centIdx + nDim))
        momVec => self % momTracks((momIdx + 1):(momIdx + matSize))
        volTrack => self % volumeTracks(cIdx)
        
        ! Update track lengths
        !self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
        volTrack = volTrack + length
        
        ! Update centroid
        !$omp simd aligned(centVec)
        do g = 1, nDim
          centVec(g) = centVec(g) + rC(g)
        end do

        ! Update spatial moment scores
        !$omp simd aligned(momVec)
        do g = 1, matSize
          momVec(g) = momVec(g) + matScore(g)
        end do
        
        call OMP_unset_lock(self % locks(cIdx))

        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
      
      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

    end do

  end subroutine transportSweep

  subroutine SphericalHarmonicCalculator(self,mu)
    !! angle: x = r sin θ cos φ, y = r sin θ sin φ, and z = r cos θ
    !! rational approximation of angles? 
    class(linearAnisotropicScattering), intent(inout) :: self
    real(defReal), dimension(self % harmonicOrder)    :: RCoeffs ! Array to store harmonic coefficients
    real(defReal)                                     :: X,Y,Z


    X = mu(1)
    Y = mu(2)
    Z = mu(3)

    ! Using a select case structure to assign coefficients based on harmonic order
    select case(self % harmonicOrder)
    case(1)  
        RCoeffs(1) = 1  

    case(4)  
        RCoeffs(1) = 1
        RCoeffs(2) = X 
        RCoeffs(3) = Y
        RCoeffs(4) = Z

    case(9) 
        RCoeffs(1) = 1
        RCoeffs(2) = X
        RCoeffs(3) = Y
        RCoeffs(4) = Z
        RCoeffs(5) = SQRT3 * Z * X
        RCoeffs(6) = SQRT3 * Y * X
        RCoeffs(7) = one_two * (3 * Y**2 - 1)
        RCoeffs(8) = SQRT3 * Z * Y
        RCoeffs(9) = (Z**2 - X**2) * SQRT3 / 2

    case(16) 
        RCoeffs(1) = 1
        RCoeffs(2) = X
        RCoeffs(3) = Y
        RCoeffs(4) = Z
        RCoeffs(5) = SQRT3 * Z * X
        RCoeffs(6) = SQRT3 * Y * X
        RCoeffs(7) = one_two * (3 * Y**2 - 1)
        RCoeffs(8) = SQRT3 * Z * Y
        RCoeffs(9) = (Z**2 - X**2) * SQRT3 / 2
        RCoeffs(10) = SQRT5_8 * X * (3 * Z**2 - X**2)
        RCoeffs(11) = SQRT15 * Y * Z * X
        RCoeffs(12) = SQRT3_8 * X * (5 * Y**2 - 1)
        RCoeffs(13) = one_two * Y * (5 * Y**2 - 3)
        RCoeffs(14) = SQRT3_8 * Z (5 * Y**2 - 1)
        RCoeffs(15) = Y * (Z**2 - X**2) * SQRT15 / 2
        RCoeffs(16) = SQRT5_8 * Z * (Z**2 - 3 * X**2)
    end select
  end subroutine SphericalHarmonicCalculator

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(linearAnisotropicScattering), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defReal), save                           :: invVol
    real(defFlt), save                            :: total, vol, NTV
    integer(shortInt)                             :: cIdx
    integer(shortInt), save                       :: g, matIdx, idx, dIdx, mIdx
    !$omp threadprivate(total, vol, idx, mIdx, dIdx, g, matIdx, invVol, NTV)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

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

        !removed - else set to zeros
      else
                ! Update centroids
        self % centroidTracks(dIdx + x) = ZERO
        self % centroidTracks(dIdx + y) = ZERO
        self % centroidTracks(dIdx + z) = ZERO
      
                ! Update spatial moments
        self % momTracks(mIdx + xx) = ZERO
        self % momTracks(mIdx + xy) = ZERO
        self % momTracks(mIdx + xz) = ZERO
        self % momTracks(mIdx + yy) = ZERO
        self % momTracks(mIdx + yz) = ZERO
        self % momTracks(mIdx + zz) = ZERO
      end if

      do g = 1, self % nG

        total = self % sigmaT((matIdx - 1) * self % nG + g)
        idx   = self % nG * (cIdx - 1) + g
        NTV = norm / (total * vol)
        
        if (vol > volume_tolerance) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * NTV
          self % scalarX(idx) = self % scalarX(idx) * NTV 
          self % scalarY(idx) = self % scalarY(idx) * NTV
          self % scalarZ(idx) = self % scalarZ(idx) * NTV
        end if

        self % scalarFlux(idx) =  (self % scalarflux(idx) + self % source(idx))

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(linearAnisotropicScattering), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, xScatter, yScatter, zScatter, &
                                                             fission, xFission, yFission, zFission, &
                                                             xSource, ySource, zSource
    real(defFlt)                                          :: invMxx, invMxy, invMxz, invMyy, invMyz, invMzz
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, &
                                                             condX, condY, condZ, inversionTest
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec, xFluxVec, yFluxVec, zFluxVec
    real(defReal), pointer, dimension(:)                  :: momVec
    real(defReal)                                         :: det, one_det 

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 

    momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize))

    ! Pre-invert the moment matrix
    ! Need to check for poor conditioning by evaluating the
    ! diagonal elements of the matrix
    condX = 0
    condY = 0
    condZ = 0

    if (momVec(xx) > 1.0E-6_defReal) condX = 1
    if (momVec(yy) > 1.0E-6_defReal) condY = 1
    if (momVec(zz) > 1.0E-6_defReal) condZ = 1

    ! Map conditions to test variable
    inversionTest = condX * 4 + condY * 2 + condZ

    select case(inversionTest)
    case(invertXYZ)
      det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
            - momVec(yy) * momVec(xz) * momVec(xz) - momVec(zz) * momVec(xy) * momVec(xy) &
            + 2 * momVec(xy) * momVec(xz) * momVec(yz)
      one_det = ONE/det
      invMxx = real(one_det * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)),defFlt)
      invMxy = real(one_det * (momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz)),defFlt)
      invMxz = real(one_det * (momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz)),defFlt)
      invMyy = real(one_det * (momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)),defFlt)
      invMyz = real(one_det * (momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz)),defFlt)
      invMzz = real(one_det * (momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)),defFlt)

    case(invertYZ)
      det = momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)
      one_det = ONE/det
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = real(one_det * momVec(zz),defFlt)
      invMyz = real(-one_det * momVec(yz),defFlt)
      invMzz = real(one_det * momVec(yy),defFlt)

    case(invertXY)
      det = momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)
      one_det = ONE/det
      invMxx = real(one_det * momVec(yy),defFlt)
      invMxy = real(-one_det * momVec(xy),defFlt)
      invMxz = 0.0_defFlt
      invMyy = real(one_det * momVec(xx),defFlt)
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt

    case(invertXZ)
      det = momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)
      one_det = ONE/det
      invMxx = real(one_det * momVec(zz),defFlt)
      invMxy = 0.0_defFlt
      invMxz = real(-one_det * momVec(xz),defFlt)
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = real(one_det * momVec(xx),defFlt)

    case(invertX)
      det = momVec(xx)
      one_det = ONE/det
      invMxx = real(one_det,defFlt)
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFLt

    case(invertY)
      det = momVec(yy)
      one_det = ONE/det
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = real(one_det,defFlt)
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt

    case(invertZ)
      det = momVec(zz)
      one_det = ONE/det
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = real(one_det,defFlt)

    case default
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt
      det = ONE
    end select

    ! Check for zero determinant
    if (abs(det) < 1E-6) then
      invMxx = 0.0_defFlt
      invMxy = 0.0_defFlt
      invMxz = 0.0_defFlt
      invMyy = 0.0_defFlt
      invMyz = 0.0_defFlt
      invMzz = 0.0_defFlt
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
        
      ! Calculate source gradients by inverting the moment matrix
      self % sourceX(baseIdx + g) = invMxx * xSource + &
              invMxy * ySource + invMxz * zSource
      self % sourceY(baseIdx + g) = invMxy * xSource + &
              invMyy * ySource + invMyz * zSource
      self % sourceZ(baseIdx + g) = invMxz * xSource + &
           invMyz * ySource + invMzz * zSource

    end do

  end subroutine sourceUpdateKernel

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(linearAnisotropicScattering), intent(inout) :: self
    real(defFlt)                                  :: fissionRate, prevFissionRate
    real(defFlt), save                            :: fissLocal, prevFissLocal, vol
    integer(shortInt), save                       :: matIdx, g, idx, mIdx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, matIdx, g, idx, mIdx, vol)

    fissionRate     = 0.0_defFlt
    prevFissionRate = 0.0_defFlt
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      if (matIdx > 100000) cycle
      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      if (.not. mat % isFissile()) cycle

      vol = real(self % volume(cIdx), defFlt)

      if (vol <= volume_tolerance) cycle

      fissLocal = 0.0_defFlt
      prevFissLocal = 0.0_defFlt
      mIdx = (matIdx - 1) * self % nG
      do g = 1, self % nG
        
        ! Source index
        idx = self % nG * (cIdx - 1) + g
        fissLocal     = fissLocal     + self % scalarFlux(idx) * self % nuSigmaF(mIdx + g)
        prevFissLocal = prevFissLocal + self % prevFlux(idx) * self % nuSigmaF(mIdx + g)

      end do

      fissionRate     = fissionRate     + fissLocal * vol
      prevFissionRate = prevFissionRate + prevFissLocal * vol

    end do
    !$omp end parallel do

    ! Update k
    self % keff = self % keff * fissionRate / prevFissionRate

  end subroutine calculateKeff
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! Same for flux moments
  !!
  subroutine resetFluxes(self)
    class(linearAnisotropicScattering), intent(inout) :: self
    integer(shortInt)                                  :: idx

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
    end do
    !$omp end parallel do
    
  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(linearAnisotropicScattering), intent(inout) :: self
    real(defReal), save                            :: flux
    integer(shortInt)                              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff

  end subroutine accumulateFluxAndKeffScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it)
    class(linearAnisotropicScattering), intent(inout) :: self
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
    self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(linearAnisotropicScattering), intent(inout) :: self
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
    class(linearAnisotropicScattering), intent(in) :: self

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
    class(linearAnisotropicScattering), intent(inout) :: self
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
      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if
    if(allocated(self % fluxMap)) then
      call self % resultsMap % kill()
      deallocate(self % fluxMap)
    end if

  end subroutine kill

end module linearAnisotropicScattering_class
