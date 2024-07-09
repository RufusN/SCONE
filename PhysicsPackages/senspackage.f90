module adjointTRRMPhysicsPackage_class !currently backward

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential, F1
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
  real(defReal), parameter :: volume_tolerance = 1.0E-12

  !!
  !! Physics package to perform The Random Ray Method (TRRM) eigenvalue calculations
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
  !!     type adjointTRRMPhysicsPackage;
  !!     dead 10;              // Dead length where rays do not score to scalar fluxes
  !!     termination 100;      // Length a ray travels before it is terminated
  !!     rays 1000;            // Number of rays to sample per iteration
  !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;           // Number of scoring cycles (would use eps otherwise)
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #fluxMap {<map>}#     // Optionally output one-group fluxes according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!     #rho 0;#              // Optional stabilisation for negative in-group scattering XSs
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
  !!   rho         -> Stabilisation factor: 0 is no stabilisation, 1 is aggressive stabilisation
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
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!   volume      -> Array of stochastically estimated cell volumes of length = nCells
  !!   cellHit     -> Array tracking whether given cells have been hit during tracking
  !!   cellFound   -> Array tracking whether a cell was ever found
  !!   cellPos     -> Array of cell positions, populated once they are found
  !!
  !!   locks       -> OpenMP locks used when writing to scalar flux during transport
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: adjointTRRMPhysicsPackage
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

    !Pertubation
    integer(shortInt), dimension(:), allocatable   :: energyId 
    integer(shortInt)  :: XStype      = 0
    integer(shortInt)  :: matPert     = 0
    real(defReal)      :: XSchange    = ZERO
    real(defReal)      :: numSum      = 0
    real(defReal)      :: denSum      = 0

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
    real(defFlt), dimension(:), allocatable     :: fission
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi
    real(defFlt), dimension(:), allocatable     :: sigmaC
    real(defFlt), dimension(:), allocatable     :: adjNuSigmaF
    real(defFlt), dimension(:), allocatable     :: adjSigmaS
    real(defFlt), dimension(:), allocatable     :: adjChi

    ! Results space
    real(defFlt)                               :: keff
    real(defFlt)                               :: adjkeff
    real(defReal)                              :: sensitivity 
    real(defReal)                              :: deltaKeff 
    real(defReal), dimension(2)                :: keffScore
    real(defReal), dimension(2)                :: adjkeffScore
    real(defReal), dimension(2)                :: sensitivityScore
    real(defReal), dimension(2)                :: deltaKeffScore
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defFlt), dimension(:), allocatable    :: adjScalarFlux
    real(defFlt), dimension(:), allocatable    :: adjPrevFlux
    real(defReal), dimension(:), allocatable   :: angularIP
    real(defReal), dimension(:,:), allocatable :: fluxScores
    real(defFlt), dimension(:), allocatable    :: source
    real(defFlt), dimension(:), allocatable    :: adjSource
    real(defReal), dimension(:), allocatable   :: volume
    real(defReal), dimension(:), allocatable   :: volumeTracks

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
    procedure, private :: initPertubation
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sourceUpdateKernel
    procedure, private :: adjointSourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: adjointNormaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: sensitivityCalculation
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings

  end type adjointTRRMPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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
    character(100), parameter :: Here = 'init (adjointTRRMPhysicsPackage_class.f90)'

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


    !!!!!!!!!!!!!!!!!!!!!!
    tempDict => dict % getDictPtr("Pertubation")
    call self % initPertubation(tempDict)
    
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
    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))

    allocate(self % adjScalarFlux(self % nCells * self % nG))
    allocate(self % adjPrevFlux(self % nCells * self % nG))
    allocate(self % adjSource(self % nCells * self % nG))

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

  subroutine initPertubation(self,dict)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(in)                   :: dict
    class(dictionary), pointer                      :: tempDict
    integer(shortInt)                               :: i, g
    character(100), parameter :: Here = 'initPertubation (adjointTRRMPhysicsPackage_class.f90)'

    call dict % getOrDefault(self % XStype, 'XStype', 0)

    call dict % getOrDefault(self % matPert, 'material', 0)

    call dict % getOrDefault(self % XSchange, 'XSchange', ZERO)

    call dict % get(self % energyId, 'energyGroups')

    if (self % XStype == 3 .and. mod(size(self % energyId), 2) /= 0) then
      call fatalError(Here, 'energyGroups for scattering XS pertubation must be given in pairs.')
    end if 

  
  end subroutine initPertubation

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self

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
  subroutine cycles(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF, adj_KEFF
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: keepRunning, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % keff       = 1.0_defFlt
    self % adjkeff    = 1.0_defFlt

    self % scalarFlux    = 0.0_defFlt
    self % prevFlux      = 1.0_defFlt
    self % source        = 0.0_defFlt
    self % adjScalarFlux = 0.0_defFlt
    self % adjPrevFlux   = 1.0_defFlt
    self % adjSource     = 0.0_defFlt
    self % angularIP     = 0.0_defReal

    self % keffScore        = ZERO
    self % adjKeffScore     = ZERO
    self % fluxScores       = ZERO
    self % deltaKeffScore   = ZERO
    self % sensitivity      = ZERO
    self % sensitivityScore = ZERO

    ! Initialise other results
    self % cellHit      = 0
    self % volume       = ZERO
    self % volumeTracks = ZERO
    self % intersectionsTotal  = 0

    ! Initialise cell information
    self % cellFound = .false.
    self % cellPos = -INFINITY


    ! Stopping criterion is on number of convergence iterations.
    ! TODO: Make this on, e.g., entropy during inactive, followed by stochastic noise during active!
    itInac = 0
    itAct  = 0
    isActive = .false.
    keepRunning = .true.
    
    ! Power iteration
    do while( keepRunning )
      
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
        call self % sourceUpdateKernel(i, ONE_KEFF)
        call self % adjointSourceUpdateKernel(i, adj_KEFF)
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
        call self % transportSweep(r,ints)
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

      self % numSum = ZERO
      self % denSum = ZERO
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
       call self % sensitivityCalculation(i ,ONE_KEFF, it)
      end do
      self % deltaKeff  = self % keff * self % keff * (self % numSum / self % denSum) 
      self % sensitivity = (self % deltaKeff / ( self % keff * self % XSchange ) ) 

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxAndKeffScores()

      ! Calculate proportion of cells that were hit
      hitRate = real(sum(self % cellHit),defFlt) / self % nCells
      self % cellHit = 0

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        keepRunning = (itAct < self % active)
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
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (adjointTRRMPhysicsPackage_class.f90)'

    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      x = self % bottom + (self % top - self % bottom) * rand3

      ! Exit if point is inside the geometry
      call self % geom % whatIsAt(matIdx, cIdx, x, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, 1, ONE)
    call self % geom % placeCoord(r % coords)

    if (.not. self % cellFound(cIdx)) then
      !$omp critical 
      self % cellFound(cIdx) = .true.
      self % cellPos(cIdx,:) = x
      !$omp end critical
    end if

  end subroutine initialiseRay

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep(self, r, ints)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, matIdx0, baseIdx, &
                                                               segCount, i, segCountCrit, segIdx, gIn, pIdx
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum, tally
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec, avgFluxVec, tau
    real(defFlt), allocatable, dimension(:)               :: tauBack, tauBackBuffer
    real(shortInt), allocatable, dimension(:)             :: cIdxBack, cIdxBackBuffer
    logical(defBool), allocatable, dimension(:)           :: vacBack, vacBackBuffer 
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec
    real(defReal), pointer, dimension(:)                  :: angularProd
    real(defFlt), allocatable, dimension(:)               :: fluxRecord, fluxRecordBuffer 
    ! real(defFlt), allocatable, dimension(:)               :: adjointRecord
    real(defReal), dimension(3)                           :: r0, mu0
        
    ! Set initial angular flux to angle average of cell source
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

    allocate(tauBack(self % nG * 50))         ! could add some kind of termination  / average cell length to get an approx length
    allocate(fluxRecord(self % nG * 50))
    allocate(cIdxBack(50))
    allocate(vacBack(50))

    do while (totalLength < self % termination + self % dead)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID

      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        ! Cache total cross section
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
      end if

      ! Remember co-ordinates to set new cell's position
      if (.not. self % cellFound(cIdx)) then
        r0 = r % rGlobal()
        mu0 = r % dirGlobal()
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
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = r0 + length * HALF * mu0
        !$omp end critical
      end if

      ints = ints + 1
      lenFlt = real(length,defFlt)
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % source((baseIdx + 1):(baseIdx + self % nG))

      !$omp simd aligned(totVec)
      do g = 1, self % nG
        tau(g) = totVec(g) * lenFlt
      end do
        
      !$omp simd 
      do g = 1, self % nG
        attenuate(g) = F1(tau(g))
      end do

      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
      end do

      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g) * tau(g)
      end do

      ! Accumulate to scalar flux
      if (activeRay) then
        segCount = segCount + 1

        if (tally) then
            segCountCrit = segCount
        end if
 
        if (segCount > size(cIdxBack) - 1) then ! doubles length if needed

          allocate(tauBackBuffer(size(tauBack) * 2))   !record expoentials 
          tauBackBuffer(1:size(tauBack)) = tauBack
          call move_alloc(tauBackBuffer, tauBack)

          allocate(cIdxBackBuffer(size(cIdxBack) * 2)) !add cell id for scalar flux update
          cIdxBackBuffer(1:size(cIdxBack)) = cIdxBack
          call move_alloc(cIdxBackBuffer, cIdxBack)

          allocate(vacBackBuffer(size(vacBack) * 2))   !check vacuum
          vacBackBuffer(1:size(vacBack)) = vacBack
          call move_alloc(vacBackBuffer, vacBack)

          allocate(fluxRecordBuffer(size(fluxRecord) * 2))   ! recording of average flux
          fluxRecordBuffer(1:size(fluxRecord)) = fluxRecord
          call move_alloc(fluxRecordBuffer, fluxRecord)
        end if

        !$omp simd
        do g = 1, self % nG
          avgFluxVec(g) = (delta(g) + sourceVec(g))
        end do
        
        !$omp simd
        do g = 1, self % nG
          tauBack((segCount - 1) * self % nG + g) = tau(g)
        end do

        cIdxBack(segCount) = cIdx

        if (tally) then

        !$omp simd
        do g = 1, self % nG
          fluxRecord((segCount - 1) * self % nG + g)  = avgFluxVec(g)  
          !fluxRecord((segCount) * self % nG - g)  = avgFluxVec(g) !???
        end do

          call OMP_set_lock(self % locks(cIdx))
          scalarVec => self % scalarFlux((baseIdx + 1):(baseIdx + self % nG))
          !$omp simd aligned(scalarVec)
          do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + delta(g) * tau(g)
          end do
          self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length 
          call OMP_unset_lock(self % locks(cIdx))
        end if

        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1

        vacBack(segCount + 1) = hitVacuum
      
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

    ! Could trim arrays here, (tauBack...)
    allocate(tauBackBuffer(segCount))
    tauBackBuffer = tauBack(1 : (segCount * self % nG))
    call move_alloc(tauBackBuffer, tauBack)

    allocate(cIdxBackBuffer(segCount))
    cIdxBackBuffer = cIdxBack(1:segCount)
    call move_alloc(cIdxBackBuffer, cIdxBack)

    allocate(vacBackBuffer(segCount))
    vacBackBuffer = vacBack(1:segCount)
    call move_alloc(vacBackBuffer, vacBack)

    allocate(fluxRecordBuffer(segCount))
    fluxRecordBuffer = fluxRecord(1 : (segCount * self % nG))
    call move_alloc(fluxRecordBuffer, fluxRecord)

    !iterate over segments
    do i = segCount, 1, -1
      
      cIdx = cIdxBack(i)
      baseIdx = (cIdx - 1) * self % nG
      segIdx =  (i - 1) * self % nG
      sourceVec => self % adjSource((baseIdx + 1):(baseIdx + self % nG))

      !$omp simd
      do g = 1, self % nG
        attenuate(g) = F1(tauBack(segIdx + g))
      end do

      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
      end do

      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g) * tauBack(segIdx + g)
      end do
      
      if (i <= segCountCrit) then

        !$omp simd
        do g = 1, self % nG
          avgFluxVec(g) = (sourceVec(g) + delta(g))
        end do

        call OMP_set_lock(self % locks(cIdx))

        scalarVec => self % adjScalarFlux((baseIdx + 1):(baseIdx + self % nG))
        angularProd => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG )
        !$omp simd 
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + delta(g) * tauBack(segIdx + g)
            do gIn = 1, self % nG
              pIdx = self % nG * (g - 1) + gIn
              angularProd(pIdx) = angularProd(pIdx) + avgFluxVec(g) * fluxRecord(segIdx + gIn)
            end do
        end do
        call OMP_unset_lock(self % locks(cIdx))

      end if  

      if (vacBack(i)) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

     end do

     !use records to computer inner product etc - seperate subroutine? 
     
  end subroutine transportSweep

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defFlt), save                            :: total, vol, sigGG, D
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx),defFlt)

      do g = 1, self % nG

        total = self % sigmaT((matIdx - 1) * self % nG + g)
        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total * vol)
        end if

        ! self % scalarFlux(idx) = self % scalarflux(idx) + self % source(idx) 
        
        ! Apply stabilisation for negative XSs
        if (matIdx < UNDEF_MAT) then
          sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)

          ! Presumes non-zero total XS
          if ((sigGG < 0) .and. (total > 0)) then
            D = -real(self % rho, defFlt) * sigGG / total
          else
            D = 0.0_defFlt
          end if
        else
          D = 0.0_defFlt
        end if

        self % scalarFlux(idx) =  (self % scalarflux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume

  subroutine adjointNormaliseFluxAndVolume(self, it)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defFlt), save                            :: total, vol, sigGG, D
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx),defFlt)

      do g = 1, self % nG

        total = self % sigmaT((matIdx - 1) * self % nG + g)
        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          self % adjScalarFlux(idx) = self % adjScalarFlux(idx) * norm / ( total * vol)
        end if

        ! self % adjScalarFlux(idx) = self % adjScalarFlux(idx) + self % adjSource(idx)
        
        ! Apply stabilisation for negative XSs
        if (matIdx < UNDEF_MAT) then

          sigGG = self % adjSigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)

          ! Presumes non-zero total XS
          if ((sigGG < 0) .and. (total > 0)) then
            D = -real(self % rho, defFlt) * sigGG / total
          else
            D = 0.0_defFlt
          end if
        else
          D = 0.0_defFlt
        end if

        self % adjScalarFlux(idx) =  (self % adjScalarflux(idx) + self % adjSource(idx) + D * self % adjPrevFlux(idx) ) / (1 + D)

      end do

    end do
    !$omp end parallel do

  end subroutine adjointNormaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, fission
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= UNDEF_MAT) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx+1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission) aligned(fluxVec)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG))

      ! Calculate scattering source
      scatter = 0.0_defFlt

      ! Sum contributions from all energies
      !$omp simd reduction(+:scatter) aligned(fluxVec)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      self % source(idx) = self % source(idx) / total(g)

    end do

  end subroutine sourceUpdateKernel

  subroutine adjointSourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, fission
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= UNDEF_MAT) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % adjSource(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % adjSigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
    nuFission => self % adjNuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % adjChi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % adjPrevFlux((baseIdx+1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission) aligned(fluxVec)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG))

      ! Calculate scattering source
      scatter = 0.0_defFlt

      ! Sum contributions from all energies
      !$omp simd reduction(+:scatter) aligned(fluxVec)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      self % adjSource(idx) = chi(g) * fission + scatter
      self % adjSource(idx) = self % adjSource(idx) / total(g)

    end do

  end subroutine adjointSourceUpdateKernel

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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
  !!
  subroutine resetFluxes(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt)                             :: idx, g

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
      self % adjPrevFlux(idx) = self % adjScalarFlux(idx)
      self % adjScalarFlux(idx) = 0.0_defFlt
      do g = 1, self % nG
        self % angularIP((idx - 1) * self % nG + g) = 0.0_defFlt
      end do
    end do
    !$omp end parallel do

  end subroutine resetFluxes

  subroutine sensitivityCalculation(self, cIdx, ONE_KEFF, it)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
    real(defFlt), intent(in)                      :: ONE_KEFF
    integer(shortInt), intent(in)                 :: it
    integer(shortInt), intent(in)                 :: cIdx
    real(defReal)                                 :: delta, fission, fission_pert, scatter_pert
    integer(shortInt)                             :: baseIdx, idx, matIdx, g, gIn, mat, g1Pert, g2pert, i, matPert
    real(defFlt), dimension(:), pointer           :: nuFission, total, chi, capture, scatterXS, &
                                                      fissVec, scatterVec
    real(defReal), dimension(:), pointer          :: IPVec

    !removed normalisation 
    
    mat  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    baseIdx = (cIdx - 1) * self % nG 
    matIdx = (mat - 1) * self % nG
    matPert = self % matPert

    IPVec => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG)
    ! total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    fissVec => self % fission((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))
    capture => self % sigmaC((matIdx + 1):(matIdx + self % nG)) 

    do g = 1, self % nG
      idx = baseIdx + g
      fission = 0.0_defFlt
      !$omp simd
      do gIn = 1, self % nG 
        fission = fission + IPVec((g - 1) * self % nG + gIn) * nuFission(gIn)
      end do
      fission  = fission * chi(g)
      self % denSum = self % denSum + fission
    end do

      if (self % XStype == 1) then ! capture - complete 
        do i = 1, size(self % energyId)
          g1Pert = self % energyId(i)
          
          !$omp simd
          do g = 1, self % nG 
            delta = 0.0_defFlt
            ! Change in XS
              if ( mat == matPert .and. g == g1Pert ) then  !add energy/mat group check if needed
                delta = self % XSchange * capture(g) * IPVec((g - 1) * self % nG + g)
              end if
            idx = baseIdx + g
            self % numSum = self % numSum - delta
          end do
        end do

      elseif (self % XStype == 2) then ! fission - complete
        do i = 1, size(self % energyId)
          g1Pert = self % energyId(i)
          
          do g = 1, self % nG 
            delta = 0.0_defFlt
            fission_pert = 0.0_defFlt
            !$omp simd
            do gIn = 1, self % nG
              if ( gIn == g1Pert ) then
                fission_pert = fission_pert + IPVec(self % nG * (g - 1) + gIn) * nuFission(gIn) * self % XSchange
              end if
            end do
            if ( g == g1Pert ) then 
              delta = IPVec((g - 1) * self % nG + g) * fissVec(g) * self % XSchange
            end if
            idx = baseIdx + g
            self % numSum = self % numSum - delta + fission_pert * chi(g) * ONE_KEFF
          end do
        end do

      elseif (self % XStype == 3) then !scatter - complete
        do i = 1, size(self % energyId), 2
          scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
          g1Pert = self % energyId(i)
          g2Pert = self % energyId(i+1)

          do g = 1, self % nG 
            delta = 0.0_defFlt
            scatter_pert = 0.0_defFlt
            if (g == g1Pert) then
              delta = IPVec((g1Pert - 1) * self % nG + g1Pert) * scatterXS((g2Pert - 1) * self % nG + g1Pert) * &
                           self % XSchange
              do gIn = 1, self % nG
                if (gIn == g2Pert) then
                  scatter_pert = scatter_pert + scatterXS( (g2Pert - 1) * self % nG + g1Pert ) * &
                                              IPVec((g2Pert - 1 ) * self % nG + g1Pert) * self % XSchange
                end if
              end do
            end if

            idx = baseIdx + g
            self % numSum = self % numSum -delta + scatter_pert
          end do
        end do

      end if



  end subroutine sensitivityCalculation

  ! subroutine sensitivityCalculation(self, ONE_KEFF, it)
  !   class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
  !   real(defFlt), intent(in)                      :: ONE_KEFF
  !   integer(shortInt), intent(in)                 :: it
  !   real(defReal)                                 :: norm
  !   real(defFlt)                                  :: XSchange
  !   real(defReal)                                 :: delta, fission, fission_pert, scatter_pert, numSum, denSum
  !   integer(shortInt)                             :: cIdx, baseIdx, idx, matIdx, g, gIn, mat, XScase, g1Pert, g2Pert
  !   real(defFlt), dimension(:), pointer           :: nuFission, total, chi, capture, scatterXS, &
  !                                                     fissVec, scatterVec
  !   real(defReal), dimension(:), pointer          :: IPVec
  !   real(defReal), dimension(self % nG * self % nCells) :: num, den

  !   !change in XS
  !   XSchange = 0.01 
  !   XScase   = 2
  !   norm     = ONE / (self % lengthPerIt) 

  !   do cIdx = 1, self % nCells
  !     do g = 1, self % nG * self % nG
  !       idx = (cIdx - 1) * self % nG * self % nG + g
  !       if (self % volume(cIdx) > volume_tolerance) then
  !         self % angularIP(idx) = self % angularIP(idx) * norm / (self % volume(cIdx))
  !       else
  !         self % angularIP(idx) = 0.0_defFlt
  !       end if
  !     end do
  !   end do

  !   do idx = 1, self % nCells * self % nG
  !     num(idx) =  0.0_defFlt
  !     den(idx) = 0.0_defFlt
  !   end do

  !   do cIdx = 1, self % nCells

  !     ! Identify material
  !     mat  =  self % geom % geom % graph % getMatFromUID(cIdx) 
  !     baseIdx = (cIdx - 1) * self % nG 
  !     matIdx = (mat - 1) * self % nG

  !     IPVec => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG)
  !     total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
  !     nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
  !     fissVec => self % fission((matIdx + 1):(matIdx + self % nG))
  !     chi => self % chi((matIdx + 1):(matIdx + self % nG))
  !     capture => self % sigmaC((matIdx + 1):(matIdx + self % nG)) 
  !     scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))


  !     do g = 1, self % nG
  !       idx = baseIdx + g
  !       fission = 0.0_defFlt
  !       !$omp simd
  !       do gIn = 1, self % nG 
  !         fission = fission + IPVec((g - 1) * self % nG + gIn) * nuFission(gIn)
  !       end do
  !       fission  = fission * chi(g)
  !       den(idx) = fission
  !     end do

  !     if (XScase == 1) then ! capture - complete 

  !       g1Pert = 1

  !       do g = 1, self % nG 

  !         delta = 0.0_defFlt

  !         ! Change in XS
  !          if ( g == g1Pert ) then  !mat == 2 .and.  !add energy/mat group check if needed
  !            delta = XSchange * capture(g) * IPVec((g - 1) * self % nG + g)
  !          end if

  !         idx = baseIdx + g
  !         num(idx) = -delta

  !       end do

  !     elseif (XScase == 2) then ! fission - complete

  !       g1Pert = 1
        
  !       do g = 1, self % nG 

  !         delta = 0.0_defFlt

  !         fission_pert = 0.0_defFlt
  !         !$omp simd
  !         do gIn = 1, self % nG
  !           if ( gIn == g1Pert ) then
  !             fission_pert = fission_pert + IPVec(self % nG * (g - 1) + gIn) * nuFission(gIn) * XSchange
  !           end if
  !         end do

  !         if ( g == g1Pert ) then 
  !           delta = IPVec((g - 1) * self % nG + g) * fissVec(g) * XSchange
  !         end if

  !         idx = baseIdx + g
  !         num(idx) = - delta + fission_pert * chi(g) * ONE_KEFF

  !       end do

  !     elseif (XScase == 3) then !scatter - complete
        
  !       ! Assume input of two numbers. 
  !       g1Pert = 2
  !       g2Pert = 3

  !       do g = 1, self % nG 

  !         delta = 0.0_defFlt
  !         scatter_pert = 0.0_defFlt
          
  !         if (g == g1Pert) then
  !           delta = IPVec((g1Pert - 1) * self % nG + g1Pert) * scatterXS((g2Pert - 1) * self % nG + g1Pert) * XSchange

  !           do gIn = 1, self % nG
  !             if (gIn == g2Pert) then
  !               scatter_pert = scatter_pert + scatterXS( (g2Pert - 1) * self % nG + g1Pert ) * &
  !                                           IPVec((g2Pert - 1 ) * self % nG + g1Pert) * XSchange
  !             end if
  !           end do
  !         end if

  !         ! delta = IPVec((Sg1 - 1) * self % nG + Sg1) * scatterXS((g2Pert - 1) * self % nG + Sg1) * XSchange
  !         ! scatter_pert = scatter_pert + scatterXS((g2Pert - 1) * self % nG + Sg1 ) * IPVec((g2Pert - 1) * self % nG + Sg1) * XSchange
  !         ! for S12, formula is = IP(1) * scatter(3) + IP(2) * scatter(3), contributions from total and scatter operator respectively. 

  !         idx = baseIdx + g
  !         num(idx) = -delta + scatter_pert
  !       end do

  !     end if

  !   end do

  !   ! sum to get IP
  !   numSum = 0.0_defFlt
  !   denSum = 0.0_defFlt
  !   do idx = 1, self % nG * self % nCells
  !     numSum = numSum + num(idx)
  !     denSum = denSum + den(idx)
  !   end do

  !   self % deltaKeff  = self % keff * self % keff * (numSum / denSum) 
  !   !self % deltaKeff = self % keff * self % keff * (numSum / denSum) / (1 - (numSum / denSum)) ! this sometimes seems closer?

  !   self % sensitivity = (self % deltaKeff / ( self % keff * XSchange ) ) 

  ! end subroutine sensitivityCalculation

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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

    !self % deltaKeffScore(1) = self % deltaKeffScore(1) * self % keffScore(1) * self % keffScore(1)
    !self % deltaKeffScore(2) = self % deltaKeffScore(2) * self % keffScore(1) !* self % keffScore(1)

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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
    call out % printResult(real(self % sensitivityScore(1),defReal), real(self % sensitivityScore(2),defReal), name)
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
    class(adjointTRRMPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY EIGENVALUE CALCULATION /\/\"
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
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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
    if(allocated(self % nuSigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % chi)) deallocate(self % chi)
    if(allocated(self % sigmaC)) deallocate(self % sigmaC)

    ! if(allocated(self % adjSigmaT)) deallocate(self % adjSigmaT)
    if(allocated(self % adjSigmaS)) deallocate(self % adjSigmaS)
    if(allocated(self % adjNuSigmaF)) deallocate(self % adjNuSigmaF)
    if(allocated(self % adjChi)) deallocate(self % adjChi)
    if(allocated(self % fission)) deallocate(self % fission)

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    self % cache       = .false.
    self % rho         = ZERO
    self % mapFission  = .false.
    self % mapFlux     = .false.
    self % plotResults = .false.
    self % printFlux   = .false.
    self % printVolume = .false.
    self % printCells  = .false.

    self % XStype   = 0
    self % matPert  = 0
    self % XSchange = 0.0_defFlt
    self % numSum   = ZERO
    self % denSum   = ZERO

    self % sensitivity      = ZERO
    self % deltaKeff        = ZERO
    self % keff             = ZERO
    self % adjkeff          = ZERO
    self % keffScore        = ZERO
    self % adjkeffScore     = ZERO
    self % sensitivityScore = ZERO
    self % deltaKeffScore   = ZERO

    if(allocated(self % energyId)) deallocate(self % energyId)

    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)

    if(allocated(self % adjScalarFlux)) deallocate(self % adjScalarFlux)
    if(allocated(self % adjPrevFlux)) deallocate(self % adjPrevFlux)
    if(allocated(self % adjSource)) deallocate(self % adjSource)
    if (allocated(self % angularIP)) deallocate(self % angularIP)

    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
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

end module adjointTRRMPhysicsPackage_class