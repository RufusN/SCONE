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
  !!     #innerProduct scalar;# // Sensitivity inner product, 'scalar' (default) or 'angular':
  !!                            //  'scalar'  - isotropic gain terms (fission production, in-scatter and
  !!                            //              the <psi+, F psi> denominator) from cell scalar-flux
  !!                            //              products; removal terms from the diagonal angular product
  !!                            //              psi+_g * psi_g, track-length weighted (nG per cell).
  !!                            //              Outputs the full sensitivity profile (every material;
  !!                            //              capture and fission per group; scatter per group pair).
  !!                            //  'angular' - legacy nG x nG angular inner product, one perturbation.
  !!     #sensBatches 20;#      // Batches of active cycles for profile uncertainties ('scalar' only)
  !!     #Perturbation {material 5; XStype 2; XSchange 0.01; energyGroups (6);}#
  !!                            // XStype 1 capture, 2 fission, 3 scatter (groups as from-to pairs).
  !!                            // Required for 'angular'. Optional for 'scalar': the requested entries
  !!                            // are summed (first order) and printed as 'sensitivity'/'delta_keff',
  !!                            // with absolute standard errors of the mean.
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


    integer(shortInt)  :: NSegMax     = 100
    integer(shortInt)  :: NSegMaxTemp = 100
    integer(shortInt)  :: XStype      = 0
    integer(shortInt)  :: matPert     = 0
    real(defReal)      :: XSchange    = ZERO

    ! Sensitivity inner product method
    logical(defBool)   :: legacyIP    = .false. ! .true. -> legacy nG x nG angular inner product
    logical(defBool)   :: hasPert     = .false. ! Was a Perturbation dictionary given?
    integer(shortInt)  :: nPM         = 0       ! Profile entries per material: 2*nG + nG*nG
    integer(shortInt)  :: nBatches    = 20      ! Batches of active cycles for uncertainties
    integer(shortInt)  :: batchLen    = 0
    integer(shortInt)  :: batchIt     = 0
    integer(shortInt)  :: batchCount  = 0

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
    integer(shortInt), dimension(:), allocatable       :: energyId 

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
    real(defReal), dimension(:), allocatable   :: diagIP       ! Cycle sum(length*(psi+ psi - q+ q)), nCells*nG
    real(defReal), dimension(:), allocatable   :: diagIPScore  ! Active sum of removal inner products
    real(defReal), dimension(:), allocatable   :: adjFluxScore ! Active sum of adjoint scalar flux
    real(defReal), dimension(:), allocatable   :: profBatch    ! Sum of cycle profiles in current batch
    real(defReal), dimension(:,:), allocatable :: profScore    ! Batch means: (:,1) sum, (:,2) sum sq -> SEM
    real(defReal), dimension(:), allocatable   :: profMean     ! Final sensitivity profile
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
    procedure, private :: accumulateSensitivity
    procedure, private :: sensitivityProfile
    procedure, private :: finaliseSensitivity
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
    character(nameLen)                            :: geomName, graphType, nucData, ipMethod
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

    ! Sensitivity inner product method
    call dict % getOrDefault(ipMethod, 'innerProduct', 'scalar')
    select case(ipMethod)
      case('scalar')
        self % legacyIP = .false.
      case('angular')
        self % legacyIP = .true.
      case default
        call fatalError(Here, "Unrecognised innerProduct: "//trim(ipMethod)//". Use 'scalar' or 'angular'.")
    end select
    call dict % getOrDefault(self % nBatches, 'sensBatches', 20)

    ! Perturbation is required by the legacy method and optional otherwise
    self % hasPert = dict % isPresent('Perturbation')
    if (self % hasPert) then
      tempDict => dict % getDictPtr("Perturbation")
      call self % initPerturbation(tempDict)
    else if (self % legacyIP) then
      call fatalError(Here, "innerProduct 'angular' requires a Perturbation dictionary.")
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

    if (self % legacyIP) then
      allocate(self % angularIP(self % nCells * self % nG * self % nG))
    else
      allocate(self % diagIP(self % nCells * self % nG))
      allocate(self % diagIPScore(self % nCells * self % nG))
      allocate(self % adjFluxScore(self % nCells * self % nG))
    end if
    
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

    ! Sensitivity profile storage and checks on the requested perturbation ('scalar' only)
    if (.not. self % legacyIP) then
      if (self % active < 1) call fatalError(Here, 'At least one active cycle is needed for sensitivities.')
      self % nPM = 2 * self % nG + self % nG * self % nG
      allocate(self % profBatch(self % nMat * self % nPM + 1))
      allocate(self % profScore(self % nMat * self % nPM + 1, 2))
      allocate(self % profMean(self % nMat * self % nPM + 1))
      self % nBatches = max(1, min(self % nBatches, self % active))
      self % batchLen = self % active / self % nBatches

      if (self % hasPert) then
        if (self % matPert < 1 .or. self % matPert > self % nMat) &
          call fatalError(Here, 'Perturbation material index out of range: '//numToChar(self % matPert))
        if (self % XStype < 1 .or. self % XStype > 3) &
          call fatalError(Here, 'Perturbation XStype must be 1 (capture), 2 (fission) or 3 (scatter).')
        if (size(self % energyId) < 1) call fatalError(Here, 'Perturbation energyGroups is empty.')
        if (any(self % energyId < 1) .or. any(self % energyId > self % nG)) &
          call fatalError(Here, 'Perturbation energyGroups must lie in 1 to '//numToChar(self % nG))
      end if
    end if

  end subroutine init

  subroutine initPerturbation(self,dict)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(in)                   :: dict
    class(dictionary), pointer                      :: tempDict
    integer(shortInt)                               :: i, g
    character(100), parameter :: Here = 'initPerturbation (adjointTRRMPhysicsPackage_class.f90)'

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
    real(defReal)                                 :: numTotal, denTotal
    real(defReal),save                            :: numSum, denSum
    real(defFlt)                                  :: hitRate, ONE_KEFF, adj_KEFF
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: keepRunning, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints, numSum, denSum)

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
    if (allocated(self % angularIP)) self % angularIP = 0.0_defReal
    if (.not. self % legacyIP) then
      self % diagIP       = ZERO
      self % diagIPScore  = ZERO
      self % adjFluxScore = ZERO
      self % profBatch    = ZERO
      self % profScore    = ZERO
      self % profMean     = ZERO
      self % batchIt      = 0
      self % batchCount   = 0
    end if
    self % sensitivity = ZERO
    self % deltaKeff   = ZERO

    self % keffScore        = ZERO
    self % adjKeffScore     = ZERO
    self % fluxScores       = ZERO
    self % deltaKeffScore   = ZERO
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
        call self % transportSweep(r,ints,isActive)
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
        if (self % legacyIP) then
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
        else
          call self % accumulateSensitivity(itAct)
        end if
        call self % accumulateFluxAndKeffScores()
      end if

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
    if (.not. self % legacyIP) call self % finaliseSensitivity(itAct)

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
  subroutine transportSweep(self, r, ints, isActive)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    logical(defBool), intent(in)                          :: isActive
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, matIdx0, baseIdx, &
                                                               segCount, i, segCountCrit, segIdx, gIn, pIdx
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum, tally,  oversized
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec, avgFluxVec, tau
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec, scalarVecFW, sourceVecFW, deltaFW
    real(defReal), pointer, dimension(:)                  :: angularProd, diagVec
    real(defReal), dimension(3)                           :: r0, mu0
    real(defFlt), dimension(self % NSegMax*2)             :: lenBack
    integer(shortInt), dimension(self % NSegMax*2)        :: cIdxBack
    logical(defBool), dimension(self % NSegMax*2)         :: vacBack
    real(defFlt), target, dimension(self % nG * self % NSegMax*2) :: deltaRecord
    real(defReal), allocatable, dimension(:)              :: lenBackOversized, lenBackBuffer
    integer(shortInt), allocatable, dimension(:)          :: cIdxBackOversized, cIdxBackBuffer
    logical(defBool), allocatable, dimension(:)           :: vacBackOversized, vacBackBuffer 
    real(defFlt), target, allocatable, dimension(:)       :: deltaRecordOversized
    real(defFlt), allocatable, dimension(:)               :: deltaRecordBuffer
        
    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % source(idx)
    end do

    ints = 0
    matIdx0 = 0
    segCount = 0
    vacBack(1) = .false.
    totalLength = ZERO
    activeRay = .false.
    tally = .true.
    oversized = .false.

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
 
        if ( (segCount >= (self % NSegMax - 1)) .and. .not. oversized) then

          oversized = .true.
          allocate(lenBackOversized(size(lenBack) * 2)) !add seg length for scalar flux update
          lenBackOversized(1:size(lenBack)) = lenBack
          allocate(cIdxBackOversized(size(cIdxBack) * 2)) !add cell id for scalar flux update
          cIdxBackOversized(1:size(cIdxBack)) = cIdxBack
          allocate(vacBackOversized(size(vacBack) * 2))   !check vacuum
          vacBackOversized(1:size(vacBack)) = vacBack
          allocate(deltaRecordOversized(size(deltaRecord) * 2))   ! recording of delta for scalar flux update
          deltaRecordOversized(1:size(deltaRecord)) = deltaRecord

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
          end if

        end if

        if (oversized) then
          cIdxBackOversized(segCount) = cIdx
          lenBackOversized(segCount) = length
          vacBackOversized(segCount + 1) = hitVacuum
        else
          cIdxBack(segCount) = cIdx
          lenBack(segCount) = length
          vacBack(segCount + 1) = hitVacuum
        end if


        if (tally) then
          if (oversized) then
            !$omp simd
            do g = 1, self % nG
              deltaRecordOversized((segCount - 1) * self % nG + g) = delta(g)  
            end do
          else
            !$omp simd
            do g = 1, self % nG
              deltaRecord((segCount - 1) * self % nG + g) = delta(g)  
            end do
          end if
        end if

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

    ! Flux guess for new adjoint deadlength
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % adjSource(idx)
    end do

    !iterate over segments
    do i = segCount, 1, -1
      segIdx =  (i - 1) * self % nG

      if (oversized) then
        lenFlt = real(lenBackOversized(i),defFlt)
        length = lenBackOversized(i)
        hitVacuum = vacBackOversized(i)
        deltaFW => deltaRecordOversized((segIdx + 1) : (segIdx + self % nG))
        cIdx = cIdxBackOversized(i)
      else
        lenFlt = real(lenBack(i),defFlt)
        length = lenBack(i)
        hitVacuum = vacBack(i)
        deltaFW => deltaRecord((segIdx + 1) : (segIdx + self % nG))
        cIdx = cIdxBack(i)
      end if

      matIdx = self % geom % geom % graph % getMatFromUID(cIdx)
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % adjSource((baseIdx + 1):(baseIdx + self % nG))
      sourceVecFW => self % source((baseIdx + 1):(baseIdx + self % nG))
      
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        ! Cache total cross section
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
      end if

      !$omp simd
      do g = 1, self % nG
        tau(g) = lenFlt * totVec(g)
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
      
      if (i <= segCountCrit) then

        !$omp simd
        do g = 1, self % nG
          avgFluxVec(g) = (sourceVec(g) + delta(g))
        end do

        call OMP_set_lock(self % locks(cIdx))

        scalarVec => self % adjScalarFlux((baseIdx + 1):(baseIdx + self % nG))
        scalarVecFW => self % scalarFlux((baseIdx + 1):(baseIdx + self % nG))

        !$omp simd aligned(scalarVec, scalarVecFW)
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + delta(g) * tau(g)
            scalarVecFW(g) = scalarVecFW(g) + deltaFW(g) * tau(g)
        end do
        if (isActive) then
          if (self % legacyIP) then
            ! Legacy: full nG x nG angular inner product, no length weighting
            angularProd => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG )
            do g = 1, self % nG
              associate(angProd => angularProd(self % nG * (g - 1) + 1 : self % nG * g))
                !$omp simd 
                do gIn = 1, self % nG
                  angProd(gIn) = angProd(gIn) + avgFluxVec(g) * (deltaFW(gIn) + sourceVecFW(gIn))
                end do
              end associate
            end do
          else
            ! Diagonal angular product only, track-length weighted. With psi = q + delta on each
            ! segment, accumulate length * (psi+ psi - q+ q); the flat-source part q+ q is added
            ! per cell with its expected track length in accumulateSensitivity.
            diagVec => self % diagIP((baseIdx + 1):(baseIdx + self % nG))
            !$omp simd
            do g = 1, self % nG
              diagVec(g) = diagVec(g) + length * &
                  (real(delta(g), defReal) * real(deltaFW(g), defReal) &
                 + real(sourceVec(g), defReal) * real(deltaFW(g), defReal) &
                 + real(delta(g), defReal) * real(sourceVecFW(g), defReal))
            end do
          end if
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
    self % NSegMax = self % NSegMaxTemp

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

        self % scalarFlux(idx) = (self % scalarflux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)

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

  subroutine normaliseInnerProduct(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
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

  subroutine sensitivityCalculation(self, cIdx, ONE_KEFF, it, numSum, denSum)
    class(adjointTRRMPhysicsPackage), target, intent(inout) :: self
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
      !$omp simd reduction(+:fission)
      do gIn = 1, self % nG 
        fission = fission + IPVec((g - 1) * self % nG + gIn) * nuFission(gIn)
      end do
      fission  = fission * chi(g)
      denSum = denSum + fission
    end do

    if (mat == self % matPert .and. self % XStype == 1) then ! capture - complete 
      g1Pert = self % energyId(1)
      delta = 0.0_defFlt
      delta = self % XSchange * capture(g1Pert) * IPVec((g1Pert - 1) * self % nG + g1Pert)
      numSum = numSum - delta

    elseif ( mat == self % matPert .and. self % XStype == 2) then ! fission - complete
      g1Pert = self % energyId(1)
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

    elseif (mat == self % matPert .and. self % XStype == 3) then !scatter - complete
      g1Pert = self % energyId(1)
      g2Pert = self % energyId(2)
      delta = IPVec((g1Pert - 1) * self % nG + g1Pert) * scatterXS((g2Pert - 1) * self % nG + g1Pert) * &
                    self % XSchange   
          
      scatter_pert = scatterXS( (g2Pert - 1) * self % nG + g1Pert ) * &
                                  IPVec((g2Pert - 1 ) * self % nG + g1Pert) * self % XSchange
      delta = scatter_pert - delta
      numSum = numSum + delta
    end if


  end subroutine sensitivityCalculation

  !!
  !! Accumulate 'scalar' inner-product sensitivity scores for one active cycle
  !!
  !! Adds this cycle's adjoint flux and removal inner product to the run totals, which give the
  !! final point estimate, and adds this cycle's full profile to the current batch, which is used
  !! only for the uncertainty. The forward flux total is already kept in fluxScores.
  !!
  subroutine accumulateSensitivity(self, itAct)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                   :: itAct
    real(defReal), dimension(:), allocatable        :: prof
    real(defReal)                                   :: norm, vol
    integer(shortInt)                               :: idx, cIdx, g

    norm = ONE / self % lengthPerIt

    ! Removal inner product per cell: sum(length * (psi+ psi - q+ q)) / lengthPerIt + V * q+ q.
    ! Using the expected track length V * lengthPerIt for the flat-source part, instead of this
    ! cycle's sampled track length, removes volume-sampling noise without changing the mean
    ! (the same reasoning as the simulation-averaged volume estimator). Converted in place: diagIP
    ! is reset every cycle.
    !$omp parallel do schedule(static) private(vol, idx, g)
    do cIdx = 1, self % nCells
      vol = self % volume(cIdx)
      do g = 1, self % nG
        idx = (cIdx - 1) * self % nG + g
        self % diagIP(idx) = self % diagIP(idx) * norm + &
                             vol * real(self % adjSource(idx), defReal) * real(self % source(idx), defReal)
        self % diagIPScore(idx)  = self % diagIPScore(idx) + self % diagIP(idx)
        self % adjFluxScore(idx) = self % adjFluxScore(idx) + real(self % adjScalarFlux(idx), defReal)
      end do
    end do
    !$omp end parallel do

    ! This cycle's profile from this cycle's fluxes
    allocate(prof(size(self % profMean)))
    call self % sensitivityProfile(prof, self % scalarFlux, self % adjScalarFlux, self % diagIP, &
                                   real(self % keff, defReal))
    self % profBatch = self % profBatch + prof
    self % batchIt   = self % batchIt + 1

    ! Close the batch; the last batch absorbs any remainder of active / nBatches
    if ((self % batchIt == self % batchLen .and. self % batchCount < self % nBatches - 1) &
        .or. itAct == self % active) then
      prof = self % profBatch / self % batchIt
      self % profScore(:,1) = self % profScore(:,1) + prof
      self % profScore(:,2) = self % profScore(:,2) + prof * prof
      self % batchCount = self % batchCount + 1
      self % profBatch  = ZERO
      self % batchIt    = 0
    end if

  end subroutine accumulateSensitivity

  !!
  !! Sensitivity profile S = (alpha/k) dk/dalpha for every material and reaction
  !!
  !! Gain terms (fission production, in-scatter and the <psi+, F psi> denominator) involve isotropic
  !! emission, so their inner products reduce exactly to scalar-flux products V * phi+_g' * phi_g.
  !! Only removal needs the angular product, and only its diagonal:
  !! removal(cIdx, g) ~ V * <psi+_g psi_g>, see accumulateSensitivity.
  !! Same first-order formulae as sensitivityCalculation, e.g. dk = k^2 num / den.
  !!
  !! Layout of prof, per material m at offset (m - 1) * nPM:
  !!   capture(g) at g, fission(g) at nG + g, scatter g1 -> g2 at 2*nG + nG*(g2 - 1) + g1.
  !!   prof(nMat*nPM + 1) is the sum of the entries requested in the Perturbation dictionary (or 0).
  !!
  subroutine sensitivityProfile(self, prof, phi, adjPhi, removal, k)
    class(adjointTRRMPhysicsPackage), intent(in) :: self
    real(defReal), dimension(:), intent(out)     :: prof
    real(defFlt), dimension(:), intent(in)       :: phi, adjPhi
    real(defReal), dimension(:), intent(in)      :: removal
    real(defReal), intent(in)                    :: k
    real(defReal), dimension(:), allocatable     :: remM, prodM, scatM
    real(defReal)                                :: vol, chiAdj, den, sigS
    integer(shortInt)                            :: cIdx, m, g, g1, g2, baseIdx, mIdx, sIdx, off, i, nG

    nG = self % nG
    allocate(remM(self % nMat * nG), prodM(self % nMat * nG), scatM(self % nMat * nG * nG))
    remM  = ZERO
    prodM = ZERO
    scatM = ZERO

    ! Reduce cell contributions to per-material sums
    !$omp parallel do schedule(static) private(m, vol, chiAdj, g, g1, g2, baseIdx, mIdx, sIdx) reduction(+: remM, prodM, scatM)
    do cIdx = 1, self % nCells
      m = self % geom % geom % graph % getMatFromUID(cIdx)
      if (m < 1 .or. m > self % nMat) cycle
      vol = self % volume(cIdx)
      if (vol <= volume_tolerance) cycle
      baseIdx = (cIdx - 1) * nG
      mIdx    = (m - 1) * nG

      ! Fission-spectrum weighted adjoint flux, chi . phi+
      chiAdj = ZERO
      do g = 1, nG
        chiAdj = chiAdj + self % chi(mIdx + g) * real(adjPhi(baseIdx + g), defReal)
      end do

      do g = 1, nG
        remM(mIdx + g)  = remM(mIdx + g)  + removal(baseIdx + g)
        prodM(mIdx + g) = prodM(mIdx + g) + vol * chiAdj * real(phi(baseIdx + g), defReal)
      end do

      ! In-scatter g1 -> g2, only where the scattering matrix is non-zero
      do g2 = 1, nG
        do g1 = 1, nG
          sIdx = nG * nG * (m - 1) + nG * (g2 - 1) + g1
          if (self % sigmaS(sIdx) == 0.0_defFlt) cycle
          scatM(sIdx) = scatM(sIdx) + vol * real(adjPhi(baseIdx + g2), defReal) * real(phi(baseIdx + g1), defReal)
        end do
      end do
    end do
    !$omp end parallel do

    ! Denominator <psi+, F psi> = sum over materials and groups of nuSigmaF_g * (chi . phi+) * phi_g
    den = ZERO
    do i = 1, self % nMat * nG
      den = den + self % nuSigmaF(i) * prodM(i)
    end do

    prof = ZERO
    if (den <= ZERO) return

    do m = 1, self % nMat
      mIdx = (m - 1) * nG
      off  = (m - 1) * self % nPM
      do g = 1, nG
        ! Capture: removal only
        prof(off + g) = -k * self % sigmaC(mIdx + g) * remM(mIdx + g) / den
        ! Fission: isotropic production minus removal
        prof(off + nG + g) = (self % nuSigmaF(mIdx + g) * prodM(mIdx + g) &
                              - k * self % fission(mIdx + g) * remM(mIdx + g)) / den
      end do
      ! Scatter g1 -> g2: in-scatter into g2 minus removal from g1
      do g2 = 1, nG
        do g1 = 1, nG
          sIdx = nG * nG * (m - 1) + nG * (g2 - 1) + g1
          sigS = real(self % sigmaS(sIdx), defReal)
          prof(off + 2 * nG + nG * (g2 - 1) + g1) = k * sigS * (scatM(sIdx) - remM(mIdx + g1)) / den
        end do
      end do
    end do

    ! Entries requested in the Perturbation dictionary, summed to first order
    if (self % hasPert) then
      off = (self % matPert - 1) * self % nPM
      select case(self % XStype)
        case(1)
          do i = 1, size(self % energyId)
            prof(size(prof)) = prof(size(prof)) + prof(off + self % energyId(i))
          end do
        case(2)
          do i = 1, size(self % energyId)
            prof(size(prof)) = prof(size(prof)) + prof(off + nG + self % energyId(i))
          end do
        case(3)
          do i = 1, size(self % energyId) - 1, 2
            prof(size(prof)) = prof(size(prof)) + &
                    prof(off + 2 * nG + nG * (self % energyId(i + 1) - 1) + self % energyId(i))
          end do
      end select
    end if

  end subroutine sensitivityProfile

  !!
  !! Final 'scalar' sensitivity profile from run-averaged fluxes and removal products, with the
  !! standard error of the mean estimated from the spread of the batch means.
  !! Must follow finaliseFluxAndKeffScores, which turns fluxScores and keffScore into means.
  !!
  subroutine finaliseSensitivity(self, itAct)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                   :: itAct
    real(defReal)                                   :: N1, B
    integer(shortInt)                               :: j

    N1 = ONE / itAct
    call self % sensitivityProfile(self % profMean, real(self % fluxScores(:,1), defFlt), &
                                   real(self % adjFluxScore * N1, defFlt), self % diagIPScore * N1, &
                                   self % keffScore(1))

    B = real(self % batchCount, defReal)
    do j = 1, size(self % profMean)
      if (self % batchCount > 1) then
        self % profScore(j,2) = (self % profScore(j,2) - self % profScore(j,1)**2 / B) / (B * (B - ONE))
        self % profScore(j,2) = sqrt(max(self % profScore(j,2), ZERO))
      else
        self % profScore(j,2) = ZERO
      end if
      self % profScore(j,1) = self % profScore(j,1) / max(B, ONE)
    end do

  end subroutine finaliseSensitivity

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(adjointTRRMPhysicsPackage), intent(inout) :: self
    real(defFlt)                                  :: fissionRate, prevFissionRate, adjfissionRate, adjprevFissionRate
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
    end do
    !$omp end parallel do

    if (allocated(self % angularIP)) then
      !$omp parallel do schedule(static)
      do g = 1, size(self % angularIP)
        self % angularIP(g) = 0.0_defFlt
      end do
      !$omp end parallel do
    end if

    if (allocated(self % diagIP)) then
      !$omp parallel do schedule(static)
      do g = 1, size(self % diagIP)
        self % diagIP(g) = ZERO
      end do
      !$omp end parallel do
    end if


  end subroutine resetFluxes

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
    integer(shortInt)                             :: m, off, j, nProf
    real(defReal)                                 :: dk, dkSTD
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

    name = 'innerProduct'
    if (self % legacyIP) then
      call out % printValue('angular', name)
    else
      call out % printValue('scalar', name)
    end if

    if (self % legacyIP) then
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

    else
      ! Requested perturbation; uncertainties are absolute standard errors of the mean
      if (self % hasPert) then
        nProf = size(self % profMean)
        dk    = self % profMean(nProf) * self % keffScore(1) * self % XSchange
        dkSTD = self % profScore(nProf,2) * self % keffScore(1) * abs(self % XSchange)

        name = 'delta_keff'
        call out % startBlock(name)
        call out % printResult(dk, dkSTD, name)
        call out % endBlock()

        name = 'pert_keff'
        call out % startBlock(name)
        call out % printResult(self % keffScore(1) + dk, self % keffScore(1) - dk, name)
        call out % endBlock()

        name = 'sensitivity'
        call out % startBlock(name)
        call out % printResult(self % profMean(nProf), self % profScore(nProf,2), name)
        call out % endBlock()
      end if

      name = 'sensitivity_batches'
      call out % printValue(self % batchCount, name)

      ! Full profile: per material, capture(g), fission(g) and scatter(g_from, g_to)
      name = 'sensitivity_profile'
      call out % startBlock(name)
      do m = 1, self % nMat
        off  = (m - 1) * self % nPM
        name = mm_matName(m)
        call out % startBlock(name)

        resArrayShape = [self % nG]
        name = 'capture'
        call out % startArray(name, resArrayShape)
        do j = 1, self % nG
          call out % addResult(self % profMean(off + j), self % profScore(off + j, 2))
        end do
        call out % endArray()

        name = 'fission'
        call out % startArray(name, resArrayShape)
        do j = 1, self % nG
          call out % addResult(self % profMean(off + self % nG + j), self % profScore(off + self % nG + j, 2))
        end do
        call out % endArray()

        ! Written with g_to fastest, so JSON reads scatter[g_from][g_to] (as in the XS file P0 matrix)
        resArrayShape = [self % nG, self % nG]
        name = 'scatter'
        call out % startArray(name, resArrayShape)
        do g1 = 1, self % nG
          do j = 1, self % nG
            call out % addResult(self % profMean(off + 2 * self % nG + self % nG * (j - 1) + g1), &
                                 self % profScore(off + 2 * self % nG + self % nG * (j - 1) + g1, 2))
          end do
        end do
        call out % endArray()

        call out % endBlock()
      end do
      call out % endBlock()
    end if

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
    if (self % legacyIP) then
      print *, "Sensitivity inner product: angular (legacy nG x nG)"
    else
      print *, "Sensitivity inner product: scalar ("//numToChar(self % nBatches)//" batches)"
    end if
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
    if(allocated(self % fission)) deallocate(self % fission)

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

    self % NSegMax  = 0 
    self % XStype   = 0
    self % matPert  = 0
    self % XSchange = 0.0_defFlt

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
    if (allocated(self % diagIP)) deallocate(self % diagIP)
    if (allocated(self % diagIPScore)) deallocate(self % diagIPScore)
    if (allocated(self % adjFluxScore)) deallocate(self % adjFluxScore)
    if (allocated(self % profBatch)) deallocate(self % profBatch)
    if (allocated(self % profScore)) deallocate(self % profScore)
    if (allocated(self % profMean)) deallocate(self % profMean)
    self % legacyIP   = .false.
    self % hasPert    = .false.
    self % nPM        = 0
    self % nBatches   = 20
    self % batchLen   = 0
    self % batchIt    = 0
    self % batchCount = 0

    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % resultsMap)) then
      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if
    if(allocated(self % fluxMap)) then
      call self % fluxMap % kill()
      deallocate(self % fluxMap)
    end if

  end subroutine kill

end module adjointTRRMPhysicsPackage_class
