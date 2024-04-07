module anisotropicRRPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential, expTau
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
                                SQRT5 = sqrt(5._defReal), &
                                SQRT15 = sqrt(15._defReal), & 
                                !SQRT3_8 = sqrt(3._defReal/8._defReal), &
                                !SQRT5_8 = sqrt(5._defReal/8._defReal), &
                                SQRT70_4 = sqrt(70._defReal)/4._defReal, &
                                SQRT105 = sqrt(105._defReal), &
                                SQRT42_4 = sqrt(42._defReal)/4._defReal, &
                                SQRT7_2 = sqrt(7._defReal)/2._defReal


  ! Convenient arithmetic parameters
  real(defFlt), parameter :: one_two = real(HALF,defFlt), &
                                two_three = real(2.0_defFlt/3.0_defFlt,defFlt)

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
  !!     type anisotropicRRPhysicsPackage;
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
  !!     #SHOrder 0;#    // Optional spherical harmonic order for higher order scattering.
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
  !!   SHOrder     -> Order of higher order scattering, isotropic 0 by default, maximum P3 / 3rd order
  !!   SHLength    -> Number of spherical harmonics for given SHOrder.
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
  !!   moments     -> Array of angular moments, dimension = [nG * nCells, SHLength], [:,1] element are scalar fluxs.
  !!   prevMoments -> Array of previous angular moments, dimension = [nG * nCells, SHLength]
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source moment values, dimension = [nG * nCells, SHLength]
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
  type, public, extends(physicsPackage) :: anisotropicRRPhysicsPackage
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
    integer(shortInt)                     :: nMatVOID    = 0
    real(defReal)                         :: lengthPerIt = ZERO

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
    integer(shortInt)  :: SHOrder  = 0
    integer(shortInt)  :: SHLength = 1
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

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:,:), allocatable   :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi

    ! Results space
    real(defFlt)                               :: keff
    real(defReal), dimension(2)                :: keffScore
    real(defFlt), dimension(:,:), allocatable  :: moments
    real(defFlt), dimension(:,:), allocatable  :: prevMoments
    real(defFlt), dimension(:,:), allocatable  :: source
    real(defReal), dimension(:,:), allocatable :: fluxScores
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
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sphericalHarmonicCalculator
    procedure, private :: sourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings


  end type anisotropicRRPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    integer(shortInt)                             :: seed_temp, i, g, g1, m, SH
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
    character(100), parameter :: Here = 'init (anisotropicRRPhysicsPackage_class.f90)'

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

    ! Higher order scattering order
    call dict % getOrDefault(self % SHOrder, 'SHOrder', 0)

    ! Number of spherical harmonic terms for given SHOrder
    self % SHLength = (self % SHOrder + 1 )**2

    ! Stabilisation factor for negative in-group scattering for P0
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

    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % moments(self % nCells * self % nG, self % SHLength))
    allocate(self % prevMoments(self % nCells * self % nG, self % SHLength))
    allocate(self % source(self % nCells * self % nG, self % SHLength))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    
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
    self % nMatVOID = self % nMat + 1
    allocate(self % sigmaT(self % nMatVOID * self % nG))
    allocate(self % nuSigmaF(self % nMatVOID * self % nG))
    allocate(self % chi(self % nMatVOID * self % nG))
    allocate(self % sigmaS(self % nMatVOID * self % nG * self % nG, self % SHOrder + 1 ))

    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, self % rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, self % rand),defFlt)
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, self % rand),defFlt)
        ! Include scattering multiplicity
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1, 1)  = &
                real(mat % getScatterXS(g1, g, self % rand) * mat % scatter % prod(g1, g) , defFlt)
          if (self % SHOrder > 0) then
            do SH = 1, self % SHOrder
              self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1, SH + 1)  = &
                    real(mat % scatter % getPnScatter(g1, g, SH ) * mat % scatter % prod(g1, g) , defFlt)
            end do
          end if
        end do
      end do
    end do

    do g = 1, self % nG
        self % sigmaT(self % nG * (self % nMatVOID - 1) + g)   = 0.0_defFlt
        self % nuSigmaF(self % nG * (self % nMatVOID - 1) + g) = 0.0_defFlt
        self % chi(self % nG * (self % nMatVOID - 1) + g)      = 0.0_defFlt
        
        do g1 = 1, self % nG
            do SH = 1, self % SHOrder
                self % sigmaS(self % nG * self % nG * (self % nMatVOID - 1) &
                + self % nG * (g - 1) + g1, SH)  = 0.0_defFlt
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
    class(anisotropicRRPhysicsPackage), intent(inout) :: self

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
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: keepRunning, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise moments
    self % keff        = 1.0_defFlt
    self % moments     = 0.0_defFlt
    self % prevMoments = 1.0_defFlt

    self % fluxScores = ZERO
    self % keffScore  = ZERO
    self % source     = 0.0_defFlt

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

      ! Calculate new k
      call self % calculateKeff()

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
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (anisotropicRRPhysicsPackage_class.f90)'

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
    class(anisotropicRRPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, matIdx0, baseIdx, &
                                                             SH
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum, newRay
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec
    real(defFlt), dimension(self % nG)                    :: currentSource
    real(defFlt), pointer, dimension(:)                   :: totVec
    real(defReal), dimension(3)                           :: r0, mu0
    real(defFlt), dimension(self % SHLength)              :: RCoeffs   
    real(defFlt), pointer, dimension(:,:)                 :: sourceVec, angularMomVec

    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % source(idx,1)
    end do

    ints = 0
    matIdx0 = 0
    totalLength = ZERO
    activeRay = .false.
    newRay = .true.

    do while (totalLength < self % termination)

        ! Get material and cell the ray is moving through
        matIdx  = r % coords % matIdx
        cIdx    = r % coords % uniqueID

        if (matIdx >= VOID_MAT) then
            matIdx = self % nMatVOID
        end if

        if (matIdx0 /= matIdx) then 
          matIdx0 = matIdx
          ! Cache total cross section
          totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
        end if

        ! Caluclate r0 and mu0, mu0 used in SHarmonics
        r0 = r % rGlobal()
        mu0 = r % dirGlobal()

        if (newRay .or. event == BOUNDARY_EV) then
            call self % sphericalHarmonicCalculator(mu0, RCoeffs) 
            newRay = .false.
        end if
            
        ! Set maximum flight distance and ensure ray is active
        if (totalLength >= self % dead) then
            length = self % termination - totalLength 
            activeRay = .true.
        else
            length = self % dead - totalLength
        end if

        ! Include for limited maximum optical length

        ! maxtot = 0.0_defFlt
        ! !$omp simd
        ! do g = 1, self % nG
        !   if (maxtot < totVec(g)) then
        !     maxtot = totVec(g)
        !   end if
        ! end do
  
        ! if ((30.0_defFlt/maxtot) < length) then
        !   length = real(30.0_defFlt/maxtot)
        ! end if  

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
        sourceVec => self % source((baseIdx + 1):(baseIdx + self % nG), :)

        !$omp simd
        do g = 1, self % nG
            currentSource(g) = 0.0_defFlt
            do SH = 1, self % SHLength
              currentSource(g) = currentSource(g) + sourceVec(g,SH) * RCoeffs(SH)
            end do 
        end do

        !$omp simd aligned(totVec)
        do g = 1, self % nG
          tau(g) = totVec(g) * lenFlt
          if (tau(g) < 1E-8) then
            tau(g) = 0.0_defFlt
          end if
        end do

        !$omp simd 
        do g = 1, self % nG
            attenuate(g) = expTau(tau(g)) * lenFlt
            delta(g) = (fluxVec(g) - currentSource(g)) * attenuate(g)     
            fluxVec(g) = fluxVec(g) - delta(g) * totVec(g)
        end do


      ! Accumulate to scalar flux
      if (activeRay) then

        angularMomVec => self % moments((baseIdx + 1):(baseIdx + self % nG), :)
      
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd aligned(angularMomVec)
        do g = 1, self % nG
          do SH = 1, self % SHLength
            angularMomVec(g,SH)  = angularMomVec(g,SH) + delta(g) * RCoeffs(SH)
          end do  
        end do
        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
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

  subroutine SphericalHarmonicCalculator(self, mu, RCoeffs)
    ! angle: x = r sin θ cos φ, y = r sin θ sin φ, and z = r cos θ
    class(anisotropicRRPhysicsPackage), intent(inout)       :: self
    real(defFlt), dimension(self % SHLength), intent(out)  :: RCoeffs ! Array to store harmonic coefficients
    real(defReal)                                           :: dirX,dirY,dirZ
    real(defReal), dimension(3), intent(in)                 :: mu


    dirX = mu(1)
    dirY = mu(2)
    dirZ = mu(3)

    ! Assign coefficients based on SHOrder
    select case(self % SHLength)
    case(1)  
        RCoeffs(1) = 1.0_defFlt 

    case(4)  
        RCoeffs(1) = 1.0_defFlt 

        RCoeffs(2) = real(SQRT3 * dirY,defFlt)
        RCoeffs(3) = real(SQRT3 * dirZ,defFlt)
        RCoeffs(4) = real(SQRT3 * dirX,defFlt)    

    case(9) 
        RCoeffs(1) = 1.0_defFlt 

        RCoeffs(2) = real(SQRT3 * dirY,defFlt)
        RCoeffs(3) = real(SQRT3 * dirZ,defFlt)
        RCoeffs(4) = real(SQRT3 * dirX,defFlt) 
        RCoeffs(5) = real(SQRT15 * HALF * dirX * dirY,defFlt) 
        RCoeffs(6) = real(SQRT15 * HALF * dirZ * dirY,defFlt) 
        RCoeffs(7) = real(SQRT5 * HALF * (3 * dirZ**2 - 1),defFlt)
        RCoeffs(8) = real(SQRT15 * HALF * dirX * dirZ,defFlt) 
        RCoeffs(9) = real(SQRT15 * HALF * (dirX**2 - dirY**2),defFlt) 

    case(16) 
        RCoeffs(1) = 1.0_defFlt 

        RCoeffs(2) = real(SQRT3 * dirY,defFlt)
        RCoeffs(3) = real(SQRT3 * dirZ,defFlt)
        RCoeffs(4) = real(SQRT3 * dirX,defFlt) 
        RCoeffs(5) = real(SQRT15 * HALF * dirX * dirY,defFlt)
        RCoeffs(6) = real(SQRT15 * HALF * dirZ * dirY,defFlt) 
        RCoeffs(7) = real(SQRT5 * HALF * (3 * dirZ**2 - 1),defFlt)
        RCoeffs(8) = real(SQRT15 * HALF * dirX * dirZ,defFlt) 
        RCoeffs(9) = real(SQRT15 * HALF * (dirX**2 - dirY**2),defFlt) 

        RCoeffs(10) = real(SQRT70_4 * dirY * (3 * dirX**2 - dirY**2),defFlt)
        RCoeffs(11) = real(SQRT105 * dirZ * dirX * dirY,defFlt)
        RCoeffs(12) = real(SQRT42_4 * dirY * (5 * dirZ**2 - 1),defFlt)
        RCoeffs(13) = real(SQRT7_2 * dirZ * (5 * dirZ**2 - 3),defFlt)
        RCoeffs(14) = real(SQRT42_4 * dirX * (5 * dirZ**2 - 1),defFlt)
        RCoeffs(15) = real(SQRT105 * dirZ * (dirX**2 - dirY**2),defFlt)
        RCoeffs(16) = real(SQRT70_4 * dirX * (dirX**2 - 3 * dirY**2),defFlt)
    end select
  end subroutine SphericalHarmonicCalculator


  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defFlt), save                            :: total, vol, D, SigGG
    integer(shortInt), save                       :: g, idx, SH, matIdx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, vol, idx, g, SH, matIdx, sigGG, D)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
        matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 

        if (matIdx >= VOID_MAT - 1) then
            matIdx = self % nMatVOID
        end if

        ! Update volume due to additional rays
        self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
        vol = real(self % volume(cIdx),defFlt)

        do g = 1, self % nG
            idx   = self % nG * (cIdx - 1) + g

            if (self % SHOrder > 0) then
              
              do SH = 1, self % SHLength
                if (vol > volume_tolerance) then
                    self % moments(idx,SH) =  self % moments(idx,SH) * norm / vol
                end if
                self % moments(idx,SH) =  self % moments(idx,SH) + self % source(idx,SH) 
              end do 

            else 

              total = self % sigmaT((matIdx - 1) * self % nG + g)
              sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g, 1)
    
              ! Presumes non-zero total XS
              if ((sigGG < 0) .and. (total > 0)) then
                D = -real(self % rho, defFlt) * sigGG / total
              else
                D = 0.0_defFlt
              end if

              if (vol > volume_tolerance) then
                self % moments(idx,1) =  self % moments(idx,1) * norm / vol
              end if

              self % moments(idx,1) =  (self % moments(idx,1) + self % source(idx,1) + D * self % prevMoments(idx,1) ) / (1 + D)
            
            end if

        end do


    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(anisotropicRRPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, fission
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, SH, SHidx
    real(defFlt), pointer, dimension(:,:)                 :: scatterVec, scatterXS
    real(defFlt), pointer, dimension(:,:)                :: angularMomVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= VOID_MAT - 1) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        do SH = 1, self % SHLength
          self % source(idx,SH) = 0.0_defFlt
        end do
      end do
      return
    end if

    ! Obtain XSs
    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG), :)
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % nG * (cIdx - 1)

    angularMomVec => self % prevMoments((baseIdx + 1):(baseIdx + self % nG), :)

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission) aligned(angularMomVec)
    do gIn = 1, self % nG
      fission = fission + angularMomVec(gIn,1) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG), :)
      idx = baseIdx + g

      do SH = 1, self % SHLength
        ! Calculate scattering source for higher order scattering
        SHidx = ceiling(sqrt(real(SH, defReal)) - 1) + 1

        scatter = 0.0_defFlt

        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter) aligned(angularMomVec)
        do gIn = 1, self % nG
          scatter = scatter + angularMomVec(gIn, SH) * scatterVec(gIn, SHidx)
        end do

        self % source(idx,SH) = scatter / total(g)

      end do

      ! Calculate scattering source for isotropic scattering / flat source

      self % source(idx,1) = self % source(idx,1) + chi(g) * fission / total(g)

  end do

  end subroutine sourceUpdateKernel

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
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
      if (matIdx >= VOID_MAT - 1) cycle

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
        fissLocal     = fissLocal     + self % moments(idx,1) * self % nuSigmaF(mIdx + g)
        prevFissLocal = prevFissLocal + self % prevMoments(idx,1) * self % nuSigmaF(mIdx + g)

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
  !!
  subroutine resetFluxes(self)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    integer(shortInt)                                 :: idx, SH

   !$omp parallel do schedule(static)
    do idx = 1, (self % nG * self % nCells)
      !$omp simd
      do SH = 1, self % SHLength
        self % prevMoments(idx,SH) = self % moments(idx,SH) 
        self % moments(idx,SH) = 0.0_defFlt
      end do

    end do
   !$omp end parallel do

  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
    real(defReal), save                            :: flux
    integer(shortInt)                              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, (self % nG * self % nCells)
      flux = real(self % moments(idx,1),defReal)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux * flux
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff

  end subroutine accumulateFluxAndKeffScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it)
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
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
    do idx = 1, (self % nG * self % nCells)
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
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
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
        if (matIdx >= VOID_MAT - 1) cycle
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
    class(anisotropicRRPhysicsPackage), intent(in) :: self

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
    class(anisotropicRRPhysicsPackage), intent(inout) :: self
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
    self % nMatVOID  = 0
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nusigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % chi)) deallocate(self % chi)

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    self % cache       = .false.
    self % SHLength    = 1
    self % SHOrder     = 0
    self % rho         = ZERO
    self % mapFission  = .false.
    self % mapFlux     = .false.
    self % plotResults = .false.
    self % printFlux   = .false.
    self % printVolume = .false.
    self % printCells  = .false.

    self % keff        = ZERO
    self % keffScore   = ZERO

    if(allocated(self % moments)) deallocate(self % moments)
    if(allocated(self % prevMoments)) deallocate(self % prevMoments)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
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

end module anisotropicRRPhysicsPackage_class

