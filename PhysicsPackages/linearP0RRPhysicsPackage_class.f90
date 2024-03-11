module linearP0RRPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential, expTau, expH, expG2, expG, expF2
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
                                SQRT70_4 = sqrt(70._defReal)/4._defReal, &
                                SQRT105 = sqrt(105._defReal), &
                                SQRT42_4 = sqrt(42._defReal)/4._defReal, &
                                SQRT7_2 = sqrt(7._defReal)/2._defReal, &
                                SQRT1_PI = sqrt(1.0_defReal/PI) * HALF

  ! Convenient arithmetic parameters
  real(defFlt), parameter :: one_two = real(HALF,defFlt), &
                                two_three = real(2.0_defFlt/3.0_defFlt,defFlt)
  
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
  !!     type linearP0RRPhysicsPackage;
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
  type, public, extends(physicsPackage) :: linearP0RRPhysicsPackage
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

    real(shortInt)     :: count       = 0 

    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
    integer(shortInt)  :: SHOrder  = 0
    integer(shortInt)  :: SHLength = 1
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

    ! Results space - Linear source terms
    real(defFlt), dimension(:), allocatable    :: scalarX
    real(defFlt), dimension(:), allocatable    :: scalarY
    real(defFlt), dimension(:), allocatable    :: scalarZ
    real(defFlt), dimension(:), allocatable    :: prevX
    real(defFlt), dimension(:), allocatable    :: prevY
    real(defFlt), dimension(:), allocatable    :: prevZ
    real(defFlt), dimension(:), allocatable    :: sourceX
    real(defFlt), dimension(:), allocatable    :: sourceY
    real(defFlt), dimension(:), allocatable    :: sourceZ
    real(defReal), dimension(:), allocatable   :: momMat
    real(defReal), dimension(:), allocatable   :: momTracks
    real(defReal), dimension(:), allocatable   :: centroid
    real(defReal), dimension(:), allocatable   :: centroidTracks

    real(defReal), dimension(:), allocatable :: tauVec
    real(defReal), dimension(:), allocatable :: errorVec 

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


  end type linearP0RRPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
    character(100), parameter :: Here = 'init (linearP0RRPhysicsPackage_class.f90)'

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
    self % SHLength = ( self % SHOrder + 1 )**2

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

    allocate(self % scalarX(self % nCells * self % nG))
    allocate(self % scalarY(self % nCells * self % nG))
    allocate(self % scalarZ(self % nCells * self % nG))
    allocate(self % prevX(self % nCells * self % nG))
    allocate(self % prevY(self % nCells * self % nG))
    allocate(self % prevZ(self % nCells * self % nG))
    allocate(self % sourceX(self % nCells * self % nG))
    allocate(self % sourceY(self % nCells * self % nG))
    allocate(self % sourceZ(self % nCells * self % nG))

    allocate(self % momMat(self % nCells * matSize))
    allocate(self % momTracks(self % nCells * matSize))
    allocate(self % centroid(self % nCells * nDim))
    allocate(self % centroidTracks(self % nCells * nDim))


    allocate(self % tauVec(1000))
    allocate(self % errorVec(1000))
    
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
                    real(mat % scatter % getPnScatter(g1, g, SH) * mat % scatter % prod(g1, g) , defFlt)
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
    class(linearP0RRPhysicsPackage), intent(inout) :: self

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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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

    self % tauVec = 0.0_defReal
    self % errorVec = 0.0_defReal

    self % count = 0_shortInt
      
    self % scalarX  = 0.0_defFlt
    self % scalarY  = 0.0_defFlt
    self % scalarZ  = 0.0_defFlt
    self % prevX    = 0.0_defFlt
    self % prevY    = 0.0_defFlt
    self % prevZ    = 0.0_defFlt
    self % sourceX  = 0.0_defFlt
    self % sourceY  = 0.0_defFlt
    self % sourceZ  = 0.0_defFlt

    self % fluxScores = ZERO
    self % keffScore  = ZERO
    self % source     = 0.0_defFlt

    ! Initialise other results
    self % cellHit      = 0
    self % volume       = ZERO
    self % volumeTracks = ZERO
    self % intersectionsTotal  = 0
 
    self % momMat         = 0.0_defReal
    self % momTracks      = 0.0_defReal
    self % centroid       = 0.0_defReal
    self % centroidTracks = 0.0_defReal

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
        call self % sourceUpdateKernel(i, ONE_KEFF, it)
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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (linearP0RRPhysicsPackage_class.f90)'

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
    class(linearP0RRPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r  
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, &
                                                             matIdx0, baseIdx, centIdx, momIdx, SH
    real(defReal)                                         :: totalLength, length, len2_12
    logical(defBool)                                      :: activeRay, hitVacuum, newRay
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt, lenFlt2_2, maxtot, vol, norm
    real(defFlt), dimension(self % nG)                    :: F1, F2, G1, G2, Gn, H, tau, delta, fluxVec, &
                                                             flatQ, gradQ, xInc, yInc, zInc, fluxVec0, &
                                                             currentSource, H1
    real(defReal), dimension(matSize)                     :: matScore
    real(defFlt),  pointer, dimension(:)                  :: totVec, xGradVec, yGradVec, zGradVec, &
                                                             xMomVec, yMomVec, zMomVec
    real(defFlt),  pointer, dimension(:,:)                :: sourceVec
    real(defFlt),  pointer, dimension(:,:)                :: angularMomVec
    real(defReal), pointer, dimension(:)                  :: mid, momVec, centVec
    real(defReal), pointer                                :: volTrack
    real(defReal), dimension(3)                           :: r0, mu0, rC, r0Norm, rNorm
    real(defFlt), dimension(self % SHLength)              :: RCoeffs   
    real(defFlt),  dimension(3)                           :: muFlt, rNormFlt, r0NormFlt

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

      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID

      if (matIdx >= VOID_MAT - 1) then
        matIdx = self % nMatVOID
      end if

      ! Get material and cell the ray is moving through
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        ! Cache total cross section (if not void)
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):(matIdx * self % nG))
      end if

      !always caluclate r0 and mu0, mu0 used in SHarmonics
      r0 = r % rGlobal()
      mu0 = r % dirGlobal()

      ! if (newRay .or. event == BOUNDARY_EV) then
        call self % sphericalHarmonicCalculator(mu0, RCoeffs) !- seems to work 
      !   newRay = .false.
      ! end if
          
      ! Set maximum flight distance and ensure ray is active
      if (totalLength >= self % dead) then
        length = self % termination - totalLength 
        activeRay = .true.
      else
        length = self % dead - totalLength
      end if

      maxtot = 0.0_defFlt
      !$omp simd
      do g = 1, self % nG
        if (maxtot < totVec(g)) then
          maxtot = totVec(g)
        end if
      end do

      if ((30.0_defFlt/maxtot) < length) then
        length = real(30.0_defFlt/maxtot)
      end if
      
      ! Move ray
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
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = rC
        !$omp end critical
      end if

      ints    = ints + 1
      baseIdx = (cIdx - 1) * self % nG

      sourceVec => self % source((baseIdx + 1):(baseIdx + self % nG), :)
      xGradVec => self % sourceX((baseIdx + 1):(baseIdx + self % nG))
      yGradVec => self % sourceY((baseIdx + 1):(baseIdx + self % nG))
      zGradVec => self % sourceZ((baseIdx + 1):(baseIdx + self % nG))
      mid => self % centroid(((cIdx - 1) * nDim + 1):(cIdx * nDim))

      ! Check cell has been visited 
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
      lenFlt  = real(length,defFlt)

      ! Calculates source for higher order moments.
      !$omp simd
      do g = 1, self % nG
        currentSource(g) = 0.0_defFlt
        do SH = 1, self % SHLength
          currentSource(g) = currentSource(g) + sourceVec(g,SH) * RCoeffs(SH) 
        end do 
      end do
      
      ! Calculate source terms
      ! aligned(xGradVec, yGradVec, zGradVec)
      !$omp simd
      do g = 1, self % nG
        flatQ(g) = rNormFlt(x) * xGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(y) * yGradVec(g)
        flatQ(g) = flatQ(g) + rNormFlt(z) * zGradVec(g)
        flatQ(g) = flatQ(g) + currentSource(g)

        gradQ(g) = muFlt(x) * xGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(y) * yGradVec(g)
        gradQ(g) = gradQ(g) + muFlt(z) * zGradVec(g)
      end do

      !$omp simd
      do g = 1, self % nG
        tau(g) = totVec(g) * lenFlt
        if (tau(g) < 1E-8) then 
          tau(g) = 0.0_defFlt
        end if
      end do

      ! Compute exponentials necessary for angular flux update
      !$omp simd
      do g = 1, self % nG
        Gn(g) = expG(tau(g))
      end do

      ! !$omp simd
      ! do g = 1, self % nG
      !   F1(g) = expTau(tau(g)) !* lenFlt !(1 - tau(g) * Gn(g)) * lenFlt ! (1.0_defFlt - tau(g) * Gn(g)) !* lenFlt !
      ! end do

      ! !$omp simd
      ! do g = 1, self % nG
      !   F2(g) = (2.0_defFlt * Gn(g) - F1(g))
      ! end do


     !$omp simd
      do g = 1, self % nG
        F1(g)  = expTau(tau(g)) * lenFlt
      end do


      !$omp simd
      do g = 1, self % nG
        F2(g) = (2 * Gn(g) * lenFlt * lenFlt - F1(g) * lenFlt )  
      end do

      !$omp simd
      do g = 1, self % nG
        delta(g) = (fluxVec(g) - flatQ(g)) * F1(g) - &
                   one_two * gradQ(g) * F2(g) 
      end do

      ! Intermediate flux variable creation or update
      !$omp simd
      do g = 1, self % nG
        fluxVec0(g) = fluxVec(g)
      end do

      ! Flux vector update
      !$omp simd
      do g = 1, self % nG
        fluxVec(g) = fluxVec(g) - delta(g) * totVec(g) 
      end do

      ! Accumulate to scalar flux
      if (activeRay) then

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
        lenFlt2_2 = lenFlt * lenFlt * one_two


        !$omp simd
        do g = 1, self % nG
          H(g) = expH(tau(g)) !* lenFlt!( F1(g) - Gn(g) ) * lenFlt !expH(tau(g)) * lenFlt !
          !G1(g) = ( one_two - H(g) ) !expH(tau(g)) * lenFlt !
       end do

        !$omp simd
        do g = 1, self % nG
          G2(g) = expG2(tau(g)) !* lenFlt
        end do
        

        !for LS-PN
        ! !$omp simd 
        ! do g = 1, self % nG
        !       G1(g) = G1(g) * flatQ(g) * lenFlt
        !       G2(g) = G2(g) * gradQ(g) * lenFlt2_2
        !       H(g)  = H(g) * fluxVec0(g) * lenFlt
        !       H(g)  = (G1(g) + G2(g) + H(g)) * lenFlt
        !       flatQ(g) = flatQ(g) * lenFlt + delta(g)
        ! end do

        ! !$omp simd
        ! do g = 1, self % nG
        !       xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) 
        !       yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) 
        !       zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g)
        ! end do

        !$omp simd 
        do g = 1, self % nG
          G2(g) = G2(g) * gradQ(g) * lenFlt2_2 !* lenFlt
          H1(g) = H(g) 
          H(g)  = H(g) * fluxVec0(g) * lenFlt !* lenFlt !* FOUR_PI
        end do

        !$omp simd
        do g = 1, self % nG
          xInc(g) = r0NormFlt(x) * delta(g) + muFlt(x) * lenFlt * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) &
                                * rNormFlt(x) * lenFlt - H1(g) * muFlt(x) * lenFlt * currentSource(g)

          yInc(g) = r0NormFlt(y) * delta(g) + muFlt(y) * lenFlt  * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) &
                                * rNormFlt(y) * lenFlt - H1(g) * lenFlt * muFlt(y) * lenFlt * currentSource(g)

          zInc(g) = r0NormFlt(z) * delta(g) + muFlt(z) * lenFlt  * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) &
                                * rNormFlt(z) * lenFlt -  H1(g) * lenFlt * muFlt(z) * lenFlt * currentSource(g)
        end do


        !for LS-P0, no higher order scattering
        !!$omp simd 
        ! do g = 1, self % nG
        ! !   !G1(g) = G1(g) * flatQ(g) * lenFlt
        ! !   !G2(g) = G2(g) * gradQ(g) * lenFlt2_2
        ! !   !H(g)  = H(g) * fluxVec0(g) * lenFlt
        ! !   !H(g) = G1(g) + G2(g) + H(g)
        ! !   !flatQ(g) = flatQ(g) * lenFlt + delta(g)
        ! !flatQ(g) = flatQ(g) - currentSource(g) + sourceVec(g,1)
        !   G2(g) = G2(g) * gradQ(g) * lenFlt2_2 * lenFlt
        !   H(g)  = H(g) * fluxVec0(g) * lenFlt * lenFlt
        ! end do

        ! !$omp simd
        ! do g = 1, self % nG
        !   !xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) * lenFlt
        !   !yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) * lenFlt
        !   !zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g) * lenFlt
        !   xInc(g) = r0NormFlt(x) * delta(g) + muFlt(x) * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) * rNormFlt(x) &
        !                         * lenFlt !* totVec(g)
        !   yInc(g) = r0NormFlt(y) * delta(g) + muFlt(y) * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) * rNormFlt(y) &
        !                         * lenFlt !* totVec(g)
        !   zInc(g) = r0NormFlt(z) * delta(g) + muFlt(z) * (H(g) + G2(g)) + (flatQ(g) - currentSource(g)) * rNormFlt(z) &
        !                         * lenFlt !* totVec(g)
        ! end do
                
        call OMP_set_lock(self % locks(cIdx))

        angularMomVec => self % moments((baseIdx + 1):(baseIdx + self % nG), :)
        xMomVec => self % scalarX((baseIdx + 1):(baseIdx + self % nG))
        yMomVec => self % scalarY((baseIdx + 1):(baseIdx + self % nG))
        zMomVec => self % scalarZ((baseIdx + 1):(baseIdx + self % nG))

        

        ! Update flux moments !aligned(angularMomVec, xMomVec, yMomVec, zMomVec)
        !$omp simd 
        do g = 1, self % nG
          xMomVec(g) = xMomVec(g) + xInc(g) 
          yMomVec(g) = yMomVec(g) + yInc(g) 
          zMomVec(g) = zMomVec(g) + zInc(g) 

          do SH = 1, self % SHLength
              angularMomVec(g, SH) = angularMomVec(g, SH) + delta(g) * RCoeffs(SH) 
          end do

        end do
        
        centVec => self % centroidTracks((centIdx + 1):(centIdx + nDim))
        momVec => self % momTracks((momIdx + 1):(momIdx + matSize))
        volTrack => self % volumeTracks(cIdx)

        volTrack = volTrack + length

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

  subroutine SphericalHarmonicCalculator(self, mu, RCoeffs)
    ! angle: x = sin θ cos φ, y = sin θ sin φ, and z = cos θ
    class(linearP0RRPhysicsPackage), intent(inout)          :: self
    real(defFlt), dimension(self % SHLength), intent(out)   :: RCoeffs ! Array to store harmonic coefficients
    real(defReal)                                           :: dirX,dirY,dirZ
    real(defReal), dimension(3), intent(in)                 :: mu

    dirX = mu(1)
    dirY = mu(2)
    dirZ = mu(3)

    select case(self % SHLength)
    case(1)  
        RCoeffs(1) = 1.0_defFlt !* SQRT1_PI

    case(4)  
        RCoeffs(1) = 1.0_defFlt !* SQRT1_PI

        RCoeffs(2) = real(SQRT3 * dirY,defFlt) !* SQRT1_PI
        RCoeffs(3) = real(SQRT3 * dirZ,defFlt) !* SQRT1_PI
        RCoeffs(4) = real(SQRT3 * dirX,defFlt) !* SQRT1_PI 
        ! RCoeffs(2) = real(dirY,defFlt) !* SQRT1_PI
        ! RCoeffs(3) = real(dirZ,defFlt) !* SQRT1_PI
        ! RCoeffs(4) = real(dirX,defFlt) !* SQRT1_PI 

    end select

  end subroutine SphericalHarmonicCalculator

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(linearP0RRPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defReal), save                           :: invVol
    real(defFlt), save                            :: vol, NTV
    integer(shortInt), save                       :: g, idx, SH, matIdx,  dIdx, mIdx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(vol, idx, g, SH, matIdx, mIdx, dIdx, invVol, NTV)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)

      if (matIdx >= VOID_MAT - 1) then
        matIdx = self % nMatVOID
      end if 

      dIdx = (cIdx - 1) * nDim
      mIdx = (cIdx - 1) * matSize

      ! Update volume due to additional rays
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

        do SH = 1, self % SHLength

          if (vol > volume_tolerance) then
              self % moments(idx,SH) = self % moments(idx,SH) * NTV
          end if
          self % moments(idx,SH) =  self % moments(idx,SH) + self % source(idx,SH) 

        end do

        if (vol > volume_tolerance) then
          self % scalarX(idx) = self % scalarX(idx) * NTV 
          self % scalarY(idx) = self % scalarY(idx) * NTV 
          self % scalarZ(idx) = self % scalarZ(idx) * NTV 
        end if

      end do

      ! if (self % moments(idx,1) < 0.0_defFlt) then
      !   !self % moments(idx,1) = 0.0_defFlt
      !   do SH = 1, self % SHLength
      !     self % moments(idx,SH) = 0.0_defFlt
      !   end do
      ! end if
      ! if (self % scalarX(idx) < 0.0_defFlt) then 
      !   self % scalarX(idx) = 0.0_defFlt
      ! end if
      ! if (self % scalarY(idx) < 0.0_defFlt) then 
      !   self % scalarY(idx) = 0.0_defFlt
      ! end if
      ! if (self % scalarZ(idx) < 0.0_defFlt) then 
      !   self % scalarZ(idx) = 0.0_defFlt
      ! end if

    end do
    !$omp end parallel do


  end subroutine normaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF, it)
    class(linearP0RRPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    integer(shortInt), intent(in)                         :: it
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, xScatter, yScatter, zScatter, &
                                                             fission, xFission, yFission, zFission, &
                                                             xSource, ySource, zSource
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, SH, SHidx
    real(defFlt), pointer, dimension(:,:)                 :: scatterVec, scatterXS
    real(defFlt), pointer, dimension(:,:)                 :: angularMomVec
    real(defFlt)                                          :: invMxx, invMxy, invMxz, invMyy, invMyz, invMzz
    integer(shortInt)                                     :: condX, condY, condZ, inversionTest
    real(defFlt), pointer, dimension(:)                   :: xFluxVec, yFluxVec, zFluxVec
    real(defReal), pointer, dimension(:)                  :: momVec
    real(defReal)                                         :: det, one_det 

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx)

    ! Guard against void cells
    if (matIdx >= VOID_MAT - 1) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        do SH = 1, self % SHLength
          self % source(idx,SH) = 0.0_defFlt
        end do
        self % sourceX(idx) = 0.0_defFlt
        self % sourceY(idx) = 0.0_defFlt
        self % sourceZ(idx) = 0.0_defFlt
      end do
      return
    end if

    ! condX = 0
    ! condY = 0
    ! condZ = 0

    momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize))

    ! if (momVec(xx) > 1.0E-6_defReal) condX = 1
    ! if (momVec(yy) > 1.0E-6_defReal) condY = 1
    ! if (momVec(zz) > 1.0E-6_defReal) condZ = 1

    ! inversionTest = condX * 4 + condY * 2 + condZ
    
    ! select case(inversionTest)
    ! case(invertXYZ)
    !   det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
    !         - momVec(yy) * momVec(xz) * momVec(xz) - momVec(zz) * momVec(xy) * momVec(xy) &
    !         + 2 * momVec(xy) * momVec(xz) * momVec(yz)
    !   one_det = ONE/det
    !   invMxx = real(one_det * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)),defFlt)
    !   invMxy = real(one_det * (momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz)),defFlt)
    !   invMxz = real(one_det * (momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz)),defFlt)
    !   invMyy = real(one_det * (momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)),defFlt)
    !   invMyz = real(one_det * (momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz)),defFlt)
    !   invMzz = real(one_det * (momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)),defFlt)

    ! case(invertYZ)
    !   det = momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)
    !   one_det = ONE/det
    !   invMxx = 0.0_defFlt
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = real(one_det * momVec(zz),defFlt)
    !   invMyz = real(-one_det * momVec(yz),defFlt)
    !   invMzz = real(one_det * momVec(yy),defFlt)

    ! case(invertXY)
    !   det = momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)
    !   one_det = ONE/det
    !   invMxx = real(one_det * momVec(yy),defFlt)
    !   invMxy = real(-one_det * momVec(xy),defFlt)
    !   invMxz = 0.0_defFlt
    !   invMyy = real(one_det * momVec(xx),defFlt)
    !   invMyz = 0.0_defFlt
    !   invMzz = 0.0_defFlt

    ! case(invertXZ)
    !   det = momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)
    !   one_det = ONE/det
    !   invMxx = real(one_det * momVec(zz),defFlt)
    !   invMxy = 0.0_defFlt
    !   invMxz = real(-one_det * momVec(xz),defFlt)
    !   invMyy = 0.0_defFlt
    !   invMyz = 0.0_defFlt
    !   invMzz = real(one_det * momVec(xx),defFlt)

    ! case(invertX)
    !   det = momVec(xx)
    !   one_det = ONE/det
    !   invMxx = real(one_det,defFlt)
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = 0.0_defFlt
    !   invMyz = 0.0_defFlt
    !   invMzz = 0.0_defFLt

    ! case(invertY)
    !   det = momVec(yy)
    !   one_det = ONE/det
    !   invMxx = 0.0_defFlt
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = real(one_det,defFlt)
    !   invMyz = 0.0_defFlt
    !   invMzz = 0.0_defFlt

    ! case(invertZ)
    !   det = momVec(zz)
    !   one_det = ONE/det
    !   invMxx = 0.0_defFlt
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = 0.0_defFlt
    !   invMyz = 0.0_defFlt
    !   invMzz = real(one_det,defFlt)

    ! case default
    !   invMxx = 0.0_defFlt
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = 0.0_defFlt
    !   invMyz = 0.0_defFlt
    !   invMzz = 0.0_defFlt
    !   det = ONE
    ! end select

    ! ! ! Check for zero determinant
    ! if (abs(det) < 1E-6) then !was 1E-6
    ! !if (one_det /= one_det) then !was 1E-6
    !   invMxx = 0.0_defFlt
    !   invMxy = 0.0_defFlt
    !   invMxz = 0.0_defFlt
    !   invMyy = 0.0_defFlt
    !   invMyz = 0.0_defFlt
    !   invMzz = 0.0_defFlt
    ! end if

    det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
    - momVec(yy) * momVec(xz) * momVec(xz) - momVec(zz) * momVec(xy) * momVec(xy) &
    + 2 * momVec(xy) * momVec(xz) * momVec(yz)

    if ((abs(det) < 1E-10) .and. (self % volume(cIdx) > 1E-6)) then
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


    matIdx = (matIdx - 1) * self % nG
    total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
    scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG), :)
    nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
    chi => self % chi((matIdx + 1):(matIdx + self % nG))

    baseIdx = self % nG * (cIdx - 1)
    angularMomVec => self % prevMoments((baseIdx + 1):(baseIdx + self % nG), :)
    xFluxVec => self % prevX((baseIdx + 1):(baseIdx + self % nG))
    yFluxVec => self % prevY((baseIdx + 1):(baseIdx + self % nG))
    zFluxVec => self % prevZ((baseIdx + 1):(baseIdx + self % nG))


    ! Calculate fission source
    fission  = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !reduction(+:fission) aligned(angularMomVec)
    !$omp simd 
    do gIn = 1, self % nG
      fission  = fission  + angularMomVec(gIn,1) * nuFission(gIn)
      xFission = xFission + xFluxVec(gIn) * nuFission(gIn)
      yFission = yFission + yFluxVec(gIn) * nuFission(gIn)
      zFission = zFission + zFluxVec(gIn) * nuFission(gIn)
    end do

    fission  = fission  * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    do g = 1, self % nG

      scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG), :)
      idx = baseIdx + g

      xScatter = 0.0_defFlt
      yScatter = 0.0_defFlt
      zScatter = 0.0_defFlt 

      ! Calculate scattering source for higher order scattering

      !reduction(+:scatter) aligned(angularMomVec)
      !$omp simd 
      do gIn = 1, self % nG
        xScatter = xScatter + xFluxVec(gIn) * scatterVec(gIn,1)
        yScatter = yScatter + yFluxVec(gIn) * scatterVec(gIn,1)
        zScatter = zScatter + zFluxVec(gIn) * scatterVec(gIn,1)
      end do
    
      do SH = 1, self % SHLength

        SHidx = ceiling(sqrt(real(SH, defReal)) - 1) + 1
        scatter = 0.0_defFlt

        !$omp simd 
        do gIn = 1, self % nG
          scatter = scatter + angularMomVec(gIn, SH) * scatterVec(gIn, SHidx)
        end do

        self % source(idx,SH) = scatter / total(g) !* ONE_FOUR_PI

      end do

      ! Calculate scattering source for isotropic scattering / flat source
      self % source(idx,1) = self % source(idx,1) + (chi(g) * fission) / total(g) !* ONE_FOUR_PI

      !if (it > 29) then !.and. self % source(idx,1) > 1E-25
        xSource = chi(g) * xFission + xScatter
        xSource = xSource / total(g) !* ONE_FOUR_PI
        ySource = chi(g) * yFission + yScatter
        ySource = ySource / total(g) !* ONE_FOUR_PI
        zSource = chi(g) * zFission + zScatter
        zSource = zSource / total(g) !* ONE_FOUR_PI

        ! Calculate source gradients by inverting the moment matrix
        self % sourceX(idx) = invMxx * xSource + &
                invMxy * ySource + invMxz * zSource 
        self % sourceY(idx) = invMxy * xSource + & 
                invMyy * ySource + invMyz * zSource 
        self % sourceZ(idx) = invMxz * xSource + &
                invMyz * ySource + invMzz * zSource 
      ! else
      !   self % sourceX(idx) = 0.0_defFlt
      !   self % sourceY(idx) = 0.0_defFlt
      !   self % sourceZ(idx) = 0.0_defFlt
      ! end if

  end do

  end subroutine sourceUpdateKernel

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
    integer(shortInt)                              :: idx, SH

   !$omp parallel do schedule(static)
    do idx = 1, (self % nG * self % nCells)
      do SH = 1, self % SHLength
        self % prevMoments(idx,SH) = self % moments(idx,SH) 
        self % moments(idx,SH) = 0.0_defFlt
        ! self % source(idx,SH) = 0.0_defFlt
      end do
      self % prevX(idx) = self % scalarX(idx)
      self % scalarX(idx) = 0.0_defFlt
      self % prevY(idx) = self % scalarY(idx)
      self % scalarY(idx) = 0.0_defFlt
      self % prevZ(idx) = self % scalarZ(idx)
      self % scalarZ(idx) = 0.0_defFlt
      ! self % sourceX(idx) = 0.0_defFlt
      ! self % sourceY(idx) = 0.0_defFlt
      ! self % sourceZ(idx) = 0.0_defFlt
      
    end do
    !$omp end parallel do

  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
        if (matIdx >= VOID_MAT - 1) cycle
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
    class(linearP0RRPhysicsPackage), intent(in) :: self

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
    class(linearP0RRPhysicsPackage), intent(inout) :: self
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
    self % SHLength  = 1
    self % SHOrder   = 0
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
    if(allocated(self % resultsMap)) then
      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if
    if(allocated(self % fluxMap)) then
      call self % resultsMap % kill()
      deallocate(self % fluxMap)
    end if

  end subroutine kill

end module linearP0RRPhysicsPackage_class
