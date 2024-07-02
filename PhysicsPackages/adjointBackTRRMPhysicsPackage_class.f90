module adjointBackTRRMPhysicsPackage_class

    use numPrecision
    use universalVariables
    use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
    use hashFunctions_func,             only : FNV_1
    use charMap_class,                  only : charMap
    use exponentialRA_func,             only : exponential, f1
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
    use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                               gr_geomIdx  => geomIdx, gr_kill    => kill
  
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
    use particle_class,                 only : ray => particle, particleState
    
    ! For locks
    use omp_lib
    
    implicit none
    private
  
    ! Parameter for when to skip a tiny volume
    real(defReal), parameter :: volume_tolerance = 1.0E-12
  
    ! Parameter for uncollided flux calculations
    integer(shortInt), parameter :: NO_UC = 0, POINT = 1, VOLUME = 2
  
    !!
    !! Physics package to perform The Random Ray Method (TRRM) fixed source calculations
    !!
    !! Modified to find which cells are actually present in the geometry as a precalculation
    !! and establish a map. Does this in a modified volume sweep.
    !!
    !! Tracks rays across the geometry, attenuating their flux. After some dead length,
    !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
    !! uniform termination length, after which it is stopped and the next ray is tracked.
    !! Once all rays have been tracked, a cycle concludes and fluxes and sources are updated.
    !!
    !! Both inactive and active cycles occur, as in eigenvalue calculations. These can be terminated
    !! after a specified number of iterations or on reaching some chosen convergence
    !! criterion (though the latter hasn't been implemented yet).
    !!
    !! Calculates relative volume of diffrent materials in the problem by performing
    !! random ray tracing in the geometry. The volume is normalised such that the total domain
    !! volume is 1.0.
    !!
    !! Can be performed with uncollided flux calculations for small or point sources.
    !! Requires an additional input and, for correct normalisation, specifying the volume of the problem.
    !!
    !! IMPORTANT N.B.: Geometry type must be extended! If shrunk, results may be dubious!
    !! This is because spatial discretisation is determined by the number of unique cells in the
    !! geometry.
    !! Also, this is obviously for multi-group calculations only.
    !!
    !! Sample Input Dictionary:
    !!   PP {
    !!     type randomRayPhysicsPackage;
    !!     dead 10;              // Dead length where rays do not score to scalar fluxes
    !!     termination 100;      // Length a ray travels before it is terminated
    !!     rays 1000;            // Number of rays to sample per iteration
    !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
    !!     active 200;           // Number of scoring cycles (would use eps otherwise)
    !!     #seed 86868;#         // Optional RNG seed
    !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
    !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
    !!
    !!     #source {             // Fixed sources for named materials and their intensities n/cm3/s
    !!                           // Intensities are in each energy group, from 1 to G
    !!         material_name1 ( s_g1 s_g2 ... s_gG );
    !!         ...
    !!      } #
    !!     #integrate ( names of materials );#
    !!     geometry {<Geometry definition>}
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
    !!   nMat        -> Number of materials in the geometry, kept for convenience.
    !!   lengthPerIt -> Distance all rays travel in a single iteration - for convenience.
    !!   normVolume  -> Volume of the geometry - used for normalisation with uncollided calculations.
    !!
    !!   termination   -> Distance a ray can travel before it is terminated
    !!   dead          -> Distance a ray must travel before it becomes active
    !!   skipLength    -> Optional distance over which to skip, i.e., not perform any operations.
    !!                    Can be used to avoid deleterious effects from tiny volumes.
    !!   pop           -> Number of rays to track per cycle
    !!   inactive      -> Number of inactive cycles to perform
    !!   active        -> Number of active cycles to perform
    !!   rho           -> Stabilisation factor for negative XSs (see Gunow)
    !!   cache         -> Perform distance caching?
    !!   itVol         -> Use the volume calculated during a given iteration? (rather than accumulated)
    !!   zeroNeg       -> Zero negative fluxes?
    !!   volCorr       -> Use a volume correction on the q/Sigma_T term?
    !!   passive       -> Volume correct only negative fluxes? 
    !!   outputFile    -> Output file name
    !!   outputFormat  -> Output file format
    !!   plotResults   -> Plot results?
    !!   printFluxes   -> Print fluxes?
    !!   printVolume   -> Print volumes?
    !!   printCells    -> Print cell positions?
    !!   viz           -> Output visualiser
    !!   samplePoints  -> Set of points from which to sample fluxes in output
    !!   sampleNames   -> Names of the corresponding sampled points
    !!   sourceIdx     -> material indices of source materials, kept for convenience
    !!   intMatIdx     -> material indices of materials over which to integrate, for convenience
    !!   intMatName    -> Name of materials over which to integrate
    !!   mapFlux       -> Integrate flux over elements of a map?
    !!   fluxMap       -> Corresponding map over which to integrate fluxes
    !!
    !!   nVolRays      -> Number of rays to use for volume precalculation
    !!   volLength     -> Length of each volume ray to trace
    !!
    !!   uncollidedType   -> Type of uncollided calculation to perform (point, volumetric, none)
    !!   uncollidedPop    -> Number of uncollided rays to trace
    !!   uncollidedCycles -> Number of cycles over which to batch uncollided stats
    !!   uncollidedLength -> Length of each uncollided ray to trace
    !!   UCNorm           -> Normalisation factor for uncollided fluxes (calculated using physical volume)
    !!   sourcePoint      -> Point source locations
    !!   sourceTop        -> Bounding box top for uncollided source
    !!   sourceBottom     -> Bounding box bottom for uncollided source
    !!   sourceStrength   -> Strength of initial source in all groups
    !!   uncollidedScores -> Stats for calculating average and std of uncollided flux
    !!
    !!   scalarFlux   -> Array of scalar flux values of length = nG * nCells
    !!   prevFlux     -> Array of previous scalar flux values of length = nG * nCells
    !!   fluxScores   -> Array of scalar flux values and squared values to be reported 
    !!                   in results, dimension =  [nG * nCells, 2]
    !!   source       -> Array of neutron source values of length = nG * nCells
    !!   fixedSource  -> Array of fixed source values of length = nG * nCells
    !!   volume       -> Array of stochastically estimated cell volumes of length = nCells
    !!   cellHit      -> Array tracking whether given cells have been hit during tracking
    !!   cellFound    -> Array tracking whether a cell was ever found
    !!   cellPos      -> Array of cell positions, populated once they are found
    !!   IDToCell     -> Map between uniqueIDs and cells present in the geometry (when using cell remapping)
    !!   CellToID     -> Map between cells present and uniqueIDs in the geometry (when using cell remapping)
    !!
    !!   sigmaT       -> Total XS vector in all materials
    !!   nuSigmaF     -> NuFission XS vector in all materials
    !!   sigmaS       -> Scatter matrix XS vector in all materials
    !!   chi          -> Fission spectrum vector in all materials
    !!
    !! Interface:
    !!   physicsPackage interface
    !!
    type, public, extends(physicsPackage) :: adjointBackTRRMPhysicsPackage
      private
      ! Components
      class(geometryStd), pointer           :: geom
      integer(shortInt)                     :: geomIdx     = 0
      real(defReal), dimension(3)           :: top         = ZERO
      real(defReal), dimension(3)           :: bottom      = ZERO
      class(baseMgNeutronDatabase), pointer :: mgData      => null()
      type(RNG)                             :: rand
      integer(shortInt)                     :: nG          = 0
      integer(shortInt)                     :: nCells      = 0
      integer(shortInt)                     :: nMat        = 0
      real(defReal)                         :: lengthPerIt = ZERO
      real(defReal)                         :: normVolume  = ONE
  
      ! Settings
      real(defReal)      :: termination = ZERO
      real(defReal)      :: dead        = ZERO
      real(defReal)      :: skipLength  = ZERO
      integer(shortInt)  :: pop         = 0
      integer(shortInt)  :: inactive    = 0
      integer(shortInt)  :: active      = 0
      real(defReal)      :: rho         = ZERO

      !Pertubation
      integer(shortInt), dimension(:), allocatable   :: energyId 
      integer(shortInt)  :: XStype      = 0
      integer(shortInt)  :: matPert     = 0
      real(defReal)      :: XSchange    = ZERO
      real(defReal)      :: numSum      = 0
      real(defReal)      :: denSum      = 0

      real(defReal)                 :: sensitivity = ZERO
      real(defReal), dimension(2)   :: sensitivityScore = ZERO
      
      logical(defBool)   :: cache       = .false.
      logical(defBool)   :: itVol       = .false.
      logical(defBool)   :: zeroNeg     = .false.
      logical(defBool)   :: volCorr     = .false.
      logical(defBool)   :: passive     = .false.
      character(pathLen) :: outputFile
      character(nameLen) :: outputFormat
      logical(defBool)   :: plotResults = .false.
      logical(defBool)   :: printFlux   = .false.
      logical(defBool)   :: printVolume = .false.
      logical(defBool)   :: printCells  = .false.
      type(visualiser)   :: viz
      integer(shortInt)   :: mapResponse  = 0
      class(tallyMap), allocatable :: resultsMap
      real(defReal), dimension(:), allocatable      :: samplePoints
      character(nameLen), dimension(:), allocatable :: sampleNames
      integer(shortInt), dimension(:), allocatable  :: sourceIdx
      integer(shortInt) ,dimension(:), allocatable  :: intMatIdx
      character(nameLen),dimension(:), allocatable  :: intMatName
      integer(shortInt) ,dimension(:), allocatable  :: detectorInt
      character(nameLen),dimension(:), allocatable  :: detectorName
      logical(defBool)   :: mapFlux     = .false.
      class(tallyMap), allocatable :: fluxMap
  
      ! Volume calculation settings
      integer(shortInt)  :: nVolRays = 0
      real(defReal)      :: volLength = ZERO
  
      ! Uncollided flux calculation settings
      integer(shortInt)           :: uncollidedType = NO_UC
      integer(shortInt)           :: uncollidedPop = 0
      integer(shortInt)           :: uncollidedCycles = 0
      real(defReal)               :: uncollidedLength = ZERO
      real(defReal)               :: UCNorm = ONE
      real(defReal), dimension(3) :: sourcePoint  = ZERO
      real(defReal), dimension(3) :: sourceTop    = ZERO
      real(defReal), dimension(3) :: sourceBottom = ZERO
      real(defReal), dimension(:), allocatable   :: sourceStrength
      real(defReal), dimension(:,:), allocatable :: uncollidedScores
  
      ! Data space - absorb all nuclear data for speed
      real(defFlt), dimension(:), allocatable     :: sigmaT
      real(defFlt), dimension(:), allocatable     :: nuSigmaF
      real(defFlt), dimension(:), allocatable     :: sigmaS
      real(defFlt), dimension(:), allocatable     :: chi
      real(defFlt), dimension(:), allocatable     :: adjNuSigmaF
      real(defFlt), dimension(:), allocatable     :: adjSigmaS
      real(defFlt), dimension(:), allocatable     :: adjChi
  
      ! Results space
      real(defFlt), dimension(:), allocatable     :: scalarFlux
      real(defFlt), dimension(:), allocatable     :: prevFlux
      real(defReal), dimension(:,:), allocatable  :: fluxScores
      real(defFlt), dimension(:), allocatable     :: source
      real(defFlt), dimension(:), allocatable     :: fixedSource
      real(defReal), dimension(:), allocatable    :: volume
      real(defReal), dimension(:), allocatable    :: volumeTracks

      real(defFlt), dimension(:), allocatable    :: adjScalarFlux
      real(defFlt), dimension(:), allocatable    :: adjPrevFlux
      real(defFlt), dimension(:), allocatable    :: adjSource

      real(defReal), dimension(:), allocatable   :: angularIP
  
      ! Tracking cell properites
      integer(shortInt), dimension(:), allocatable :: IDToCell
      integer(shortInt), dimension(:), allocatable :: CellToID
      integer(shortInt), dimension(:), allocatable :: cellHit
      logical(defBool), dimension(:), allocatable  :: cellFound
      real(defReal), dimension(:,:), allocatable   :: cellPos
      
      ! OMP locks
      integer(kind=omp_lock_kind), dimension(:), allocatable :: locks
  
      ! Timer bins
      integer(shortInt) :: timerMain
      integer(shortInt) :: timerTransport
      integer(shortInt) :: timerUC
      real (defReal)    :: time_transport   = ZERO
      real (defReal)    :: time_volume      = ZERO
      real (defReal)    :: time_UC          = ZERO
      real (defReal)    :: time_transportUC = ZERO
      real (defReal)    :: time_cell = ZERO
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
      procedure, private :: initialiseSource
      procedure, private :: initialiseRay
      procedure, private :: transportSweep
      procedure, private :: sourceUpdateKernel
      procedure, private :: adjointSourceUpdateKernel
      procedure, private :: firstCollidedSourceKernel
      procedure, private :: normaliseFluxAndVolume
      procedure, private :: adjointNormaliseFluxAndVolume
      procedure, private :: normaliseFluxUncollided
      procedure, private :: resetFluxes
      procedure, private :: sensitivityCalculation
      procedure, private :: accumulateFluxScores
      procedure, private :: finaliseFluxScores
      procedure, private :: printResults
      procedure, private :: printSettings
      procedure, private :: volumeSweep
      procedure, private :: volumeCalculation
      procedure, private :: cellMapCalculation
      procedure, private :: uncollidedSweep
      procedure, private :: uncollidedCalculation
  
    end type adjointBackTRRMPhysicsPackage
  
  contains
  
    !!
    !! Initialise Physics Package from dictionary
    !!
    !! See physicsPackage_inter for details
    !!
    subroutine init(self,dict)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      class(dictionary), intent(inout)                    :: dict
      integer(shortInt)                                   :: seed_temp, n, nPoints, i, m, g, g1
      integer(longInt)                                    :: seed
      character(10)                                       :: time
      character(8)                                        :: date
      character(:),allocatable                            :: string
      class(dictionary),pointer                           :: tempDict, graphDict, sourceDict
      real(defReal), dimension(:), allocatable            :: tempArray
      class(mgNeutronDatabase),pointer                    :: db
      character(nameLen)                                  :: geomName, graphType, nucData
      character(nameLen),dimension(:), allocatable        :: names
      class(charMap), pointer                             :: matMap
      class(geometry), pointer                            :: geom
      type(outputFile)                                    :: test_out
      class(baseMgNeutronMaterial), pointer               :: mat
      class(materialHandle), pointer                      :: matPtr
      logical(defBool)                                    :: cellCheck
      character(100), parameter :: Here = 'init (adjointBackTRRMPhysicsPackage_class.f90)'
  
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
      
      ! Perform distance caching?
      call dict % getOrDefault(self % skipLength, 'skip', ZERO)
      
      ! Use iteration-wise volumes?
      call dict % getOrDefault(self % itVol, 'itVol', .false.)
  
      ! Apply volume correction?
      call dict % getOrDefault(self % volCorr, 'volCorr', .false.)
      if (self % volCorr .and. self % itVol) call fatalError(Here,&
              'Can only apply either volume correction or iteration-wise volumes')
  
      ! Volume correct only negative fluxes?
      call dict % getOrDefault(self % passive, 'passive', .false.)
      
      ! Zero negative fluxes?
      call dict % getOrDefault(self % zeroNeg, 'zeroNeg', .false.)
  
      ! Stabilisation factor for negative in-group scattering
      call dict % getOrDefault(self % rho, 'rho', ZERO)
      ! Shouldn't need this... but need to do a bit more work to allow these to be done together
      if (self % volCorr) self % rho = ZERO
  
      ! Print fluxes?
      call dict % getOrDefault(self % printFlux, 'printFlux', .false.)
  
      ! Print volumes?
      call dict % getOrDefault(self % printVolume, 'printVolume', .false.)
  
      ! Print cell positions?
      call dict % getOrDefault(self % printCells, 'printCells', .false.)
  
      ! Read normalised volume by which to scale dimensionless volume estimates
      call dict % getOrDefault(self % normVolume, 'volume', ONE)

      ! Check response map
      call dict % getOrDefault(self % mapResponse, 'responseType', 0)

      !!!!!!!!!!!!!!!!!!!!!!
      tempDict => dict % getDictPtr("Pertubation")
      call self % initPertubation(tempDict)
      
      ! Check whether there is a map for outputting mapped fluxes
      ! If so, read and initialise the map to be used
      if (dict % isPresent('fluxMap')) then
        self % mapFlux = .true.
        tempDict => dict % getDictPtr('fluxMap')
        call new_tallyMap(self % fluxMap, tempDict)
      else
        self % mapFlux = .false.
      end if

      ! Check whether there is a map for outputting reaction rates
      ! If so, read and initialise the map to be used
      if (dict % isPresent('responseType')) then
        tempDict => dict % getDictPtr('responseMap')
        call new_tallyMap(self % resultsMap, tempDict)
      end if

      ! Return flux values at sample points?
      ! Store a set of points to return values at on concluding the simulation
      if (dict % isPresent('samplePoints')) then
  
        tempDict => dict % getDictPtr('samplePoints')
        call tempDict % keys(self % sampleNames)
        nPoints = size(self % sampleNames)
        allocate(self % samplePoints(3*nPoints))
        do n = 1, nPoints
  
          call tempDict % get(tempArray, self % sampleNames(n))
          if (size(tempArray) /= 3) call fatalError(Here,&
                 'Sample points must be 3 dimensional')
          self % samplePoints((1+3*(n-1)):(3*n)) = tempArray
  
        end do
  
      end if
  
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
      if (self % termination <= self % dead) call fatalError(Here,'Ray termination length must be greater than ray dead length')
  
      ! Register timer
      self % timerMain = registerTimer('simulationTime')
      self % timerTransport = registerTimer('transportTime')
      self % timerUC = registerTimer('uncollidedTime')
  
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
      call gr_addGeom(geomName, tempDict)
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
  
      ! Activate nuclear data
      call ndReg_activate(P_NEUTRON_MG, nucData, self % geom % activeMats())
  
      ! Ensure that nuclear data is multi-group
      db => ndReg_getNeutronMG()
      if (.NOT. associated(db)) call fatalError(Here,&
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
  
      ! Check for materials to integrate over
      if (dict % isPresent('integrate')) then
        call dict % get(names,'integrate')
        
        allocate(self % intMatIdx(size(names)))
        allocate(self % intMatName(size(names)))
        self % intMatName = names
        
        matMap => self % mgData % matNamesMap()
  
        ! Check that materials exist in the geometry
        ! and remember the matIdx
        do i = 1,size(names)
          self % intMatIdx(i) = matMap % get(names(i))
        end do
      end if

      if (dict % isPresent('detector')) then
        call dict % get(names,'detector')
        
        allocate(self % detectorInt(size(names)))
        allocate(self % detectorName(size(names)))
        self % detectorName = names
        
        matMap => self % mgData % matNamesMap()
  
        ! Check that materials exist in the geometry
        ! and remember the matIdx
        do i = 1,size(names)
          self % detectorInt(i) = matMap % get(names(i))
        end do
      end if
      
      ! Store number of cells in geometry for convenience
      self % nCells = self % geom % numberOfCells()
      allocate(self % IDToCell(self % nCells))
      allocate(self % CellToID(self % nCells))
      do i = 1, self % nCells
        self % IDToCell(i) = i
        self % CellToID(i) = i
      end do
  
      ! Check how many cells are actually present?
      call dict % getOrDefault(cellCheck, 'cellCheck', .false.)
      
      ! Perform a calculation to see which cells are present in the geometry
      ! Establishes a map between uniqueIDs in the geometry graph and the
      ! cell IDs of cells which are present
      if (cellCheck) then
        print *,'Performing cell remapping'
        ! Allocate these for convenience... Delete them after
        allocate(self % cellHit(self % nCells))
        allocate(self % cellFound(self % nCells))
        allocate(self % cellPos(self % nCells,3))
        call self % cellMapCalculation()
        deallocate(self % cellHit)
        deallocate(self % cellFound)
        deallocate(self % cellPos)
      end if
      
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
  
      ! Allocate results space
      allocate(self % scalarFlux(self % nCells * self % nG))
      allocate(self % prevFlux(self % nCells * self % nG))
      allocate(self % fluxScores(self % nCells * self % nG, 2))
      allocate(self % source(self % nCells * self % nG))
      allocate(self % fixedSource(self % nCells * self % nG))
      allocate(self % volume(self % nCells))
      allocate(self % volumeTracks(self % nCells))
      allocate(self % cellHit(self % nCells))
      allocate(self % cellFound(self % nCells))
      allocate(self % cellPos(self % nCells, 3))

      allocate(self % adjScalarFlux(self % nCells * self % nG))
      allocate(self % adjPrevFlux(self % nCells * self % nG))
      allocate(self % adjSource(self % nCells * self % nG))
      
      self % scalarFlux = 0.0_defFlt
      self % prevFlux = 0.0_defFlt
      self % fluxScores = ZERO
      self % source = 0.0_defFlt
      self % fixedSource = 0.0_defFlt
      self % volume = ZERO
      self % volumeTracks = ZERO
      self % cellHit = 0
      self % cellFound = .false.
      self % cellPos = -INFINITY
      
      self % adjScalarFlux = 0.0_defFlt
      self % adjPrevFlux   = 1.0_defFlt
      self % adjSource     = 0.0_defFlt
  
      ! Check whether to precompute volumes
      call dict % getOrDefault(self % nVolRays,'volRays',0)
      if (self % nVolRays > 0) call dict % get(self % volLength, 'volLength')
      
      ! Perform uncollided flux treatment?
      if (dict % isPresent('uncollided')) then
      
        if (self % nVolRays <= 0) call fatalError(Here,'Cannot perform uncollided '//&
                'calculations without precomputing cell volumes')
  
        allocate(self % uncollidedScores(self % nCells * self % nG,2))
        self % uncollidedScores = ZERO
        tempDict => dict % getDictPtr('uncollided')
        
        ! Note 1 = point, 2 = volume
        call tempDict % get(self % uncollidedType, 'type' )
  
        call tempDict % get(self % uncollidedPop, 'pop')
        call tempDict % get(self % uncollidedLength, 'length')
        call tempDict % get(self % uncollidedCycles, 'cycles')
  
        select case(self % uncollidedType)
        
        ! Read a point from which to sample
        case(POINT) 
          call tempDict % get(tempArray,'point')
          if (size(tempArray) /= 3) call fatalError(Here,&
                 'Uncollided point must be 3 dimensional')
          self % sourcePoint = tempArray
          allocate(self % sourceStrength(self % nG))
          call tempDict % get(tempArray, 'strength')
          if (size(tempArray) /= self % nG) call fatalError(Here,&
                  'Point source strength must have as many entries as there are energy groups')
          self % sourceStrength = tempArray
          self % UCNorm = ONE / self % normVolume
  
        ! Read a bounding box from which to sample
        case(VOLUME)
          call tempDict % get(tempArray,'top')
          if (size(tempArray) /= 3) call fatalError(Here,&
                 'Uncollided bounding box top must be 3 dimensional')
          self % sourceTop = tempArray
          call tempDict % get(tempArray,'bottom')
          if (size(tempArray) /= 3) call fatalError(Here,&
                 'Uncollided bounding box bottom must be 3 dimensional')
          self % sourceBottom = tempArray
         
          ! Must provide the RELATIVE volume of the source for normalisation 
          call tempDict % get(self % UCNorm, 'relVolume')
  
        case default
          call fatalError(Here,'Unrecognised uncollided source type: should be point or volume')
        end select
  
      end if
  
      ! Point sources cannot have cells with fixed sources in them
      ! For volume sources with uncollided, this will have to be overwritten subsequently
      if (self % uncollidedType /= POINT) then
        sourceDict => dict % getDictPtr('source')
        call self % initialiseSource(sourceDict)
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
          do g1 = 1, self % nG
            self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                    real(mat % getScatterXS(g1, g, self % rand) * mat % scatter % prod(g1, g) , defFlt)
  
          end do
        end do
      end do

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
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
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
    !! Initialises the fixed source to be used in the simulation
    !! Takes a dictionary containing names of materials in the geometry and
    !! source strengths in each energy group and places these in the appropriate
    !! elements of the fixed source vector
    !!
    !! Also sets options for uncollided flux calculations
    !!
    subroutine initialiseSource(self, dict)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      class(dictionary), intent(inout)                    :: dict
      character(nameLen),dimension(:), allocatable        :: names
      real(defReal), dimension(:), allocatable            :: sourceStrength
      integer(shortInt)                                   :: i, nSource, cIdx
      integer(shortInt), save                             :: g, matIdx, idx, id
      logical(defBool)                                    :: found
      character(nameLen)                                  :: sourceName 
      character(nameLen), save                            :: localName
      character(100), parameter :: Here = 'initialiseSource (adjointBackTRRMPhysicsPackage_class.f90)'
      !$omp threadprivate(matIdx, localName, idx, g, id)
  
      call dict % keys(names)
  
      nSource = size(names)
      allocate(self % sourceIdx(nSource))
  
      ! Cycle through entries of the dictionary
      do i = 1, nSource
  
        sourceName = names(i)
        call dict % get(sourceStrength, sourceName)
  
        ! Ensure correct number of energy groups
        if (size(sourceStrength) /= self % nG) call fatalError(Here,'Source '//sourceName//&
                ' has '//numToChar(size(sourceStrength))//' groups rather than '//numToChar(self % nG))
        
        ! Make sure that the source corresponds to a material present in the geometry
        found = .false.
        !$omp parallel do schedule(static) 
        do cIdx = 1, self % nCells
  
          id        = self % CellToID(cIdx)
          matIdx    = self % geom % geom % graph % getMatFromUID(id)
          localName = mm_matName(matIdx)
  
          if (localName == sourceName) then
  
            if (.not. found) then
              !$omp critical
              self % sourceIdx(i) = matIdx
              !$omp end critical
            end if
  
            found = .true.
            do g = 1, self % nG
              idx = (cIdx - 1) * self % nG + g
              self % fixedSource(idx) = real(sourceStrength(g),defFlt)
            end do
  
          end if
  
        end do
        !$omp end parallel do
  
        !print*, sum(self % fixedSource), 'fixed source cells'
  
        if (.not. found) call fatalError(Here,'The source '//trim(sourceName)//' does not correspond to '//&
                'any material found in the geometry.')
  
      end do
  
    end subroutine initialiseSource

  
    !!
    !! Run calculation
    !!
    !! See physicsPackage_inter for details
    !!
    subroutine run(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
  
      call self % printSettings()
      if (self % nVolRays > 0) call self % volumeCalculation()
      if (self % uncollidedType > NO_UC) call self % uncollidedCalculation()
      call self % cycles()
      call self % printResults()
  
    end subroutine run
    
    !!
    !! Calculates which cells are present in the geometry by ray tracing
    !!
    !! Randomly places the ray starting point and direction uniformly.
    !! Rays are tracked until they reach some specified termination length.
    !!
    subroutine cellMapCalculation(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), save                                     :: r
      type(RNG), target, save                             :: pRNG
      real(defReal)                                       :: hitRate
      integer(shortInt)                                   :: i, id, cIdx
      logical(defBool)                                    :: doVolume
      !$omp threadprivate(pRNG, r)
  
      ! Reset and start timer for volume calculation
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
  
      self % cellHit = 0
      self % cellFound = .false.
      self % cellPos = -INFINITY
      doVolume = .false.
        
      !$omp parallel do schedule(dynamic)
      do i = 1, self % pop * 20 
        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 
          
        call self % initialiseRay(r)
        call self % volumeSweep(r, self % termination, doVolume)
  
      end do
      !$omp end parallel do
   
      hitRate = real(sum(self % cellHit),defReal) / self % nCells
  
      ! Go through cell hit and determine which cells are present in the geometry
      deallocate(self % CellToID)
      allocate(self % CellToID(count(self % cellFound)))
      self % CellToID(:) = -1
      self % IDToCell(:) = -1
      cIdx = 0
      !$omp parallel do
      do id = 1, self % nCells
        if (self % cellFound(id)) then
          !$omp critical
          cIdx = cIdx + 1
          self % IDToCell(id) = cIdx
          self % CellToID(cIdx) = id
          !$omp end critical
        end if
      end do
      !$omp end parallel do
      self % nCells = cIdx
  
      call timerStop(self % timerTransport)
      self % time_cell = timerTime(self % timerTransport)
        
      ! Print/save volume calculation time
      call printFishLineR(1)
      print *, 'Cell map calculation complete'
      print *, 'Cell hit rate: ', trim(numToChar(hitRate))
      print *, 'Cell map calculation time: ', trim(secToChar(self % time_cell))
  
      ! Update RNG 
      call self % rand % stride(self % pop * 200 + 1)
      
    end subroutine cellMapCalculation
    
    !!
    !! Calculates volumes in the geometry by ray tracing
    !!
    !! Randomly places the ray starting point and direction uniformly.
    !! Rays are tracked until they reach some specified termination length.
    !! During tracking, fluxes are attenuated (and adjusted according to BCs),
    !! scoring to volume estimates.
    !!
    subroutine volumeCalculation(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), save                                     :: r
      type(RNG), target, save                             :: pRNG
      real(defReal)                                       :: hitRate
      integer(shortInt)                                   :: i
      logical(defBool)                                    :: doVolume
      !$omp threadprivate(pRNG, r)
  
      ! Reset and start timer for volume calculation
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
  
      doVolume = .true.
        
      !$omp parallel do schedule(dynamic)
      do i = 1, self % nVolRays 
        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 
          
        call self % initialiseRay(r)
        call self % volumeSweep(r, self % volLength, doVolume)
  
      end do
      !$omp end parallel do
   
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        self % volume(i) = self % volumeTracks(i) /(self % nVolRays * self % volLength)
      end do
      !$omp end parallel do
      call timerStop(self % timerTransport)
  
      hitRate = real(sum(self % cellHit),defReal) / self % nCells
      self % cellHit = 0
      self % time_volume = timerTime(self % timerTransport)
        
      ! Print/save volume calculation time
      call printFishLineR(1)
      print *, 'Volume calculation complete'
      print *, 'Cell hit rate: ', trim(numToChar(hitRate))
      print *, 'Volume calculation time: ', trim(secToChar(self % time_volume))
  
      ! Update RNG 
      call self % rand % stride(self % nVolRays + 1)
      
    end subroutine volumeCalculation
  
    !!
    !! Calculates uncollided flux in the geometry by MoC
    !!
    !! Simulates many rays traversing the geometry starting from particular
    !! locations. 
    !! Rays are tracked until they reach some specified termination length.
    !! During tracking, fluxes are attenuated (and adjusted according to BCs).
    !!
    subroutine uncollidedCalculation(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), save                                     :: r
      type(RNG), target, save                             :: pRNG
      real(defReal)                                       :: hitRate
      real(defReal)                                       :: elapsed_T, end_T, T_toEnd, transport_T
      integer(shortInt)                                   :: i, it
      integer(longInt), save                              :: ints
      integer(longInt)                                    :: intersections
      !$omp threadprivate(pRNG, r, ints)
  
      ! Reset and start timer for uncollided flux
      call timerReset(self % timerUC)
      call timerStart(self % timerUC)
          
      do it = 1, self % uncollidedCycles
  
        ! Reset and start timer for uncollided flux
        call timerReset(self % timerTransport)
        call timerStart(self % timerTransport)
        intersections = 0
        
        !$omp parallel do schedule(dynamic) reduction(+: intersections)
        do i = 1, self % uncollidedPop 
            
          ! Set seed
          pRNG = self % rand
          call pRNG % stride(i)
          r % pRNG => pRNG 
      
          call self % uncollidedSweep(r, ints)
          intersections = intersections + ints
  
        end do
        !$omp end parallel do
        call timerStop(self % timerTransport)
  
        ! Normalise flux estimate 
        ! Assumes volume has already been calculated!
        call self % normaliseFluxUncollided(self % UCNorm / self % uncollidedPop)
     
        call self % accumulateFluxScores()
        call self % resetFluxes()
  
        hitRate = real(sum(self % cellHit),defReal) / self % nCells
        self % cellHit = 0
        
        ! Calculate times
        call timerStop(self % timerUC)
        elapsed_T = timerTime(self % timerUC)
        transport_T = timerTime(self % timerTransport)
        self % time_transportUC = self % time_transportUC + transport_T
  
        ! Predict time to end
        end_T = real(self % uncollidedCycles, defReal) * elapsed_T / it
        T_toEnd = max(ZERO, end_T - elapsed_T)
          
        ! Print/save uncollided flux calculation time
        call printFishLineR(it)
        print *
        print *, 'Uncollided flux calculation iteration '//numToChar(it)//' of '//numToChar(self % uncollidedCycles)
        print *, 'Uncollided flux calculation time: ', trim(secToChar(elapsed_T))
        print *, 'End time:     ', trim(secToChar(end_T))
        print *, 'Time to end:  ', trim(secToChar(T_toEnd))
        print *, 'Time per integration (ns): ', &
            trim(numToChar(transport_T*10**9/(self % nG * intersections)))
        print *, 'Cell hit rate: ', trim(numToChar(hitRate))
          
        ! Update RNG
        call self % rand % stride(self % uncollidedPop + 1)
  
      end do
  
      self % time_UC = elapsed_T
  
      print *, 'Uncollided flux calculation complete '
  
      ! Finalise flux scores
      call self % finaliseFluxScores(self % uncollidedCycles)
      self % uncollidedScores = self % fluxScores
      
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        call self % firstCollidedSourceKernel(i)
      end do
      !$omp end parallel do
      
      self % fluxScores = ZERO
      self % prevFlux = 0.0_defFlt
  
    end subroutine uncollidedCalculation
  
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
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), save                                     :: r
      type(RNG), target, save                             :: pRNG
      real(defReal)                                       :: hitRate
      real(defReal)                                       :: elapsed_T, end_T, T_toEnd, transport_T
      logical(defBool)                                    :: stoppingCriterion, isActive
      integer(shortInt)                                   :: itInac, itAct, it, p, nCells, i, idx
      integer(longInt), save                              :: ints
      integer(longInt)                                    :: intersections
      integer(shortInt), dimension(:), allocatable        :: cellOnes
      !$omp threadprivate(pRNG, r, ints)
  
      ! Reinitialise volumes
      ! Perhaps it is worth doing something more clever if volumes have been precomputed...
      self % volumeTracks = ZERO
      self % volume = ZERO
  
      ! Update the cell number after several iterations
      ! Allows for better diagnostics on ray coverage
      nCells = self % nCells
      allocate(cellOnes(self % nCells))
      cellOnes = 1
  
      ! TODO: Come up with a better stopping criterion. Currently just number of iterations
      itInac = 0
      itAct  = 0
      isActive = (self % inactive < 1)
      stoppingCriterion = .true.
  
      print *,'Starting calculation'
  
      call timerReset(self % timerMain)
      call timerStart(self % timerMain)
      
      ! Source iteration
      do while( stoppingCriterion )
        
        if (isActive) then
          itAct = itAct + 1
        else
          itInac = itInac + 1
        end if
        it = itInac + itAct
  
        if (it == 20) nCells = sum(cellOnes,MASK=self % cellFound)
  
        !$omp parallel do schedule(static)
        do i = 1, self % nCells
          call self % sourceUpdateKernel(i)
          call self % adjointSourceUpdateKernel(i)
        end do
        !$omp end parallel do
       
        ! Reset and start transport timer
        call timerReset(self % timerTransport)
        call timerStart(self % timerTransport)
        intersections = 0
        
        !$omp parallel do schedule(dynamic) reduction(+: intersections)
        do p = 1, self % pop
  
          ! Set seed
          pRNG = self % rand
          call pRNG % stride(p)
          r % pRNG => pRNG 
  
          ! Set ray attributes
          call self % initialiseRay(r)
          
          ! Transport ray until termination criterion met
          call self % transportSweep(r,ints)
          intersections = intersections + ints
        
        end do
        !$omp end parallel do
        
        call timerStop(self % timerTransport)
  
        ! Update RNG 
        call self % rand % stride(self % pop + 1)
  
        ! Normalise flux estimate and combines with source
        call self % normaliseFluxAndVolume(self % lengthPerIt, it)
        call self % adjointNormaliseFluxAndVolume(self % lengthPerIt, it)
        
        ! Accumulate flux scores
        if (isActive) then
          self % numSum = ZERO
          self % denSum = ZERO
          !$omp parallel do schedule(static)
          do i = 1, self % nCells
           call self % sensitivityCalculation(i, it)
          end do
          self % sensitivity = (self % numSum / self % denSum) 

          call self % accumulateFluxScores()
        end if
  
        ! Calculate proportion of cells that were hit
        hitRate = real(sum(self % cellHit),defReal) / nCells
        self % cellHit = 0
  
        ! Evaluate stopping criterion for active or inactive iterations
        if (isActive) then
          stoppingCriterion = (itAct < self % active)
        else
          isActive = (itInac >= self % inactive)
        end if
  
        ! Set previous iteration flux to scalar flux and zero scalar flux
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
        print *, 'Cell hit rate: ', trim(numToChar(hitRate))
        print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
        print *, 'End time:     ', trim(secToChar(end_T))
        print *, 'Time to end:  ', trim(secToChar(T_toEnd))
        print *, 'Time per integration (ns): ', &
                trim(numToChar(transport_T*10**9/(self % nG * intersections)))
  
      end do
  
      ! Finalise flux scores
      call self % finaliseFluxScores(itAct)
  
      ! Add collided and uncollided results
      if (self % uncollidedType > NO_UC) then
        !$omp parallel do schedule(static)
        do idx = 1, size(self % adjscalarFlux)
          self % fluxScores(idx,2) = sqrt(self % fluxScores(idx,2)**2 * self % fluxScores(idx,1)**2 + &
                self % uncollidedScores(idx,2)**2 * self % uncollidedScores(idx,1)**2)
          self % fluxScores(idx,1) = self % fluxScores(idx,1) + self % uncollidedScores(idx,1)
          if (self % fluxScores(idx,1) > ZERO) then
            self % fluxScores(idx,2) = self % fluxScores(idx,2) / self % fluxScores(idx,1)
          else
            self % fluxScores(idx,2) = ZERO
          end if
        end do
        !$omp end parallel do
      end if
  
    end subroutine cycles
  
    !!
    !! Initialises rays: samples initial position and direction,
    !! and performs the build operation
    !!
    subroutine initialiseRay(self, r)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), intent(inout)                            :: r
      real(defReal)                                       :: mu, phi
      real(defReal), dimension(3)                         :: u, rand3, x
      integer(shortInt)                                   :: i, matIdx, id, cIdx
      character(100), parameter :: Here = 'initialiseRay (adjointBackTRRMPhysicsPackage_class.f90)'
  
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
        call self % geom % whatIsAt(matIdx, id, x, u)
        
        cIdx = self % IDToCell(id)
        if (matIdx /= OUTSIDE_MAT .and. cIdx > 0) exit rejection
  
        i = i + 1
        if (i > 5000) then
          call fatalError(Here, 'Infinite loop when searching for ray start in the geometry.')
        end if
      end do rejection
  
      ! Place in the geometry & process the ray
      call r % build(x, u, 1, ONE)
      call self % geom % placeCoord(r % coords)
  
      if (.NOT. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = x
        !$omp end critical
      end if
  
    end subroutine initialiseRay
    
    !!
    !! Moves ray through geometry, scoring volume
    !! Also used for constructing the cell map
    !!
    subroutine volumeSweep(self, r, maxLength, doVolume)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(ray), intent(inout)                            :: r
      real(defReal), intent(in)                           :: maxLength
      logical(defBool), intent(in)                        :: doVolume
      integer(shortInt)                                   :: event, cIdx
      integer(longInt)                                    :: ints
      real(defReal)                                       :: totalLength, length
      type(distCache)                                     :: cache
      real(defReal), dimension(3)                         :: r0, mu0
      logical(defBool)                                    :: hitVacuum
      
      totalLength = ZERO
      ints = 0_longInt
      
      do while (totalLength < maxLength)
  
        ! Get cell the ray is moving through
        cIdx = self % IDToCell(r % coords % uniqueID)
        !cIdx = r % coords % uniqueID
          
        r0 = r % rGlobal()
        mu0 = r % dirGlobal()
  
        ! Set maximum flight distance 
        length = maxLength - totalLength 
  
        ! Move ray
        ! Use distance caching or standard ray tracing
        ! Distance caching seems a little bit more unstable
        ! due to FP error accumulation, but is faster.
        ! This can be fixed by resetting the cache after X number
        ! of distance calculations.
        if (self % cache) then
          if (mod(ints,20_longInt) == 0)  cache % lvl = 0
          call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
        else
          call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)
        end if
        totalLength = totalLength + length
  
        if (doVolume) then
          !$omp atomic
          self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
        end if
          
        ! Set new cell's position. Use half distance across cell
        ! to try and avoid FP error
        if (.not. self % cellFound(cIdx)) then
          !$omp critical 
          self % cellFound(cIdx) = .true.
          self % cellPos(cIdx,:) = r0 + length * HALF * mu0
          !$omp end critical
        end if
        
        if (self % cellHit(cIdx) /= 1) then
          !$omp critical
          self % cellHit(cIdx) = 1
          !$omp end critical
        end if
  
      end do
  
    end subroutine volumeSweep
  
    !!
    !! Sweep to perform first collided flux calculations.
    !!
    !! Can take one of two approaches
    !!
    !! For point source: samples ray isotropically at the source point,
    !! with flux given by source intensity and traces until some termination length
    !!
    !! For volumetric source: samples ray inside source region, transports them
    !! and deposits flux as they move -- not quite like standard MoC
    !!
    !! To get good/meaningful results, volumes must be precomputed already, 
    !! otherwise the volume estimation across the geometry would be biased by the
    !! uneven ray distribution to begin
    !!
    !!
    subroutine uncollidedSweep(self, r, ints)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
      type(ray), intent(inout)                              :: r
      integer(longInt), intent(out)                         :: ints
      integer(shortInt)                                     :: matIdx, g, event, matIdx0, i, cIdx, baseIdx
      real(defReal)                                         :: totalLength, length, mu, phi
      real(defFlt)                                          :: lenFlt
      logical(defBool)                                      :: hitVacuum
      type(distCache)                                       :: cache
      real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec
      real(defFlt), pointer, dimension(:)                   :: scalarVec, totVec
      real(defReal), dimension(3)                           :: r0, mu0, u, x0, rand3
      character(100), parameter :: Here = 'uncollidedSweep (adjointBackTRRMPhysicsPackage_class.f90)'
      
      ! If point source, position and direction sample is straightforward
      ! Flux is determined by source
      if (self % uncollidedType == POINT) then
  
        mu = TWO * r % pRNG % get() - ONE
        phi = TWO_PI * r % pRNG % get()
        u = rotateVector([ONE, ZERO, ZERO], mu, phi)
        call r % build(self % sourcePoint, u, 1, ONE)
  
        !$omp simd 
        do g = 1, self % nG
          fluxVec(g) = real(self % sourceStrength(g), defFlt)
        end do
      
      ! If volumetric source, rejection sample in region source occupies
      elseif (self % uncollidedType == VOLUME) then
  
        ! Rejection sample position in source volume
        u = [ONE, ZERO, ZERO]
        i = 0
        rejectSource : do
          rand3(1) = r % pRNG % get()
          rand3(2) = r % pRNG % get()
          rand3(3) = r % pRNG % get()
          x0 = self % sourceBottom + &
                  (self % sourceTop - self % sourceBottom) * rand3
  
          ! Exit if point is inside the source
          call self % geom % whatIsAt(matIdx, cIdx, x0, u)
          if (any(matIdx == self % sourceIdx)) exit rejectSource
  
          i = i + 1
          if (i > 5000) then
            call fatalError(Here, 'Could not find source position.')
          end if
  
        end do rejectSource
  
        mu = TWO * r % pRNG % get() - ONE
        phi = TWO_PI * r % pRNG % get()
        u = rotateVector([ONE, ZERO, ZERO], mu, phi)
        call r % build(x0, u, 1, ONE)
  
        ! Set flux after for this type - need to identify cell
  
      else
        call fatalError(Here,'No compatible source type specified')
      end if
      call self % geom % placeCoord(r % coords)
      
      cIdx = self % IDToCell(r % coords % uniqueID)
      
      matIdx0 = 0
      ints = 0
      totalLength = ZERO
      
      ! Set initial condition for ray
      if (self % uncollidedType == VOLUME) then
        matIdx  = r % coords % matIdx
        matIdx0 = matIdx
        baseIdx = (cIdx - 1) * self % nG
        totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):((matIdx - 1) * self % nG + self % nG))
        !$omp simd 
        do g = 1, self % nG
          fluxVec(g) = self % fixedSource(baseIdx + g) 
        end do
      end if
  
      do while (totalLength < self % uncollidedLength)
  
        ! Get material and cell the ray is moving through
        matIdx  = r % coords % matIdx
        cIdx    = self % IDToCell(r % coords % uniqueID)
        if (matIdx0 /= matIdx) then
          matIdx0 = matIdx
          
          ! Cache total cross section
          totVec => self % sigmaT(((matIdx - 1) * self % nG + 1):((matIdx - 1) * self % nG + self % nG))
        end if
  
        ! Remember co-ordinates to set new cell's position
        if (.not. self % cellFound(cIdx)) then
          r0 = r % rGlobal()
          mu0 = r % dirGlobal()
        end if
            
        ! Set maximum flight distance
        length = self % uncollidedLength - totalLength 
  
        ! Move ray
        ! Use distance caching or standard ray tracing
        ! Distance caching seems a little bit more unstable
        ! due to FP error accumulation, but is faster.
        ! This can be fixed by resetting the cache after X number
        ! of distance calculations.
        if (self % cache) then
          if (mod(ints,20_longInt) == 0)  cache % lvl = 0
          call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
        else
          call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)
        end if
        totalLength = totalLength + length
        lenFlt = real(length,defFlt)
        
        ! Set new cell's position. Use half distance across cell
        ! to try and avoid FP error
        if (.not. self % cellFound(cIdx)) then
          !$omp critical 
          self % cellFound(cIdx) = .true.
          self % cellPos(cIdx,:) = r0 + length * HALF * mu0
          !$omp end critical
        end if
  
        ints = ints + 1
  
        baseIdx = (cIdx - 1) * self % nG
        ! No need for sourceVec when only depositing source
        scalarVec => self % scalarFlux((baseIdx + 1):(baseIdx + self % nG))
  
        !$omp simd
        do g = 1, self % nG
          attenuate(g) = exponential(totVec(g) * lenFlt)
          delta(g) = fluxVec(g) * attenuate(g)
          fluxVec(g) = fluxVec(g) - delta(g)
        end do
  
        ! Accumulate to scalar flux
        ! Assume no volume scoring due to non-uniform sampling!
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
          scalarVec(g) = scalarVec(g) + delta(g)
        end do
        call OMP_unset_lock(self % locks(cIdx))
  
        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
        
        ! Check for a vacuum hit
        ! Exit if vacuum 
        ! (can't acquire more source and UCF rays needn't travel same length)
        if (hitVacuum) then
          return
        end if
  
      end do
  
    end subroutine uncollidedSweep

    subroutine transportSweep(self, r, ints)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
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
    !! Normalise flux from uncollided calculation
    !!
    subroutine normaliseFluxUncollided(self, norm)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      real(defReal), intent(in)                           :: norm
      real(defFlt)                                        :: normFlt
      real(defFlt), save                                  :: total
      integer(shortInt), save                             :: g, matIdx, idx
      integer(shortInt)                                   :: cIdx
      !$omp threadprivate(total, idx, g, matIdx)
  
      normFlt = real(norm, defFlt)
      
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
      
        ! Presume that volumes are known otherwise this may go badly!
        if (self % volume(cIdx) > volume_tolerance) then
          matIdx =  self % geom % geom % graph % getMatFromUID(self % CellToID(cIdx))
          
          ! Guard against void cells
          if (matIdx > 2000000) then
            do g = 1, self % nG
              idx   = self % nG * (cIdx - 1) + g
              self % scalarFlux(idx) = 0.0_defFlt
            end do
            cycle
          end if 
        
          do g = 1, self % nG
            total = self % sigmaT((matIdx - 1) * self % nG + g)
            idx   = self % nG * (cIdx - 1) + g
            self % scalarFlux(idx) = self % scalarFlux(idx) * normFlt / (total * real(self % volume(cIdx),defFlt))
          end do
  
        end if
  
      end do
      !$omp end parallel do
  
    end subroutine normaliseFluxUncollided
    
    !!
    !! Normalise flux and volume by total track length and increments
    !! the flux by the neutron source
    !!
    subroutine normaliseFluxAndVolume(self, lengthPerIt, it)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      real(defReal), intent(in)                           :: lengthPerIt
      integer(shortInt), intent(in)                       :: it
      real(defReal)                                       :: normVol
      real(defFlt)                                        :: norm
      real(defFlt), save                                  :: total, sigGG, D
      real(defReal), save                                 :: vol, corr
      integer(shortInt), save                             :: g, matIdx, idx
      integer(shortInt)                                   :: cIdx
      !$omp threadprivate(total, idx, g, matIdx, vol, corr, sigGG, D)
  
      norm = real(ONE / lengthPerIt, defFlt)
      normVol = ONE / ( lengthPerIt * it)
  
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        matIdx =  self % geom % geom % graph % getMatFromUID(self % CellToID(cIdx)) 
          
        ! Guard against void cells
        if (matIdx > 2000000) then
          do g = 1, self % nG
            idx   = self % nG * (cIdx - 1) + g
            self % scalarFlux(idx) = 0.0_defFlt
          end do
          cycle
        end if 
        
        ! Update volume due to additional rays unless volume was precomputed
        !if (self % nVolRays <= 0) then 
        ! Forget the above - use precomputed volumes only for first collided
        
        if (self % itVol) then
          ! Iteration wise approach
          self % volume(cIdx) = self % volumeTracks(cIdx) * norm
          self % volumeTracks(cIdx) = ZERO
        else if (self % volCorr) then
          ! Correct the cell volume
          corr = self % volumeTracks(cIdx) * norm 
          self % volume(cIdx) = self % volume(cIdx) * self % lengthPerIt * (it - 1) + self % volumeTracks(cIdx)
          self % volume(cIdx) = self % volume(cIdx) * normVol
          corr = corr / self % volume(cIdx) 
          self % volumeTracks(cIdx) = ZERO
          if (corr /= corr) corr = ONE
        else
          ! Standard volume approach
          self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
        end if
  
        vol = self % volume(cIdx)
  
        do g = 1, self % nG
  
          idx   = self % nG * (cIdx - 1) + g
          total = self % sigmaT((matIdx - 1) * self % nG + g)
          
          if (vol > volume_tolerance) then
            self % scalarFlux(idx) = self % scalarFlux(idx) * norm / (total * real(vol,defFlt))
          else
            corr = ONE
          end if

          self % scalarFlux(idx) =  self % scalarflux(idx) + self % source(idx) 
          
          ! ! TODO: I think this may not work when volume correction is applied... 
          ! if (matIdx < UNDEF_MAT) then
          !   sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)
  
          !   ! Assumes non-zero total XS
          !   if ((sigGG < 0) .and. (total > 0)) then
          !     D = -real(self % rho, defFlt) * sigGG / total
          !   else
          !     D = 0.0_defFlt
          !   end if
          ! else
          !   D = 0.0_defFlt
          ! end if
  
          ! self % scalarFlux(idx) =  real((self % scalarflux(idx) + self % source(idx) &
          !         / total + D * self % prevFlux(idx) ) / (1 + D), defFlt)
          
          ! Apply volume correction only to negative flux cells
          if (self % volCorr .and. self % passive) then
            if (self % scalarFlux(idx) < 0) self % scalarFlux(idx) = real(self % scalarFlux(idx) + &
                    (corr - 1.0_defFlt) * self % source(idx) / total, defFlt)
          ! Apply volume correction to all cells
          elseif (self % volCorr) then
            self % scalarFlux(idx) = real(self % scalarFlux(idx) + (corr - 1.0_defFlt) * self % source(idx) / total, defFlt)
          end if
  
          ! This will probably affect things like neutron conservation...
          if ((self % scalarFlux(idx) < 0) .and. self % zeroNeg) self % scalarFlux(idx) = 0.0_defFlt
          ! DELETE THIS MAYBE?
          !if ((self % scalarFlux(idx) < 0) .and. self % zeroNeg) self % scalarFlux(idx) = self % source(idx) / total
  
          ! NaN check - kill calculation
          if (self % scalarFlux(idx) /= self % scalarFlux(idx)) then
            print *, self % scalarFlux(idx), total, self % source(idx), matIdx
                  call fatalError('normaliseFluxAndVolume','NaNs appeared in group '//numToChar(g))
          end if
  
        end do
  
      end do
      !$omp end parallel do
  
    end subroutine normaliseFluxAndVolume

    subroutine adjointNormaliseFluxAndVolume(self, lengthPerIt, it)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      integer(shortInt), intent(in)                 :: it
      real(defReal), intent(in)                     :: lengthPerIt
      real(defFlt)                                  :: norm
      real(defReal)                                 :: normVol
      real(defFlt), save                            :: total, vol ! sigGG, D
      integer(shortInt), save                       :: g, matIdx, idx
      integer(shortInt)                             :: cIdx
      !$omp threadprivate(total, vol, idx, g, matIdx)
  
      norm = real(ONE / lengthPerIt)
      normVol = ONE / (lengthPerIt * it)
  
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        matIdx =  self % geom % geom % graph % getMatFromUID(self % CellToID(cIdx)) 

        ! Guard against void cells
        if (matIdx > 2000000) then
          do g = 1, self % nG
            idx   = self % nG * (cIdx - 1) + g
            self % adjscalarFlux(idx) = 0.0_defFlt
          end do
          cycle
        end if 
        
        ! Update volume due to additional rays
        self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
        vol = real(self % volume(cIdx),defFlt)
  
        do g = 1, self % nG
  
          total = self % sigmaT((matIdx - 1) * self % nG + g)
          idx   = self % nG * (cIdx - 1) + g
  
          if (vol > volume_tolerance) then
            self % adjScalarFlux(idx) = self % adjScalarFlux(idx) * norm / ( total * vol)
          end if
  
          self % adjScalarFlux(idx) = self % adjScalarFlux(idx) + self % adjSource(idx)
          
          ! ! Apply stabilisation for negative XSs
          ! if (matIdx < VOID_MAT) then
          !   sigGG = self % sigmaS(self % nG * self % nG * (matIdx - 1) + self % nG * (g - 1) + g)
  
          !   ! Presumes non-zero total XS
          !   if ((sigGG < 0) .and. (total > 0)) then
          !     D = -real(self % rho, defFlt) * sigGG / total
          !   else
          !     D = 0.0_defFlt
          !   end if
          ! else
          !   D = 0.0_defFlt
          ! end if
  
          ! self % scalarFlux(idx) =  (self % scalarflux(idx) + self % source(idx) + D * self % prevFlux(idx) ) / (1 + D)
  
        end do
  
      end do
      !$omp end parallel do

      !print *, self % scalarFlux(1:5), self % adjScalarFlux(1:5)
  
    end subroutine adjointNormaliseFluxAndVolume
    
  
    !!
    !! Kernel to update sources given a cell index
    !!
    subroutine sourceUpdateKernel(self, cIdx)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
      integer(shortInt), intent(in)                         :: cIdx
      real(defFlt)                                          :: scatter, fission
      real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
      integer(shortInt)                                     :: matIdx, g, gIn, id, baseIdx, idx
      real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec
  
      ! Identify material
      id      =  self % CellToID(cIdx)
      matIdx  =  self % geom % geom % graph % getMatFromUID(id) 
  
      ! Hack to guard against non-material cells
      if (matIdx >= UNDEF_MAT) return
      
      ! Obtain XSs
      matIdx = (matIdx - 1) * self % nG
      total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
      scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
      nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
      chi => self % chi((matIdx + 1):(matIdx + self % nG))
  
      baseIdx = self % nG * (cIdx - 1)
      fluxVec => self % prevFlux((baseIdx+1):(baseIdx + self % nG))
  
      ! Calculate fission source
      fission = 0.0_defFlt
      !$omp simd reduction(+:fission)
      do gIn = 1, self % nG
        fission = fission + fluxVec(gIn) * nuFission(gIn)
      end do
  
      do g = 1, self % nG
  
        scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * (g - 1) + self % nG))
  
        ! Calculate scattering source
        scatter = 0.0_defFlt
  
        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
          scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        end do
  
        ! Output index
        idx = baseIdx + g
  
        self % source(idx) = chi(g) * fission + scatter + self % fixedSource(idx)
        self % source(idx) = self % source(idx) / total(g)
      end do
  
    end subroutine sourceUpdateKernel
    
    
    subroutine adjointSourceUpdateKernel(self, cIdx)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
      integer(shortInt), intent(in)                         :: cIdx
      class(baseMgNeutronMaterial), pointer                 :: mat
      class(materialHandle), pointer                        :: matPtr
      real(defFlt)                                          :: scatter, fission
      real(defFlt)                                          :: Sigma
      real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
      integer(shortInt)                                     :: matIdx, g, gIn, id, baseIdx, idx, material
      real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec
      logical(defBool)                                      :: detMat
  
      ! Identify material
      id       =  self % CellToID(cIdx)
      matIdx   =  self % geom % geom % graph % getMatFromUID(id) 
      material =  self % geom % geom % graph % getMatFromUID(id) 
      
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

      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)

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
      fission = fission
  
      do g = 1, self % nG
  
        scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * self % nG))
  
        ! Calculate scattering source
        scatter = 0.0_defFlt
  
        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter) aligned(fluxVec)
        do gIn = 1, self % nG
          scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        end do

        if (allocated(self % detectorInt)) then
          detMat = any(self % detectorInt == material)
        else
          detMat = .true.
        end if

        if (detMat) then
          select case(self % mapResponse)
          case(1)
            Sigma = real(mat % getFissionXS(g, self % rand),defFlt)
          case(2)
            Sigma = real(mat % getTotalXS(g, self % rand),defFlt)
          case(3)
            Sigma = real(mat % getCaptureXS(g, self % rand),defFlt)
          end select
        else
          Sigma = 0.0_defFlt
        end if
        ! Sigma = real(mat % getFissionXS(g, self % rand),defFlt)

        ! Output index
        idx = baseIdx + g
        
        self % adjSource(idx) = chi(g) * fission + scatter + Sigma
        self % adjSource(idx) = self % adjSource(idx) / total(g)
  
      end do
  
    end subroutine adjointSourceUpdateKernel
  
    !!
    !! Kernel to calculate the first collided source
    !! Overwrites any existing fixed source
    !!
    subroutine firstCollidedSourceKernel(self, cIdx)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
      integer(shortInt), intent(in)                         :: cIdx
      real(defFlt)                                          :: scatter, fission
      real(defFlt), dimension(:), pointer                   :: nuFission, chi, scatterXS 
      integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx
      real(defFlt), pointer, dimension(:)                   :: scatterVec
      real(defReal), pointer, dimension(:)                  :: fluxVec
  
      ! Identify material
      matIdx  =  self % geom % geom % graph % getMatFromUID(self % CellToID(cIdx)) 
      
      ! Hack to guard against non-material cells
      if (matIdx >= UNDEF_MAT) return
  
      ! Obtain XSs
      matIdx = (matIdx - 1) * self % nG
      scatterXS => self % sigmaS((matIdx * self % nG + 1):(matIdx * self % nG + self % nG*self % nG))
      nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
      chi => self % chi((matIdx + 1):(matIdx + self % nG))
  
      baseIdx = self % nG * (cIdx - 1)
      fluxVec => self % uncollidedScores((baseIdx+1):(baseIdx + self % nG),1)
  
      ! Calculate fission source
      fission = 0.0_defFlt
      !$omp simd reduction(+:fission)
      do gIn = 1, self % nG
        fission = fission + real(fluxVec(gIn),defFlt) * nuFission(gIn)
      end do
  
      do g = 1, self % nG
  
        scatterVec => scatterXS((self % nG * (g - 1) + 1):(self % nG * (g - 1) + self % nG))
  
        ! Calculate scattering source
        scatter = 0.0_defFlt
  
        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
          scatter = scatter + real(fluxVec(gIn),defFlt) * scatterVec(gIn)
        end do
  
        ! Output index
        idx = baseIdx + g
  
        ! Don't scale by 1/SigmaT - that occurs in the sourceUpdateKernel
        self % fixedSource(idx) = chi(g) * fission + scatter
  
      end do
  
    end subroutine firstCollidedSourceKernel
  
    !!
    !! Sets prevFlux to scalarFlux and zero's scalarFlux
    !!
    subroutine resetFluxes(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      integer(shortInt)                                   :: idx
  
      !$omp parallel do schedule(static)
      do idx = 1, size(self % scalarFlux)
        self % prevFlux(idx) = self % scalarFlux(idx)
        self % scalarFlux(idx) = 0.0_defFlt

        self % adjPrevFlux(idx) = self % adjScalarFlux(idx)
        self % adjScalarFlux(idx) = 0.0_defFlt
      end do
      !$omp end parallel do
  
    end subroutine resetFluxes

    subroutine sensitivityCalculation(self, cIdx, it)
      class(adjointBackTRRMPhysicsPackage), target, intent(inout) :: self
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
      mat = self % matPert
  
      IPVec => self % angularIP((baseIdx * self % nG + 1):(baseIdx + self % nG) * self % nG)
      ! total => self % sigmaT((matIdx + 1):(matIdx + self % nG))
      nuFission => self % nuSigmaF((matIdx + 1):(matIdx + self % nG))
      chi => self % chi((matIdx + 1):(matIdx + self % nG))
  
      ! do g = 1, self % nG
      !   idx = baseIdx + g
      !   fission = 0.0_defFlt
      !   !$omp simd
      !   do gIn = 1, self % nG 
      !     fission = fission + IPVec((g - 1) * self % nG + gIn) * nuFission(gIn)
      !   end do
      !   fission  = fission * chi(g)
      ! end do

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
            self % numSum = self % numSum - delta + fission_pert * chi(g) !* ONE_KEFF
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
  
    !!
    !! Accumulate flux scores for stats
    !!
    subroutine accumulateFluxScores(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      real(defReal), save                                 :: flux
      integer(shortInt)                                   :: idx
      !$omp threadprivate(flux)
  
      !$omp parallel do schedule(static)
      do idx = 1, size(self % adjscalarFlux)
        flux = real(self % adjScalarFlux(idx),defReal)
        self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
        self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
      end do
      !$omp end parallel do

      self % sensitivityScore(1) = self % sensitivityScore(1) + self % sensitivity
      self % sensitivityScore(2) = self % sensitivityScore(2) + self % sensitivity * self % sensitivity
  
    end subroutine accumulateFluxScores
    
    !!
    !! Finalise flux scores for stats
    !!
    subroutine finaliseFluxScores(self,it)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      integer(shortInt), intent(in)                       :: it
      integer(shortInt)                                   :: idx
      real(defReal)                                       :: N1, Nm1
  
      if (it > 1) then
        Nm1 = ONE/(it - 1)
      else
        Nm1 = ONE
      end if
      N1 = ONE/it
  
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
          if (abs(self % fluxScores(idx,1)) > ZERO) then
            self % fluxScores(idx,2) = self % fluxScores(idx,2) / abs(self % fluxScores(idx,1))
          end if
        end if
      end do
      !$omp end parallel do

      self % sensitivityScore(1) = self % sensitivityScore(1) * N1
      self % sensitivityScore(2) = self % sensitivityScore(2) * N1
      self % sensitivityScore(2) = (Nm1*(self % sensitivityScore(2) - &
              self % sensitivityScore(1) * self % sensitivityScore(1))) 
      self % sensitivityScore(2) = sqrt(abs(self % sensitivityScore(2)))
   
    end subroutine finaliseFluxScores
    
    !!
    !! Output calculation results to a file
    !!
    !! Args:
    !!   None
    !!
    subroutine printResults(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      type(outputFile)                                    :: out
      character(nameLen)                                  :: name
      integer(shortInt)                                   :: g1, cIdx
      integer(shortInt), save                             :: matIdx, g, idx, i
      real(defReal), save                                 :: vol, Sigma
      type(particleState), save                           :: s
      type(ray), save                                     :: pointRay
      real(defReal)                                       :: res, std, totalVol
      integer(shortInt),dimension(:),allocatable          :: resArrayShape
      class(baseMgNeutronMaterial), pointer, save         :: mat
      class(materialHandle), pointer, save                :: matPtr
      real(defReal), dimension(:), allocatable            :: groupFlux, flxOut, flxOutSTD, Response, ResponseSTD
      !$omp threadprivate(idx, matIdx, i, vol, s, g, mat, matptr, sigma, pointRay)
  
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
      
      if (self % nVolRays > 0) then
        name = 'Total_Volume_Time'
        call out % printValue(self % time_volume,name)
      end if
  
      if (self % uncollidedType > NO_UC) then
        name = 'Total_UC_Time'
        call out % printValue(self % time_UC,name)
        name = 'Total_UC_Transport_Time'
        call out % printValue(self % time_transportUC,name)
      end if
  
      name = 'Clock_Time'
      call out % printValue(timerTime(self % timerMain),name)
      
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
          call out % addResult(self % volume(cIdx) * self % normVolume, ZERO)
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
            call out % addResult(self % fluxScores(idx,1),&
                    self % fluxScores(idx,2))
          end do
          call out % endArray()
          call out % endBlock()
        end do
      end if

      ! Send fission rates to map output
      if (self % mapResponse > 0) then
        resArrayShape = self % resultsMap % binArrayShape()
        allocate(Response(self % resultsMap % bins(0)))
        allocate(ResponseSTD(self % resultsMap % bins(0)))
        Response   = ZERO
        ResponseSTD = ZERO

        ! Find whether cells are in map and sum their contributions
        !$omp parallel do reduction(+: Response, ResponseSTD)
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

              Sigma = real(mat % getFissionXS(g, self % rand),defFlt)

              
              idx = (cIdx - 1) * self % nG + g
              Response(i) = Response(i) + vol * self % fluxScores(idx,1) * self % fixedSource(idx)!* Sigma
              ! Is this correct? Also neglects uncertainty in volume - assumed small.
              ResponseSTD(i) = ResponseSTD(i) + &
                      vol * vol * self % fluxScores(idx,2) *self % fluxScores(idx,2) * &
                            self % fixedSource(idx) * self % fixedSource(idx)
            end do
          end if

        end do
        !$omp end parallel do

        !$omp parallel do
        do cidx = 1,size(ResponseSTD)
          if (Response(i) > 0) then
            ResponseSTD(i) = sqrt(ResponseSTD(i)) / Response(i)
          else
            ResponseSTD(cIdx) = ZERO
          end if
        end do
        !$omp end parallel do

        name = 'responseRate'
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        ! Add all map elements to results
        do idx = 1, self % resultsMap % bins(0)
          call out % addResult(real(Response(idx),defReal), real(ResponseSTD(idx),defReal))
        end do
        call out % endArray()
        ! Output tally map
        call self % resultsMap % print(out)
        call out % endBlock()
        
        deallocate(Response)
        deallocate(ResponseSTD)
      end if
      
      ! Send fluxes to map output
      if (self % mapFlux) then
        resArrayShape = [self % nG, self % fluxMap % binArrayShape()]
        allocate(flxOut(self % fluxMap % bins(0)*self % nG))
        allocate(flxOutSTD(self % fluxMap % bins(0)*self % nG))
        flxOut    = ZERO
        flxOutSTD = ZERO
  
        ! Find whether cells are in map and sum their contributions
        !$omp parallel do reduction(+: flxOut, flxOutSTD)
        do cIdx = 1, self % nCells
          
          vol    =  self % volume(cIdx)
          if (vol < volume_tolerance) cycle
  
          ! Fudge a particle state to search tally map
          s % r = self % cellPos(cIdx,:)
          i = self % fluxMap % map(s)
  
          if (i > 0) then
            do g = 1, self % nG
              idx = (cIdx - 1)* self % nG + g
              flxOut((i-1) * self % nG + g ) = flxOut((i-1) * self % nG + g) + vol * self % fluxScores(idx,1)
              ! Is this correct? Also neglects uncertainty in volume - assumed small.
              flxOutSTD((i-1) * self % nG + g) = flxOutSTD((i-1)* self % nG + g) + &
                      self % fluxScores(idx,2)*self % fluxScores(idx,2) * vol * vol * &
                      self % fluxScores(idx,1)*self % fluxScores(idx,1)
            end do
          end if
        end do
        !$omp end parallel do
  
        !$omp parallel do
        do cIdx = 1,size(flxOutSTD)
          if (flxOut(cIdx) > 0) then
            flxOutSTD(cIdx) = sqrt(flxOutSTD(cIdx))/flxOut(cIdx)
          else
            flxOutSTD(cIdx) = ZERO
          end if
        end do
        !$omp end parallel do
  
        name = 'fluxMap'
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        ! Add all map elements to results
        do idx = 1, size(flxOut)
          call out % addResult(real(flxOut(idx),defReal), real(flxOutSTD(idx),defReal))
        end do
        call out % endArray()
        ! Output tally map
        call self % fluxMap % print(out)
        call out % endBlock()
        
        deallocate(flxOut)
        deallocate(flxOutSTD)
      end if
  
      ! Print sample point values if requested
      if (allocated(self % sampleNames)) then
        resArrayShape = [self % nG]
        do i = 1, size(self % sampleNames)
          name = self % sampleNames(i)
          call out % startBlock(name)
          call out % startArray(name, resArrayShape)
          s % r = self % samplePoints(1+3*(i-1):3*i)
          pointRay = s
          call self % geom % placeCoord(pointRay % coords)
          cIdx = self % IDToCell(pointRay % coords % uniqueID)
          do g = 1, self % nG
            idx = (cIdx - 1)* self % nG + g
            res = self % fluxScores(idx,1)
            std = self % fluxScores(idx,2)
            call out % addResult(res, std)
          end do
          call out % endArray()
          call out % endBlock()
        end do
      end if
  
      ! Print material integrated fluxes if requested
      if (allocated(self % intMatIdx)) then
        name = 'integral reaction rate'
        resArrayShape = [1]
        call out % startBlock(name)
        do i = 1, size(self % intMatIdx)
          call out % startArray(self % intMatName(i), resArrayShape)
          res = ZERO
          std = ZERO
          totalVol = ZERO
  
          do cIdx = 1, self % nCells
            matIdx  =  self % geom % geom % graph % getMatFromUID(self % CellToID(cIdx)) 
  
            if (self % intMatIdx(i) == matIdx) then
              vol = self % normVolume * self % volume(cIdx)
              if (vol < volume_tolerance) continue
              totalVol = totalVol + real(vol,defReal)
              do g = 1, self % nG
                idx = (cIdx - 1)* self % nG + g
                res = res + self % fluxScores(idx,1) * vol * self % fixedSource(idx)
                std = std + self % fluxScores(idx,2)**2 * self % fluxScores(idx,1)**2 &
                  * vol * vol * self % fixedSource(idx) * self % fixedSource(idx)
              end do
            end if
          end do
          if (res > ZERO) then
            std = sqrt(std)/res
          else
            std = ZERO
          end if
          call out % addResult(res, std)
          call out % endArray()
        end do
        call out % endBlock()
      end if
  
      call out % writeToFile(self % outputFile)
  
      ! Send all fluxes and stds to VTK
      if (self % plotResults) then
        allocate(groupFlux(self % nCells))
        groupFlux = ZERO
        do g1 = 1, self % nG
          name = 'flux_g'//numToChar(g1)
          !$omp parallel do schedule(static)
          do cIdx = 1, self % nCells
            idx = (cIdx - 1) * self % nG + g1
            groupFlux(cIdx) = self % fluxScores(idx,1)
          end do
          !$omp end parallel do
          call self % viz % addVTKData(groupFlux,name)
        end do
        groupFlux = ZERO
        do g1 = 1, self % nG
          name = 'std_g'//numToChar(g1)
          !$omp parallel do schedule(static)
          do cIdx = 1, self % nCells
            idx = (cIdx - 1) * self % nG + g1
            groupFlux(cIdx) = self % fluxScores(idx,2)
          end do
          !$omp end parallel do
          call self % viz % addVTKData(groupFlux,name)
        end do
        groupFlux = ZERO
        do g1 = 1, self % nG
          name = 'source_g'//numToChar(g1)
          !$omp parallel do schedule(static)
          do cIdx = 1, self % nCells
            idx = (cIdx - 1) * self % nG + g1
            groupFlux(cIdx) = self % source(idx)
          end do
          !$omp end parallel do
          call self % viz % addVTKData(groupFlux,name)
        end do
        groupFlux = ZERO
        name = 'volume'
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          groupFlux(cIdx) = self % volume(cIdx) * self % normVolume
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
        if (self % uncollidedType > NO_UC) then
          groupFlux = ZERO
          do g1 = 1, self % nG
            name = 'uncollided_g'//numToChar(g1)
            !$omp parallel do schedule(static)
            do cIdx = 1, self % nCells
              idx = (cIdx - 1)* self % nG + g1
              groupFlux(cIdx) = self % uncollidedScores(idx,1)
            end do
            !$omp end parallel do
            call self % viz % addVTKData(groupFlux,name)
          end do
        end if
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
      class(adjointBackTRRMPhysicsPackage), intent(in) :: self
  
      print *, repeat("<>", MAX_COL/2)
      print *, "/\/\ RANDOM RAY FIXED SOURCE CALCULATION /\/\"
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
      if (self % nVolRays > 0) print *, "Precomputes volume"
      if (self % uncollidedType > NO_UC) print *, "Performs uncollided flux calculation"
      print *, repeat("<>", MAX_COL/2)
  
    end subroutine printSettings
  
    !!
    !! Return to uninitialised state
    !!
    subroutine kill(self)
      class(adjointBackTRRMPhysicsPackage), intent(inout) :: self
      integer(shortInt) :: i
  
      ! Clean Nuclear Data, Geometry and visualisation
      call gr_kill()
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
      self % nCells    = 0
      self % nMat      = 0
      
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

      if(allocated(self % adjSigmaS)) deallocate(self % adjSigmaS)
      if(allocated(self % adjNuSigmaF)) deallocate(self % adjNuSigmaF)
      if(allocated(self % adjChi)) deallocate(self % adjChi)
  
  
      self % termination = ZERO
      self % dead        = ZERO
      self % lengthPerIt = ZERO
      self % skipLength  = ZERO
      self % pop         = 0
      self % inactive    = 0
      self % active      = 0
      self % rho         = ZERO
      self % cache       = .false.
      self % plotResults = .false.
      self % printFlux   = .false.
      self % printVolume = .false.
      self % printCells  = .false.
      self % nVolRays    = 0
      self % volLength   = ZERO
      self % normVolume  = ONE
      self % itVol       = .false.
      self % zeroNeg     = .false.
      self % volCorr     = .false.
      self % passive     = .false.
      self % mapFlux     = .false.
      self % mapResponse = 0
  
      self % uncollidedType   = NO_UC
      self % uncollidedPop    = 0
      self % uncollidedCycles = 0
      self % uncollidedLength = ZERO
      self % UCNorm           = ZERO
      self % sourcePoint      = ZERO
      self % sourceTop        = ZERO
      self % sourceBottom     = ZERO

      self % XStype   = 0
      self % matPert  = 0
      self % XSchange = 0.0_defFlt
      self % numSum   = ZERO
      self % denSum   = ZERO
  
      self % sensitivity      = ZERO
      self % sensitivityScore = ZERO


      if(allocated(self % energyId)) deallocate(self % energyId)

      if(allocated(self % adjScalarFlux)) deallocate(self % adjScalarFlux)
      if(allocated(self % adjPrevFlux)) deallocate(self % adjPrevFlux)
      if(allocated(self % adjSource)) deallocate(self % adjSource)

      if(allocated(self % angularIP)) deallocate(self % angularIP)
  
      if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
      if(allocated(self % prevFlux)) deallocate(self % prevFlux)
      if(allocated(self % fluxScores)) deallocate(self % fluxScores)
      if(allocated(self % source)) deallocate(self % source)
      if(allocated(self % fixedSource)) deallocate(self % fixedSource)
      if(allocated(self % volume)) deallocate(self % volume)
      if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
      if(allocated(self % cellHit)) deallocate(self % cellHit)
      if(allocated(self % cellFound)) deallocate(self % cellFound)
      if(allocated(self % cellPos)) deallocate(self % cellPos)
      if(allocated(self % cellToID)) deallocate(self % cellToID)
      if(allocated(self % IDToCell)) deallocate(self % IDToCell)
      if(allocated(self % sampleNames)) deallocate(self % sampleNames)
      if(allocated(self % samplePoints)) deallocate(self % samplePoints)
      if(allocated(self % sourceIdx)) deallocate(self % sourceIdx)
      if(allocated(self % intMatIdx)) deallocate(self % intMatIdx)
      if(allocated(self % intMatName)) deallocate(self % intMatName)
      if(allocated(self % detectorInt)) deallocate(self % detectorInt)
      if(allocated(self % detectorName)) deallocate(self % detectorName)
      if(allocated(self % fluxMap)) then
        call self % fluxMap % kill()
        deallocate(self % fluxMap)
      end if
      if(allocated(self % sourceStrength)) deallocate(self % sourceStrength)
      if(allocated(self % uncollidedScores)) deallocate(self % uncollidedScores)
  
    end subroutine kill
  
  end module adjointBackTRRMPhysicsPackage_class
  