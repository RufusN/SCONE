module aceNoMT_class

  use numPrecision
  use endfConstants
  use genericProcedures ,      only : openToRead, binarySearch, interpolate, searchError, fatalError
  use aceCard_class,           only : aceCard
  use RNG_class,               only : RNG
  use emissionFromACE_func,    only : emissionFromACE
  use emissionENDF_class,      only : emissionENDF, emissionENDF_ptr
  use xsEnergyPointNoMT_class, only : xsEnergyPointNoMT
  use materialDataNoMT_class,  only : materialDataNoMT


  implicit none
  private

  integer(shortInt), parameter :: captureIdx  = 1, &
                                  escatterIdx = 2, &
                                  fissionIdx  = 3
  !!
  !! Class that stores nuclide reaction cross-sections read from ACE data card.
  !!
  type, public :: aceNoMT
    real(defReal),public      :: atomWeight    !! Atomic weight ratio A/m_neutron
    real(defReal),public      :: temp          !! Temperature of nucleide [MeV]
    character(zzIdLen),public :: zzId          !! Nuclide ZZid ie. ZZZAAA.nnC
    character(10)             :: date          !! Date the data were processed
    character(70)             :: comment       !! Quick! Guess what is it!
    character(10)             :: MATid         !! MAT indentifier (see p. F-9 Appendic F MCNP 4 manual)
!    integer(shortInt),public  :: nReact        !! number of reaction channels

    logical(defBool)   :: isFissile = .false.

!    integer(shortInt),dimension(:), pointer         :: xsMT  => null()  !! MT numbers of cross-sections in xs
    type(emissionENDF_ptr),dimension(3)             :: emissionData
    real(defReal),dimension(:),allocatable          :: energyGrid         !! Energy grid for xs data
    type(xsEnergyPointNoMT),dimension(:),pointer    :: xsData => null()


  contains
    procedure :: init
    procedure :: energyIdxFor
    procedure :: sampleMuEout
    procedure :: releaseAt
    procedure :: isInCMframe

    procedure,private :: readXS

  end type aceNoMT

contains

  !!
  !! Searches for the given energy and returns index of the bottom of the interval
  !!
  function energyIdxFor(self,E) result(idx)
    class(aceNoMT), intent(in)  :: self
    real(defReal), intent(in)   :: E
    integer(shortInt)           :: idx
    character(100), parameter   :: Here = 'searchFor (aceNoMT_class.f90)'

    idx = binarySearch(self % energyGrid,E)
    call searchError(idx,Here)

  end function energyIdxFor


  !!
  !! Sample deflection angle and emission energy for a given MT number
  !!
  subroutine sampleMuEout(self,mu,E_out,E_in,rand,MT)
    class(aceNoMT), intent(in)    :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(out)    :: E_out
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand
    integer(shortInt), intent(in) :: MT
    character(100), parameter     :: Here = 'sampleMuEout (aceNoMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        call self % emissionData(escatterIdx) % sampleAngleEnergy(mu,E_out,E_in,rand)

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        call self % emissionData(captureIdx)  % sampleAngleEnergy(mu,E_out,E_in,rand)

      case (anyFission, N_fission)
        call self % emissionData(fissionIdx)  % sampleAngleEnergy(mu,E_out,E_in,rand)

      case default
        call fatalError(Here,'Unknown MT number')

    end select

  end subroutine sampleMuEout

  !!
  !! Function that returns .true. if emission data for a given MT is stored in Center-Of-Mass frame
  !!
  function isInCMframe(self,MT) result (isIt)
    class(aceNoMT), intent(in)    :: self
    integer(shortInt), intent(in) :: MT
    logical(defBool)              :: isIt
    character(100), parameter     :: Here = 'isInCMframe (aceNoMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        isIt = self % emissionData(escatterIdx) % isInCMframe()

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        isIt = self % emissionData(captureIdx)  % isInCMframe()

      case (anyFission, N_fission)
        isIt = self % emissionData(fissionIdx)  % isInCMframe()

      case default
        call fatalError(Here,'Unknown MT number')

    end select


  end function isInCMframe


  !!
  !! Sample give value of average neutron emission for given energy and MT number
  !!
  function releaseAt(self,E_in,MT) result (nu)
    class(aceNoMT), intent(in)    :: self
    real(defReal), intent(in)     :: E_in
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: nu
    character(100), parameter     :: Here = 'releaseAt (aceNoMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        nu = self % emissionData(escatterIdx) % releaseAt(E_in)

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        nu = self % emissionData(captureIdx)  % releaseAt(E_in)

      case (anyFission, N_fission)
        nu = self % emissionData(fissionIdx)  % releaseAt(E_in)

      case default
        call fatalError(Here,'Unknown MT number')

    end select

  end function releaseAt

  !!
  !! Load reaction cross-section data from ACE file.
  !!
  subroutine init(self,filePath,line)
    class(aceNoMT), intent(inout)           :: self
    character(*), intent(in)                :: filePath  !! Path to file with ACE data card
    integer(shortInt), intent(in)           :: line      !! Line at which the ACE data card begins in the file
    type(aceCard)                           :: ACE

    ! Load ACE card
    call ACE % readFromFile(filePath,line)

    ! Copy nuclide header
    self % ZZid       = ACE % ZAID
    self % atomWeight = ACE % AW
    self % temp       = ACE % TZ
    self % date       = ACE % HD
    self % comment    = ACE % HK
    self % MATid      = ACE % HM

    call self % readXS(ACE)

  end subroutine init

  !!
  !! Read cross-sections into type variables from NXS,JXS and XSS tabels
  !! This implementation ignores reaction other then total,elastic scattering, total capture and
  !! fission.
  !!
  subroutine readXS(self,ACE)
    class(aceNoMT), intent(inout)               :: self
    type(aceCard), intent(inout)                :: ACE
    real(defReal), dimension(:), allocatable    :: xsEScatter
    real(defReal), dimension(:), allocatable    :: xsCapture
    real(defReal), dimension(:), allocatable    :: xsFission
    integer(shortInt)                           :: reactionNum
    integer(shortInt)                           :: i,j, Ngrid
    integer(shortInt)                           :: firstIdx, numIdx

    if (allocated(self % energyGrid)) deallocate(self % energyGrid)
    if (associated(self % xsData)) deallocate(self % xsData)

    ! Check if nuclide is fissile
    self % isFissile = ACE % isFissile()

    ! Load Energy Grid
    call ACE % ESZblock(self % energyGrid,'energyGrid')

    ! Load cross-sections from ESZ block of ACE data card.
    call ACE % ESZblock(xsCapture,'absorptionXS')  ! Total capture XS
    call ACE % ESZblock(xsEScatter,'elasticXS')    ! Elastic Scattering XS

    ! Find size of energy grid
    Ngrid = size( self % energyGrid)

    ! Load fission data if present
    if (self % isFissile) then
      ! First index on energy grid for which fission data is present
      i = ACE % firstIdxFiss()

      ! Number of consecutive eneteries for data
      j = ACE % numXsPointsFiss()

      ! Allocate adn load fission XSs
      allocate(xsFission(Ngrid))
      xsFission(i:i+j-1) = ACE % xsFiss()

    end if

    ! Allocate xsData
    allocate(self % xsData(size(self % energyGrid)))

    ! Move cross section data into xsEnergyPointType
    call self % xsData % load (xsEScatter,xsCapture,xsFission)

    !**** DEBUG pass NXS, JXS and XSS directly
    associate( NXS => ACE % NXS, JXS => ACE % JXS, XSS => ACE % XSS)

    ! Read emission data
    self % emissionData(captureIdx)   = emissionFromACE(ACE,NXS,JXS,XSS,N_gamma)
    self % emissionData(escatterIdx)  = emissionFromACE(ACE,NXS,JXS,XSS,N_N_elastic)

    if (self % isFissile) then
      self % emissionData(fissionIdx) = emissionFromACE(ACE,NXS,JXS,XSS,N_fission)

    else
      ! Put placeholder emission type equivalent to capture for non=present fission
      self % emissionData(fissionIdx) = emissionFromACE(ACE,NXS,JXS,XSS,N_gamma)

    end if

    end associate
    !***** DEBUG ZONE END
  end subroutine readXS


end module aceNoMT_class


