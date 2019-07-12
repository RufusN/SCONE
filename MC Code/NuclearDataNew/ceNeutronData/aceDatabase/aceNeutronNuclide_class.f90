module aceNeutronNuclide_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG
  use aceCard_class,     only : aceCard
  use stack_class,       only : stackInt

  ! Nuclear Data Interfaces
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase
  use ceNeutronNuclide_inter,       only : ceNeutronNuclide, kill_super => kill
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE
  use elasticNeutronScatter_class,  only : elasticNeutronScatter
  use fissionCE_class,              only : fissionCE
  use neutronScatter_class,         only : neutronScatter
  use pureCapture_class,            only : pureCapture

  implicit none
  private

  ! Grid location parameters
  integer(shortInt), parameter :: ENERGY_GRID   = 1
  integer(shortInt), parameter :: TOTAL_XS      = 2
  integer(shortInt), parameter :: ESCATTER_XS   = 3
  integer(shortInt), parameter :: IESCATTER_XS  = 4
  integer(shortInt), parameter :: CAPTURE_XS    = 5
  integer(shortInt), parameter :: FISSION_XS    = 6
  integer(shortInt), parameter :: NU_FISSION    = 7


  !!
  !!
  !!
  type, public :: reactionMT
    integer(shortInt)                         :: MT       = 0
    integer(shortInt)                         :: firstIdx = 0
    real(defReal),dimension(:),allocatable    :: xs
    class(uncorrelatedReactionCE),allocatable :: kinematics
  end type reactionMT

  !!
  !! Doc Doc Docstring
  !! IGNORES MT=5 (N_ANYTHING) in its current implementation !
  !! IN JEF 3.1.1 It will Reduce accuracy of Tc99 collision processing
  !!
  !!
  type, public, extends(ceNeutronNuclide) :: aceNeutronNuclide
    character(nameLen)                          :: ZAID    = ''
    real(defReal), dimension(:,:), allocatable  :: mainData
    type(reactionMT), dimension(:), allocatable :: MTdata
    integer(shortInt)                           :: nMT     = 0

    type(elasticNeutronScatter) :: elasticScatter
    type(fissionCE)             :: fission

  contains
    ! Superclass Interface
    procedure :: invertInelastic
    procedure :: xsOf

    ! Local interface
    procedure :: init
    procedure :: kill
    procedure :: display

  end type aceNeutronNuclide

contains

  !!
  !! Invert PDF of inelastic stattering
  !!
  !! See ceNeutronNuclide documentation
  !!
  function invertInelastic(self, E, rand) result(MT)
    class(aceNeutronNuclide), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: MT

    MT = 1

  end function invertInelastic

  !!
  !! Return Cross-Section of reaction MT at energy E
  !!
  !! See ceNeutronNuclide documentation
  !!
  function xsOf(self, MT, E) result(xs)
    class(aceNeutronNuclide), intent(in) :: self
    integer(shortInt), intent(in)        :: MT
    real(defReal), intent(in)            :: E
    real(defReal)                        :: xs

    xs = ONE

  end function xsOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aceNeutronNuclide), intent(inout) :: self

    ! Call superclass
    call kill_super(self)


  end subroutine kill

  !!
  !! Initialise from an ACE Card
  !!
  !! Args:
  !!   ACE [inout]   -> ACE card
  !!   nucIdx [in]   -> nucIdx of the nuclide
  !!   database [in] -> pointer to ceNeutronDatabase to set XSs to cache
  !!
  !! Errors:
  !!
  !!
  subroutine init(self, ACE, nucIdx, database)
    class(aceNeutronNuclide), intent(inout)       :: self
    class(aceCard), intent(inout)                 :: ACE
    integer(shortInt), intent(in)                 :: nucIdx
    class(ceNeutronDatabase), pointer, intent(in) :: database
    integer(shortInt)                             :: Ngrid, N, K, i, j, MT, bottom, top
    type(stackInt)                                :: scatterMT, absMT

    ! Reset nuclide just in case
    call self % kill()

    ! Load ACE ZAID
    self % ZAID = ACE % ZAID

    ! Read key data into the superclass
    call self % set( fissile  = ACE % isFissile(), &
                     mass     = ACE % AW,          &
                     kT       = ACE % TZ,          &
                     nucIdx   = nucIdx,            &
                     database = database )

    ! Get size of the grid
    Ngrid = ACE % gridSize()

    ! Allocate space for main XSs
    if(self % isFissile()) then
      N = 7
    else
      N = 5
    end if
    allocate(self % mainData(Ngrid, N))

    self % mainData = ZERO

    ! Load Main XSs
    self % mainData(:,ENERGY_GRID)  = ACE % ESZ_XS('energyGrid')
    self % mainData(:,TOTAL_XS)     = ACE % ESZ_XS('totalXS')
    self % mainData(:,ESCATTER_XS)  = ACE % ESZ_XS('elasticXS')
    self % mainData(:,CAPTURE_XS)   = ACE % ESZ_XS('absorptionXS')

    ! Get elastic kinematics
    call self % elasticScatter % init(ACE, N_N_ELASTIC)

    if(self % isFissile()) then
      call self % fission % init(ACE, N_FISSION)
      N = ACE % firstIdxFiss()
      K = ACE % numXSPointsFiss()
      self % mainData(N:N+K-1,FISSION_XS) = ACE % xsFiss()

      ! Calculate nuFission
      do i = N, K
        self % mainData(i,NU_FISSION) = self % mainData(i,FISSION_XS) * &
                                        self % fission % release(self % mainData(i,ENERGY_GRID))
      end do
    end if

    ! Read data for MT reaction

    ! Create a stack of MT reactions, devide them into ones that produce 2nd-ary
    ! particlues and pure absorbtion
    associate (MTs => ACE % getScatterMTs())
      do i=1,size(MTs)
        if (MTs(i) == N_ANYTHING) cycle
        call scatterMT % push(MTs(i))
      end do
    end associate

    associate (MTs => [ACE % getFissionMTs(), ACE % getCaptureMTs()])
      do i=1,size(MTs)
        if(MTs(i) == N_FISSION) cycle ! MT=18 is already included with FIS block
        call absMT % push(MTs(i))
      end do
    end associate

    ! Allocate space
    allocate(self % MTdata(scatterMT % size() + absMT % size()))

    ! Load scattering reactions
    N = scatterMT % size()
    self % nMT = N
    do i =1,N
      call scatterMT % pop(MT)
      self % MTdata(i) % MT       = MT
      self % MTdata(i) % firstIdx = ACE % firstIdxMT(MT)
      self % MTdata(i) % xs       = ACE % xsMT(MT)

      allocate(neutronScatter :: self % MTdata(i) % kinematics)
      call self % MTdata(i) % kinematics % init(ACE, MT)
    end do

    ! Load capture reactions
    K = absMT % size()
    do i = N+1,N+K
      call absMT % pop(MT)
      self % MTdata(i) % MT       = MT
      self % MTdata(i) % firstIdx = ACE % firstIdxMT(MT)
      self % MTdata(i) % xs       = ACE % xsMT(MT)

      allocate(pureCapture :: self % MTdata(i) % kinematics)
      call self % MTdata(i) % kinematics % init(ACE, MT)
    end do

    ! Calculate Inelastic scattering XS
    do i=1,self % nMT
      do j=1,size(self % mainData, 1)
        ! Find bottom and Top of the grid
        bottom = self % MTdata(i) % firstIdx
        top    = size(self % MTdata(i) % xs)
        if( j>= bottom .and. j <= top + bottom) then
          self % mainData(j, IESCATTER_XS) = self % mainData(j, IESCATTER_XS) + &
                                             self % MTdata(i) % xs(j-bottom)
        end if
      end do
    end do

    ! Recalculate totalXS
    if(self % isFissile()) then
      K = FISSION_XS
    else
      K = CAPTURE_XS
    end if
    self % mainData(:, TOTAL_XS) = sum(self % mainData(:,ESCATTER_XS:K))

  end subroutine init

  !!
  !! A Procedure that displays information about the nuclide to the screen
  !!
  !! Prints:
  !!   nucIdx; Size of Energy grid; Maximum and Minumum energy; Mass; Temperature; ACE ZAID;
  !!   MT numbers used in tranposrt (active MTs) and not used in transport (inactiveMTs)
  !!
  !! NOTE: For Now Formatting is horrible. Can be used only for debug
  !!       MT = 18 may appear as inactive MT but is included from FIS block !
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  subroutine display(self)
    class(aceNeutronNuclide), intent(in)  :: self
    integer(shortInt)                     :: N, sMT, allMT
    real(defReal)                         :: E_min, E_max, M, kT

    ! Print Most relevant information
    print *, repeat('#',40)
    print '(A)', "Nuclide: " // trim(self % ZAID)

    N = size(self % mainData,1)
    E_min = self % mainData(1, ENERGY_GRID)
    E_max = self % mainDATA(N, ENERGY_GRID)
    M = self % getMass()
    kT = self % getkT()
    print '(A)', "Mass: " //numToChar(M) //" Temperature: "//numToChar(kT) //" [MeV]"
    print '(A)', "Energy grid: " //numToChar(N)// " points from "// &
                  numToChar(E_min)//" to "// numToChar(E_max) // " [MeV] "

    ! Print MT information
    sMT = self % nMT
    allMT = size(self % MTdata)
    print '(A)', "Active MTs: "  // numToChar(self % MTdata(1:sMT) % MT)
    print '(A)', "Inactive MTs: "// numToChar(self % MTdata(sMT+1:allMT) % MT)
    print '(A)', "This implementation ignores MT=5 (N,anything) !"

  end subroutine display
end module aceNeutronNuclide_class