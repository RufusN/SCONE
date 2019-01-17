module tallyClerk_inter

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  implicit none
  private

  !!
  !! Abstract interface for a single tallyClerk.
  !! It recives reports from the admin and processed them into scores and estimates.
  !! Its responsibilites are as follows:
  !! 1) Score some result by accepting a subset of all avalible reports
  !! 2) Display implementation determined measure of convergance (usually some variance)
  !! 3) Can return information about reports it requires
  !!
  type, public,abstract :: tallyClerk
    private
    ! Location of the tallyClerk on scoreMemory
    integer(longInt) :: memAdress = -1

  contains
    !! File reports and check status
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: isConverged

    !! Assign memory and check memory id requested
    procedure                  :: setMemAddress
    procedure, non_overridable :: getMemAddress


    procedure(validReports), deferred :: validReports
    procedure(display), deferred      :: display
    procedure(init),deferred          :: init
    procedure(print),deferred         :: print

  end type tallyClerk

  abstract interface
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    function validReports(self) result(validCodes)
      import :: tallyClerk ,&
                shortInt
      class(tallyClerk),intent(in)               :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
    end function validReports

    !!
    !! Display convergance progress on the console
    !!
    subroutine display(self)
      import :: tallyClerk
      class(tallyClerk), intent(in)  :: self
    end subroutine display

    !!
    !! Initialise tally clerk from a dictionary
    !!
    subroutine init(self,dict,name)
      import :: tallyClerk, &
                dictionary, &
                nameLen
      class(tallyClerk), intent(inout) :: self
      class(dictionary), intent(in)    :: dict
      character(nameLen), intent(in)   :: name
    end subroutine init

    !!
    !! Write contents of the clerk to output file
    !!
    subroutine print(self,outFile)
      import :: tallyClerk, &
                outputFile
      class(tallyClerk), intent(in)    :: self
      class(outputFile), intent(inout) :: outFile
    end subroutine print


  end interface

contains

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self,p)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    character(100),parameter    :: Here = 'reportInColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self,p,MT,muL)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    character(100),parameter  :: Here = 'reportOutColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,p,L)
    class(tallyClerk), intent(inout)     :: self
    class(particle), intent(in)          :: p
    real(defReal), intent(in)            :: L
    character(100),parameter  :: Here = 'reportPath (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,p)
    class(tallyClerk), intent(inout)     :: self
    class(particle), intent(in)          :: p
    character(100),parameter  :: Here = 'reportTrans (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! Particle is associated with one of the fate codes in tallyCodes
  !!
  subroutine reportHist(self,p)
    class(tallyClerk), intent(inout) :: self
    class(particle), intent(in)      :: p
    character(100),parameter  :: Here = 'reportHist (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start
    character(100),parameter  :: Here = 'reportCycleStart (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    character(100),parameter  :: Here = 'reportCycleEnd (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleEnd

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(tallyClerk), intent(in) :: self
    logical(defBool)              :: isIt
    character(100),parameter  :: Here = 'isConverged (tallyClerk_inter.f90)'

    call fatalError(Here,'Convergence check is not implemented in the instance')

    ! Avoid warning
    isIt =.false.

  end function isConverged

  !!
  !! Set memory adress for the clerk
  !!
  subroutine setMemAddress(self, addr)
    class(tallyClerk), intent(inout) :: self
    integer(longInt), intent(in)     :: addr

    self % memAdress = addr

  end subroutine setMemAddress

  !!
  !! Return memory address of the clerk
  !!
  elemental function getMemAddress(self) result(addr)
    class(tallyClerk), intent(in) :: self
    integer(longInt)              :: addr

    addr = self % memAdress

  end function getMemAddress


end module tallyClerk_inter
