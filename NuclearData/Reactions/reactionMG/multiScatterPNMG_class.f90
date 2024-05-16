module multiScatterPNMG_class

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, numToChar
  use legendrePoly_func,     only : sampleLegendre
  use RNG_class,             only : RNG
  use reactionHandle_inter,  only : reactionHandle
  use dataDeck_inter,        only : dataDeck
  use dictDeck_class,        only : dictDeck
  use dictionary_class,      only : dictionary
  use multiScatterMG_class,  only : multiScatterMG, kill_super => kill, &
                                    buildFromDict_super => buildFromDict



  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: multiScatterPNMG_TptrCast

  !!
  !! MG multiplicative scattering
  !!
  !! Public Members:
  !!   order -> scattering order
  !!   PN    -> contains scattering matrices from P1 up to P7 max
  !!
  !! Interface
  !!   multiScatterMG interface
  !!
  !! NOTE:
  !!   It is sufficient to extend concreate build procedure. We can still use init in the
  !!   superclass!
  !!
  type, public, extends(multiScatterMG) :: multiScatterPNMG
    integer(shortInt) :: order
    real(defReal), dimension(:,:,:), allocatable :: PN

  contains
    ! Override some procedures
    procedure :: kill
    procedure :: sampleOut
    procedure :: getPnScatter
    procedure :: buildFromDict
  end type multiScatterPNMG

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multiScatterPNMG), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Clean own memory
    if(allocated(self % PN)) deallocate(self % PN)

  end subroutine kill

  !!
  !! Sample outgoing particle
  !!
  !! See reactionMG documentation for details
  !!
  subroutine sampleOut(self, mu, phi, G_out, G_in, rand)
    class(multiScatterPNMG), intent(in)   :: self
    real(defReal), intent(out)     :: mu
    real(defReal), intent(out)     :: phi
    integer(shortInt), intent(out) :: G_out
    integer(shortInt), intent(in)  :: G_in
    class(RNG), intent(inout)      :: rand
    character(100),parameter :: Here = 'sampleOut (multiScatterPNMG_class.f90)'

    ! Sample G_out
    G_out = self % sampleGout(G_in, rand)

    ! Sample deflection
    phi = TWO_PI * rand % get()
    mu  = TWO * rand % get() - ONE

  end subroutine sampleOut

  !!
  !! Public function to read scattering matrices
  !!
  !! Args:
  !!  G_out [in] -> outgoing energy group
  !!  G_in [in]  -> ingoing energy group
  !!  N [in]     -> scattering order
  !!
  !! Errors:
  !!   FatalError if size of scattering order requested is too high
  !!
  function getPnScatter(self, G_in, G_out, N) result(xs)
    class(multiScatterPNMG), intent(in) :: self
    integer(shortInt), intent(in)       :: G_in
    integer(shortInt), intent(in)       :: G_out
    integer(shortInt), intent(in)       :: N
    real(defReal)                       :: xs
    character(100),parameter :: Here = 'getPnScatter (multiScatterPNMG_class.f90)'

    if (N > self % order) call fatalError(Here, 'Scattering order requested is too high')

    xs = self % PN(G_out,G_in,N)

  end function getPnScatter

  !!
  !! Builds multiScatterPNMG from SCONE dictionary
  !!
  !! Extends multiScatterMG procedure!
  !! See its documentation for extra details.
  !!
  !! Errors:
  !!   FatalError if size of scattering matrices does not match numer of group
  !!
  subroutine buildFromDict(self, dict)
    class(multiScatterPNMG), intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    real(defReal),dimension(:),allocatable :: temp
    integer(shortInt)                      :: nG, i
    character(nameLen)                     :: name
    character(100),parameter :: Here = 'buildFromDict (multiScatterPNMG_class.f90)'

    ! Call superclass procedure
    call buildFromDict_super(self, dict)

    ! Re-read number of groups
    call dict % get(nG,'numberOfGroups')

    ! Determine the scattering order - only odd options
    if (dict % isPresent('P7')) then
      self % order = 7
    elseif (dict % isPresent('P5')) then
      self % order = 5
    elseif (dict % isPresent('P3')) then
      self % order = 3
    elseif (dict % isPresent('P2')) then
      self % order = 2
    elseif (dict % isPresent('P1')) then
      self % order = 1
    else
      call fatalError(Here, 'Scatteting order must be at least P1. No data in the xss file. ')
    end if

    ! Allocate scattering matrix
    allocate(self % PN(nG, nG, self % order))

    do i = 1,self % order

      ! Read right order
      name = 'P'//numToChar(i)

      ! Read scattering matrix
      call dict % get(temp, trim(name))

      ! Check size
      if( size(temp) /= nG*nG) then
        call fatalError(Here,'Invalid size of matrix P'//numToChar(i)//'. Expected: '//numToChar(nG**2)//&
                             ' got: '//numToChar(size(temp)))
      end if

      ! Assign values
      self % PN(:,:,i) = reshape(temp,[nG, nG])

    end do

  end subroutine buildFromDict

  !!
  !! Cast reactionHandle pointer to multiScatterP1MG pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of multiScatterP1MG type
  !!   Target points to source if source is multiScatterP1MG type
  !!
  pure function multiScatterPNMG_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(multiScatterPNMG), pointer              :: ptr

    select type(source)
    type is(multiScatterPNMG)
        ptr => source

      class default
        ptr => null()
    end select

  end function multiScatterPNMG_TptrCast


end module multiScatterPNMG_class
