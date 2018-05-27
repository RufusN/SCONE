module tallyResponse_inter

  use numPrecision
  use particle_class, only : particle, phaseCoord

  implicit none
  private

  !!
  !! Abstract interface for collision tally response definition
  !!
  type, public,abstract :: tallyResponse
    integer(shortInt) :: responseCode = 0
  contains
    procedure(getScore),deferred :: getScore
  end type tallyResponse

  abstract interface

    !!
    !! Given collision parameters return score function (f(x)) value for collision
    !! * f needs to be multiplied by flux estimate before scoring
    !!
    function getScore(self,pre,post,MT,muL) result (f)
      import :: tallyResponse, &
                particle, &
                phaseCoord, &
                shortInt, &
                defReal
      class(tallyResponse), intent(inout) :: self
      class(phaseCoord), intent(in)  :: pre
      class(particle), intent(in)    :: post
      integer(shortInt), intent(in)  :: MT
      real(defReal), intent(in)      :: muL
      real(defReal)                  :: f
    end function getScore

  end interface
    
end module tallyResponse_inter