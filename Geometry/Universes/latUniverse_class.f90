module latUniverse_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use box_class,          only : box
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill

  implicit none
  private

  ! Parameters
  ! Note X/Y/Z MIN/MAX MUST have the value they have (or generated by a function)
  integer(shortInt), parameter :: X_MIN = -1, X_MAX = -2, Y_MIN = -3, Y_MAX = -4, Z_MIN = -5, &
                                  Z_MAX = -6, OUTLINE_SURF = -7

  !!
  !! 2D or 3D Cartesian lattice with constant pitch
  !!
  !! Universe consists of a lattice of fixed, finite size (e.g 17x17x2). Centre of the
  !! lattice is placed at the origin. An additional cell is placed beyond the lattice
  !! called background (or out) cell.
  !!
  !! Local ID is 1 in bottom X, Y & Z corner. It increases first with X then Y and lastly Z.
  !! Cells inside the lattice can only be filled with a universe (given as integer ID).
  !! Background cell can have any filling given by keyword (material or universe)
  !!
  !! Every lattice cell has an offset to its centere (so the centre of the nested universe
  !! is in the center of the lattice cell).
  !!
  !! Minimum lattice pitch is set to 10 * SURF_TOL
  !!
  !! Sample Input Dictionary (3D):
  !!   latt { id 7;
  !!          type latUniverse;
  !!          #origin (0.0 0.0 0.0); #
  !!          #rotation (30.0 0.0 0.0); #
  !!          shape (3 2 2);
  !!          pitch (1.0 1.0 1.0);
  !!          padMat <u13>;
  !!          map (  1  2  3    // Top layer
  !!                 4  5  6    // Lower Y row
  !!                 7  8  9    // Bottom layer
  !!                10 11 12  )
  !!   }
  !!
  !! Sample Input Dictionary (2D):
  !!   latt2D { id 8;
  !!            shape (2 2 0);       // 0 indicates infinite extent in z-axis
  !!            pitch (1.0 1.0 0.0); // Any pitch is allowed in z, use 0.0 for clarity
  !!            padMat void;
  !!            map ( 1 2
  !!                  2 1);
  !!   }
  !!
  !! NOTE: Input in MAP for a single layer is WYSIWYG. Lower row in map is lower row in
  !!   geometry. There is no inversion like in other (e.g. Serpent) MC codes. Basically
  !!   input vector is flipped in Y and Z direction.
  !!
  !! Private Members:
  !!  pitch      -> Values of lattice pitch in x, y & z directions
  !!  sizeN      -> Number of lattice cells in x, y & z directions
  !!  corner     -> Location of the minumum corner
  !!  a_bar      -> Halfwidth of lattice cell reduced by surface tolerance
  !!  outline    -> Box type surface that is a boundary between lattice & background
  !!  outLocalID -> LocalID of the background cell
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: latUniverse
    private
    real(defReal), dimension(3)     :: pitch = ZERO
    integer(shortInt), dimension(3) :: sizeN = 0
    real(defReal), dimension(3)     :: corner = ZERO
    real(defReal), dimension(3)     :: a_bar  = ZERO
    type(box)                       :: outline
    integer(shortInt)               :: outLocalID = 0
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
  end type latUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if input is invalid.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(latUniverse), intent(inout)                         :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    real(defReal), dimension(:), allocatable       :: temp
    integer(shortInt), dimension(:), allocatable   :: tempI
    integer(shortInt)                              :: id, N, i, j, outFill
    type(dictionary)                               :: tempDict
    integer(shortInt), dimension(:,:), allocatable :: tempMap
    character(nameLen)                             :: name
    character(100), parameter :: Here = 'init (latUniverse_class.f90)'

    ! Load basic data
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: '//numToChar(id))
    call self % setId(id)

    ! Load origin
    if (dict % isPresent('origin')) then
      call dict % get(temp, 'origin')

      if (size(temp) /= 3) then
        call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(size(temp)))
      end if
      call self % setTransform(origin=temp)

    end if

    ! Load rotation
    if (dict % isPresent('rotation')) then
      call dict % get(temp, 'rotation')

      if (size(temp) /= 3) then
        call fatalError(Here, '3 rotation angles must be given. Has only: '//numToChar(size(temp)))
      end if
      call self % setTransform(rotation=temp)
    end if

    ! Load pitch
    call dict % get(temp, 'pitch')
    N = size(temp)

    if (N /= 3) then
      call fatalError(Here, 'Pitch must have size 3. Has: '//numToChar(N))
    end if
    self % pitch = temp

    ! Load Size
    call dict % get(tempI, 'shape')
    N = size(tempI)

    if (N /= 3) then
      call fatalError(Here, 'Shape must have size 3. Has: '//numToChar(N))
    else if (any(tempI < 0)) then
      call fatalError(Here, 'Shape contains -ve entries')
    end if
    self % sizeN = tempI

    ! Detect reduced Z dimension
    if (self % sizeN(3) == 0) then
      self % sizeN(3) = 1
      self % pitch(3) = TWO * INF
    end if

    ! Check X & Y for 0 size
    if (any( self % sizeN == 0)) call fatalError(Here, 'Shape in X and Y axis cannot be 0.')

    ! Check for invalid pitch
    if (any(self % pitch < 10 * SURF_TOL)) then
     call fatalError(Here, 'Pitch size must be larger than: '//numToChar( 10 * SURF_TOL))
   end if

    ! Calculate halfwidth and corner
    self % a_bar = self % pitch * HALF - SURF_TOL
    self % corner = -(self % sizeN * HALF * self % pitch)

    ! Calculate local ID of the background
    self % outLocalID = product(self % sizeN) + 1

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', [ZERO, ZERO, ZERO])
    call tempDict % store('halfwidth', abs(self % corner))
    call self % outline % init(tempDict)

    ! Construct fill array
    call dict % get(tempI, 'map')

    ! Flip array up-down for more natural input
    ! Reshape into rank 2 array
    tempMap = reshape(tempI, [self % sizeN(1), self % sizeN(2) * self % sizeN(3)])
    N = size(tempMap, 2)
    do i = 1, N/2
      call swap(tempMap(:,i), tempMap(:,N - i + 1))
    end do

    ! Find background fill and change to tempMap to uniID
    tempMap = -tempMap
    call dict % get(name, 'padMat')
    outFill = charToFill(name, mats, Here)

    ! Build fill array
    allocate(fill( self % outLocalID))
    N = size(tempMap, 1)
    do j = 1, size(tempMap, 2)
      do i = 1, N
        fill(i + (j-1) * N) = tempMap(i, j)
      end do
    end do
    fill(self % outLocalID) = outFill

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(latUniverse), intent(inout)       :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), dimension(3)         :: ijk
    integer(shortInt)                       :: i, inc
    real(defReal), dimension(3)             :: r_bar

    ! Find lattice location in x,y&z
    ijk = floor((r - self % corner) / self % pitch) + 1

    ! Get position wrt middle of the lattice cell
    r_bar = r - self % corner - ijk * self % pitch + HALF * self % pitch

    ! Check if position is within surface tolerance
    ! If it is, push it to next cell
    do i = 1, 3
      if (abs(r_bar(i)) > self % a_bar(i) .and. r_bar(i)*u(i) > ZERO) then

        ! Select increment. Ternary expression
        if (u(i) < ZERO) then
          inc = -1
        else
          inc = 1
        end if

        ijk(i) = ijk(i) + inc
      end if
    end do

    ! Set localID & cellIdx
    if (any(ijk <= 0 .or. ijk > self % sizeN)) then ! Point is outside lattice
      localID = self % outLocalID

    else
      localID = ijk(1) + self % sizeN(1) * (ijk(2)-1 + self % sizeN(2) * (ijk(3)-1))

    end if
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(latUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    real(defReal), dimension(3)        :: r_bar, u, bounds
    real(defReal)                      :: test_d
    integer(shortInt)                  :: i, ax

    ! Catch case if particle is outside the lattice
    if (coords % localID == self % outLocalID) then
      surfIdx = OUTLINE_SURF
      d = self % outline % distance(coords % r, coords % dir)
      return

    end if

    ! Find position wrt lattice cell centre
    ! Need to use localID to properly handle under and overshoots
    u = coords % dir
    r_bar = coords % r - self % corner
    r_bar = r_bar - (get_ijk(coords % localID, self % sizeN) - HALF) * self % pitch

    ! Select surfaces in the direction of the particle
    bounds = sign(self % pitch * HALF, u)

    ! Find minimum distance
    ! Relay on IEEE 754 standard (for floating point numbers)
    ! 0.0/0.0 = NaN and (NaN < A = false; for every A)
    ! A/0.0 = Infinity (if A > 0.0)
    !
    ! Provide default axis to ensure no out of bounds array access if
    ! all distances happen to be infinite
    d = INF
    ax = 1 
    do i = 1, 3
      ! Nominator and denominator will have the same sign (by ealier bounds selection)
      test_d = (bounds(i) - r_bar(i)) / u(i)

      if (test_d < d) then
        d = test_d
        ax = i
      end if
    end do

    ! Cap distance value
    d = max(ZERO, d)
    d = min(INF, d)

    ! Generate surface memento
    surfIdx = ax * 2
    if (u(ax) < ZERO) surfIdx = surfIdx - 1
    surfIdx = -surfIdx

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  subroutine cross(self, coords, surfIdx)
    class(latUniverse), intent(inout) :: self
    type(coord), intent(inout)        :: coords
    integer(shortInt), intent(in)     :: surfIdx

    call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(latUniverse), intent(in)  :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    if (coords % localID == self % outLocalID) then
      offset = ZERO

    else
      offset = (get_ijk(coords % localID, self % sizeN) - HALF) * self % pitch + self % corner

    end if

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(latUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    self % pitch = ZERO
    self % sizeN = 0
    self % corner = ZERO
    self % a_bar  = ZERO
    call self % outline % kill()
    self % outLocalID = 0

  end subroutine kill

  !!
  !! Generate ijk from localID and shape
  !!
  !! Args:
  !!   localID [in] -> Local id of the cell between 1 and product(sizeN)
  !!   sizeN [in]   -> Number of cells in each cardinal direction x,y&z
  !!
  !! Result:
  !!   Array ijk which has integer position in each cardinal direction
  !!
  pure function get_ijk(localID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: localID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = localID - 1

    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = temp - sizeN(2) * base + 1

    ijk(3) = base + 1

  end function get_ijk

end module latUniverse_class
