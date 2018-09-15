module dictionary_test
  use numPrecision
  use dictionary_class, only : dictionary
  use pfunit_mod
  implicit none

  public :: test_dictionary

@TestCase
  type, extends(TestCase) :: test_dictionary
    type(dictionary) :: dict
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_dictionary

  !! Parameters
  real(defReal),parameter                  :: realVal     = 3.3_defReal
  integer(shortInt), parameter             :: intVal      = 1_shortInt
  character(nameLen),parameter             :: charNameLen = 'GoFortran_DownWithCpp'
  character(pathLen), parameter            :: charPathLen ='/home/KyloRen/VaderFanFic'
  real(defReal), dimension(2), parameter   :: realArray = [-1.0E-17_defReal, 14.7_defReal]
  integer(shortInt),dimension(3),parameter :: intArray =[-6475_shortInt, 13_shortInt, 14_shortInt]
  character(nameLen), dimension(1), parameter :: charNameLenArray = ['TK-421']
  character(pathLen), dimension(2), parameter :: charPathLenArray = ['C:\User\Tarkin\DeathStarPlans              ', &
                                                                     '/home/Dodonna/Articles/whyRebelsUseUNIX.odt']


contains

  !!
  !! Sets up test_dictionary object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary) :: tempDict

    call this % dict % init(1)
    call this % dict % store('myReal', realVal)
    call this % dict % store('myInt', intVal )
    call this % dict % store('myCharNameLen', charNameLen)
    call this % dict % store('myCharPathLen', charPathLen)
    call this % dict % store('realArray', realArray)
    call this % dict % store('intArray', intArray)
    call this % dict % store('charNameLenArray', charNameLenArray)
    call this % dict % store('charPathLenArray', charPathLenArray)

    tempDict = this % dict
    call this % dict % store('nestedDict', tempDict)

  end subroutine setUp

  !!
  !! SKills test_dictionary object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_dictionary), intent(inout) :: this

    call this % dict % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!  TESTS PROPER BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test extracting Real value
!!
@test
  subroutine testGettingReal(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal)                            :: tempReal

    call this % dict % get(tempReal,'myReal')
    @assertEqual(realVal, tempReal, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(tempReal,'myReal',7.0_defReal)
    @assertEqual(realVal, tempReal, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(tempReal,'invalid',7.0_defReal)
    @assertEqual(7.0_defReal, tempReal, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingReal

!!
!! Test extracting Real Array
!!
@test
  subroutine testGettingRealArray(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal),dimension(:),allocatable   :: tempReal

    call this % dict % get(tempReal,'realArray')
    @assertEqual(realArray, tempReal, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(tempReal,'realArray', [7.0_defReal])
    @assertEqual(realArray, tempReal, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(tempReal,'myRealNon',[7.0_defReal])
    @assertEqual([7.0_defReal], tempReal, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingRealArray

!!
!! Test extracting Integer Value
!!
@test
  subroutine testGettingInt(this)
    class(test_dictionary), intent(inout)    :: this
    integer(shortInt)                        :: temp

    call this % dict % get(temp,'myInt')
    @assertEqual(intVal, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myInt',7_shortInt)
    @assertEqual(intVal, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', 7_shortInt)
    @assertEqual(7_shortInt, temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingInt

!!
!! Test extracting Integer Array
!!
@test
  subroutine testGettingIntArray(this)
    class(test_dictionary), intent(inout)      :: this
    integer(shortInt),dimension(:),allocatable :: temp

    call this % dict % get(temp,'intArray')
    @assertEqual(intArray, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'intArray',[7_shortInt])
    @assertEqual(intArray, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', [7_shortInt])
    @assertEqual([7_shortInt], temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingIntArray

!!
!! Test extracting nameLen long character
!!
@test
  subroutine testGettingNameLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(nameLen)                    :: temp
    character(nameLen),parameter :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharNameLen')
    @assertEqual(charNameLen, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myCharNameLen',default)
    @assertEqual(charNameLen, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    @assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingNameLenChar

!!
!! Test extracting pathLen long character
!!
@test
  subroutine testGettingPathLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(pathLen)                    :: temp
    character(pathLen),parameter  :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharPathLen')
    @assertEqual(charPathLen, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myCharPathLen', default)
    @assertEqual(charPathLen, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    @assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingPathLenChar

!!
!! Test extracting nameLen long array of Chars
!!
@test
  subroutine testGettingNameLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: temp
    character(nameLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charNameLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charNameLenArray == temp)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'charNameLenArray', default)
    isSame = all(charNameLenArray == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

  end subroutine testGettingNameLenCharArray

!!
!! Test extracting pathLen long array of Chars
!!
@test
  subroutine testGettingPathLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(pathLen),dimension(:),allocatable  :: temp
    character(pathLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charPathLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charPathLenArray == temp)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'charPathLenArray', default)
    isSame = all(charPathLenArray == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

  end subroutine testGettingPathLenCharArray

@test
  subroutine testGettingNestedDictionary(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary)                      :: temp
    real(defReal)                         :: tempReal
    integer(shortInt)                     :: tempInt
    character(nameLen)                    :: tempCharNameLen
    character(pathLen)                    :: tempCharPathLen
    real(defReal),dimension(:),allocatable        :: tempRealArray
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    character(nameLen), dimension(:), allocatable :: tempCharArrayNameLen
    character(pathLen), dimension(:), allocatable :: tempCharArrayPathLen
    logical(defBool)                      :: isSame

    call this % dict % get(temp,'nestedDict')

    ! Get all contents of the dictionary
    call temp % get(tempReal, 'myReal')
    call temp % get(tempint, 'myInt')
    call temp % get(tempCharNameLen, 'myCharNameLen')
    call temp % get(tempCharPathLen, 'myCharPathLen')
    call temp % get(tempRealArray, 'realArray')
    call temp % get(tempIntArray, 'intArray')
    call temp % get(tempCharArrayNameLen, 'charNameLenArray')
    call temp % get(tempCharArrayPathLen, 'charPathLenArray')

    ! Verify that content was not deformed
    isSame = tempReal == realVal
    isSame = isSame .and. tempInt == intVal
    isSame = isSame .and. tempCharNameLen == charNameLen
    isSame = isSame .and. tempCharPathLen == charPathLen
    isSame = isSame .and. all(tempIntArray == intArray)
    isSame = isSame .and. all(tempRealArray == realArray)
    isSame = isSame .and. all(tempCharArrayNameLen == charNameLenArray)
    isSame = isSame .and. all(tempCharArrayPathLen == charPathLenArray)

    @assertTrue(isSame, 'Contents of nested dictionary were changed')

  end subroutine testGettingNestedDictionary

!!
!! Test keys retrival
!!
@test
  subroutine testKeys(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: tempAll
    character(nameLen),dimension(:),allocatable  :: tempReal
    character(nameLen),dimension(:),allocatable  :: tempInt
    character(nameLen),dimension(:),allocatable  :: tempChar
    character(nameLen),dimension(:),allocatable  :: tempRealArray
    character(nameLen),dimension(:),allocatable  :: tempIntArray
    character(nameLen),dimension(:),allocatable  :: tempCharArray
    character(nameLen),dimension(:),allocatable  :: tempDict
    character(nameLen)                           :: keyword
    logical(defBool) :: isValid

    ! Obtain all keys arrays
    call this % dict % keys(tempAll)
    call this % dict % keysReal(tempReal)
    call this % dict % keysInt(tempInt)
    call this % dict % keysChar(tempChar)
    call this % dict % keysDict(tempDict)

    call this % dict % keysRealArray(tempRealArray)
    call this % dict % keysIntArray(tempIntArray)
    call this % dict % keysCharArray(tempCharArray)

    ! Verify keys for real
    keyword = 'myReal'
    isValid = contains(tempAll, keyword) .and. contains(tempReal, keyword)
    @assertTrue(isValid, 'Keywords failed for real')

    ! Verify keys for int
    keyword = 'myInt'
    isValid = contains(tempAll, keyword) .and. contains(tempInt, keyword)
    @assertTrue(isValid,' Keywords failed for int')

    ! Verify keys for char
    keyword = 'myCharNameLen'
    isValid = contains(tempAll, keyword) .and. contains(tempChar, keyword)
    keyword = 'myCharPathLen'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempChar, keyword)
    @assertTrue(isValid,' Keywords failed for char')

    ! Verify keys for dict
    keyword = 'nestedDict'
    isValid = contains(tempAll, keyword) .and. contains(tempDict, keyword)
    @assertTrue(isValid,' Keywords failed for nested dictionary')

    ! Verify keys for realArray
    keyword = 'realArray'
    isValid = contains(tempAll, keyword) .and. contains(tempRealArray, keyword)
    @assertTrue(isValid, 'Keywords failed for real array ')

    ! Verify keys for int
    keyword = 'intArray'
    isValid = contains(tempAll, keyword) .and. contains(tempIntArray, keyword)
    @assertTrue(isValid,' Keywords failed for int array')

    ! Verify keys for char
    keyword = 'charNameLenArray'
    isValid = contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
    keyword = 'charPathLenArray'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
    @assertTrue(isValid,' Keywords failed for char array')


  contains
    function contains(array,keyword) result(doesIt)
      character(nameLen), dimension(:) :: array
      character(*)                     :: keyword
      logical(defBool)                 :: doesIt

      doesIt = count(array == keyword) == 1

    end function contains
  end subroutine testKeys

  !!
  !! Test isPresent function of a dictionary
  !!
@test
  subroutine testIsPresent(this)
    class(test_dictionary), intent(inout) :: this
    logical(defBool)                      :: isPresent

    isPresent = this % dict % isPresent('nestedDict')
    @assertTrue(isPresent)

    isPresent = this % dict % isPresent('invalid')
    @assertFalse(isPresent)

  end subroutine testIsPresent

end module dictionary_test
