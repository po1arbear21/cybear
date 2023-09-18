module test_util_m

  use string_m
  use test_case_m, only: test_case
  use util_m,      only: c2fstring, cstrlen, f2cstring, int2str, split_string, select_int

  implicit none

  private
  public test_util

contains

  subroutine test_util()
    type(test_case) :: tc

    call tc%init("util")

    ! int2str
    block
      character(:), allocatable :: cval
      integer                   :: i
      type(string)              :: exp, val

      allocate (character(0) :: cval)   ! removes gfortran warning

      ! standard call, no min_len specified
      i     = 123
      cval  = int2str(i)
      val%s = cval
      exp%s = "123"
      call tc%assert_eq(exp, val, "int2str")

      i     = -4567
      cval  = int2str(i)
      val%s = cval
      exp%s = "-4567"
      call tc%assert_eq(exp, val, "int2str")

      ! min_len smaller than actual length
      i     = 89012
      cval  = int2str(i, min_len=2)
      val%s = cval
      exp%s = "89012"
      call tc%assert_eq(exp, val, "int2str")

      i     = -123
      cval  = int2str(i, min_len=2)
      val%s = cval
      exp%s = "-123"
      call tc%assert_eq(exp, val, "int2str")

      ! min_len equal actual length
      i     = 1234
      cval  = int2str(i, min_len=4)
      val%s = cval
      exp%s = "1234"
      call tc%assert_eq(exp, val, "int2str")

      i     = -1234
      cval  = int2str(i, min_len=4)
      val%s = cval
      exp%s = "-1234"
      call tc%assert_eq(exp, val, "int2str")

      ! min_len larger than actual length
      i     = 56
      cval  = int2str(i, min_len=6)
      val%s = cval
      exp%s = "000056"
      call tc%assert_eq(exp, val, "int2str")

      i     = -56
      cval  = int2str(i, min_len=6)
      val%s = cval
      exp%s = "-00056"
      call tc%assert_eq(exp, val, "int2str")
    end block

    ! cstrlen (length of a string)
    block
      character(1), allocatable :: c(:)
      integer                   :: length(3), i
      type(string)              :: str(3)

      str(1)%s = "Hello World"
      str(2)%s = "Ki.D.d!ng"
      str(3)%s = ""

      length = [11, 9, 0]

      do i = 1, 3
        c = f2cstring(str(i)%s)
        call tc%assert_eq(length(i), cstrlen(c), "cstrlen")
      end do
    end block

    ! c2fstring (converting c to f strings)
    block
      character(1), allocatable :: c(:)
      integer                   :: i
      type(string)              :: str, str_exp(3)

      str_exp(1)%s = "Hello World"
      str_exp(2)%s = "Ki.D.d!ng"
      str_exp(3)%s = ""

      do i = 1, 3
        c     = f2cstring(str_exp(i)%s)
        str%s = c2fstring(c)
        call tc%assert_eq(str_exp(i), str, "c2f")
      end do
    end block

    ! split_string
    block
      ! function split_string(str, delim) result(res)

      character(:), allocatable :: str
      type(string), allocatable :: r(:), r_exp(:)

      str = " ;;  aaa bb  c;ddddd;;,:; e,:ffffffff   gg  ,,, "

      allocate (r_exp(7))
      r_exp(1)%s = "aaa"
      r_exp(2)%s = "bb"
      r_exp(3)%s = "c"
      r_exp(4)%s = "ddddd"
      r_exp(5)%s = "e"
      r_exp(6)%s = "ffffffff"
      r_exp(7)%s = "gg"

      r = split_string(str, [" ", ";", ",", ":"])
      call tc%assert_eq(r_exp, r, "split_string")
    end block

    ! select_int
    block
      integer :: ints(5)
      logical :: mask(5)

      ! exactly one hit
      ints = [     -1,     500,      3,       4,       5]
      mask = [.false., .false., .true., .false., .false.]
      call tc%assert_eq(3, select_int(mask, ints), "select_int 1")

      ! two hits
      ints = [    -3,     500,     -3,       4,       5]
      mask = [.true., .false., .true., .false., .false.]
      call tc%assert_eq(-1, select_int(mask, ints), "select_int 2")

      ! no hits
      ints = [     -1,     500,       3,       4,       5]
      mask = [.false., .false., .false., .false., .false.]
      call tc%assert_eq(0, select_int(mask, ints), "select_int 3")
    end block

    call tc%finish()
  end subroutine

end module
