module test_util_m

  use string_m
  use test_case_m, only: test_case
  use util_m,      only: c2fstring, cstrlen, f2cstring, int2str, real2str, select_int

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
      call tc%assert_eq(exp, val, "int2str 1")

      i     = -4567
      cval  = int2str(i)
      val%s = cval
      exp%s = "-4567"
      call tc%assert_eq(exp, val, "int2str 2")

      ! format length equal actual length
      i     = 1234
      cval  = int2str(i, fmt="(I4)")
      val%s = cval
      exp%s = "1234"
      call tc%assert_eq(exp, val, "int2str 3")

      i     = -1234
      cval  = int2str(i, fmt="(I5)")
      val%s = cval
      exp%s = "-1234"
      call tc%assert_eq(exp, val, "int2str 4")

      ! format length longer than actual length
      i     = 56
      cval  = int2str(i, fmt="(I6.6)")
      val%s = cval
      exp%s = "000056"
      call tc%assert_eq(exp, val, "int2str 5")

      i     = -56
      cval  = int2str(i, fmt="(I6.5)")
      val%s = cval
      exp%s = "-00056"
      call tc%assert_eq(exp, val, "int2str 6")
    end block

    ! real2str
    block
      character(:), allocatable :: cval
      real                      :: r
      type(string)              :: exp, val

      r     = 123.27858499357
      cval  = real2str(r)
      val%s = cval
      exp%s = "1.2327858499357001E+002"
      call tc%assert_eq(exp, val, "real2str 1")

      r     = -4567.45655
      cval  = real2str(r)
      val%s = cval
      exp%s = "-4.5674565499999999E+003"
      call tc%assert_eq(exp, val, "real2str 2")

      r     = 3.141592653589793238462643
      cval  = real2str(r, fmt="(F10.2)")
      val%s = cval
      exp%s = "3.14"
      call tc%assert_eq(exp, val, "real2str 3")
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
