module test_util_m

  use test_case_m
  use util_m
  use string_m

  implicit none

  private
  public test_util

contains

  subroutine test_util()
    type(test_case) :: tc

    call tc%init("util")

    ! int2str
    block
      integer                   :: i
      type(string)              :: exp, val
      character(:), allocatable :: cval

      allocate (character(0) :: cval)   ! removes gfortran warning

      ! standard call, no min_len specified
      i     = 123
      cval  = int2str(i)
      val%s = cval
      exp%s = '123'
      call tc%assert_eq(exp, val, "int2str")

      i     = -4567
      cval  = int2str(i)
      val%s = cval
      exp%s = '-4567'
      call tc%assert_eq(exp, val, "int2str")

      ! min_len smaller than actual length
      i     = 89012
      cval  = int2str(i, min_len=2)
      val%s = cval
      exp%s = '89012'
      call tc%assert_eq(exp, val, "int2str")

      i     = -123
      cval  = int2str(i, min_len=2)
      val%s = cval
      exp%s = '-123'
      call tc%assert_eq(exp, val, "int2str")

      ! min_len equal actual length
      i     = 1234
      cval  = int2str(i, min_len=4)
      val%s = cval
      exp%s = '1234'
      call tc%assert_eq(exp, val, "int2str")

      i     = -1234
      cval  = int2str(i, min_len=4)
      val%s = cval
      exp%s = '-1234'
      call tc%assert_eq(exp, val, "int2str")

      ! min_len larger than actual length
      i     = 56
      cval  = int2str(i, min_len=6)
      val%s = cval
      exp%s = '000056'
      call tc%assert_eq(exp, val, "int2str")

      i     = -56
      cval  = int2str(i, min_len=6)
      val%s = cval
      exp%s = '-00056'
      call tc%assert_eq(exp, val, "int2str")
    end block

    call tc%finish
  end subroutine

end module
