module test_vector_m
  use test_case_m
  use vector_m
  implicit none

  private
  public :: test_vector

  type(test_case) :: tc

contains

  subroutine test_vector
    type(vector_int) :: vec
    integer, allocatable :: arr_exp(:), arr(:)

    print "(1A)", "test_vector"
    call tc%init("vector")

    ! test1: set/get array
    arr_exp = [2, 4, 5, 2, 1]
    call vec%init(size(arr_exp))
    call vec%set(arr_exp)
    arr = vec%get()
    call tc%assert_eq(arr_exp, arr, "get/set array")

    ! test2: first/last
    call tc%assert_eq(arr_exp(1), vec%first(), "first")
    call tc%assert_eq(arr_exp(5), vec%last(), "last")

    ! push arrays
    arr = [2, 4, 5, 2, 1]
    call vec%reset
    call vec%push(arr)
    call vec%push(arr)

    deallocate (arr_exp)
    allocate (arr_exp(10))
    arr_exp(1:5)  = arr
    arr_exp(6:10) = arr
    arr = vec%get()
    call tc%assert_eq(arr_exp, arr, "push elems")

    call tc%finish
  end subroutine

end module
