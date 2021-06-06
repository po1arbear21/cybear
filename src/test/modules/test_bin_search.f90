module test_bin_search_m

  use bin_search_m
  use test_case_m, only: test_case

  implicit none

  private
  public test_bin_search

contains

  subroutine test_bin_search()
    integer         :: a(7)
    type(test_case) :: tc

    call tc%init("bin_search")

    a = [1, 4, 6, 10, 20, 21, 22]

    ! BS_NEAR
    call tc%assert_eq(2, bin_search(a,  4, mode = BS_NEAR), "BS_NEAR 1")
    call tc%assert_eq(2, bin_search(a,  3, mode = BS_NEAR), "BS_NEAR 2")
    call tc%assert_eq(1, bin_search(a,  0, mode = BS_NEAR), "BS_NEAR 3")
    call tc%assert_eq(7, bin_search(a, 25, mode = BS_NEAR), "BS_NEAR 4")

    ! BS_LESS
    call tc%assert_eq(4, bin_search(a, 10, mode = BS_LESS), "BS_LESS 1")
    call tc%assert_eq(1, bin_search(a,  3, mode = BS_LESS), "BS_LESS 2")
    call tc%assert_eq(7, bin_search(a, 25, mode = BS_LESS), "BS_LESS 3")
    call tc%assert_eq(1, bin_search(a,  0, mode = BS_LESS), "BS_LESS 4")

    ! BS_GREAT
    call tc%assert_eq(6, bin_search(a, 21, mode = BS_GREAT), "BS_GREAT 1")
    call tc%assert_eq(4, bin_search(a,  7, mode = BS_GREAT), "BS_GREAT 2")
    call tc%assert_eq(1, bin_search(a,  0, mode = BS_GREAT), "BS_GREAT 3")
    call tc%assert_eq(7, bin_search(a, 25, mode = BS_GREAT), "BS_GREAT 4")

    call tc%finish()
  end subroutine

end module
