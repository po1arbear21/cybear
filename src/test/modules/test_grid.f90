#include "../../util/macro.f90.inc"

module test_grid_m
  use test_case_m
  use grid_m
  use grid1D_m
  use tensor_grid_m
  use math_m
  implicit none

contains

  subroutine test_grid()
    type(test_case) :: tc

    print "(A)", "test_grid"
    call tc%init("grid")

    block
      type(grid1D)      :: g
      real              :: p(1)
      real, allocatable :: x(:)

      x = linspace(0.0, 10.0, 101)

      call g%init(x)

      call g%get_vertex([  1], p)
      call tc%assert_eq(x(1), p(1), 0.0, "grid1D: get_vertex 1")
      call g%get_vertex([  5], p)
      call tc%assert_eq(x(5), p(1), 0.0, "grid1D: get_vertex 2")
      call g%get_vertex([ 43], p)
      call tc%assert_eq(x(43), p(1), 0.0, "grid1D: get_vertex 3")
      call g%get_vertex([101], p)
      call tc%assert_eq(x(101), p(1), 0.0, "grid1D: get_vertex 4")
    end block

    call tc%finish()
  end subroutine

end module