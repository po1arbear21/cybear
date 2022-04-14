m4_include(../../util/macro.f90.inc)

module test_grid_table_m

  use grid_m,       only: IDX_EDGE, IDX_CELL, IDX_FACE, IDX_VERTEX
  use grid_table_m, only: grid_table
  use grid1D_m,     only: grid1D
  use math_m,       only: logspace
  use test_case_m,  only: test_case
  use util_m,       only: int2str

  implicit none

  private
  public test_grid_table

contains

  subroutine test_grid_table()
    type(test_case)      :: tc
    type(grid1D), target :: g
    type(grid_table)     :: g_cell

    call tc%init("grid_table")

    ! init grid table by grid
    ! choose all cells except one
    block
      integer            :: ncell, ic, idx_bnd(1)
      integer, parameter :: nx=101, ic0(1)=[3]

      call g%init("x", logspace(1.0, 10.0, nx))
      call g_cell%init('cells', g, IDX_CELL, 0)

      ! enable each cell
      call g%get_idx_bnd(IDX_CELL, 0, idx_bnd)
      ncell = idx_bnd(1)
      call g_cell%flags%set([(.true., ic=1,ncell)])

      ! disable one cell to make testing more interesting
      call g_cell%flags%set(ic0, .false.)

      ! finalize init
      call g_cell%init_final()
    end block

    ! test get_flat, get_idx
    block
      integer            :: idx(1), flat, i
      integer, parameter ::  idx_exp(6) = [1, 2, 3, 4, 5, 6]
      integer, parameter :: flat_exp(6) = [1, 2, 0, 3, 4, 5]

      do i = 1, size(flat_exp)
        flat = g_cell%get_flat(idx_exp(i:i))
        call tc%assert_eq(flat_exp(i), flat, "grid_table: get flat. i:"//int2str(i))

        if (flat_exp(i) == 0) cycle
        idx = g_cell%get_idx(flat_exp(i))
        call tc%assert_eq(idx_exp(i:i), idx, "grid_table: get idx. i:"//int2str(i))
      end do
    end block

    call tc%finish()
  end subroutine
end module
