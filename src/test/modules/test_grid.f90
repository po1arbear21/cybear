#include "../../util/macro.f90.inc"

module test_grid_m
  use test_case_m
  use equation_m
  use res_equation_m
  use jacobian_chain_m
  use stencil_m
  use grid_m
  use grid0D_m
  use grid1D_m
  use triang_grid_m, only: triang_grid
  use tensor_grid_m
  use variable_m
  use vselector_m
  use jacobian_matrix_m
  use esystem_m
  use math_m, only: logspace
  use util_m, only: int2str
  use qsort_m, only: qsort

  implicit none

contains

  subroutine test_grid()
    type(test_case) :: tc

    print "(A)", "test_grid"
    call tc%init("grid")

    call test_grid1D()
    call test_triang_grid()

    call tc%finish()

  contains

    subroutine test_grid1D()
      integer, parameter :: nx=101
      type(grid1D)       :: g
      real, allocatable  :: x(:)

      allocate (x(nx))                ! remove gfortran warning
      x = logspace(1.0, 10.0, nx)     ! logspacing -> non-equidistant!

      call g%init(x)

      ! testing attributes
      call tc%assert_eq( 1,  g%dim,      "grid1D: dim")
      call tc%assert_eq( 1,  g%idx_dim,  "grid1D: idx_dim")
      call tc%assert_eq([1], g%face_dim, "grid1D: face_dim")
      call tc%assert_eq( 2,  g%cell_dim, "grid1D: cell_dim")
      call tc%assert_eq( x,  g%x, 0.0,   "grid1D: x")

      ! testing get_idx_bnd
      block
        integer :: bnd(1)

        call g%get_idx_bnd(IDX_VERTEX, 0, bnd)
        call tc%assert_eq(nx, bnd(1), "grid1D: idx_bnd: IDX_VERTEX")

        call g%get_idx_bnd(IDX_EDGE, 1, bnd)
        call tc%assert_eq(nx-1, bnd(1), "grid1D: idx_bnd: IDX_EDGE")

        call g%get_idx_bnd(IDX_FACE, 1, bnd)
        call tc%assert_eq(nx, bnd(1), "grid1D: idx_bnd: IDX_FACE")

        call g%get_idx_bnd(IDX_CELL, 0, bnd)
        call tc%assert_eq(nx-1, bnd(1), "grid1D: idx_bnd: IDX_CELL")
      end block

      ! get_vertex
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx]
        integer            :: i, i0
        real               :: p(1)

        do i = 1, size(i_arr)
          i0 = i_arr(i)
          call g%get_vertex([ i0], p)
          call tc%assert_eq(x(i0), p(1), 0.0, "grid1D: get_vertex "//int2str(i0))
        end do
      end block

      ! get_edge
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx-1]
        integer            :: i, i0
        real               :: p(1,2)

        do i = 1, size(i_arr)
          i0 = i_arr(i)
          call g%get_edge([i0], 1, p)
          call tc%assert_eq(x(i0:i0+1), p(1,:), 0.0, "grid1D: get_edge "//int2str(i0))
        end do
      end block

      ! get_face
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx]
        integer            :: i, i0
        real               :: p(1,1)

        do i = 1, size(i_arr)
          i0 = i_arr(i)
          call g%get_face([i0], 1, p)
          call tc%assert_eq(x(i0), p(1,1), 0.0, "grid1D: get_face "//int2str(i0))
        end do
      end block

      ! get_cell
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx-1]
        integer            :: i, i0
        real               :: p(1,2)

        do i = 1, size(i_arr)
          i0 = i_arr(i)
          call g%get_cell([i0], p)
          call tc%assert_eq(x(i0:i0+1), p(1,:), 0.0, "grid1D: get_cell "//int2str(i0))
        end do
      end block

      ! get_vol
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx-1]
        integer            :: i, i0
        real               :: vol

        do i = 1, size(i_arr)
          i0  = i_arr(i)
          vol = g%get_vol([i0])
          call tc%assert_eq(x(i0+1)-x(i0), vol, 5e-16, "grid1D: get_vol "//int2str(i0))
        end do
      end block

      ! get_max_neighb
      call tc%assert_eq(3, g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "grid1D: get_max_neighb V V")
      call tc%assert_eq(2, g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   1), "grid1D: get_max_neighb V E")
      call tc%assert_eq(1, g%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   1), "grid1D: get_max_neighb V F")
      call tc%assert_eq(2, g%get_max_neighb(IDX_VERTEX, 0, IDX_CELL,   0), "grid1D: get_max_neighb V C")
      call tc%assert_eq(3, g%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   1), "grid1D: get_max_neighb E E")
      call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   1), "grid1D: get_max_neighb E F")
      call tc%assert_eq(1, g%get_max_neighb(IDX_EDGE,   1, IDX_CELL,   0), "grid1D: get_max_neighb E C")
      call tc%assert_eq(3, g%get_max_neighb(IDX_FACE,   1, IDX_FACE,   1), "grid1D: get_max_neighb F F")
      call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_CELL,   0), "grid1D: get_max_neighb F C")
      call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_CELL,   0), "grid1D: get_max_neighb C C")

      ! get_neighb
      block
        integer :: nidx2, idx2(1,3)

        associate (V => IDX_VERTEX, E => IDX_EDGE, F => IDX_FACE, C => IDX_CELL)
          ! V(ix=1) -> V_neigh=?
          call g%get_neighb(V, 0, V, 0, [1], idx2(:,:g%get_max_neighb(V, 0, V, 0)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb V V: n")
          call tc%assert_eq([1, 2], idx2(1,:nidx2), "grid1D: get_neighb V V: v")

          ! V(ix=1) -> E_neigh=?
          call g%get_neighb(V, 0, E, 1, [1], idx2(:,:g%get_max_neighb(V, 0, E, 1)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb V E: n")
          call tc%assert_eq([1], idx2(1,:nidx2), "grid1D: get_neighb V E: v")

          ! V(ix=1) -> F_neigh=?
          call g%get_neighb(V, 0, F, 1, [1], idx2(:,:g%get_max_neighb(V, 0, F, 1)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb V F: n")
          call tc%assert_eq([1], idx2(1,:nidx2), "grid1D: get_neighb V F: v")

          ! V(ix=1) -> C_neigh=?
          call g%get_neighb(V, 0, C, 0, [1], idx2(:,:g%get_max_neighb(V, 0, C, 0)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb V C: n")
          call tc%assert_eq([1], idx2(1,:nidx2), "grid1D: get_neighb V C: v")

          ! E(ix=1) -> E_neigh=?
          call g%get_neighb(E, 1, E, 1, [1], idx2(:,:g%get_max_neighb(E, 1, E, 1)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb E E: n")
          call tc%assert_eq([1, 2], idx2(1,:nidx2), "grid1D: get_neighb E E: v")

          ! E(ix=1) -> F_neigh=?
          call g%get_neighb(E, 1, F, 1, [1], idx2(:,:g%get_max_neighb(E, 1, F, 1)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb E F: n")
          call tc%assert_eq([1, 2], idx2(1,:nidx2), "grid1D: get_neighb E F: v")

          ! E(ix=1) -> C_neigh=?
          call g%get_neighb(E, 1, C, 0, [1], idx2(:,:g%get_max_neighb(E, 1, C, 0)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb E C: n")
          call tc%assert_eq([1], idx2(1,:nidx2), "grid1D: get_neighb E C: v")

          ! F(ix=1) -> F_neigh=?
          call g%get_neighb(F, 1, F, 1, [1], idx2(:,:g%get_max_neighb(F, 1, F, 1)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb E F: n")
          call tc%assert_eq([1, 2], idx2(1,:nidx2), "grid1D: get_neighb E F: v")

          ! F(ix=1) -> C_neigh=?
          call g%get_neighb(F, 1, C, 0, [1], idx2(:,:g%get_max_neighb(F, 1, C, 0)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb E C: n")
          call tc%assert_eq([1], idx2(1,:nidx2), "grid1D: get_neighb E C: v")

          ! C(ix=1) -> C_neigh=?
          call g%get_neighb(C, 0, C, 0, [1], idx2(:,:g%get_max_neighb(C, 0, C, 0)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb E C: n")
          call tc%assert_eq([1, 2], idx2(1,:nidx2), "grid1D: get_neighb E C: v")

          ! V(ix=5) -> V_neigh=?
          call g%get_neighb(V, 0, V, 0, [5], idx2(:,:g%get_max_neighb(V, 0, V, 0)), nidx2)
          call tc%assert_eq(3,         nidx2,          "grid1D: get_neighb V V: n")
          call tc%assert_eq([4, 5, 6], idx2(1,:nidx2), "grid1D: get_neighb V V: v")

          ! V(ix=5) -> E_neigh=?
          call g%get_neighb(V, 0, E, 1, [5], idx2(:,:g%get_max_neighb(V, 0, E, 1)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb V E: n")
          call tc%assert_eq([4, 5], idx2(1,:nidx2), "grid1D: get_neighb V E: v")

          ! V(ix=5) -> F_neigh=?
          call g%get_neighb(V, 0, F, 1, [5], idx2(:,:g%get_max_neighb(V, 0, F, 1)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb V F: n")
          call tc%assert_eq([5], idx2(1,:nidx2), "grid1D: get_neighb V F: v")

          ! V(ix=5) -> C_neigh=?
          call g%get_neighb(V, 0, C, 0, [5], idx2(:,:g%get_max_neighb(V, 0, C, 0)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb V C: n")
          call tc%assert_eq([4, 5], idx2(1,:nidx2), "grid1D: get_neighb V C: v")

          ! E(ix=5) -> E_neigh=?
          call g%get_neighb(E, 1, E, 1, [5], idx2(:,:g%get_max_neighb(E, 1, E, 1)), nidx2)
          call tc%assert_eq(3,         nidx2,          "grid1D: get_neighb E E: n")
          call tc%assert_eq([4, 5, 6], idx2(1,:nidx2), "grid1D: get_neighb E E: v")

          ! E(ix=5) -> F_neigh=?
          call g%get_neighb(E, 1, F, 1, [5], idx2(:,:g%get_max_neighb(E, 1, F, 1)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb E F: n")
          call tc%assert_eq([5, 6], idx2(1,:nidx2), "grid1D: get_neighb E F: v")

          ! E(ix=5) -> C_neigh=?
          call g%get_neighb(E, 1, C, 0, [5], idx2(:,:g%get_max_neighb(E, 1, C, 0)), nidx2)
          call tc%assert_eq(1,   nidx2,          "grid1D: get_neighb E C: n")
          call tc%assert_eq([5], idx2(1,:nidx2), "grid1D: get_neighb E C: v")

          ! F(ix=5) -> F_neigh=?
          call g%get_neighb(F, 1, F, 1, [5], idx2(:,:g%get_max_neighb(F, 1, F, 1)), nidx2)
          call tc%assert_eq(3,         nidx2,          "grid1D: get_neighb F F: n")
          call tc%assert_eq([4, 5, 6], idx2(1,:nidx2), "grid1D: get_neighb F F: v")

          ! F(ix=5) -> C_neigh=?
          call g%get_neighb(F, 1, C, 0, [5], idx2(:,:g%get_max_neighb(F, 1, C, 0)), nidx2)
          call tc%assert_eq(2,      nidx2,          "grid1D: get_neighb F C: n")
          call tc%assert_eq([4, 5], idx2(1,:nidx2), "grid1D: get_neighb F C: v")

          ! C(ix=5) -> C_neigh=?
          call g%get_neighb(C, 0, C, 0, [5], idx2(:,:g%get_max_neighb(C, 0, C, 0)), nidx2)
          call tc%assert_eq(3,         nidx2,          "grid1D: get_neighb C C: n")
          call tc%assert_eq([4, 5, 6], idx2(1,:nidx2), "grid1D: get_neighb C C: v")
        end associate
      end block
    end subroutine

    subroutine test_triang_grid()
      !! testing triang_grid_m
      integer, parameter :: nvert=7, ncell=6
      type(triang_grid)  :: g
      real               :: vert(2,nvert)
      integer            :: icell(3,ncell)

      ! init grid
      vert(:,1) = [0,  0]
      vert(:,2) = [1,  1]
      vert(:,3) = [2,  0]
      vert(:,4) = [0, -2]
      vert(:,5) = [1,  0]
      vert(:,6) = [0,  1]
      vert(:,7) = [0,  2]

      icell(:,1) = [1, 4, 5]
      icell(:,2) = [5, 3, 2]
      icell(:,3) = [5, 4, 3]
      icell(:,4) = [1, 2, 6]
      icell(:,5) = [1, 5, 2]
      icell(:,6) = [6, 2, 7]

      call g%init(vert, icell)

      ! testing attributes
      call tc%assert_eq( 2,  g%dim,      "triang_grid: dim")
      call tc%assert_eq( 1,  g%idx_dim,  "triang_grid: idx_dim")
      call tc%assert_eq([2], g%face_dim, "triang_grid: face_dim")
      call tc%assert_eq( 3,  g%cell_dim, "triang_grid: cell_dim")

      ! ! testing get_idx_bnd
      block
        integer :: bnd(1)

        call g%get_idx_bnd(IDX_VERTEX, 0, bnd)
        call tc%assert_eq(nvert, bnd(1), "triang_grid: idx_bnd: IDX_VERTEX")

        call g%get_idx_bnd(IDX_EDGE, 0, bnd)
        call tc%assert_eq(12, bnd(1), "triang_grid: idx_bnd: IDX_EDGE")

        call g%get_idx_bnd(IDX_FACE, 0, bnd)
        call tc%assert_eq(12, bnd(1), "triang_grid: idx_bnd: IDX_FACE")

        call g%get_idx_bnd(IDX_CELL, 0, bnd)
        call tc%assert_eq(ncell, bnd(1), "triang_grid: idx_bnd: IDX_CELL")
      end block

      ! get_vertex
      block
        integer :: iv
        real    :: p(2)

        do iv = 1, nvert
          call g%get_vertex([iv], p)
          call tc%assert_eq(vert(:,iv), p, 0.0, "triang_grid: get_vertex "//int2str(iv))
        end do
      end block

      ! get_cell
      block
        integer            :: ic, iv
        real               :: p(2,3)

        do ic = 1, ncell
          call g%get_cell([ic], p)
          do iv = 1, 3
            call tc%assert_eq(vert(:,icell(iv,ic)), p(:,iv), 0.0, "triang_grid: get_cell "//int2str(ic)//' '//int2str(iv))
          end do
        end do
      end block

      ! get_vol
      block
        integer         :: ic
        real, parameter :: vol(6) = [1.0, 0.5, 1.0, 0.5, 0.5, 0.5]

        do ic = 1, ncell
          call tc%assert_eq(vol(ic), g%get_vol([ic]), 5e-16, "triang_grid: get_vol "//int2str(ic))
        end do
      end block

      ! get_max_neighb
      call tc%assert_eq(6, g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "triang_grid: get_max_neighb V V")
      call tc%assert_eq(5, g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   0), "triang_grid: get_max_neighb V E")
      call tc%assert_eq(5, g%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   0), "triang_grid: get_max_neighb V F")
      call tc%assert_eq(4, g%get_max_neighb(IDX_VERTEX, 0, IDX_CELL,   0), "triang_grid: get_max_neighb V C")

      call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   0, IDX_VERTEX, 0), "triang_grid: get_max_neighb E V")
      call tc%assert_eq(8, g%get_max_neighb(IDX_EDGE,   0, IDX_EDGE,   0), "triang_grid: get_max_neighb E E")
      call tc%assert_eq(8, g%get_max_neighb(IDX_EDGE,   0, IDX_FACE,   0), "triang_grid: get_max_neighb E F")
      call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   0, IDX_CELL,   0), "triang_grid: get_max_neighb E C")

      call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   0, IDX_VERTEX, 0), "triang_grid: get_max_neighb F V")
      call tc%assert_eq(8, g%get_max_neighb(IDX_FACE,   0, IDX_EDGE,   0), "triang_grid: get_max_neighb F E")
      call tc%assert_eq(8, g%get_max_neighb(IDX_FACE,   0, IDX_FACE,   0), "triang_grid: get_max_neighb F F")
      call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   0, IDX_CELL,   0), "triang_grid: get_max_neighb F C")

      call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_VERTEX, 0), "triang_grid: get_max_neighb C V")
      call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   0), "triang_grid: get_max_neighb C E")
      call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_FACE,   0), "triang_grid: get_max_neighb C F")
      call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_CELL,   0), "triang_grid: get_max_neighb C C")

      ! get_neighb
      block
        integer              :: idx1(1), nidx2
        integer, allocatable :: idx2(:,:), exp_idx2(:,:)

        ! V V
        allocate (idx2(1,g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0)))

        idx1     = [1]
        exp_idx2 = reshape([2, 4, 5, 6], [1, 4])
        call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V V  idx2 "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 3, 4], [1, 4])
        call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V V  idx2 "//int2str(idx1(1)))

        idx1     = [7]
        exp_idx2 = reshape([2, 6], [1, 2])
        call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V V  idx2 "//int2str(idx1(1)))

        ! V E
        deallocate (idx2)
        allocate (idx2(1,g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE, 0)))

        idx1 = [1]
        call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call tc%assert_eq(4, nidx2, "triang_grid: get_neighb V E nidx2 "//int2str(idx1(1)))

        idx1 = [3]
        call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call tc%assert_eq(3, nidx2, "triang_grid: get_neighb V E nidx2 "//int2str(idx1(1)))

        idx1 = [2]
        call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call tc%assert_eq(5, nidx2, "triang_grid: get_neighb V E nidx2 "//int2str(idx1(1)))

        ! V C
        deallocate (idx2)
        allocate (idx2(1,g%get_max_neighb(IDX_VERTEX, 0, IDX_CELL, 0)))

        idx1     = [1]
        exp_idx2 = reshape([1, 4, 5], [1, 3])
        call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V C  idx2 "//int2str(idx1(1)))

        idx1     = [4]
        exp_idx2 = reshape([1, 3], [1, 2])
        call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V C  idx2 "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 3, 5], [1, 4])
        call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb V C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb V C  idx2 "//int2str(idx1(1)))

        ! C V
        deallocate (idx2)
        allocate (idx2(1,g%get_max_neighb(IDX_CELL, 0, IDX_VERTEX, 0)))

        idx1     = [1]
        exp_idx2 = reshape([1, 4, 5], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C V  idx2 "//int2str(idx1(1)))

        idx1     = [4]
        exp_idx2 = reshape([1, 2, 6], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C V  idx2 "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 5, 2], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, idx2, nidx2)
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C V nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C V  idx2 "//int2str(idx1(1)))

        ! C E
        deallocate (idx2)
        allocate (idx2(1,g%get_max_neighb(IDX_CELL, 0, IDX_EDGE, 0)))

        idx1     = [1]
        exp_idx2 = reshape([1, 7, 9], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C E nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C E  idx2 "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([2, 3, 9], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C E nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C E  idx2 "//int2str(idx1(1)))

        idx1     = [6]
        exp_idx2 = reshape([4, 5, 12], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C E nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C E  idx2 "//int2str(idx1(1)))

        ! C C
        deallocate (idx2)
        allocate (idx2(1,g%get_max_neighb(IDX_CELL, 0, IDX_CELL, 0)))

        idx1     = [1]
        exp_idx2 = reshape([3, 5], [1, 2])
        call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C C  idx2 "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 4], [1, 3])
        call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C C  idx2 "//int2str(idx1(1)))

        idx1     = [6]
        exp_idx2 = reshape([4], [1, 1])
        call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, idx2, nidx2)
        call qsort(idx2(1,:nidx2))
        call tc%assert_eq(size(exp_idx2, dim=2), nidx2,  "triang_grid: get_neighb C C nidx2 "//int2str(idx1(1)))
        call tc%assert_eq(exp_idx2(1,:), idx2(1,:nidx2), "triang_grid: get_neighb C C  idx2 "//int2str(idx1(1)))
      end block
    end subroutine

  end subroutine
end module
