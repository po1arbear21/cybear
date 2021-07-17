module test_grid_m

  use grid_m,        only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use grid0D_m,      only: grid0D
  use grid1D_m,      only: grid1D
  use grid_data_m
  use math_m,        only: linspace, logspace
  use qsort_m,       only: qsort
  use test_case_m,   only: test_case
  use tensor_grid_m, only: tensor_grid
  use triang_grid_m, only: triang_grid
  use util_m,        only: int2str
  use vector_m,      only: vector_int

  ! no tests implemented but used s.t. modules get compiled at all
  use sum_grid_m

  implicit none

  private
  public test_grid

contains

  subroutine test_grid()
    type(test_case) :: tc

    call tc%init("grid")

    call test_grid_data()
    call test_grid0D()
    call test_grid1D()
    call test_triang_grid()
    call test_tensor_grid()

    call tc%finish()

  contains

    subroutine test_grid_data()
      integer, parameter :: nx=10, ny=10
      real, allocatable  :: x(:), y(:)
      type(grid1D)       :: gx, gy
      type(grid_data2_real) :: par
      type(tensor_grid)     :: g

      allocate (x(nx), y(ny))
      x = linspace(1.0, nx, nx)
      y = linspace(1.0, ny, ny)

      ! initialize 1D grids
      call gx%init(x)
      call gy%init(y)

      ! initialize tensor grid and grid_data
      call g%init([gx%get_ptr(), gy%get_ptr()])
      call par%init(g, IDX_CELL, 0)

      ! check attributes
      call tc%assert_eq(3,   par%idx_type,         "grid_data: idx_type")
      call tc%assert_eq(0.0, par%d0, 0.0,          "grid_data: d0")
      call tc%assert_eq(0.0, par%get([1, 1]), 0.0, "grid_data: get default")

      ! check set and get
      call par%set([2, 3], 10.0)
      call tc%assert_eq(10.0, par%get([2, 3]), 0.0, "grid_data: set")

      ! check update
      call par%update([2, 3], -5.0)
      call par%update([1, 1],  5.0)
      call tc%assert_eq(-5.0, par%get([2, 3]), 0.0, "grid_data: update idx")
      call par%update(2.0)
      call tc%assert_eq( 2.0, par%get([2, 3]), 0.0, "grid_data: update all 1")
      call tc%assert_eq( 2.0, par%get([5, 8]), 0.0, "grid_data: update all 2")
    end subroutine

    subroutine test_grid0D()
      type(grid0D) :: g
      integer      :: bnd(0), idx1(0), idx2(0)
      logical      :: status
      real         :: p(0)

      call g%init()

      ! testing attributes
      call tc%assert_eq(0,      g%dim,       "grid0D: dim")
      call tc%assert_eq(0,      g%idx_dim,   "grid0D: idx_dim")
      call tc%assert_eq(0, size(g%face_dim), "grid0D: face_dim")
      call tc%assert_eq(0,      g%cell_dim,  "grid0D: cell_dim")

      ! get_idx_bnd
      call g%get_idx_bnd(IDX_VERTEX, 0, bnd)

      ! get_vertex
      call g%get_vertex(idx1, p)

      ! get_max_neighb
      call tc%assert_eq(0, g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "grid0D: get_max_neighb")

      ! get_neighb
      call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, 1, idx2, status)
      call tc%assert(.not. status, "grid0D: get_neighb")
    end subroutine

    subroutine test_grid1D()
      integer, parameter :: nx=101
      type(grid1D)       :: g
      real, allocatable  :: x(:)

      allocate (x(nx))                ! remove gfortran warning
      x = logspace(1.0, 10.0, nx)     ! logspacing -> non-equidistant!

      call g%init(x)

      ! check attributes
      block
        call tc%assert_eq( 1,  g%dim,      "grid1D: dim")
        call tc%assert_eq( 1,  g%idx_dim,  "grid1D: idx_dim")
        call tc%assert_eq([1], g%face_dim, "grid1D: face_dim")
        call tc%assert_eq( 2,  g%cell_dim, "grid1D: cell_dim")
        call tc%assert_eq( x,  g%x, 0.0,   "grid1D: x")
      end block

      ! get_idx_bnd
      block
        integer :: bnd_n(1), bnd

        call g%get_idx_bnd(IDX_VERTEX, 0, bnd_n)
        call tc%assert_eq([nx], bnd_n, "grid1D: idx_bnd_n: IDX_VERTEX")
        call g%get_idx_bnd(IDX_VERTEX, 0, bnd)
        call tc%assert_eq( nx,  bnd,   "grid1D: idx_bnd_1: IDX_VERTEX")

        call g%get_idx_bnd(IDX_EDGE, 1, bnd_n)
        call tc%assert_eq([nx-1], bnd_n, "grid1D: idx_bnd_n: IDX_EDGE")
        call g%get_idx_bnd(IDX_EDGE, 1, bnd)
        call tc%assert_eq( nx-1,  bnd,   "grid1D: idx_bnd_1: IDX_EDGE")

        call g%get_idx_bnd(IDX_FACE, 1, bnd_n)
        call tc%assert_eq([nx], bnd_n, "grid1D: idx_bnd_n: IDX_FACE")
        call g%get_idx_bnd(IDX_FACE, 1, bnd)
        call tc%assert_eq( nx,  bnd,   "grid1D: idx_bnd_1: IDX_FACE")

        call g%get_idx_bnd(IDX_CELL, 0, bnd_n)
        call tc%assert_eq([nx-1], bnd_n, "grid1D: idx_bnd_n: IDX_CELL")
        call g%get_idx_bnd(IDX_CELL, 0, bnd)
        call tc%assert_eq( nx-1,  bnd,   "grid1D: idx_bnd_1: IDX_CELL")
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

      ! get_len
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx-1]
        integer            :: i, i0
        real               :: len

        do i = 1, size(i_arr)
          i0  = i_arr(i)
          len = g%get_len([i0], 1)
          call tc%assert_eq(x(i0+1)-x(i0), len, 5e-16, "grid1D: get_len "//int2str(i0))
        end do
      end block

      ! get_surf
      block
        integer, parameter :: i_arr(4) = [1, 5, 43, nx]
        integer            :: i, i0
        real               :: surf

        do i = 1, size(i_arr)
          i0  = i_arr(i)
          surf = g%get_surf([i0], 1)
          call tc%assert_eq(1.0, surf, 0.0, "grid1D: get_surf "//int2str(i0))
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
      block
        call tc%assert_eq(2, g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "grid1D: get_max_neighb V V")
        call tc%assert_eq(2, g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   1), "grid1D: get_max_neighb V E")
        call tc%assert_eq(1, g%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   1), "grid1D: get_max_neighb V F")
        call tc%assert_eq(2, g%get_max_neighb(IDX_VERTEX, 0, IDX_CELL,   0), "grid1D: get_max_neighb V C")
        call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_VERTEX, 0), "grid1D: get_max_neighb E V")
        call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   1), "grid1D: get_max_neighb E E")
        call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   1), "grid1D: get_max_neighb E F")
        call tc%assert_eq(1, g%get_max_neighb(IDX_EDGE,   1, IDX_CELL,   0), "grid1D: get_max_neighb E C")
        call tc%assert_eq(1, g%get_max_neighb(IDX_FACE,   1, IDX_VERTEX, 0), "grid1D: get_max_neighb F V")
        call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_EDGE,   1), "grid1D: get_max_neighb F E")
        call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_FACE,   1), "grid1D: get_max_neighb F F")
        call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_CELL,   0), "grid1D: get_max_neighb F C")
        call tc%assert_eq(2, g%get_max_neighb(IDX_CELL,   0, IDX_VERTEX, 0), "grid1D: get_max_neighb C V")
        call tc%assert_eq(1, g%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   1), "grid1D: get_max_neighb C E")
        call tc%assert_eq(2, g%get_max_neighb(IDX_CELL,   0, IDX_FACE,   1), "grid1D: get_max_neighb C F")
        call tc%assert_eq(2, g%get_max_neighb(IDX_CELL,   0, IDX_CELL,   0), "grid1D: get_max_neighb C C")
      end block

      ! get_neighb
      block
        integer              :: idx2(1), j
        integer, allocatable :: exp_idx2(:)
        logical              :: status

        associate (V => IDX_VERTEX, E => IDX_EDGE, F => IDX_FACE, C => IDX_CELL)
          ! V(ix=1) -> V_neigh=?
          exp_idx2 = [2]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, V, 0, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V V: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V V: idx2"  )
          end do

          ! V(ix=1) -> E_neigh=?
          exp_idx2 = [1]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, E, 1, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V E: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V E: idx2"  )
          end do

          ! V(ix=1) -> F_neigh=?
          exp_idx2 = [1]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, F, 1, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V F: idx2"  )
          end do

          ! V(ix=1) -> C_neigh=?
          exp_idx2 = [1]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, C, 0, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V C: idx2"  )
          end do

          ! E(ix=1) -> E_neigh=?
          exp_idx2 = [2]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, E, 1, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E E: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E E: idx2"  )
          end do

          ! E(ix=1) -> F_neigh=?
          exp_idx2 = [1, 2]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, F, 1, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E F: idx2"  )
          end do

          ! E(ix=1) -> C_neigh=?
          exp_idx2 = [1]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, C, 0, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E C: idx2"  )
          end do

          ! F(ix=1) -> F_neigh=?
          exp_idx2 = [2]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(F, 1, F, 1, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb F F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb F F: idx2"  )
          end do

          ! F(ix=1) -> C_neigh=?
          exp_idx2 = [1]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(F, 1, C, 0, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb F C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb F C: idx2"  )
          end do

          ! C(ix=1) -> C_neigh=?
          exp_idx2 = [2]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(C, 0, C, 0, [1], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb C C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb C C: idx2"  )
          end do

          ! V(ix=5) -> V_neigh=?
          exp_idx2 = [4, 6]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, V, 0, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V V: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V V: idx2"  )
          end do

          ! V(ix=5) -> E_neigh=?
          exp_idx2 = [4, 5]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, E, 1, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V E: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V E: idx2"  )
          end do

          ! V(ix=5) -> F_neigh=?
          exp_idx2 = [5]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, F, 1, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V F: idx2"  )
          end do

          ! V(ix=5) -> C_neigh=?
          exp_idx2 = [4, 5]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(V, 0, C, 0, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb V C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb V C: idx2"  )
          end do

          ! E(ix=5) -> E_neigh=?
          exp_idx2 = [4, 6]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, E, 1, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E E: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E E: idx2"  )
          end do

          ! E(ix=5) -> F_neigh=?
          exp_idx2 = [5, 6]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, F, 1, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E F: idx2"  )
          end do

          ! E(ix=5) -> C_neigh=?
          exp_idx2 = [5]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(E, 1, C, 0, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb E C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb E C: idx2"  )
          end do

          ! F(ix=5) -> F_neigh=?
          exp_idx2 = [4, 6]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(F, 1, F, 1, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb F F: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb F F: idx2"  )
          end do

          ! F(ix=5) -> C_neigh=?
          exp_idx2 = [4, 5]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(F, 1, C, 0, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb F C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb F C: idx2"  )
          end do

          ! C(ix=5) -> C_neigh=?
          exp_idx2 = [4, 6]
          do j = 1, size(exp_idx2)+1
            call g%get_neighb(C, 0, C, 0, [5], j, idx2, status)
            call tc%assert((status .eqv. (j <= size(exp_idx2))), "grid1D: get_neighb C C: status")
            if (status) call tc%assert_eq(exp_idx2(j:j), idx2,   "grid1D: get_neighb C C: idx2"  )
          end do
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

      ! testing get_idx_bnd
      block
        integer :: bnd_n(1), bnd

        call g%get_idx_bnd(IDX_VERTEX, 0, bnd_n)
        call tc%assert_eq([nvert], bnd_n, "triang_grid: idx_bnd_n: IDX_VERTEX")
        call g%get_idx_bnd(IDX_VERTEX, 0, bnd)
        call tc%assert_eq( nvert,  bnd,   "triang_grid: idx_bnd_1: IDX_VERTEX")

        call g%get_idx_bnd(IDX_EDGE, 1, bnd_n)
        call tc%assert_eq([12], bnd_n, "triang_grid: idx_bnd_n: IDX_EDGE")
        call g%get_idx_bnd(IDX_EDGE, 1, bnd)
        call tc%assert_eq( 12,  bnd,   "triang_grid: idx_bnd_1: IDX_EDGE")

        call g%get_idx_bnd(IDX_FACE, 1, bnd_n)
        call tc%assert_eq([12], bnd_n, "triang_grid: idx_bnd_n: IDX_FACE")
        call g%get_idx_bnd(IDX_FACE, 1, bnd)
        call tc%assert_eq( 12,  bnd,   "triang_grid: idx_bnd_1: IDX_FACE")

        call g%get_idx_bnd(IDX_CELL, 0, bnd_n)
        call tc%assert_eq([ncell], bnd_n, "triang_grid: idx_bnd_n: IDX_CELL")
        call g%get_idx_bnd(IDX_CELL, 0, bnd)
        call tc%assert_eq( ncell,  bnd,   "triang_grid: idx_bnd_1: IDX_CELL")
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
      block
        call tc%assert_eq(5, g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "triang_grid: get_max_neighb V V")
        call tc%assert_eq(5, g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   1), "triang_grid: get_max_neighb V E")
        call tc%assert_eq(5, g%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   1), "triang_grid: get_max_neighb V F")
        call tc%assert_eq(4, g%get_max_neighb(IDX_VERTEX, 0, IDX_CELL,   0), "triang_grid: get_max_neighb V C")

        call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_VERTEX, 0), "triang_grid: get_max_neighb E V")
        call tc%assert_eq(7, g%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   1), "triang_grid: get_max_neighb E E")
        call tc%assert_eq(7, g%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   1), "triang_grid: get_max_neighb E F")
        call tc%assert_eq(2, g%get_max_neighb(IDX_EDGE,   1, IDX_CELL,   0), "triang_grid: get_max_neighb E C")

        call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_VERTEX, 0), "triang_grid: get_max_neighb F V")
        call tc%assert_eq(7, g%get_max_neighb(IDX_FACE,   1, IDX_EDGE,   1), "triang_grid: get_max_neighb F E")
        call tc%assert_eq(7, g%get_max_neighb(IDX_FACE,   1, IDX_FACE,   1), "triang_grid: get_max_neighb F F")
        call tc%assert_eq(2, g%get_max_neighb(IDX_FACE,   1, IDX_CELL,   0), "triang_grid: get_max_neighb F C")

        call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_VERTEX, 0), "triang_grid: get_max_neighb C V")
        call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   1), "triang_grid: get_max_neighb C E")
        call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_FACE,   1), "triang_grid: get_max_neighb C F")
        call tc%assert_eq(3, g%get_max_neighb(IDX_CELL,   0, IDX_CELL,   0), "triang_grid: get_max_neighb C C")
      end block

      ! get_neighb
      block
        integer              :: idx1(1), nidx2, idx2(1), j
        integer, allocatable :: idx2_arr(:,:), exp_idx2(:,:)
        logical :: status
        type(vector_int)     :: idx2_vec

        call idx2_vec%init(0, c=10)
        allocate (exp_idx2(0,0))

        ! V V
        idx1     = [1]
        exp_idx2 = reshape([2, 4, 5, 6], [1, 4])
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V V, idx1: "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 3, 4], [1, 4])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V V, idx1: "//int2str(idx1(1)))

        idx1     = [7]
        exp_idx2 = reshape([2, 6], [1, 2])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V V, idx1: "//int2str(idx1(1)))

        ! V E
        idx1  = [1]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(4, nidx2, "triang_grid: get_neighb V E, nidx2: "//int2str(idx1(1)))

        idx1  = [2]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(5, nidx2, "triang_grid: get_neighb V E, nidx2: "//int2str(idx1(1)))

        idx1  = [3]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(3, nidx2, "triang_grid: get_neighb V E, nidx2: "//int2str(idx1(1)))

        ! V C
        idx1     = [1]
        exp_idx2 = reshape([1, 4, 5], [1, 3])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V C, idx1: "//int2str(idx1(1)))

        idx1     = [4]
        exp_idx2 = reshape([1, 3], [1, 2])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V C, idx1: "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 3, 5], [1, 4])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb V C, idx1: "//int2str(idx1(1)))

        ! edge to XXX tests
        block
          integer :: ie(nvert,nvert), ie_tmp(1), iv1, iv2, k, iv2_tmp(1)

          ie = -1
          ! determine edge indices
          do iv1 = 1, nvert
            do iv2 = 1, nvert
              NEIGHB: do j = 1, 100
                call g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, [iv1], j, ie_tmp, status)
                if (.not. status) exit
                do k = 1, 2
                  call g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, ie_tmp, k, iv2_tmp, status)
                  if (iv2_tmp(1) == iv2) then
                    ie(iv1,iv2) = ie_tmp(1)
                    exit NEIGHB
                  end if
                end do
              end do NEIGHB
            end do
          end do

          ! E V
          idx1     = [ie(1,2)]
          exp_idx2 = reshape([1, 2], [1, 2])
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call qsort(idx2_arr(1,:))
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E V, idx1: "//int2str(idx1(1)))

          idx1     = [ie(6,7)]
          exp_idx2 = reshape([6, 7], [1, 2])
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call qsort(idx2_arr(1,:))
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E V, idx1: "//int2str(idx1(1)))

          ! E E
          idx1     = [ie(1,2)]
          exp_idx2 = reshape([ie(1,6), ie(1,5), ie(1,4), ie(2,7), ie(2,6), ie(2,5), ie(2,3)], [1, 7])
          call qsort(exp_idx2(1,:))
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_EDGE, 1, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call qsort(idx2_arr(1,:))
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E E, idx1: "//int2str(idx1(1)))

          idx1     = [ie(6,7)]
          exp_idx2 = reshape([ie(6,1), ie(6,2), ie(7,2)], [1, 3])
          call qsort(exp_idx2(1,:))
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_EDGE, 1, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call qsort(idx2_arr(1,:))
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E E, idx1: "//int2str(idx1(1)))

          ! E V
          idx1     = [ie(1,2)]
          exp_idx2 = reshape([4, 5], [1, 2])
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_CELL, 0, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call qsort(idx2_arr(1,:))
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E C, idx1: "//int2str(idx1(1)))

          idx1     = [ie(6,7)]
          exp_idx2 = reshape([6], [1, 1])
          call idx2_vec%reset()
          do j = 1, 100
            call g%get_neighb(IDX_EDGE, 1, IDX_CELL, 0, idx1, j, idx2, status)
            if (.not. status) exit
            call idx2_vec%push(idx2(1))
          end do
          idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
          call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb E C, idx1: "//int2str(idx1(1)))
        end block

        ! C V
        idx1     = [1]
        exp_idx2 = reshape([1, 4, 5], [1, 3])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C V, idx1: "//int2str(idx1(1)))

        idx1     = [4]
        exp_idx2 = reshape([1, 2, 6], [1, 3])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C V, idx1: "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 5], [1, 3])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C V, idx1: "//int2str(idx1(1)))

        ! C E
        idx1  = [1]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(3, nidx2, "triang_grid: get_neighb V C, nidx2: "//int2str(idx1(1)))

        idx1  = [5]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(3, nidx2, "triang_grid: get_neighb V C, nidx2: "//int2str(idx1(1)))

        idx1  = [6]
        nidx2 = 0
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx1, j, idx2, status)
          if (.not. status) exit
          nidx2 = nidx2+1
        end do
        call tc%assert_eq(3, nidx2, "triang_grid: get_neighb V C, nidx2: "//int2str(idx1(1)))

        ! C C
        idx1     = [1]
        exp_idx2 = reshape([3, 5], [1, 2])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C C, idx1: "//int2str(idx1(1)))

        idx1     = [5]
        exp_idx2 = reshape([1, 2, 4], [1, 3])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C C, idx1: "//int2str(idx1(1)))

        idx1     = [6]
        exp_idx2 = reshape([4], [1, 1])
        call idx2_vec%reset()
        do j = 1, 100
          call g%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx1, j, idx2, status)
          if (.not. status) exit
          call idx2_vec%push(idx2(1))
        end do
        idx2_arr = reshape(idx2_vec%to_array(), [1, idx2_vec%n])
        call qsort(idx2_arr(1,:))
        call tc%assert_eq(exp_idx2(1,:), idx2_arr(1,:), "triang_grid: get_neighb C C, idx1: "//int2str(idx1(1)))
      end block
    end subroutine

    subroutine test_tensor_grid()
      integer, parameter :: nx=31, ny=21, nz=11, nx2=(nx+1)/2, ny2=(ny+1)/2, nz2=(nz+1)/2
      real, allocatable  :: x(:), y(:), z(:)
      type(grid1D)       :: gx, gy, gz
      type(tensor_grid)  :: tg

      ! initialize 1D grids
      allocate (x(nx), y(ny), z(nz))  ! remove gfortran warning
      x = logspace(1.0, 10.0, nx)     ! logspacing -> non-equidistant!
      y = logspace(2.0,  5.0, ny)
      z = logspace(3.0, 20.0, nz)
      call gx%init(x)
      call gy%init(y)
      call gz%init(z)

      ! initialize tensor grid
      call tg%init([gx%get_ptr(), gy%get_ptr(), gz%get_ptr()])

      ! check attributes
      block
        call tc%assert_eq(      3, tg%dim,      "tensor_grid: dim")
        call tc%assert_eq(      3, tg%idx_dim,  "tensor_grid: idx_dim")
        call tc%assert_eq([4,4,4], tg%face_dim, "tensor_grid: face_dim")
        call tc%assert_eq(      8, tg%cell_dim, "tensor_grid: cell_dim")
      end block

      ! get_idx_bnd
      block
        integer :: bnd(3)

        call tg%get_idx_bnd(IDX_VERTEX, 0, bnd)
        call tc%assert_eq([nx, ny, nz], bnd, "tensor_grid: idx_bnd: IDX_VERTEX")

        call tg%get_idx_bnd(IDX_EDGE, 1, bnd)
        call tc%assert_eq([nx-1, ny, nz], bnd, "tensor_grid: idx_bnd: IDX_EDGE X")
        call tg%get_idx_bnd(IDX_EDGE, 2, bnd)
        call tc%assert_eq([nx, ny-1, nz], bnd, "tensor_grid: idx_bnd: IDX_EDGE Y")
        call tg%get_idx_bnd(IDX_EDGE, 3, bnd)
        call tc%assert_eq([nx, ny, nz-1], bnd, "tensor_grid: idx_bnd: IDX_EDGE Z")

        call tg%get_idx_bnd(IDX_FACE, 1, bnd)
        call tc%assert_eq([nx, ny-1, nz-1], bnd, "tensor_grid: idx_bnd: IDX_FACE X")
        call tg%get_idx_bnd(IDX_FACE, 2, bnd)
        call tc%assert_eq([nx-1, ny, nz-1], bnd, "tensor_grid: idx_bnd: IDX_FACE Y")
        call tg%get_idx_bnd(IDX_FACE, 3, bnd)
        call tc%assert_eq([nx-1, ny-1, nz], bnd, "tensor_grid: idx_bnd: IDX_FACE Z")

        call tg%get_idx_bnd(IDX_CELL, 0, bnd)
        call tc%assert_eq([nx-1, ny-1, nz-1], bnd, "tensor_grid: idx_bnd: IDX_CELL")
      end block

      ! get_vertex
      block
        integer, parameter :: ix(4) = [1,  5, 7, nx]
        integer, parameter :: iy(4) = [1, ny, 2, ny]
        integer, parameter :: iz(4) = [1,  1, 2, nz]
        integer            :: i, idx(3)
        real               :: p(3)

        do i = 1, size(ix)
          idx = [ix(i), iy(i), iz(i)]
          call tg%get_vertex(idx, p)
          call tc%assert_eq([x(ix(i)), y(iy(i)), z(iz(i))], p, 0.0, "tensor_grid: get_vertex "//int2str(i))
        end do
      end block

      ! get_edge, get_len
      block
        integer, parameter :: ix(4) = [1,    5, 7, nx  ]
        integer, parameter :: iy(4) = [1, ny-1, 2, ny-1]
        integer, parameter :: iz(4) = [1,    1, 2, nz  ]
        integer            :: i, idx(3)
        real               :: p(3,2), len

        do i = 1, size(ix)
          idx = [ix(i), iy(i), iz(i)]
          call tg%get_edge(idx, 2, p) ! get y edge
          call tc%assert_eq([x(ix(i)), y(iy(i)  ), z(iz(i))], p(:,1), 0.0, "tensor_grid: get_edge "//int2str(i)//" A")
          call tc%assert_eq([x(ix(i)), y(iy(i)+1), z(iz(i))], p(:,2), 0.0, "tensor_grid: get_edge "//int2str(i)//" B")
          len = tg%get_len(idx, 2)
          call tc%assert_eq(y(iy(i)+1)-y(iy(i)), len, 5e-16, "tensor_grid: get_len "//int2str(i))
        end do
      end block

      ! get_face, get_surf
      block
        integer, parameter :: ix(4) = [1,    5, 7, nx-1]
        integer, parameter :: iy(4) = [1, ny-1, 2, ny-1]
        integer, parameter :: iz(4) = [1,    1, 2, nz  ]
        integer            :: i, idx(3)
        real               :: p(3,4), surf

        do i = 1, size(ix)
          idx = [ix(i), iy(i), iz(i)]
          call tg%get_face(idx, 3, p) ! get z face
          call tc%assert_eq([x(ix(i)  ), y(iy(i)  ), z(iz(i))], p(:,1), 0.0, "tensor_grid: get_face "//int2str(i)//" A")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)  ), z(iz(i))], p(:,2), 0.0, "tensor_grid: get_face "//int2str(i)//" B")
          call tc%assert_eq([x(ix(i)  ), y(iy(i)+1), z(iz(i))], p(:,3), 0.0, "tensor_grid: get_face "//int2str(i)//" C")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)+1), z(iz(i))], p(:,4), 0.0, "tensor_grid: get_face "//int2str(i)//" D")
          surf = tg%get_surf(idx, 3)
          call tc%assert_eq((x(ix(i)+1)-x(ix(i)))*(y(iy(i)+1)-y(iy(i))), surf, 5e-16, "tensor_grid: get_surf "//int2str(i))
        end do
      end block

      ! get cell, get_vol
      block
        integer, parameter :: ix(4) = [1,    5, 7, nx-1]
        integer, parameter :: iy(4) = [1, ny-1, 2, ny-1]
        integer, parameter :: iz(4) = [1,    1, 2, nz-1]
        integer            :: i, idx(3)
        real               :: p(3,8), vol

        do i = 1, size(ix)
          idx = [ix(i), iy(i), iz(i)]
          call tg%get_cell(idx, p)
          call tc%assert_eq([x(ix(i)  ), y(iy(i)  ), z(iz(i)  )], p(:,1), 0.0, "tensor_grid: get_cell "//int2str(i)//" A")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)  ), z(iz(i)  )], p(:,2), 0.0, "tensor_grid: get_cell "//int2str(i)//" B")
          call tc%assert_eq([x(ix(i)  ), y(iy(i)+1), z(iz(i)  )], p(:,3), 0.0, "tensor_grid: get_cell "//int2str(i)//" C")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)+1), z(iz(i)  )], p(:,4), 0.0, "tensor_grid: get_cell "//int2str(i)//" D")
          call tc%assert_eq([x(ix(i)  ), y(iy(i)  ), z(iz(i)+1)], p(:,5), 0.0, "tensor_grid: get_cell "//int2str(i)//" E")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)  ), z(iz(i)+1)], p(:,6), 0.0, "tensor_grid: get_cell "//int2str(i)//" F")
          call tc%assert_eq([x(ix(i)  ), y(iy(i)+1), z(iz(i)+1)], p(:,7), 0.0, "tensor_grid: get_cell "//int2str(i)//" G")
          call tc%assert_eq([x(ix(i)+1), y(iy(i)+1), z(iz(i)+1)], p(:,8), 0.0, "tensor_grid: get_cell "//int2str(i)//" H")
          vol = tg%get_vol(idx)
          call tc%assert_eq((x(ix(i)+1)-x(ix(i)))*(y(iy(i)+1)-y(iy(i)))*(z(iz(i)+1)-z(iz(i))), vol, 5e-16, "tensor_grid: get_vol "//int2str(i))
        end do
      end block

      ! get_max_neighb
      block
        call tc%assert_eq(6, tg%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0), "tensor_grid: get_max_neighb V V")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   1), "tensor_grid: get_max_neighb V EX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   2), "tensor_grid: get_max_neighb V EY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE,   3), "tensor_grid: get_max_neighb V EZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   1), "tensor_grid: get_max_neighb V FX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   2), "tensor_grid: get_max_neighb V FY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_VERTEX, 0, IDX_FACE,   3), "tensor_grid: get_max_neighb V FZ")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_VERTEX, 0, IDX_CELL,   0), "tensor_grid: get_max_neighb V C")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   1, IDX_VERTEX, 0), "tensor_grid: get_max_neighb EX V")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   2, IDX_VERTEX, 0), "tensor_grid: get_max_neighb EY V")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   3, IDX_VERTEX, 0), "tensor_grid: get_max_neighb EZ V")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   1), "tensor_grid: get_max_neighb EX EX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   2), "tensor_grid: get_max_neighb EX EY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   1, IDX_EDGE,   3), "tensor_grid: get_max_neighb EX EZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   2, IDX_EDGE,   1), "tensor_grid: get_max_neighb EY EX")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_EDGE,   2, IDX_EDGE,   2), "tensor_grid: get_max_neighb EY EY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   2, IDX_EDGE,   3), "tensor_grid: get_max_neighb EY EZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   3, IDX_EDGE,   1), "tensor_grid: get_max_neighb EZ EX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   3, IDX_EDGE,   2), "tensor_grid: get_max_neighb EZ EY")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_EDGE,   3, IDX_EDGE,   3), "tensor_grid: get_max_neighb EZ EZ")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   1), "tensor_grid: get_max_neighb EX FX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   2), "tensor_grid: get_max_neighb EX FY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   1, IDX_FACE,   3), "tensor_grid: get_max_neighb EX FZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   2, IDX_FACE,   1), "tensor_grid: get_max_neighb EY FX")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_EDGE,   2, IDX_FACE,   2), "tensor_grid: get_max_neighb EY FY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   2, IDX_FACE,   3), "tensor_grid: get_max_neighb EY FZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   3, IDX_FACE,   1), "tensor_grid: get_max_neighb EZ FX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_EDGE,   3, IDX_FACE,   2), "tensor_grid: get_max_neighb EZ FY")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_EDGE,   3, IDX_FACE,   3), "tensor_grid: get_max_neighb EZ FZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   1, IDX_CELL,   0), "tensor_grid: get_max_neighb EX C")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   2, IDX_CELL,   0), "tensor_grid: get_max_neighb EY C")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_EDGE,   3, IDX_CELL,   0), "tensor_grid: get_max_neighb EZ C")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   1, IDX_VERTEX, 0), "tensor_grid: get_max_neighb FX V")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   2, IDX_VERTEX, 0), "tensor_grid: get_max_neighb FY V")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   3, IDX_VERTEX, 0), "tensor_grid: get_max_neighb FZ V")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_FACE,   1, IDX_EDGE,   1), "tensor_grid: get_max_neighb FX EX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   1, IDX_EDGE,   2), "tensor_grid: get_max_neighb FX EY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   1, IDX_EDGE,   3), "tensor_grid: get_max_neighb FX EZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   2, IDX_EDGE,   1), "tensor_grid: get_max_neighb FY EX")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_FACE,   2, IDX_EDGE,   2), "tensor_grid: get_max_neighb FY EY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   2, IDX_EDGE,   3), "tensor_grid: get_max_neighb FY EZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   3, IDX_EDGE,   1), "tensor_grid: get_max_neighb FZ EX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   3, IDX_EDGE,   2), "tensor_grid: get_max_neighb FZ EY")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_FACE,   3, IDX_EDGE,   3), "tensor_grid: get_max_neighb FZ EZ")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_FACE,   1, IDX_FACE,   1), "tensor_grid: get_max_neighb FX FX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   1, IDX_FACE,   2), "tensor_grid: get_max_neighb FX FY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   1, IDX_FACE,   3), "tensor_grid: get_max_neighb FX FZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   2, IDX_FACE,   1), "tensor_grid: get_max_neighb FY FX")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_FACE,   2, IDX_FACE,   2), "tensor_grid: get_max_neighb FY FY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   2, IDX_FACE,   3), "tensor_grid: get_max_neighb FY FZ")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   3, IDX_FACE,   1), "tensor_grid: get_max_neighb FZ FX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_FACE,   3, IDX_FACE,   2), "tensor_grid: get_max_neighb FZ FY")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_FACE,   3, IDX_FACE,   3), "tensor_grid: get_max_neighb FZ FZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   1, IDX_CELL,   0), "tensor_grid: get_max_neighb FX C")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   2, IDX_CELL,   0), "tensor_grid: get_max_neighb FY C")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_FACE,   3, IDX_CELL,   0), "tensor_grid: get_max_neighb FZ C")
        call tc%assert_eq(8, tg%get_max_neighb(IDX_CELL,   0, IDX_VERTEX, 0), "tensor_grid: get_max_neighb C V")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   1), "tensor_grid: get_max_neighb C EX")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   2), "tensor_grid: get_max_neighb C EY")
        call tc%assert_eq(4, tg%get_max_neighb(IDX_CELL,   0, IDX_EDGE,   3), "tensor_grid: get_max_neighb C EZ")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_CELL,   0, IDX_FACE,   1), "tensor_grid: get_max_neighb C FX")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_CELL,   0, IDX_FACE,   1), "tensor_grid: get_max_neighb C FY")
        call tc%assert_eq(2, tg%get_max_neighb(IDX_CELL,   0, IDX_FACE,   1), "tensor_grid: get_max_neighb C FZ")
        call tc%assert_eq(6, tg%get_max_neighb(IDX_CELL,   0, IDX_CELL,   0), "tensor_grid: get_max_neighb C C")
      end block

      ! get_neighb
      block
        integer :: idx(3,5), idx2(3)
        logical :: status

        idx(:,1) = [  1,  1,  1]
        idx(:,2) = [ nx, ny, nz]
        idx(:,3) = [  1,ny2,  1]
        idx(:,4) = [nx2,  1,nz2]
        idx(:,5) = [nx2,ny2,nz2]

        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([2,1,1], idx2, "tensor_grid: get_neighb V V 1 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,1), 2, idx2, status)
        call tc%assert_eq([1,2,1], idx2, "tensor_grid: get_neighb V V 1 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,1), 3, idx2, status)
        call tc%assert_eq([1,1,2], idx2, "tensor_grid: get_neighb V V 1 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,1), 4, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb V V 1 4")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny,nz], idx2, "tensor_grid: get_neighb V V 2 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,2), 2, idx2, status)
        call tc%assert_eq([nx,ny-1,nz], idx2, "tensor_grid: get_neighb V V 2 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,2), 3, idx2, status)
        call tc%assert_eq([nx,ny,nz-1], idx2, "tensor_grid: get_neighb V V 2 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,2), 4, idx2, status)
        call tc%assert(         .not. status, "tensor_grid: get_neighb V V 2 4")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,3), 1, idx2, status)
        call tc%assert_eq([2,ny2  ,1], idx2, "tensor_grid: get_neighb V V 3 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2-1,1], idx2, "tensor_grid: get_neighb V V 3 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,3), 3, idx2, status)
        call tc%assert_eq([1,ny2+1,1], idx2, "tensor_grid: get_neighb V V 3 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,3), 4, idx2, status)
        call tc%assert_eq([1,ny2  ,2], idx2, "tensor_grid: get_neighb V V 3 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,3), 5, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb V V 3 5")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2  ], idx2, "tensor_grid: get_neighb V V 4 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2+1,1,nz2  ], idx2, "tensor_grid: get_neighb V V 4 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2  ,2,nz2  ], idx2, "tensor_grid: get_neighb V V 4 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2-1], idx2, "tensor_grid: get_neighb V V 4 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 5, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2+1], idx2, "tensor_grid: get_neighb V V 4 5")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,4), 6, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb V V 4 6")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb V V 5 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2+1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb V V 5 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2  ,ny2-1,nz2  ], idx2, "tensor_grid: get_neighb V V 5 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2  ,ny2+1,nz2  ], idx2, "tensor_grid: get_neighb V V 5 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 5, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2-1], idx2, "tensor_grid: get_neighb V V 5 5")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 6, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2+1], idx2, "tensor_grid: get_neighb V V 5 6")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx(:,5), 7, idx2, status)
        call tc%assert(                .not. status, "tensor_grid: get_neighb V V 5 7")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb V EX 1 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx(:,1), 2, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb V EX 1 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 2, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx,ny-1,nz], idx2, "tensor_grid: get_neighb V EY 2 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 2, idx(:,2), 2, idx2, status)
        call tc%assert(         .not. status, "tensor_grid: get_neighb V EY 2 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 2, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2-1,1], idx2, "tensor_grid: get_neighb V EY 3 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 2, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2  ,1], idx2, "tensor_grid: get_neighb V EY 3 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 2, idx(:,3), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb V EY 3 3")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 3, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2,1,nz2-1], idx2, "tensor_grid: get_neighb V EZ 4 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 3, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2,1,nz2  ], idx2, "tensor_grid: get_neighb V EZ 4 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 3, idx(:,4), 3, idx2, status)
        call tc%assert(          .not. status, "tensor_grid: get_neighb V EZ 4 3")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2,nz2], idx2, "tensor_grid: get_neighb V EX 5 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2], idx2, "tensor_grid: get_neighb V EX 5 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, idx(:,5), 3, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb V EX 5 3")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb V FX 1 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,1), 2, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb V FX 1 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny,nz-1], idx2, "tensor_grid: get_neighb V FY 2 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,2), 2, idx2, status)
        call tc%assert(           .not. status, "tensor_grid: get_neighb V FY 2 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 3, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2-1,1], idx2, "tensor_grid: get_neighb V FZ 3 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 3, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2  ,1], idx2, "tensor_grid: get_neighb V FZ 3 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 3, idx(:,3), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb V FZ 3 3")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2-1], idx2, "tensor_grid: get_neighb V FY 4 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2-1], idx2, "tensor_grid: get_neighb V FY 4 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2  ], idx2, "tensor_grid: get_neighb V FY 4 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2  ], idx2, "tensor_grid: get_neighb V FY 4 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 2, idx(:,4), 5, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb V FY 4 5")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2,ny2-1,nz2-1], idx2, "tensor_grid: get_neighb V FX 5 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2,ny2  ,nz2-1], idx2, "tensor_grid: get_neighb V FX 5 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2,ny2-1,nz2  ], idx2, "tensor_grid: get_neighb V FX 5 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb V FX 5 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_FACE, 1, idx(:,5), 5, idx2, status)
        call tc%assert(              .not. status, "tensor_grid: get_neighb V FX 5 5")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb V C 1 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,1), 2, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb V C 1 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz-1], idx2, "tensor_grid: get_neighb V C 2 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,2), 2, idx2, status)
        call tc%assert(             .not. status, "tensor_grid: get_neighb V C 2 2")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2-1,1], idx2, "tensor_grid: get_neighb V C 3 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2  ,1], idx2, "tensor_grid: get_neighb V C 3 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,3), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb V C 3 3")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2-1], idx2, "tensor_grid: get_neighb V C 4 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2-1], idx2, "tensor_grid: get_neighb V C 4 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2  ], idx2, "tensor_grid: get_neighb V C 4 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2  ], idx2, "tensor_grid: get_neighb V C 4 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,4), 5, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb V C 4 5")

        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2-1,nz2-1], idx2, "tensor_grid: get_neighb V C 5 1")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2  ,ny2-1,nz2-1], idx2, "tensor_grid: get_neighb V C 5 2")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2-1], idx2, "tensor_grid: get_neighb V C 5 3")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2-1], idx2, "tensor_grid: get_neighb V C 5 4")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 5, idx2, status)
        call tc%assert_eq([nx2-1,ny2-1,nz2  ], idx2, "tensor_grid: get_neighb V C 5 5")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 6, idx2, status)
        call tc%assert_eq([nx2  ,ny2-1,nz2  ], idx2, "tensor_grid: get_neighb V C 5 6")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 7, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb V C 5 7")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 8, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb V C 5 8")
        call tg%get_neighb(IDX_VERTEX, 0, IDX_CELL, 0, idx(:,5), 9, idx2, status)
        call tc%assert(                .not. status, "tensor_grid: get_neighb V C 4 9")

        idx(:,1) = [  1,   1,  1]
        idx(:,2) = [ nx,ny-1, nz]
        idx(:,3) = [  1, ny2,  1]
        idx(:,4) = [nx2,   1,nz2]
        idx(:,5) = [nx2, ny2,nz2]

        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb EX V 1 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,1), 2, idx2, status)
        call tc%assert_eq([2,1,1], idx2, "tensor_grid: get_neighb EX V 1 2")
        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,1), 3, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb EX V 1 3")

        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx,ny-1,nz], idx2, "tensor_grid: get_neighb EY V 2 1")
        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,2), 2, idx2, status)
        call tc%assert_eq([nx,ny  ,nz], idx2, "tensor_grid: get_neighb EY V 2 2")
        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,2), 3, idx2, status)
        call tc%assert(         .not. status, "tensor_grid: get_neighb EY V 2 3")

        call tg%get_neighb(IDX_EDGE, 3, IDX_VERTEX, 0, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2,1], idx2, "tensor_grid: get_neighb EZ V 3 1")
        call tg%get_neighb(IDX_EDGE, 3, IDX_VERTEX, 0, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2,2], idx2, "tensor_grid: get_neighb EZ V 3 2")
        call tg%get_neighb(IDX_EDGE, 3, IDX_VERTEX, 0, idx(:,3), 3, idx2, status)
        call tc%assert(      .not. status, "tensor_grid: get_neighb EZ V 3 3")

        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2,1,nz2], idx2, "tensor_grid: get_neighb EY V 4 1")
        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2,2,nz2], idx2, "tensor_grid: get_neighb EY V 4 2")
        call tg%get_neighb(IDX_EDGE, 2, IDX_VERTEX, 0, idx(:,4), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb EY V 4 3")

        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2], idx2, "tensor_grid: get_neighb EX V 5 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2+1,ny2,nz2], idx2, "tensor_grid: get_neighb EX V 5 2")
        call tg%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(:,5), 3, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb EX V 5 3")

        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 2, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb EX EY 1 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 2, idx(:,1), 2, idx2, status)
        call tc%assert_eq([2,1,1], idx2, "tensor_grid: get_neighb EX EY 1 2")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 2, idx(:,1), 3, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb EX EY 1 3")

        call tg%get_neighb(IDX_EDGE, 2, IDX_EDGE, 2, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz  ], idx2, "tensor_grid: get_neighb EY EY 2 1")
        call tg%get_neighb(IDX_EDGE, 2, IDX_EDGE, 2, idx(:,2), 2, idx2, status)
        call tc%assert_eq([nx  ,ny-2,nz  ], idx2, "tensor_grid: get_neighb EY EY 2 2")
        call tg%get_neighb(IDX_EDGE, 2, IDX_EDGE, 2, idx(:,2), 3, idx2, status)
        call tc%assert_eq([nx  ,ny-1,nz-1], idx2, "tensor_grid: get_neighb EY EY 2 3")
        call tg%get_neighb(IDX_EDGE, 2, IDX_EDGE, 2, idx(:,2), 4, idx2, status)
        call tc%assert(             .not. status, "tensor_grid: get_neighb EY EY 2 4")

        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 3, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2-1], idx2, "tensor_grid: get_neighb EX EZ 5 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 3, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2+1,ny2,nz2-1], idx2, "tensor_grid: get_neighb EX EZ 5 2")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 3, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2  ], idx2, "tensor_grid: get_neighb EX EZ 5 3")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 3, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2+1,ny2,nz2  ], idx2, "tensor_grid: get_neighb EX EZ 5 4")
        call tg%get_neighb(IDX_EDGE, 1, IDX_EDGE, 3, idx(:,5), 5, idx2, status)
        call tc%assert(              .not. status, "tensor_grid: get_neighb EX EZ 5 5")

        call tg%get_neighb(IDX_EDGE, 1, IDX_FACE, 1, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb EX FX 1 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_FACE, 1, idx(:,1), 2, idx2, status)
        call tc%assert_eq([2,1,1], idx2, "tensor_grid: get_neighb EX FX 1 2")
        call tg%get_neighb(IDX_EDGE, 1, IDX_FACE, 1, idx(:,1), 3, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb EX FX 1 3")

        call tg%get_neighb(IDX_EDGE, 2, IDX_FACE, 3, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz], idx2, "tensor_grid: get_neighb EY FZ 2 1")
        call tg%get_neighb(IDX_EDGE, 2, IDX_FACE, 3, idx(:,2), 2, idx2, status)
        call tc%assert(           .not. status, "tensor_grid: get_neighb EY FZ 2 2")

        call tg%get_neighb(IDX_EDGE, 3, IDX_FACE, 2, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2,nz2], idx2, "tensor_grid: get_neighb EZ FY 5 1")
        call tg%get_neighb(IDX_EDGE, 3, IDX_FACE, 2, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2], idx2, "tensor_grid: get_neighb EZ FY 5 1")
        call tg%get_neighb(IDX_EDGE, 3, IDX_FACE, 2, idx(:,5), 3, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb EZ FY 5 3")

        call tg%get_neighb(IDX_EDGE, 1, IDX_CELL, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb EX C 1 1")
        call tg%get_neighb(IDX_EDGE, 1, IDX_CELL, 0, idx(:,1), 2, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb EX C 1 2")

        call tg%get_neighb(IDX_EDGE, 2, IDX_CELL, 0, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2-1], idx2, "tensor_grid: get_neighb EY C 4 1")
        call tg%get_neighb(IDX_EDGE, 2, IDX_CELL, 0, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2-1], idx2, "tensor_grid: get_neighb EY C 4 2")
        call tg%get_neighb(IDX_EDGE, 2, IDX_CELL, 0, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2  ], idx2, "tensor_grid: get_neighb EY C 4 3")
        call tg%get_neighb(IDX_EDGE, 2, IDX_CELL, 0, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2  ], idx2, "tensor_grid: get_neighb EY C 4 4")
        call tg%get_neighb(IDX_EDGE, 2, IDX_CELL, 0, idx(:,4), 5, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb EY C 4 5")

        idx(:,1) = [   1,  1,   1]
        idx(:,2) = [nx-1, ny,nz-1]
        idx(:,3) = [   1,ny2,   1]
        idx(:,4) = [ nx2,  1, nz2]
        idx(:,5) = [ nx2,ny2, nz2]

        call tg%get_neighb(IDX_FACE, 1, IDX_VERTEX, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb FX V 1 1")
        call tg%get_neighb(IDX_FACE, 1, IDX_VERTEX, 0, idx(:,1), 2, idx2, status)
        call tc%assert_eq([1,2,1], idx2, "tensor_grid: get_neighb FX V 1 2")
        call tg%get_neighb(IDX_FACE, 1, IDX_VERTEX, 0, idx(:,1), 3, idx2, status)
        call tc%assert_eq([1,1,2], idx2, "tensor_grid: get_neighb FX V 1 3")
        call tg%get_neighb(IDX_FACE, 1, IDX_VERTEX, 0, idx(:,1), 4, idx2, status)
        call tc%assert_eq([1,2,2], idx2, "tensor_grid: get_neighb FX V 1 4")
        call tg%get_neighb(IDX_FACE, 1, IDX_VERTEX, 0, idx(:,1), 5, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb FX V 1 5")

        call tg%get_neighb(IDX_FACE, 2, IDX_VERTEX, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2  ], idx2, "tensor_grid: get_neighb FY V 5 1")
        call tg%get_neighb(IDX_FACE, 2, IDX_VERTEX, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2+1,ny2,nz2  ], idx2, "tensor_grid: get_neighb FY V 5 2")
        call tg%get_neighb(IDX_FACE, 2, IDX_VERTEX, 0, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2  ,ny2,nz2+1], idx2, "tensor_grid: get_neighb FY V 5 3")
        call tg%get_neighb(IDX_FACE, 2, IDX_VERTEX, 0, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2+1,ny2,nz2+1], idx2, "tensor_grid: get_neighb FY V 5 4")
        call tg%get_neighb(IDX_FACE, 2, IDX_VERTEX, 0, idx(:,5), 5, idx2, status)
        call tc%assert(              .not. status, "tensor_grid: get_neighb FY V 5 5")

        call tg%get_neighb(IDX_FACE, 3, IDX_EDGE, 1, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2  ,1], idx2, "tensor_grid: get_neighb FZ EX 3 1")
        call tg%get_neighb(IDX_FACE, 3, IDX_EDGE, 1, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2+1,1], idx2, "tensor_grid: get_neighb FZ EX 3 2")
        call tg%get_neighb(IDX_FACE, 3, IDX_EDGE, 1, idx(:,3), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb FZ EX 3 3")

        call tg%get_neighb(IDX_FACE, 2, IDX_EDGE, 2, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2  ], idx2, "tensor_grid: get_neighb FY EY 4 1")
        call tg%get_neighb(IDX_FACE, 2, IDX_EDGE, 2, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2+1,1,nz2  ], idx2, "tensor_grid: get_neighb FY EY 4 2")
        call tg%get_neighb(IDX_FACE, 2, IDX_EDGE, 2, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2+1], idx2, "tensor_grid: get_neighb FY EY 4 3")
        call tg%get_neighb(IDX_FACE, 2, IDX_EDGE, 2, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2+1,1,nz2+1], idx2, "tensor_grid: get_neighb FY EY 4 4")
        call tg%get_neighb(IDX_FACE, 2, IDX_EDGE, 2, idx(:,4), 5, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb FY EY 4 5")

        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb FX EX 5 1")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb FX EX 5 2")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2-1,ny2+1,nz2  ], idx2, "tensor_grid: get_neighb FX EX 5 3")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2  ,ny2+1,nz2  ], idx2, "tensor_grid: get_neighb FX EX 5 4")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 5, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2+1], idx2, "tensor_grid: get_neighb FX EX 5 5")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 6, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2+1], idx2, "tensor_grid: get_neighb FX EX 5 6")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 7, idx2, status)
        call tc%assert_eq([nx2-1,ny2+1,nz2+1], idx2, "tensor_grid: get_neighb FX EX 5 7")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 8, idx2, status)
        call tc%assert_eq([nx2  ,ny2+1,nz2+1], idx2, "tensor_grid: get_neighb FX EX 5 8")
        call tg%get_neighb(IDX_FACE, 1, IDX_EDGE, 1, idx(:,5), 9, idx2, status)
        call tc%assert(                .not. status, "tensor_grid: get_neighb FX EX 5 9")

        call tg%get_neighb(IDX_FACE, 1, IDX_FACE, 1, idx(:,3), 1, idx2, status)
        call tc%assert_eq([2,ny2  ,1], idx2, "tensor_grid: get_neighb FX FX 3 1")
        call tg%get_neighb(IDX_FACE, 1, IDX_FACE, 1, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2-1,1], idx2, "tensor_grid: get_neighb FX FX 3 2")
        call tg%get_neighb(IDX_FACE, 1, IDX_FACE, 1, idx(:,3), 3, idx2, status)
        call tc%assert_eq([1,ny2+1,1], idx2, "tensor_grid: get_neighb FX FX 3 3")
        call tg%get_neighb(IDX_FACE, 1, IDX_FACE, 1, idx(:,3), 4, idx2, status)
        call tc%assert_eq([1,ny2  ,2], idx2, "tensor_grid: get_neighb FX FX 3 3")
        call tg%get_neighb(IDX_FACE, 1, IDX_FACE, 1, idx(:,3), 5, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb FX FX 3 5")

        call tg%get_neighb(IDX_FACE, 2, IDX_FACE, 3, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz-1], idx2, "tensor_grid: get_neighb FY FZ 2 1")
        call tg%get_neighb(IDX_FACE, 2, IDX_FACE, 3, idx(:,2), 2, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz  ], idx2, "tensor_grid: get_neighb FY FZ 2 2")
        call tg%get_neighb(IDX_FACE, 2, IDX_FACE, 3, idx(:,2), 3, idx2, status)
        call tc%assert(             .not. status, "tensor_grid: get_neighb FY FZ 2 3")

        call tg%get_neighb(IDX_FACE, 1, IDX_CELL, 0, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2,1], idx2, "tensor_grid: get_neighb FX C 3 1")
        call tg%get_neighb(IDX_FACE, 1, IDX_CELL, 0, idx(:,3), 2, idx2, status)
        call tc%assert(      .not. status, "tensor_grid: get_neighb FX C 3 2")

        call tg%get_neighb(IDX_FACE, 3, IDX_CELL, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2,ny2,nz2-1], idx2, "tensor_grid: get_neighb FZ C 5 1")
        call tg%get_neighb(IDX_FACE, 3, IDX_CELL, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2,ny2,nz2  ], idx2, "tensor_grid: get_neighb FZ C 5 2")
        call tg%get_neighb(IDX_FACE, 3, IDX_CELL, 0, idx(:,5), 3, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb FZ C 5 3")

        idx(:,1) = [   1,   1,   1]
        idx(:,2) = [nx-1,ny-1,nz-1]
        idx(:,3) = [   1, ny2,   1]
        idx(:,4) = [ nx2,   1, nz2]
        idx(:,5) = [ nx2, ny2, nz2]

        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 1, idx2, status)
        call tc%assert_eq([1,1,1], idx2, "tensor_grid: get_neighb C V 1 1")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 2, idx2, status)
        call tc%assert_eq([2,1,1], idx2, "tensor_grid: get_neighb C V 1 2")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 3, idx2, status)
        call tc%assert_eq([1,2,1], idx2, "tensor_grid: get_neighb C V 1 3")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 4, idx2, status)
        call tc%assert_eq([2,2,1], idx2, "tensor_grid: get_neighb C V 1 4")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 5, idx2, status)
        call tc%assert_eq([1,1,2], idx2, "tensor_grid: get_neighb C V 1 5")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 6, idx2, status)
        call tc%assert_eq([2,1,2], idx2, "tensor_grid: get_neighb C V 1 6")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 7, idx2, status)
        call tc%assert_eq([1,2,2], idx2, "tensor_grid: get_neighb C V 1 7")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 8, idx2, status)
        call tc%assert_eq([2,2,2], idx2, "tensor_grid: get_neighb C V 1 8")
        call tg%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx(:,1), 9, idx2, status)
        call tc%assert(    .not. status, "tensor_grid: get_neighb C V 1 9")

        call tg%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx(:,2), 1, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz-1], idx2, "tensor_grid: get_neighb C EX 2 1")
        call tg%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx(:,2), 2, idx2, status)
        call tc%assert_eq([nx-1,ny  ,nz-1], idx2, "tensor_grid: get_neighb C EX 2 2")
        call tg%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx(:,2), 3, idx2, status)
        call tc%assert_eq([nx-1,ny-1,nz  ], idx2, "tensor_grid: get_neighb C EX 2 3")
        call tg%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx(:,2), 4, idx2, status)
        call tc%assert_eq([nx-1,ny  ,nz  ], idx2, "tensor_grid: get_neighb C EX 2 4")
        call tg%get_neighb(IDX_CELL, 0, IDX_EDGE, 1, idx(:,2), 5, idx2, status)
        call tc%assert(             .not. status, "tensor_grid: get_neighb C EX 2 5")

        call tg%get_neighb(IDX_CELL, 0, IDX_FACE, 2, idx(:,3), 1, idx2, status)
        call tc%assert_eq([1,ny2  ,1], idx2, "tensor_grid: get_neighb C FY 3 1")
        call tg%get_neighb(IDX_CELL, 0, IDX_FACE, 2, idx(:,3), 2, idx2, status)
        call tc%assert_eq([1,ny2+1,1], idx2, "tensor_grid: get_neighb C FY 3 2")
        call tg%get_neighb(IDX_CELL, 0, IDX_FACE, 2, idx(:,3), 3, idx2, status)
        call tc%assert(        .not. status, "tensor_grid: get_neighb C FY 3 3")

        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 1, idx2, status)
        call tc%assert_eq([nx2-1,1,nz2  ], idx2, "tensor_grid: get_neighb C C 4 1")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 2, idx2, status)
        call tc%assert_eq([nx2+1,1,nz2  ], idx2, "tensor_grid: get_neighb C C 4 2")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 3, idx2, status)
        call tc%assert_eq([nx2  ,2,nz2  ], idx2, "tensor_grid: get_neighb C C 4 3")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 4, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2-1], idx2, "tensor_grid: get_neighb C C 4 4")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 5, idx2, status)
        call tc%assert_eq([nx2  ,1,nz2+1], idx2, "tensor_grid: get_neighb C C 4 5")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,4), 6, idx2, status)
        call tc%assert(            .not. status, "tensor_grid: get_neighb C C 4 6")

        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 1, idx2, status)
        call tc%assert_eq([nx2-1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb C C 5 1")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 2, idx2, status)
        call tc%assert_eq([nx2+1,ny2  ,nz2  ], idx2, "tensor_grid: get_neighb C C 5 2")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 3, idx2, status)
        call tc%assert_eq([nx2  ,ny2-1,nz2  ], idx2, "tensor_grid: get_neighb C C 5 3")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 4, idx2, status)
        call tc%assert_eq([nx2  ,ny2+1,nz2  ], idx2, "tensor_grid: get_neighb C C 5 4")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 5, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2-1], idx2, "tensor_grid: get_neighb C C 5 5")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 6, idx2, status)
        call tc%assert_eq([nx2  ,ny2  ,nz2+1], idx2, "tensor_grid: get_neighb C C 5 6")
        call tg%get_neighb(IDX_CELL, 0, IDX_CELL, 0, idx(:,5), 7, idx2, status)
        call tc%assert(                .not. status, "tensor_grid: get_neighb C C 5 7")
      end block
    end subroutine

  end subroutine

end module
