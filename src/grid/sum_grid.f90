#include "../util/macro.f90.inc"

module sum_grid_m
  use error_m, only: assert_failed, program_error
  use grid_m,  only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL

  implicit none

  private
  public sum_grid

  type, extends(grid) :: sum_grid
    !! Combination of multiple sub-grids by direct sum
    class(grid), pointer :: g(:)
      !! sub-grids
  contains
    procedure :: init           => sum_grid_init
    procedure :: get_idx_bnd    => sum_grid_get_idx_bnd
    procedure :: get_vertex     => sum_grid_get_vertex
    procedure :: get_edge       => sum_grid_get_edge
    procedure :: get_face       => sum_grid_get_face
    procedure :: get_cell       => sum_grid_get_cell
    procedure :: get_len        => sum_grid_get_len
    procedure :: get_surf       => sum_grid_get_surf
    procedure :: get_vol        => sum_grid_get_vol
    procedure :: get_max_neighb => sum_grid_get_max_neighb
    procedure :: get_neighb     => sum_grid_get_neighb
    procedure :: to_subgrid     => sum_grid_to_subgrid
    procedure :: to_sumgrid     => sum_grid_to_sumgrid
  end type
contains
  subroutine sum_grid_init(this, g)
    !! initialize sum grid
    class(sum_grid),     intent(out) :: this
    class(grid), target, intent(in)  :: g(:)
      !! sub-grids

    ! init base
    call this%grid_init(g(1)%dim, g(1)%idx_dim, g(1)%face_dim, g(1)%cell_dim)

    ! save sub-grid pointers
    this%g => g

    ! init tables
    call this%init_tab_all()
  end subroutine

  subroutine sum_grid_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(out) :: idx_bnd(:)
      !! output upper bound for each index (idx_dim)

    integer :: i

    ASSERT(this%idx_allowed(idx_type, idx_dir))
    ASSERT(size(idx_bnd) == this%idx_dim)

    idx_bnd = 0
    do i = 1, size(this%g)
      call this%g(i)%get_idx_bnd(idx_type, idx_dir, idx_bnd)
    end do
  end subroutine

  subroutine sum_grid_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx(:)
      !! vertex indices
    real,            intent(out) :: p(:)
      !! output vertex coordinates (dim)

    integer :: i, idx_sub(this%idx_dim)

    ASSERT(size(p) == this%dim)

    call this%to_subgrid(IDX_VERTEX, 0, idx, i, idx_sub)
    call this%g(i)%get_vertex(idx_sub, p)
  end subroutine

  subroutine sum_grid_get_edge(this, idx, idx_dir, p)
    !! get edge coordinates from grid indices
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx(:)
      !! edge indices
    integer,         intent(in)  :: idx_dir
      !! edge direction
    real,            intent(out) :: p(:,:)
      !! output edge coordinates (dim x 2)

    integer :: i, idx_sub(this%idx_dim)

    ASSERT(all(shape(p) == [this%dim, 2]))

    call this%to_subgrid(IDX_EDGE, idx_dir, idx, i, idx_sub)
    call this%g(i)%get_edge(idx_sub, idx_dir, p)
  end subroutine

  subroutine sum_grid_get_face(this, idx, idx_dir, p)
    !! get face coordinates from grid indices
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx(:)
      !! face indices
    integer,         intent(in)  :: idx_dir
      !! face direction
    real,            intent(out) :: p(:,:)
      !! output face coordinates (dim x face_dim(idx_dir))

    integer :: i, idx_sub(this%idx_dim)

    ASSERT(all(shape(p) == [this%dim, this%face_dim(idx_dir)]))

    call this%to_subgrid(IDX_FACE, idx_dir, idx, i, idx_sub)
    call this%g(i)%get_face(idx_sub, idx_dir, p)
  end subroutine

  subroutine sum_grid_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx(:)
      !! cell indices
    real,            intent(out) :: p(:,:)
      !! output cell coordinates (dim x cell_dim)

    integer :: i, idx_sub(this%idx_dim)

    ASSERT(all(shape(p) == [this%dim, this%cell_dim]))

    call this%to_subgrid(IDX_CELL, 0, idx, i, idx_sub)
    call this%g(i)%get_cell(idx_sub, p)
  end subroutine

  function sum_grid_get_len(this, idx, idx_dir) result(len)
    !! get edge length
    class(sum_grid), intent(in) :: this
    integer,         intent(in) :: idx(:)
      !! edge indices (idx_dim)
    integer,         intent(in) :: idx_dir
      !! edge direction
    real                        :: len
      !! return edge length

    integer :: i, idx_sub(this%idx_dim)

    call this%to_subgrid(IDX_EDGE, idx_dir, idx, i, idx_sub)
    len = this%g(i)%get_len(idx_sub, idx_dir)
  end function

  function sum_grid_get_surf(this, idx, idx_dir) result(surf)
    !! get size of face
    class(sum_grid), intent(in) :: this
    integer,         intent(in) :: idx(:)
      !! face indices
    integer,         intent(in) :: idx_dir
      !! face direction
    real                        :: surf
      !! return size of face

    integer :: i, idx_sub(this%idx_dim)

    call this%to_subgrid(IDX_FACE, idx_dir, idx, i, idx_sub)
    surf = this%g(i)%get_surf(idx_sub, idx_dir)
  end function

  function sum_grid_get_vol(this, idx) result(vol)
    !! get cell volume
    class(sum_grid), intent(in) :: this
    integer,         intent(in) :: idx(:)
      !! cell indices
    real                        :: vol
      !! return cell volume

    integer :: i, idx_sub(this%idx_dim)

    call this%to_subgrid(IDX_CELL, 0, idx, i, idx_sub)
    vol = this%g(i)%get_vol(idx_sub)
  end function

  function sum_grid_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
    !! get maximal number of nearest neighbours
    class(sum_grid), intent(in) :: this
    integer,         intent(in) :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in) :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(in) :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in) :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer                     :: max_neighb
      !! return maximal number of nearest neighbours

    integer :: i

    ASSERT(this%idx_allowed(idx1_type, idx1_dir))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))

    max_neighb = maxval([(this%g(i)%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir), i = 1, size(this%g))])
  end function

  subroutine sum_grid_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
      !! get j-th neighbour.
      !!
      !! j: we count neighbors from 1,2,...,N.
      !! N: depends on idx1 (e.g. boundary nodes might have fewer neighbors).
      !! status: indicates if j-th neighbor exists
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in)  :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(in)  :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in)  :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(in)  :: idx1(:)
      !! first indices. size: (idx_dim)
    integer,         intent(in)  :: j
      !! j-th neighbor
    integer,         intent(out) :: idx2(:)
      !! output neighbour indices. size: (idx_dim)
    logical,         intent(out) :: status
      !! does j-th neighbor exist?

    integer :: i, idx1_sub(this%idx_dim), idx2_sub(this%idx_dim)

    ASSERT(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))
    ASSERT(size(idx2) == this%idx_dim)

    call this%to_subgrid(idx1_type, idx1_dir, idx1, i, idx1_sub)

    call this%g(i)%get_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir, idx1_sub, j, idx2_sub, status)
    if (status) idx2 = this%to_sumgrid(idx2_type, idx2_dir, i, idx2_sub)
  end subroutine

  subroutine sum_grid_to_subgrid(this, idx_type, idx_dir, idx_sum, i, idx_sub)
    !! transform index of sum_grid to sub-grid index
    class(sum_grid), intent(in)  :: this
    integer,         intent(in)  :: idx_type
      !! index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(in)  :: idx_sum(:)
      !! indices wrt to sum_grid
    integer,         intent(out) :: i
      !! output sub-grid index
    integer,         intent(out) :: idx_sub(:)
      !! output indices wrt sub-grid i

    integer :: idx_bnd(this%idx_dim)

    ASSERT(this%idx_allowed(idx_type, idx_dir, idx=idx_sum))
    ASSERT(size(idx_sub) == this%idx_dim)

    idx_sub = idx_sum
    do i = 1, size(this%g)
      call this%g(i)%get_idx_bnd(idx_type, idx_dir, idx_bnd)
      if (idx_sub(1) <= idx_bnd(1)) exit
      idx_sub = idx_sub - idx_bnd
    end do
  end subroutine

  function sum_grid_to_sumgrid(this, idx_type, idx_dir, i, idx_sub) result(idx_sum)
    !! transform index of sum_grid to sub-grid index
    class(sum_grid), intent(in) :: this
    integer,         intent(in) :: idx_type
      !! index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,         intent(in) :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,         intent(in) :: i
      !! sub-grid index
    integer,         intent(in) :: idx_sub(:)
      !! indices wrt sub-grid i
    integer                     :: idx_sum(this%idx_dim)
      !! output indices wrt to sum_grid

    integer :: j, idx_bnd(this%idx_dim)

    ASSERT(this%g(i)%idx_allowed(idx_type, idx_dir, idx=idx_sub))

    idx_sum = idx_sub
    do j = 1, i-1
      call this%g(j)%get_idx_bnd(idx_type, idx_dir, idx_bnd)
      idx_sum = idx_sub + idx_bnd
    end do
  end function
end module
