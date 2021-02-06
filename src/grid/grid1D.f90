#include "../util/macro.f90.inc"

module grid1D_m
  use error_m
  use grid_m
  implicit none

  type, extends(grid) :: grid1D
    !! 1D grid
    real, allocatable :: x(:)
      !! grid points
  contains
    procedure :: init        => grid1D_init
    procedure :: get_idx_bnd => grid1D_get_idx_bnd
    procedure :: get_vertex  => grid1D_get_vertex
    procedure :: get_edge    => grid1D_get_edge
    procedure :: get_face    => grid1D_get_face
    procedure :: get_cell    => grid1D_get_cell
    procedure :: get_surf    => grid1D_get_surf
    procedure :: get_vol     => grid1D_get_vol
  end type

contains

  subroutine grid1D_init(this, x)
    !! initialize 1D grid
    class(grid1D), intent(out) :: this
    real,          intent(in)  :: x(:)
      !! grid points

    ! init base
    call this%grid_init(1, 1, [1], 2)

    ! save grid points
    this%x = x
  end subroutine

  subroutine grid1D_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:)
      !! output upper bound for each index (1)

    ASSERT(size(idx_bnd) == 1)

    select case (idx_type)
      case (IDX_VERTEX)
        ASSERT(idx_dir == 0)
        idx_bnd(1) = size(this%x)
      case (IDX_EDGE)
        ASSERT(idx_dir == 1)
        idx_bnd(1) = size(this%x) - 1
      case (IDX_FACE)
        ASSERT(idx_dir == 1)
        idx_bnd(1) = size(this%x)
      case (IDX_CELL)
        ASSERT(idx_dir == 0)
        idx_bnd(1) = size(this%x) - 1
    end select
  end subroutine

  subroutine grid1D_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex indices
    real,          intent(out) :: p(:)
      !! output vertex coordinates (dim)

    ASSERT(size(idx) == 1)
    ASSERT(size(p  ) == 1)

    p(1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_edge(this, idx, idx_dir, p)
    !! get edge coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! edge indices
    integer,       intent(in)  :: idx_dir
      !! edge direction
    real,          intent(out) :: p(:,:)
      !! output edge coordinates (dim x 2)

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 2)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  subroutine grid1D_get_face(this, idx, idx_dir, p)
    !! get face coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! face indices
    integer,       intent(in)  :: idx_dir
      !! face direction
    real,          intent(out) :: p(:,:)
      !! output face coordinates (dim x face_dim(idx_dir))

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 1)

    p(1,1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell indices
    real,          intent(out) :: p(:,:)
      !! output cell coordinates (dim x cell_dim)

    ASSERT(size(idx) == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 2)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  function grid1D_get_surf(this, idx, idx_dir) result(surf)
    !! get size of face
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face indices
    integer,       intent(in) :: idx_dir
      !! face direction
    real                      :: surf
      !! return size of face

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)

    IGNORE(this)
    IGNORE(idx)
    IGNORE(idx_dir)

    surf = 1.0
  end function

  function grid1D_get_vol(this, idx) result(vol)
    !! get cell volume
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! cell indices
    real                      :: vol
      !! return cell volume

    ASSERT(size(idx) == 1)

    vol = this%x(idx(1)+1) - this%x(idx(1))
  end function

end module