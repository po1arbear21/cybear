#include "../util/macro.f90.inc"

module grid1D_m
  use error_m
  use grid_m, only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  implicit none

  private
  public grid1D

  type, extends(grid) :: grid1D
    !! 1D grid
    real, allocatable :: x(:)
      !! grid points
  contains
    procedure :: init           => grid1D_init
    procedure :: get_idx_bnd    => grid1D_get_idx_bnd
    procedure :: get_vertex     => grid1D_get_vertex
    procedure :: get_edge       => grid1D_get_edge
    procedure :: get_face       => grid1D_get_face
    procedure :: get_cell       => grid1D_get_cell
    procedure :: get_surf       => grid1D_get_surf
    procedure :: get_vol        => grid1D_get_vol
    procedure :: get_max_neighb => grid1D_get_max_neighb
    procedure :: get_neighb     => grid1D_get_neighb
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
      !! output: upper bound for each index. size: (1)

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
    !! get single vertex: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex' grid indices. size: (idx_dim)
    real,          intent(out) :: p(:)
      !! output: vertex' coordinates. size: (dim)

    ASSERT(size(idx) == 1)
    ASSERT(size(p  ) == 1)

    p(1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_edge(this, idx, idx_dir, p)
    !! get single edge: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! edge's indices. size: (idx_dim)
    integer,       intent(in)  :: idx_dir
      !! edge's direction
    real,          intent(out) :: p(:,:)
      !! output: edge's coordinates. size: (dim, 2)

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 2)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  subroutine grid1D_get_face(this, idx, idx_dir, p)
    !! get single face: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! face's indices. size: (idx_dim)
    integer,       intent(in)  :: idx_dir
      !! face's direction
    real,          intent(out) :: p(:,:)
      !! output: face's coordinates. size: (dim, face_dim(idx_dir))

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 1)

    p(1,1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_cell(this, idx, p)
    !! get single cell: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell's indices. size: (idx_dim)
    real,          intent(out) :: p(:,:)
      !! output: cell's coordinates. size: (dim, cell_dim)

    ASSERT(size(idx) == 1)
    ASSERT(size(p,1) == 1)
    ASSERT(size(p,2) == 2)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  function grid1D_get_surf(this, idx, idx_dir) result(surf)
    !! get single face's surface
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face's indices. size: (idx_dim)
    integer,       intent(in) :: idx_dir
      !! face's direction
    real                      :: surf
      !! output: size of face

    ASSERT(size(idx) == 1)
    ASSERT(idx_dir   == 1)

    IGNORE(this)
    IGNORE(idx)
    IGNORE(idx_dir)

    surf = 1.0
  end function

  function grid1D_get_vol(this, idx) result(vol)
    !! get single cell's volume
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! cell's indices. size: (idx_dim)
    real                      :: vol
      !! return cell volume

    ASSERT(size(idx) == 1)

    vol = this%x(idx(1)+1) - this%x(idx(1))
  end function

  function grid1D_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
    !! get maximal number of nearest neighbours
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in) :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(in) :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in) :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer                   :: max_neighb
      !! return maximal number of nearest neighbours

    !   V E F C
    ! V 3 2 1 2
    ! E 2 3 2 1
    ! F 1 2 3 2
    ! C 2 1 2 3
    integer, parameter :: n(4,4) = reshape([3, 2, 1, 2, 2, 3, 2, 1, 1, 2, 3, 2, 2, 1, 2, 3], [4, 4])

    ASSERT((((idx1_type == IDX_VERTEX) .or. (idx1_type == IDX_CELL)) .and. (idx1_dir == 0)) \
      .or. (((idx1_type == IDX_EDGE  ) .or. (idx1_type == IDX_FACE)) .and. (idx1_dir == 1)))
    ASSERT((((idx2_type == IDX_VERTEX) .or. (idx2_type == IDX_CELL)) .and. (idx2_dir == 0)) \
      .or. (((idx2_type == IDX_EDGE  ) .or. (idx2_type == IDX_FACE)) .and. (idx2_dir == 1)))

    max_neighb = n(idx1_type,idx2_type)
  end function

  subroutine grid1D_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, idx2, nidx2)
    !! get nearest neighbours
    !!
    !! output might not be fully set. consider boundary nodes which have fewer neighbours.
    !!    idx2(:,nidx2+1:) will be empty.
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(in)  :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(in)  :: idx1(:)
      !! first's indices. size: (idx_dim)
    integer,       intent(out) :: idx2(:,:)
      !! output: neighbours' indices. size: (idx_dim, max_neighb)
    integer,       intent(out) :: nidx2
      !! output: actual number of neighbours.

    integer :: max_neighb

    ASSERT((((idx1_type == IDX_VERTEX) .or. (idx1_type == IDX_CELL)) .and. (idx1_dir == 0)) \
      .or. (((idx1_type == IDX_EDGE  ) .or. (idx1_type == IDX_FACE)) .and. (idx1_dir == 1)))
    ASSERT((((idx2_type == IDX_VERTEX) .or. (idx2_type == IDX_CELL)) .and. (idx2_dir == 0)) \
      .or. (((idx2_type == IDX_EDGE  ) .or. (idx2_type == IDX_FACE)) .and. (idx2_dir == 1)))

    max_neighb = this%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir)

    select case (max_neighb)
      case (1)
        nidx2 = 1
        idx2(1,1) = idx1(1)

      case (2)
        if ((idx1_type == IDX_VERTEX) .or. (idx1_type == IDX_FACE)) then
          nidx2 = 0
          if (idx1(1) > 1) then
            nidx2 = nidx2 + 1
            idx2(1,nidx2) = idx1(1) - 1
          end if
          if (idx1(1) < size(this%x)) then
            nidx2 = nidx2 + 1
            idx2(1,nidx2) = idx1(1)
          end if
        else
          nidx2 = 2
          idx2(1,1) = idx1(1)
          idx2(1,2) = idx1(1) + 1
        end if

      case (3)
        nidx2 = 0
        if (idx1(1) > 1) then
          nidx2 = nidx2 + 1
          idx2(1,nidx2) = idx1(1) - 1
        end if
        nidx2 = nidx2 + 1
        idx2(1,nidx2) = idx1(1)
        if (idx1(1) < size(this%x)) then
          nidx2 = nidx2 + 1
          idx2(1,nidx2) = idx1(1) + 1
        end if
    end select
  end subroutine

end module
