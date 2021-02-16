#include "../util/macro.f90.inc"

module grid_m
  use error_m
  implicit none

  private
  public grid
  public grid_ptr
  public IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL

  ! grid index types
  integer, parameter :: IDX_VERTEX = 1
  integer, parameter :: IDX_EDGE   = 2
  integer, parameter :: IDX_FACE   = 3
  integer, parameter :: IDX_CELL   = 4

  type, abstract :: grid
    !! Base grid

    integer              :: dim
      !! grid dimension = number of coordinates per point
    integer              :: idx_dim
      !! index dimension = number of indices per point
    integer, allocatable :: face_dim(:)
      !! number of points per face depending on direction (1:idx_dim)
    integer              :: cell_dim
      !! number of points per cell
  contains
    procedure                                :: grid_init
    procedure                                :: idx_allowed => grid_idx_allowed
    procedure(grid_get_idx_bnd),    deferred :: get_idx_bnd
    procedure(grid_get_vertex),     deferred :: get_vertex
    procedure(grid_get_edge),       deferred :: get_edge
    procedure(grid_get_face),       deferred :: get_face
    procedure(grid_get_cell),       deferred :: get_cell
    procedure(grid_get_surf),       deferred :: get_surf
    procedure(grid_get_vol),        deferred :: get_vol
    procedure(grid_get_max_neighb), deferred :: get_max_neighb
    procedure(grid_get_neighb),     deferred :: get_neighb
  end type

  type grid_ptr
    class(grid), pointer :: p => null()
  end type

  abstract interface

    subroutine grid_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
      !! get grid index bounds
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx_type
        !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in)  :: idx_dir
        !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer,     intent(out) :: idx_bnd(:)
        !! output: upper bound for each index. size: (idx_dim)
    end subroutine

    subroutine grid_get_vertex(this, idx, p)
      !! get single vertex: from grid indices to coordinates
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! vertex' grid indices. size: (idx_dim)
      real,        intent(out) :: p(:)
        !! output: vertex' coordinates. size: (dim)
    end subroutine

    subroutine grid_get_edge(this, idx, idx_dir, p)
      !! get single edge: from grid indices to coordinates
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! edge's indices. size: (idx_dim)
      integer,     intent(in)  :: idx_dir
        !! edge's direction
      real,        intent(out) :: p(:,:)
        !! output: edge's coordinates. size: (dim, 2)
    end subroutine

    subroutine grid_get_face(this, idx, idx_dir, p)
      !! get single face: from grid indices to coordinates
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! face's indices. size: (idx_dim)
      integer,     intent(in)  :: idx_dir
        !! face's direction
      real,        intent(out) :: p(:,:)
        !! output: face's coordinates. size: (dim, face_dim(idx_dir))
    end subroutine

    subroutine grid_get_cell(this, idx, p)
      !! get single cell: from grid indices to coordinates
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! cell's indices. size: (idx_dim)
      real,        intent(out) :: p(:,:)
        !! output: cell's coordinates. size: (dim, cell_dim)
    end subroutine

    function grid_get_surf(this, idx, idx_dir) result(surf)
      !! get single face's surface
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! face's indices. size: (idx_dim)
      integer,     intent(in) :: idx_dir
        !! face's direction
      real                    :: surf
        !! output: size of face
    end function

    function grid_get_vol(this, idx) result(vol)
      !! get single cell's volume
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! cell's indices. size: (idx_dim)
      real                    :: vol
        !! return cell volume
    end function

    function grid_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
      !! get maximal number of nearest neighbours
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx1_type
        !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in) :: idx1_dir
        !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer,     intent(in) :: idx2_type
        !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in) :: idx2_dir
        !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer                 :: max_neighb
        !! return maximal number of nearest neighbours
    end function

    subroutine grid_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, idx2, nidx2)
      !! get nearest neighbours
      !!
      !! output might not be fully set. consider boundary nodes which have fewer neighbours.
      !!    idx2(:,nidx2+1:) will be empty.
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx1_type
        !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in)  :: idx1_dir
        !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer,     intent(in)  :: idx2_type
        !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in)  :: idx2_dir
        !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer,     intent(in)  :: idx1(:)
        !! first indices. size: (idx_dim)
      integer,     intent(out) :: idx2(:,:)
        !! output neighbour indices. size: (idx_dim, max_neighb)
      integer,     intent(out) :: nidx2
        !! output actual number of neighburs
    end subroutine
  end interface

contains

  subroutine grid_init(this, dim, idx_dim, face_dim, cell_dim)
    !! initialize grid
    class(grid), intent(out) :: this
    integer,     intent(in)  :: dim
      !! grid dimension = number of coordinates per point
    integer,     intent(in)  :: idx_dim
      !! index dimension = number of indices per point
    integer,     intent(in)  :: face_dim(:)
      !! number of points per face depending on direction (idx_dim)
    integer,     intent(in)  :: cell_dim
      !! number of points per cell

    ASSERT(          dim  >  0      )
    ASSERT(      idx_dim  >  0      )
    ASSERT(size(face_dim) == idx_dim)
    ASSERT(all( face_dim  >= 0)     )
    ASSERT(     cell_dim  >= 0      )

    this%dim      =      dim
    this%idx_dim  =  idx_dim
    this%face_dim = face_dim
    this%cell_dim = cell_dim
  end subroutine

  logical function grid_idx_allowed(this, idx_type, idx_dir, idx) result(allowed)
    !! checks if given indices are allowed for grid.
    class(grid), intent(in)           :: this
    integer,     intent(in)           :: idx_type
      !! index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,     intent(in)           :: idx_dir
      !! index direction.
    integer,     intent(in), optional :: idx(:)
      !! index. size: (idx_dim)

    integer :: idx_bnd(this%idx_dim)

    if ((idx_type == IDX_VERTEX) .or. (idx_type == IDX_CELL)) then
      allowed = (idx_dir == 0)

    else if ((idx_type == IDX_EDGE  ) .or. (idx_type == IDX_FACE)) then
      allowed = ((idx_dir > 0) .and. (idx_dir <= this%idx_dim))
    end if

    if (present(idx)) then
      allowed = allowed .and. (size(idx) == this%idx_dim) .and. all(idx > 0)

      call this%get_idx_bnd(idx_type, idx_dir, idx_bnd)
      allowed = allowed .and. all(idx <= idx_bnd)
    end if
  end function

end module
