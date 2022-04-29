m4_include(../util/macro.f90.inc)

module grid_m

  use error_m,       only: assert_failed, program_error
  use output_file_m, only: output_file

  implicit none

  private
  public IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  public grid, grid_ptr

  ! grid index types
  integer, parameter :: IDX_VERTEX = 1
  integer, parameter :: IDX_EDGE   = 2
  integer, parameter :: IDX_FACE   = 3
  integer, parameter :: IDX_CELL   = 4

  ! grid index type names
  character(6), parameter :: IDX_NAME(4) = ["Vertex", "Edge  ", "Face  ", "Cell  "]

  type, abstract :: grid
    !! Base grid

    character(:), allocatable :: name
      !! grid name
    integer                   :: dim
      !! grid dimension = number of coordinates per point
    integer                   :: idx_dim
      !! index dimension = number of indices per point
    integer,      allocatable :: face_dim(:)
      !! number of points per face depending on direction (1:idx_dim)
    integer                   :: cell_dim
      !! number of points per cell
  contains
    procedure :: grid_init
    procedure :: get_ptr     => grid_get_ptr
    procedure :: idx_allowed => grid_idx_allowed
    generic   :: get_idx_bnd => get_idx_bnd_n

    procedure(grid_get_idx_bnd_n),  deferred :: get_idx_bnd_n
    procedure(grid_get_vertex),     deferred :: get_vertex
    procedure(grid_get_edge),       deferred :: get_edge
    procedure(grid_get_face),       deferred :: get_face
    procedure(grid_get_cell),       deferred :: get_cell
    procedure(grid_get_len),        deferred :: get_len
    procedure(grid_get_surf),       deferred :: get_surf
    procedure(grid_get_vol),        deferred :: get_vol
    procedure(grid_get_max_neighb), deferred :: get_max_neighb
    procedure(grid_get_neighb),     deferred :: get_neighb
    procedure(grid_output),         deferred :: output
  end type

  type grid_ptr
    class(grid), pointer :: p => null()
  end type

  abstract interface
    subroutine grid_get_idx_bnd_n(this, idx_type, idx_dir, idx_bnd)
      !! get grid index bounds
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx_type
        !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
      integer,     intent(in)  :: idx_dir
        !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
      integer,     intent(out) :: idx_bnd(:,:)
        !! output: lower/upper bound for each index. size: (2, idx_dim)
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

    function grid_get_len(this, idx, idx_dir) result(len)
      !! get edge length
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! edge indices (idx_dim)
      integer,     intent(in) :: idx_dir
        !! edge direction
      real                    :: len
        !! return edge length
    end function

    function grid_get_surf(this, idx, idx_dir) result(surf)
      !! get face area
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! face indices (idx_dim)
      integer,     intent(in) :: idx_dir
        !! face direction
      real                    :: surf
        !! return face area
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

    subroutine grid_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
      !! get j-th neighbour (does not include idx1!).
      !!
      !! j: we count neighbors from 1,2,...,N.
      !! N: depends on idx1 (e.g. boundary nodes might have fewer neighbors).
      !! status: indicates if j-th neighbor exists
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
      integer,     intent(in)  :: j
        !! j-th neighbor
      integer,     intent(out) :: idx2(:)
        !! output neighbour indices. size: (idx_dim)
      logical,     intent(out) :: status
        !! does j-th neighbor exist?
    end subroutine

    subroutine grid_output(this, of, unit)
      !! output grid
      import grid, output_file
      class(grid),            intent(in)    :: this
      type(output_file),      intent(inout) :: of
        !! output file handle
      character(*), optional, intent(in)    :: unit
        !! physical unit of coordinates
    end subroutine
  end interface

contains

  subroutine grid_init(this, name, dim, idx_dim, face_dim, cell_dim)
    !! initialize grid
    class(grid),  intent(out) :: this
    character(*), intent(in) :: name
      !! grid name
    integer,      intent(in)  :: dim
      !! grid dimension = number of coordinates per point
    integer,      intent(in)  :: idx_dim
      !! index dimension = number of indices per point
    integer,      intent(in)  :: face_dim(:)
      !! number of points per face depending on direction (idx_dim)
    integer,      intent(in)  :: cell_dim
      !! number of points per cell

    m4_assert(          dim  >= 0      )
    m4_assert(      idx_dim  >= 0      )
    m4_assert(size(face_dim) == idx_dim)
    m4_assert(all( face_dim  >= 0)     )
    m4_assert(     cell_dim  >= 0      )

    this%name     =     name
    this%dim      =      dim
    this%idx_dim  =  idx_dim
    this%face_dim = face_dim
    this%cell_dim = cell_dim
  end subroutine

  function grid_get_ptr(this) result(ptr)
    !! returns pointer type to this grid
    class(grid), target, intent(in) :: this
    type(grid_ptr)                  :: ptr

    ptr%p => this
  end function

  recursive function grid_idx_allowed(this, idx_type, idx_dir, idx) result(allowed)
    !! checks if given indices are allowed for grid.
    class(grid), intent(in)           :: this
    integer,     intent(in)           :: idx_type
      !! index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,     intent(in)           :: idx_dir
      !! index direction.
    integer,     intent(in), optional :: idx(:)
      !! index. size: (idx_dim)
    logical                           :: allowed

    integer :: idx_bnd(2,this%idx_dim)

    ! check idx_type, idx_dir
    select case (idx_type)
      case (IDX_VERTEX, IDX_CELL)
        allowed = (idx_dir == 0)
      case (IDX_EDGE, IDX_FACE)
        allowed = ((idx_dir > 0) .and. (idx_dir <= this%idx_dim))
      case default
        allowed = .false.
    end select

    ! check idx
    if (present(idx)) then
      allowed = allowed .and. (size(idx) == this%idx_dim)

      call this%get_idx_bnd(idx_type, idx_dir, idx_bnd)
      allowed = allowed .and. all(idx >= idx_bnd(1,:)) .and. all(idx <= idx_bnd(2,:))
    end if
  end function

end module
