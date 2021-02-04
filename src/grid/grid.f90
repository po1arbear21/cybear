module grid_m
  implicit none

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
    procedure                             :: grid_init
    procedure(grid_get_idx_bnd), deferred :: get_idx_bnd
    procedure(grid_get_vertex),  deferred :: get_vertex
    procedure(grid_get_edge),    deferred :: get_edge
    procedure(grid_get_face),    deferred :: get_face
    procedure(grid_get_cell),    deferred :: get_cell
  end type

  type grid_ptr
    class(grid), pointer :: p => null()
  end type

  abstract interface

    subroutine grid_get_idx_bnd(this, idx_type, dir, idx_bnd)
      !! get grid index bounds
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx_type
        !! grid index type (e.g. IDX_VERTEX)
      integer,     intent(in)  :: dir
        !! index of direction (only used for IDX_EDGE and IDX_FACE)
      integer,     intent(out) :: idx_bnd(:)
        !! output upper bound for each index (1:idx_dim)
    end subroutine

    subroutine grid_get_vertex(this, idx, p)
      !! get vertex coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! vertex indices
      real,        intent(out) :: p(:)
        !! output vertex coordinates (1:dim)
    end subroutine

    subroutine grid_get_edge(this, idx, dir, p)
      !! get edge coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! edge indices
      integer,     intent(in)  :: dir
        !! edge direction
      real,        intent(out) :: p(:,:)
        !! output edge coordinates (1:dim x 1:2)
    end subroutine

    subroutine grid_get_face(this, idx, dir, p)
      !! get face coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! face indices
      integer,     intent(in)  :: dir
        !! face direction
      real,        intent(out) :: p(:,:)
        !! output face coordinates (1:dim x 1:face_dim(dir))
    end subroutine

    subroutine grid_get_cell(this, idx, p)
      !! get cell coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! cell indices
      real,        intent(out) :: p(:,:)
        !! output cell coordinates (1:dim x 1:cell_dim)
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
      !! number of points per face depending on direction (1:idx_dim)
    integer,     intent(in)  :: cell_dim
      !! number of points per cell

    this%dim      = dim
    this%idx_dim  = idx_dim
    this%face_dim = face_dim
    this%cell_dim = cell_dim
  end subroutine

end module
