#include "../util/macro.f90.inc"

module grid_m
  use error_m
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
    procedure(grid_get_surf),    deferred :: get_surf
    procedure(grid_get_vol),     deferred :: get_vol
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
        !! output upper bound for each index (idx_dim)
    end subroutine

    subroutine grid_get_vertex(this, idx, p)
      !! get vertex coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! vertex indices
      real,        intent(out) :: p(:)
        !! output vertex coordinates (dim)
    end subroutine

    subroutine grid_get_edge(this, idx, idx_dir, p)
      !! get edge coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! edge indices
      integer,     intent(in)  :: idx_dir
        !! edge direction
      real,        intent(out) :: p(:,:)
        !! output edge coordinates (dim x 2)
    end subroutine

    subroutine grid_get_face(this, idx, idx_dir, p)
      !! get face coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! face indices
      integer,     intent(in)  :: idx_dir
        !! face direction
      real,        intent(out) :: p(:,:)
        !! output face coordinates (dim x face_dim(idx_dir))
    end subroutine

    subroutine grid_get_cell(this, idx, p)
      !! get cell coordinates from grid indices
      import grid
      class(grid), intent(in)  :: this
      integer,     intent(in)  :: idx(:)
        !! cell indices
      real,        intent(out) :: p(:,:)
        !! output cell coordinates (dim x cell_dim)
    end subroutine

    function grid_get_surf(this, idx, idx_dir) result(surf)
      !! get size of face
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! face indices
      integer,     intent(in) :: idx_dir
        !! face direction
      real                    :: surf
        !! return size of face
    end function

    function grid_get_vol(this, idx) result(vol)
      !! get cell volume
      import grid
      class(grid), intent(in) :: this
      integer,     intent(in) :: idx(:)
        !! cell indices
      real                    :: vol
        !! return cell volume
    end function
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

end module
