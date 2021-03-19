#include "../util/macro.f90.inc"

module grid_m

  use error_m, only: assert_failed, program_error

  implicit none

  private
  public IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  public grid, grid_ptr
  public allocate_grid_data
  ! public (see grid_data_def.f90.inc)
  public grid_table, grid_table_ptr

  ! grid index types
  integer, parameter :: IDX_VERTEX = 1
  integer, parameter :: IDX_EDGE   = 2
  integer, parameter :: IDX_FACE   = 3
  integer, parameter :: IDX_CELL   = 4

  ! grid index type names
  character(6), parameter :: IDX_NAME(4) = ["Vertex", "Edge  ", "Face  ", "Cell  "]

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

    type(grid_table), allocatable :: tab_all(:,:)
      !! select all grid indices (idx_type, idx_dir)
  contains
    procedure                                :: grid_init
    procedure                                :: init_tab_all => grid_init_tab_all
    procedure                                :: get_ptr      => grid_get_ptr
    procedure                                :: idx_allowed  => grid_idx_allowed
    procedure(grid_get_idx_bnd),    deferred :: get_idx_bnd
    procedure(grid_get_vertex),     deferred :: get_vertex
    procedure(grid_get_edge),       deferred :: get_edge
    procedure(grid_get_face),       deferred :: get_face
    procedure(grid_get_cell),       deferred :: get_cell
    procedure(grid_get_len),        deferred :: get_len
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
  end interface

#define T int
#define TT integer
#include "grid_data_def.f90.inc"

#define T log
#define TT logical
#define TLOG
#include "grid_data_def.f90.inc"

#define T real
#define TT real
#include "grid_data_def.f90.inc"

#define T cmplx
#define TT complex
#include "grid_data_def.f90.inc"

  type grid_table
    !! table to select points from grid

    character(:), allocatable :: name
      !! grid table name

    class(grid), pointer :: g => null()
      !! pointer to grid

    integer             :: idx_type
      !! index type
    integer             :: idx_dir
      !! index direction for edges and faces
    class(grid_data_log), allocatable :: flags
      !! include/exclude point (product idx_bnd)

    integer                           :: n
      !! number of entries
    integer,              allocatable :: flat2idx(:,:)
      !! flat index to grid indices (idx_dim x n)
    class(grid_data_int), allocatable :: idx2flat
      !! grid indices to flat index
  contains
    procedure :: init       => grid_table_init
    procedure :: init_final => grid_table_init_final
    procedure :: get_idx    => grid_table_get_idx
    procedure :: get_flat   => grid_table_get_flat
    procedure :: get_ptr    => grid_table_get_ptr
  end type

  type grid_table_ptr
    type(grid_table), pointer :: p => null()
  end type

  interface
    module subroutine grid_table_init(this, name, g, idx_type, idx_dir, initial_flags)
      class(grid_table),   intent(out) :: this
      character(*),        intent(in)  :: name
      class(grid), target, intent(in)  :: g
      integer,             intent(in)  :: idx_type
      integer,             intent(in)  :: idx_dir
      logical, optional,   intent(in)  :: initial_flags
    end subroutine

    module subroutine grid_table_init_final(this)
      class(grid_table), intent(inout) :: this
    end subroutine

    module function grid_table_get_idx(this, i) result(idx)
      class(grid_table), intent(in)  :: this
      integer,           intent(in)  :: i
      integer                        :: idx(this%g%idx_dim)
    end function

    module function grid_table_get_flat(this, idx) result(i)
      class(grid_table), intent(in)  :: this
      integer,           intent(in)  :: idx(:)
      integer                        :: i
    end function

    module function grid_table_get_ptr(this) result(ptr)
      class(grid_table), target, intent(in) :: this
      type(grid_table_ptr)                  :: ptr
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

    ASSERT(          dim  >= 0      )
    ASSERT(      idx_dim  >= 0      )
    ASSERT(size(face_dim) == idx_dim)
    ASSERT(all( face_dim  >= 0)     )
    ASSERT(     cell_dim  >= 0      )

    this%dim      =      dim
    this%idx_dim  =  idx_dim
    this%face_dim = face_dim
    this%cell_dim = cell_dim
  end subroutine

  subroutine grid_init_tab_all(this)
    !! initialize this%tab_all
    class(grid), intent(inout) :: this

    integer       :: idx_type, idx_dir, idir0(4), idir1(4)
    character(32) :: name

    ! initialize tables
    allocate (this%tab_all(4,0:this%idx_dim))
    idir0 = [0,       1,       1, 0]
    idir1 = [0, this%idx_dim, this%idx_dim, 0]
    do idx_type = 1, 4
      do idx_dir = idir0(idx_type), idir1(idx_type)
        if (idx_dir > 0) then
          write (name, "(A,A,I0)") "All", trim(IDX_NAME(idx_type)), idx_dir
        else
          write (name, "(A,A)") "All", trim(IDX_NAME(idx_type))
        end if
        call this%tab_all(idx_type,idx_dir)%init(trim(name), this, idx_type, idx_dir, initial_flags = .true.)
        call this%tab_all(idx_type,idx_dir)%init_final()
      end do
    end do
  end subroutine

  function grid_get_ptr(this) result(ptr)
    !! returns pointer type to this grid
    class(grid), target, intent(in) :: this
    type(grid_ptr)                  :: ptr

    ptr%p => this
  end function

  function grid_idx_allowed(this, idx_type, idx_dir, idx) result(allowed)
    !! checks if given indices are allowed for grid.
    class(grid), intent(in)           :: this
    integer,     intent(in)           :: idx_type
      !! index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,     intent(in)           :: idx_dir
      !! index direction.
    integer,     intent(in), optional :: idx(:)
      !! index. size: (idx_dim)
    logical                           :: allowed

    integer :: idx_bnd(this%idx_dim)

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
      allowed = allowed .and. (size(idx) == this%idx_dim) .and. all(idx > 0)

      call this%get_idx_bnd(idx_type, idx_dir, idx_bnd)
      allowed = allowed .and. all(idx <= idx_bnd)
    end if
  end function

end module
