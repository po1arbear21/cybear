m4_include(../util/macro.f90.inc)

module grid0D_m

  use error_m,       only: assert_failed
  use grid_m,        only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use json_m,        only: json_object
  use output_file_m, only: output_file

  implicit none

  private
  public grid0D
  public get_dummy_grid

  type, extends(grid) :: grid0D
    !! 0D pseudo grid (consists of single vertex at x=0), can be used for global scalar variables
  contains
    procedure :: init           => grid0D_init
    procedure :: get_idx_bnd_n  => grid0D_get_idx_bnd_n
    procedure :: get_vertex     => grid0D_get_vertex
    procedure :: get_edge       => grid0D_get_edge
    procedure :: get_face       => grid0D_get_face
    procedure :: get_cell       => grid0D_get_cell
    procedure :: get_len        => grid0D_get_len
    procedure :: get_surf       => grid0D_get_surf
    procedure :: get_vol        => grid0D_get_vol
    procedure :: get_max_neighb => grid0D_get_max_neighb
    procedure :: get_neighb     => grid0D_get_neighb
    procedure :: output         => grid0D_output
  end type

  type(grid0D), target :: dum_grid

contains

  function get_dummy_grid() result(ptr)
    !! return pointer to module global grid0D.
    !! used in variable(grid=optional).
    type(grid0D), pointer :: ptr

    if (.not. allocated(dum_grid%face_dim)) call dum_grid%init("dummy")
    ptr => dum_grid
  end function

  subroutine grid0D_init(this, name)
    !! initialize 0D pseudo grid
    class(grid0D), intent(out) :: this
    character(*),  intent(in)  :: name
      !! grid name

    integer :: face_dim(0)

    ! init base
    call this%grid_init(name, 0, 0, face_dim, 0)
  end subroutine

  subroutine grid0D_get_idx_bnd_n(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (only IDX_VERTEX allowed)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:,:)
      !! output lower/upper bound for each index (2, idx_dim = 0)

    m4_assert(this%idx_allowed(idx_type, idx_dir))
    m4_assert(size(idx_bnd,1) == 2)
    m4_assert(size(idx_bnd,2) == 0)

    m4_ignore(this    )
    m4_ignore(idx_type)
    m4_ignore(idx_dir )
    idx_bnd(:,:) = 0
  end subroutine

  subroutine grid0D_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex indices
    real,          intent(out) :: p(:)
      !! output vertex coordinates (dim)

    m4_assert(size(idx) == 0)
    m4_assert(size(  p) == 0)

    m4_ignore(this)
    m4_ignore(idx )
    p(:) = 0.0
  end subroutine

  subroutine grid0D_get_edge(this, idx, idx_dir, p)
    !! get edge coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! edge indices
    integer,       intent(in)  :: idx_dir
      !! edge direction
    real,          intent(out) :: p(:,:)
      !! output edge coordinates (dim x 2)

    m4_assert(size(idx) == 0)
    m4_assert(size(p,1) == 0)

    m4_ignore(this   )
    m4_ignore(idx    )
    m4_ignore(idx_dir)
    p(:,:) = 0.0
  end subroutine

  subroutine grid0D_get_face(this, idx, idx_dir, p)
    !! get face coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! face indices
    integer,       intent(in)  :: idx_dir
      !! face direction
    real,          intent(out) :: p(:,:)
      !! output face coordinates (dim x face_dim(idx_dir))

    m4_assert(size(idx) == 0)
    m4_assert(size(p,1) == 0)

    m4_ignore(this   )
    m4_ignore(idx    )
    m4_ignore(idx_dir)
    p(:,:) = 0.0
  end subroutine

  subroutine grid0D_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell indices
    real,          intent(out) :: p(:,:)
      !! output cell coordinates (dim x cell_dim)

    m4_assert(size(idx) == 0)
    m4_assert(size(p,1) == 0)

    m4_ignore(this)
    m4_ignore(idx )
    p(:,:) = 0.0
  end subroutine

  function grid0D_get_len(this, idx, idx_dir) result(len)
    !! get edge length
    class(grid0D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! edge indices (idx_dim)
    integer,       intent(in) :: idx_dir
      !! edge direction
    real                      :: len
      !! return edge length

    m4_assert(size(idx) == 0)
    m4_assert(idx_dir   == 0)

    m4_ignore(this)
    m4_ignore(idx)
    m4_ignore(idx_dir)

    len = 0
  end function

  function grid0D_get_surf(this, idx, idx_dir) result(surf)
    !! get size of face
    class(grid0D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face indices
    integer,       intent(in) :: idx_dir
      !! face direction
    real                      :: surf
      !! return size of face

    m4_assert(size(idx) == 0)
    m4_assert(idx_dir   == 0)

    m4_ignore(this)
    m4_ignore(idx)
    m4_ignore(idx_dir)

    surf = 0
  end function

  function grid0D_get_vol(this, idx) result(vol)
    !! get cell volume
    class(grid0D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! cell indices
    real                      :: vol
      !! return cell volume

    m4_assert(size(idx) == 0)

    m4_ignore(this)
    m4_ignore(idx )

    vol = 0
  end function

  function grid0D_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
    !! get maximal number of nearest neighbours
    class(grid0D), intent(in) :: this
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

    m4_assert(idx1_dir == 0)
    m4_assert(idx2_dir == 0)

    m4_ignore(this     )
    m4_ignore(idx1_type)
    m4_ignore(idx1_dir )
    m4_ignore(idx2_type)
    m4_ignore(idx2_dir )

    max_neighb = 0
  end function

  subroutine grid0D_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
    !! get j-th neighbor
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(in)  :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(in)  :: idx1(:)
      !! first indices
    integer,       intent(in)  :: j
      !! j-th neighbor
    integer,       intent(out) :: idx2(:)
      !! output neighbour indices (idx_dim)
    logical,       intent(out) :: status
      !! does j-th neighb exist?

    m4_assert(idx1_dir   == 0)
    m4_assert(idx2_dir   == 0)
    m4_assert(size(idx1) == 0)
    m4_assert(size(idx2) == 0)

    m4_ignore(this     )
    m4_ignore(idx1_type)
    m4_ignore(idx1_dir )
    m4_ignore(idx2_type)
    m4_ignore(idx2_dir )
    m4_ignore(idx1     )
    m4_ignore(j        )

    idx2(:) = 0
    status  = .false.
  end subroutine

  subroutine grid0D_output(this, of, unit)
    !! output 0D grid
    class(grid0D),          intent(in)    :: this
    type(output_file),      intent(inout) :: of
      !! output file handle
    character(*), optional, intent(in)    :: unit
      !! physical unit of coordinates (ignored)

    type(json_object), pointer :: obj

    m4_ignore(unit)

    obj => of%new_object("Grids")
    call obj%add_string("Name", this%name)
    call obj%add_string("Type", "0D")
  end subroutine

end module
