#include "../util/macro.f90.inc"

module grid0D_m
  use error_m
  use grid_m
  implicit none

  type, extends(grid) :: grid0D
    !! 0D pseudo grid (consists of single vertex at x=0), can be used for global scalar variables
  contains
    procedure :: init           => grid0D_init
    procedure :: get_idx_bnd    => grid0D_get_idx_bnd
    procedure :: get_vertex     => grid0D_get_vertex
    procedure :: get_edge       => grid0D_get_edge
    procedure :: get_face       => grid0D_get_face
    procedure :: get_cell       => grid0D_get_cell
    procedure :: get_surf       => grid0D_get_surf
    procedure :: get_vol        => grid0D_get_vol
    procedure :: get_max_neighb => grid0D_get_max_neighb
    procedure :: get_neighb     => grid0D_get_neighb
  end type

contains

  subroutine grid0D_init(this)
    !! initialize 0D pseudo grid
    class(grid0D), intent(out) :: this

    ! init base
    call this%grid_init(1, 1, [0], 0)
  end subroutine

  subroutine grid0D_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (only IDX_VERTEX allowed)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:)
      !! output upper bound for each index (1)

    ASSERT(idx_type      == IDX_VERTEX)
    ASSERT(idx_dir       == 0         )
    ASSERT(size(idx_bnd) == 1         )

    IGNORE(this    )
    IGNORE(idx_type)
    IGNORE(idx_dir )

    idx_bnd(1) = 1
  end subroutine

  subroutine grid0D_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex indices
    real,          intent(out) :: p(:)
      !! output vertex coordinates (dim)

    ASSERT(size(idx) == 1)
    ASSERT(idx(1)    == 1)
    ASSERT(size(p  ) == 1)

    IGNORE(this)
    IGNORE(idx )

    p(1) = 0
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

    IGNORE(this   )
    IGNORE(idx    )
    IGNORE(idx_dir)
    IGNORE(p      )

    call program_error("0D Grid does not have edges")
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

    IGNORE(this   )
    IGNORE(idx    )
    IGNORE(idx_dir)
    IGNORE(p      )

    call program_error("0D Grid does not have faces")
  end subroutine

  subroutine grid0D_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(grid0D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell indices
    real,          intent(out) :: p(:,:)
      !! output cell coordinates (dim x cell_dim)

    IGNORE(this)
    IGNORE(idx )
    IGNORE(p   )

    call program_error("0D Grid does not have cells")
  end subroutine

  function grid0D_get_surf(this, idx, idx_dir) result(surf)
    !! get size of face
    class(grid0D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face indices
    integer,       intent(in) :: idx_dir
      !! face direction
    real                      :: surf
      !! return size of face

    IGNORE(this)
    IGNORE(idx)
    IGNORE(idx_dir)

    call program_error("0D Grid does not have faces")

    surf = 0
  end function

  function grid0D_get_vol(this, idx) result(vol)
    !! get cell volume
    class(grid0D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! cell indices
    real                      :: vol
      !! return cell volume

    IGNORE(this)
    IGNORE(idx )

    call program_error("0D Grid does not have cells")

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

    IGNORE(this     )
    IGNORE(idx1_type)
    IGNORE(idx1_dir )
    IGNORE(idx2_type)
    IGNORE(idx2_dir )

    max_neighb = 0
  end function

  subroutine grid0D_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, idx2, nidx2)
    !! get nearest neighbours
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
    integer,       intent(out) :: idx2(:,:)
      !! output neighbour indices (idx_dim x max_neighb)
    integer,       intent(out) :: nidx2
      !! output actual number of neighburs

    IGNORE(this     )
    IGNORE(idx1_type)
    IGNORE(idx1_dir )
    IGNORE(idx2_type)
    IGNORE(idx2_dir )
    IGNORE(idx1     )

    idx2  = 0
    nidx2 = 0
  end subroutine

end module
