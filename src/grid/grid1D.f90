m4_include(../util/macro.f90.inc)

module grid1D_m

  use error_m,         only: assert_failed
  use json_m,          only: json_object
  use grid_m,          only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use normalization_m, only: denorm
  use output_file_m,   only: output_file

  implicit none

  private
  public grid1D

  type, extends(grid) :: grid1D
    !! 1D grid
    real, allocatable :: x(:)
      !! grid points
    integer           :: i0, i1
      !! lower/upper bounds of x
    integer           :: n
      !! size of x (n = i1 - i0 + 1)
  contains
    procedure :: init           => grid1D_init
    procedure :: get_idx_bnd_n  => grid1D_get_idx_bnd_n
    procedure ::                   grid1D_get_idx_bnd_1
    generic   :: get_idx_bnd    => grid1D_get_idx_bnd_1
    procedure :: get_vertex     => grid1D_get_vertex
    procedure :: get_edge       => grid1D_get_edge
    procedure :: get_face       => grid1D_get_face
    procedure :: get_cell       => grid1D_get_cell
    procedure :: get_len        => grid1D_get_len
    procedure :: get_surf       => grid1D_get_surf
    procedure :: get_vol        => grid1D_get_vol
    procedure :: get_max_neighb => grid1D_get_max_neighb
    procedure :: get_neighb     => grid1D_get_neighb
    procedure :: get_adjoint    => grid1D_get_adjoint
    procedure :: output         => grid1D_output
  end type

contains

  subroutine grid1D_init(this, name, x, i0)
    !! initialize 1D grid
    class(grid1D),     intent(out) :: this
    character(*),      intent(in)  :: name
      !! grid name
    real,              intent(in)  :: x(:)
      !! grid points
    integer, optional, intent(in)  :: i0
      !! lower bound of x

    ! init base
    call this%grid_init(name, 1, 1, [1], 2, [1])

    ! grid points
    this%i0 = 1
    if (present(i0)) this%i0 = i0
    this%n = size(x)
    this%i1 = this%n + this%i0 - 1
    allocate (this%x(this%i0:this%i1), source = x)
  end subroutine

  subroutine grid1D_get_idx_bnd_n(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:,:)
      !! output: lower/upper bound for each index. size: (2, idx_dim = 1)

    m4_assert(size(idx_bnd,1) == 2)
    m4_assert(size(idx_bnd,2) == 1)

    call this%get_idx_bnd(idx_type, idx_dir, idx_bnd(:,1))
  end subroutine

  subroutine grid1D_get_idx_bnd_1(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:)
      !! output: lower/upper bound for each index.

    m4_assert(size(idx_bnd) == 2)
    m4_assert(this%idx_allowed(idx_type, idx_dir))

    m4_ignore(idx_dir)

    idx_bnd(1) = this%i0
    select case (idx_type)
      case (IDX_VERTEX)
        idx_bnd(2) = this%i1
      case (IDX_EDGE)
        idx_bnd(2) = this%i1 - 1
      case (IDX_FACE)
        idx_bnd(2) = this%i1
      case (IDX_CELL)
        idx_bnd(2) = this%i1 - 1
    end select
  end subroutine

  subroutine grid1D_get_vertex(this, idx, p)
    !! get single vertex: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex grid indices (idx_dim)
    real,          intent(out) :: p(:)
      !! output: vertex coordinates (dim=1)

    m4_assert(this%idx_allowed(IDX_VERTEX, 0, idx=idx))
    m4_assert(size(p) == this%dim)

    p(1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_edge(this, idx, idx_dir, p)
    !! get single edge: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! edge indices (idx_dim)
    integer,       intent(in)  :: idx_dir
      !! edge direction
    real,          intent(out) :: p(:,:)
      !! output: edge coordinates (dim, 2)

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, 2]))

    m4_ignore(idx_dir)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  subroutine grid1D_get_face(this, idx, idx_dir, p)
    !! get single face: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! face indices (idx_dim)
    integer,       intent(in)  :: idx_dir
      !! face direction
    real,          intent(out) :: p(:,:)
      !! output: face coordinates (dim=1, face_nvert(idx_dir)=1)

    m4_assert(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%face_nvert(idx_dir)]))

    m4_ignore(idx_dir)

    p(1,1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_cell(this, idx, p)
    !! get single cell: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell indices (idx_dim)
    real,          intent(out) :: p(:,:)
      !! output: cell coordinates (dim=1, cell_nvert=2)

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%cell_nvert]))

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

  function grid1D_get_len(this, idx, idx_dir) result(len)
    !! get edge length
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! edge indices (idx_dim)
    integer,       intent(in) :: idx_dir
      !! edge direction
    real                      :: len
      !! return edge length

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))

    m4_ignore(idx_dir)

    len = this%x(idx(1)+1) - this%x(idx(1))
  end function

  function grid1D_get_surf(this, idx, idx_dir) result(surf)
    !! get single surface
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face indices. size: (idx_dim)
    integer,       intent(in) :: idx_dir
      !! face direction
    real                      :: surf
      !! output: size of face

    m4_assert(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))

    m4_ignore(this)
    m4_ignore(idx)
    m4_ignore(idx_dir)

    surf = 1.0
  end function

  function grid1D_get_vol(this, idx) result(vol)
    !! get single cell volume
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! cell indices. size: (idx_dim)
    real                      :: vol
      !! return cell volume

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))

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
    ! V 2 2 1 2
    ! E 2 2 2 1
    ! F 1 2 2 2
    ! C 2 1 2 2
    integer, parameter :: n(4,4) = reshape([2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 2, 2], [4, 4])

    m4_assert(this%idx_allowed(idx1_type, idx1_dir))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))

    m4_ignore(this)
    m4_ignore(idx1_dir)
    m4_ignore(idx2_dir)

    max_neighb = n(idx1_type,idx2_type)
  end function

  subroutine grid1D_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
    !! get j-th neighbour.
    !!
    !! j: we count neighbors from 1,2,...,N.
    !! N: depends on idx1 (e.g. boundary nodes might have fewer neighbors).
    !! status: indicates if j-th neighbor exists
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
      !! first indices. size: (idx_dim)
    integer,       intent(in)  :: j
      !! j-th neighbor
    integer,       intent(out) :: idx2(:)
      !! output neighbour indices. size: (idx_dim)
    logical,       intent(out) :: status
      !! does j-th neighbor exist?

    integer :: idx_bnd(2), shift

    m4_assert(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))
    m4_assert(size(idx2) == this%idx_dim)

    m4_ignore(idx1_dir)
    m4_ignore(idx2_dir)

    idx2  = this%i0 - 1
    shift = 0

    if (idx1_type == idx2_type) then
      ! idx1 is not a neighbour of itself => either subtract or add 1
      if (idx1(1) == this%i0) shift = 1
      if (j + shift == 1) then
        idx2(1) = idx1(1) - 1
      else if (j + shift == 2) then
        idx2(1) = idx1(1) + 1
      end if
    else if (((idx1_type == IDX_VERTEX) .and. (idx2_type == IDX_FACE  )) .or. &
      &      ((idx1_type == IDX_EDGE  ) .and. (idx2_type == IDX_CELL  )) .or. &
      &      ((idx1_type == IDX_FACE  ) .and. (idx2_type == IDX_VERTEX)) .or. &
      &      ((idx1_type == IDX_CELL  ) .and. (idx2_type == IDX_EDGE  ))) then
      ! vertices == faces and edges == cells
      if (j == 1) then
        idx2(1) = idx1(1)
      end if
    else if ((idx1_type == IDX_VERTEX) .or. (idx1_type == IDX_FACE)) then
      if (idx1(1) == this%i0) shift = 1
      if (j + shift == 1) then
        idx2(1) = idx1(1) - 1
      else if (j + shift == 2) then
        idx2(1) = idx1(1)
      end if
    else ! if ((idx1_type == IDX_EDGE) .or. (idx1_type == IDX_CELL)) then
      if (j == 1) then
        idx2(1) = idx1(1)
      else if (j == 2) then
        idx2(1) = idx1(1) + 1
      end if
    end if

    ! make sure idx2 is valid
    call this%get_idx_bnd(idx2_type, idx2_dir, idx_bnd)
    status = ((idx2(1) >= idx_bnd(1)) .and. (idx2(1) <= idx_bnd(2)))
  end subroutine

  subroutine grid1D_get_adjoint(this, idx, len, surf, vol)
    !! get adjoint grid information per cell
    class(grid1D),  intent(in)  :: this
    integer,        intent(in)  :: idx(:)
      !! cell indices. size: (idx_dim)
    real, optional, intent(out) :: len(:,:)
      !! edge lengths (max_cell_nedge, idx_dim)
    real, optional, intent(out) :: surf(:,:)
      !! adjoint surface parts per edge (max_cell_nedge, idx_dim)
    real, optional, intent(out) :: vol(:)
      !! adjoint volume parts per vertex (cell_nvert)

    m4_assert(size(idx) == 1)

    if (present(len)) then
      m4_assert(size(len,1) == 1)
      m4_assert(size(len,2) == 1)

      len(1,1) = this%get_len(idx, 1)
    end if

    if (present(surf)) then
      m4_assert(size(surf,1) == 1)
      m4_assert(size(surf,2) == 1)

      surf(1,1) = 1.0
    end if

    if (present(vol)) then
      m4_assert(size(vol) == 2)

      vol = 0.5 * this%get_len(idx, 1)
    end if
  end subroutine

  subroutine grid1D_output(this, of, unit)
    !! output 1D grid
    class(grid1D),          intent(in)    :: this
    type(output_file),      intent(inout) :: of
      !! output file handle
    character(*), optional, intent(in)    :: unit
      !! physical unit of coordinates; default = "um"

    character(:), allocatable  :: unit_
    type(json_object), pointer :: obj

    unit_ = "um"
    if (present(unit)) unit_ = unit

    obj => of%new_object("Grids")
    call obj%add_string("Name", this%name)
    call obj%add_string("Type", "1D")
    call of%write(obj, "Vertices", this%x, unit_)
  end subroutine

end module
