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
    procedure :: get_len        => grid1D_get_len
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

    ! init tables
    call this%init_tab_all()
  end subroutine

  subroutine grid1D_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,       intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,       intent(out) :: idx_bnd(:)
      !! output: upper bound for each index. size: (idx_dim=1)

    ASSERT(this%idx_allowed(idx_type, idx_dir))
    ASSERT(size(idx_bnd) == 1)

    IGNORE(idx_dir)

    select case (idx_type)
      case (IDX_VERTEX)
        idx_bnd(1) = size(this%x)
      case (IDX_EDGE)
        idx_bnd(1) = size(this%x) - 1
      case (IDX_FACE)
        idx_bnd(1) = size(this%x)
      case (IDX_CELL)
        idx_bnd(1) = size(this%x) - 1
    end select
  end subroutine

  subroutine grid1D_get_vertex(this, idx, p)
    !! get single vertex: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex' grid indices. size: (idx_dim)
    real,          intent(out) :: p(:)
      !! output: vertex' coordinates. size: (dim=1)

    ASSERT(this%idx_allowed(IDX_VERTEX, 0, idx=idx))
    ASSERT(size(p) == this%dim)

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

    ASSERT(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))
    ASSERT(all(shape(p) == [this%dim, 2]))

    IGNORE(idx_dir)

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
      !! output: face's coordinates. size: (dim=1, face_dim(idx_dir)=1)

    ASSERT(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))
    ASSERT(all(shape(p) == [this%dim, this%face_dim(idx_dir)]))

    IGNORE(idx_dir)

    p(1,1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_cell(this, idx, p)
    !! get single cell: from grid indices to coordinates
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell's indices. size: (idx_dim)
    real,          intent(out) :: p(:,:)
      !! output: cell's coordinates. size: (dim=1, cell_dim=2)

    ASSERT(this%idx_allowed(IDX_CELL, 0, idx=idx))
    ASSERT(all(shape(p) == [this%dim, this%cell_dim]))

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

    ASSERT(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))

    IGNORE(idx_dir)

    len = this%x(idx(1)+1) - this%x(idx(1))
  end function

  function grid1D_get_surf(this, idx, idx_dir) result(surf)
    !! get single face's surface
    class(grid1D), intent(in) :: this
    integer,       intent(in) :: idx(:)
      !! face's indices. size: (idx_dim)
    integer,       intent(in) :: idx_dir
      !! face's direction
    real                      :: surf
      !! output: size of face

    ASSERT(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))

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

    ASSERT(this%idx_allowed(IDX_CELL, 0, idx=idx))

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

    ASSERT(this%idx_allowed(idx1_type, idx1_dir))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))

    IGNORE(this)
    IGNORE(idx1_dir)
    IGNORE(idx2_dir)

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

    integer :: idx_bnd(1), shift

    ASSERT(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    ASSERT(this%idx_allowed(idx2_type, idx2_dir))
    ASSERT(size(idx2) == this%idx_dim)

    IGNORE(idx1_dir)
    IGNORE(idx2_dir)

    idx2  = 0
    shift = 0

    if (idx1_type == idx2_type) then
      ! idx1 is not a neighbour of itself => either subtract or add 1
      if (idx1(1) == 1) shift = 1
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
      if (idx1(1) == 1) shift = 1
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
    status = ((idx2(1) >= 1) .and. (idx2(1) <= idx_bnd(1)))
  end subroutine

end module
