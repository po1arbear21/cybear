module grid1D_m
  use grid_m
  implicit none

  type, extends(grid) :: grid1D
    !! 1D grid
    real, allocatable :: x(:)
      !! grid points
  contains
    procedure :: init        => grid1D_init
    procedure :: get_idx_bnd => grid1D_get_idx_bnd
    procedure :: get_vertex  => grid1D_get_vertex
    procedure :: get_edge    => grid1D_get_edge
    procedure :: get_face    => grid1D_get_face
    procedure :: get_cell    => grid1D_get_cell
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

  subroutine grid1D_get_idx_bnd(this, idx_type, dir, idx_bnd)
    !! get grid index bounds
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx_type
      !! grid index type (e.g. IDX_VERTEX)
    integer,       intent(in)  :: dir
      !! index of direction (only used for IDX_EDGE and IDX_FACE; range = 1:idx_dim)
    integer,       intent(out) :: idx_bnd(:)
      !! output upper bound for each index (1)

    select case (idx_type)
      case (IDX_VERTEX)
        idx_bnd(1) = size(this%x)
      case (IDX_EDGE)
        if (dir == 1) then
          idx_bnd(1) = size(this%x) - 1
        else
          idx_bnd(1) = size(this%x)
        end if
      case (IDX_FACE)
        if (dir == 1) then
          idx_bnd(1) = size(this%x)
        else
          idx_bnd(1) = size(this%x) - 1
        end if
      case (IDX_CELL)
        idx_bnd(1) = size(this%x) - 1
    end select
  end subroutine

  subroutine grid1D_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! vertex indices
    real,          intent(out) :: p(:)
      !! output vertex coordinates (1:dim)

    p(1) = this%x(idx(1))
  end subroutine

  subroutine grid1D_get_edge(this, idx, dir, p)
    !! get edge coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! edge indices
    integer,       intent(in)  :: dir
      !! edge direction
    real,          intent(out) :: p(:,:)
      !! output edge coordinates (1:dim x 1:2)

    if (dir == 1) then
      p(1,1:2) = this%x(idx(1):idx(1)+1)
    else
      p(1,1:2) = this%x(idx(1))
    end if
  end subroutine

  subroutine grid1D_get_face(this, idx, dir, p)
    !! get face coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! face indices
    integer,       intent(in)  :: dir
      !! face direction
    real,          intent(out) :: p(:,:)
      !! output face coordinates (1:dim x 1:face_dim(dir))

    if (dir == 1) then
      p(1,1) = this%x(idx(1))
    else
      p(1,1:2) = this%x(idx(1):idx(1)+1)
    end if
  end subroutine

  subroutine grid1D_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(grid1D), intent(in)  :: this
    integer,       intent(in)  :: idx(:)
      !! cell indices
    real,          intent(out) :: p(:,:)
      !! output cell coordinates (1:dim x 1:cell_dim)

    p(1,1:2) = this%x(idx(1):idx(1)+1)
  end subroutine

end module