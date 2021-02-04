module tensor_grid_m
  use grid_m
  implicit none

  type, extends(grid) :: tensor_grid
    !! Combination of multiple sub-grids by a tensor product
    type(grid_ptr), allocatable :: g(:)
      !! sub-grids
  contains
    procedure :: init        => tensor_grid_init
    procedure :: get_idx_bnd => tensor_grid_get_idx_bnd
    procedure :: get_vertex  => tensor_grid_get_vertex
    procedure :: get_edge    => tensor_grid_get_edge
    procedure :: get_face    => tensor_grid_get_face
    procedure :: get_cell    => tensor_grid_get_cell
  end type

contains

  subroutine tensor_grid_init(this, g)
    !! initialize tensor grid
    class(tensor_grid), intent(out) :: this
    type(grid_ptr),     intent(in)  :: g(:)
      !! sub-grids

    integer              :: i, j, j0, j1, k, dim, idx_dim, cell_dim
    integer, allocatable :: face_dim(:)

    ! get dimension and index dimension
    dim = 0
    idx_dim = 0
    do i = 1, size(g)
      dim     = dim     + g(i)%p%dim
      idx_dim = idx_dim + g(i)%p%idx_dim
    end do

    ! get number of points per face
    allocate (face_dim(idx_dim))
    j1 = 0
    do i = 1, size(g)
      j0 = j1 + 1
      j1 = j1 + g(i)%p%idx_dim
      do j = j0, j1
        face_dim(j) = g(i)%p%face_dim(j)
        do k = 1, size(g)
          if (k == i) cycle
          face_dim(j) = face_dim(j) * g(k)%p%cell_dim
        end do
      end do
    end do

    ! get number of points per cell
    cell_dim = 1
    do i = 1, size(g)
      cell_dim = cell_dim * g(i)%p%cell_dim
    end do

    ! init base
    call this%grid_init(dim, idx_dim, face_dim, cell_dim)

    ! save sub-grid pointers
    this%g = g
  end subroutine

  subroutine tensor_grid_get_idx_bnd(this, idx_type, dir, idx_bnd)
    !! get grid index bounds
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (e.g. IDX_VERTEX)
    integer,            intent(in)  :: dir
      !! index of direction (only used for IDX_EDGE and IDX_FACE; range = 1:idx_dim)
    integer,            intent(out) :: idx_bnd(:)
      !! output upper bound for each index (1:idx_dim)

    integer :: i, j0, j1

    j1 = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim
      call this%g(i)%p%get_idx_bnd(idx_type, dir-j0+1, idx_bnd(j0:j1))
    end do
  end subroutine

  subroutine tensor_grid_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! vertex indices
    real,               intent(out) :: p(:)
      !! output vertex coordinates (1:dim)

    integer :: i, i0, i1, j0, j1

    i1 = 0
    j1 = 0
    do i = 1, size(this%g)
      i0 = i1 + 1
      i1 = i1 + this%g(i)%p%dim
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim
      call this%g(i)%p%get_vertex(idx(j0:j1), p(i0:i1))
    end do
  end subroutine

  subroutine tensor_grid_get_edge(this, idx, dir, p)
    !! get edge coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! edge indices
    integer,            intent(in)  :: dir
      !! edge direction
    real,               intent(out) :: p(:,:)
      !! output edge coordinates (1:dim x 1:2)

    integer :: i, i0, i1, j0, j1

    i1 = 0
    j1 = 0
    do i = 1, size(this%g)
      i0 = i1 + 1
      i1 = i1 + this%g(i)%p%dim
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim
      call this%g(i)%p%get_edge(idx(j0:j1), dir-j0+1, p(i0:i1,1:2))
    end do
  end subroutine

  subroutine tensor_grid_get_face(this, idx, dir, p)
    !! get face coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! face indices
    integer,            intent(in)  :: dir
      !! face direction
    real,               intent(out) :: p(:,:)
      !! output face coordinates (1:dim x 1:face_dim(dir))

    ! FIXME
  end subroutine

  subroutine tensor_grid_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell indices
    real,               intent(out) :: p(:,:)
      !! output cell coordinates (1:dim x 1:cell_dim)

    ! FIXME
  end subroutine

end module