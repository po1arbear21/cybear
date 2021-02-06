#include "../util/macro.f90.inc"

module tensor_grid_m
  use error_m
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
    procedure :: get_surf    => tensor_grid_get_surf
    procedure :: get_vol     => tensor_grid_get_vol
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

  subroutine tensor_grid_get_idx_bnd(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(out) :: idx_bnd(:)
      !! output upper bound for each index (idx_dim)

    integer :: i, j0, j1, rdir

    ASSERT((((idx_type == IDX_VERTEX) .or. (idx_type == IDX_CELL)) .and. (idx_dir == 0)) \
      .or. (((idx_type == IDX_EDGE) .or. (idx_type == IDX_FACE)) .and. ((idx_dir >= 1) .and. (idx_dir <= this%idx_dim))))
    ASSERT(size(idx_bnd) == this%idx_dim)

    j1 = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      if ((idx_type == IDX_EDGE) .or. (idx_type == IDX_FACE)) then
        rdir = idx_dir - j0 + 1 ! relative direction
        if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
          if (idx_type == IDX_EDGE) then
            call this%g(i)%p%get_idx_bnd(IDX_VERTEX, 0, idx_bnd(j0:j1))
          else ! idx_type == IDX_FACE
            call this%g(i)%p%get_idx_bnd(IDX_CELL, 0, idx_bnd(j0:j1))
          end if
        else
          call this%g(i)%p%get_idx_bnd(idx_type, rdir, idx_bnd(j0:j1))
        end if
      else
        call this%g(i)%p%get_idx_bnd(idx_type, 0, idx_bnd(j0:j1))
      end if
    end do
  end subroutine

  subroutine tensor_grid_get_vertex(this, idx, p)
    !! get vertex coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! vertex indices
    real,               intent(out) :: p(:)
      !! output vertex coordinates (dim)

    integer :: i, i0, i1, j0, j1

    ASSERT(size(idx) == this%idx_dim)
    ASSERT(size(p  ) == this%dim    )

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

  subroutine tensor_grid_get_edge(this, idx, idx_dir, p)
    !! get edge coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! edge indices
    integer,            intent(in)  :: idx_dir
      !! edge direction
    real,               intent(out) :: p(:,:)
      !! output edge coordinates (dim x 2)

    integer :: i, i0, i1, j0, j1, rdir

    ASSERT(size(idx) == this%idx_dim)
    ASSERT((idx_dir >= 1) .and. (idx_dir <= this%idx_dim))
    ASSERT(size(p  ) == this%dim    )

    i1 = 0
    j1 = 0
    do i = 1, size(this%g)
      i0 = i1 + 1
      i1 = i1 + this%g(i)%p%dim
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! relative direction for i-th grid
      rdir = idx_dir - j0 + 1

      ! get coordinates
      if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
        call this%g(i)%p%get_vertex(idx(j0:j1), p(i0:i1,1))
        p(i0:i1,2) = p(i0:i1,1)
      else
        call this%g(i)%p%get_edge(idx(j0:j1), rdir, p(i0:i1,1:2))
      end if
    end do
  end subroutine

  subroutine tensor_grid_get_face(this, idx, idx_dir, p)
    !! get face coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! face indices
    integer,            intent(in)  :: idx_dir
      !! face direction
    real,               intent(out) :: p(:,:)
      !! output face coordinates (dim x face_dim(idx_dir))

    integer :: i, j, k, l, m, n, c, i0, i1, j0, j1, rdir
    real    :: tmp(this%dim,this%face_dim(idx_dir))

    ASSERT(size(idx) == this%idx_dim)
    ASSERT((idx_dir >= 1) .and. (idx_dir <= this%idx_dim))
    ASSERT(size(p,1) == this%dim)
    ASSERT(size(p,2) == this%face_dim(idx_dir))

    ! tensor product point counter
    c = 1

    i1 = 0
    j1 = 0
    do i = 1, size(this%g)
      i0 = i1 + 1
      i1 = i1 + this%g(i)%p%dim
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! relative direction for i-th grid
      rdir = idx_dir - j0 + 1

      ! get face or cell points for i-th grid
      if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
        n = this%g(i)%p%cell_dim
        call this%g(i)%p%get_cell(idx(j0:j1), tmp(i0:i1,1:n))
      else
        n = this%g(i)%p%face_dim(rdir)
        call this%g(i)%p%get_face(idx(j0:j1), rdir, tmp(i0:i1,1:n))
      end if

      ! tensor product of points (combine points from this grid with all points from other grids)
      m = 0
      do j = 1, this%face_dim(idx_dir)/(n*c) ! repeat until result points are filled
        do k = 1, n                          ! loop over all n points from this grid
          do l = 1, c                        ! replicate point c times
            m = m + 1
            p(i0:i1,m) = tmp(i0:i1,k)
          end do
        end do
      end do

      ! update counter
      c = c * n
    end do
  end subroutine

  subroutine tensor_grid_get_cell(this, idx, p)
    !! get cell coordinates from grid indices
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell indices
    real,               intent(out) :: p(:,:)
      !! output cell coordinates (dim x cell_dim)

    integer :: i, j, k, l, m, n, c, i0, i1, j0, j1
    real    :: tmp(this%dim,this%cell_dim)

    ASSERT(size(idx) == this%idx_dim)
    ASSERT(size(p,1) == this%dim)
    ASSERT(size(p,2) == this%cell_dim)

    ! tensor product point counter
    c = 1

    i1 = 0
    j1 = 0
    do i = 1, size(this%g)
      i0 = i1 + 1
      i1 = i1 + this%g(i)%p%dim
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! get cell for i-th grid
      n = this%g(i)%p%cell_dim
      call this%g(i)%p%get_cell(idx(j0:j1), tmp(i0:i1,1:n))

      ! tensor product of points (combine points from this grid with all points from other grids)
      m = 0
      do j = 1, this%cell_dim/(n*c) ! repeat until result points are filled
        do k = 1, n                 ! loop over all n points from this grid
          do l = 1, c               ! replicate point c times
            m = m + 1
            p(i0:i1,m) = tmp(i0:i1,k)
          end do
        end do
      end do

      ! update counter
      c = c * n
    end do
  end subroutine

  function tensor_grid_get_surf(this, idx, idx_dir) result(surf)
    !! get size of face
    class(tensor_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! face indices
    integer,            intent(in) :: idx_dir
      !! face direction
    real                           :: surf
      !! return size of face

    integer :: i, j0, j1, rdir

    ASSERT(size(idx) == this%idx_dim)
    ASSERT((idx_dir >= 1) .and. (idx_dir <= this%idx_dim))

    surf = 1.0
    j1   = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! relative direction for i-th grid
      rdir = idx_dir - j0 + 1

      ! get surf or vol for i-th grid
      if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
        surf = surf * this%g(i)%p%get_vol(idx(j0:j1))
      else
        surf = surf * this%g(i)%p%get_surf(idx(j0:j1), rdir)
      end if
    end do
  end function

  function tensor_grid_get_vol(this, idx) result(vol)
    !! get cell volume
    class(tensor_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! cell indices
    real                           :: vol
      !! return cell volume

    integer :: i, j0, j1

    ASSERT(size(idx) == this%idx_dim)

    vol = 1.0
    j1  = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim
      vol = vol * this%g(i)%p%get_vol(idx(j0:j1))
    end do
  end function

end module