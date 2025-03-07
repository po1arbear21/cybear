m4_include(../util/macro.f90.inc)

module tensor_grid_m

  use array_m,       only: array2_int
  use error_m,       only: assert_failed
  use grid_m,        only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, grid, grid_ptr
  use grid_data_m,   only: allocate_grid_data, grid_data_int
  use hashmap_m,     only: hashmap_int
  use json_m,        only: json_array, json_object
  use output_file_m, only: output_file
  use vector_m,      only: vector_int

  implicit none

  private
  public tensor_grid

  type, extends(grid) :: tensor_grid
    !! Combination of multiple sub-grids by a tensor product
    type(grid_ptr), allocatable :: g(:)
      !! sub-grids

    integer,              allocatable :: max_neighb(:,:,:,:)
      !! maximal number of neighbours (idx1_type, idx1_dir, idx2_type, idx2_dir)
    class(grid_data_int), allocatable :: neighb_i0(:,:,:,:)
      !! neighbour start indices (idx1_type, idx1_dir, idx2_type, idx2_dir)
    class(grid_data_int), allocatable :: neighb_i1(:,:,:,:)
      !! neighbour end indices (idx1_type, idx1_dir, idx2_type, idx2_dir)
    type(array2_int),     allocatable :: neighb(:,:,:,:)
      !! neighbour data (delta_idx, i) x (idx1_type, idx1_dir, idx2_type, idx2_dir)
  contains
    procedure :: init           => tensor_grid_init
    procedure :: get_idx_bnd_n  => tensor_grid_get_idx_bnd_n
    procedure :: get_vertex     => tensor_grid_get_vertex
    procedure :: get_edge       => tensor_grid_get_edge
    procedure :: get_face       => tensor_grid_get_face
    procedure :: get_cell       => tensor_grid_get_cell
    procedure :: get_len        => tensor_grid_get_len
    procedure :: get_surf       => tensor_grid_get_surf
    procedure :: get_vol        => tensor_grid_get_vol
    procedure :: get_max_neighb => tensor_grid_get_max_neighb
    procedure :: get_neighb     => tensor_grid_get_neighb
    procedure :: get_adjoint    => tensor_grid_get_adjoint
    procedure :: output         => tensor_grid_output

    procedure, private :: init_neighb      => tensor_grid_init_neighb
    procedure, private :: init_neighb_loop => tensor_grid_init_neighb_loop
  end type

contains

  subroutine tensor_grid_init(this, name, g, no_neighbours)
    !! initialize tensor grid
    class(tensor_grid), intent(out) :: this
    character(*),       intent(in)  :: name
      !! grid name
    type(grid_ptr),     intent(in)  :: g(:)
      !! sub-grids
    logical, optional,  intent(in)  :: no_neighbours
      !! do not compute neighbour table (grid can not be used in conjunction with jacobians)

    integer              :: i, j, j0, j1, k, dim, idx_dim, cell_nvert
    integer, allocatable :: face_nvert(:), cell_nedge(:)
    logical              :: no_neighbours_

    ! get dimension and index dimension
    dim = 0
    idx_dim = 0
    do i = 1, size(g)
      dim     = dim     + g(i)%p%dim
      idx_dim = idx_dim + g(i)%p%idx_dim
    end do

    ! get number of vertices per face
    allocate (face_nvert(idx_dim))
    j1 = 0
    do i = 1, size(g)
      j0 = j1 + 1
      j1 = j1 + g(i)%p%idx_dim
      do j = j0, j1
        face_nvert(j) = g(i)%p%face_nvert(j-j0+1)
        do k = 1, size(g)
          if (k == i) cycle
          face_nvert(j) = face_nvert(j) * g(k)%p%cell_nvert
        end do
      end do
    end do

    ! get number of vertices per cell
    cell_nvert = 1
    do i = 1, size(g)
      cell_nvert = cell_nvert * g(i)%p%cell_nvert
    end do

    ! get number of edges per cell
    allocate (cell_nedge(idx_dim))
    j1 = 0
    do i = 1, size(g)
      j0 = j1 + 1
      j1 = j1 + g(i)%p%idx_dim
      do j = j0, j1
        cell_nedge(j) = g(i)%p%cell_nedge(j-j0+1)
        do k = 1, size(g)
          if (k == i) cycle
          cell_nedge(j) = cell_nedge(j) * g(k)%p%cell_nvert
        end do
      end do
    end do

    ! init base
    call this%grid_init(name, dim, idx_dim, face_nvert, cell_nvert, cell_nedge)

    ! unit
    allocate (this%unit(dim))
    this%unit = [(g(i)%p%unit, i = 1, size(g))]

    ! save sub-grid pointers
    this%g = g

    ! init neighbour data
    no_neighbours_ = .false.
    if (present(no_neighbours)) no_neighbours_ = no_neighbours
    if (.not. no_neighbours_) call this%init_neighb()
  end subroutine

  subroutine tensor_grid_get_idx_bnd_n(this, idx_type, idx_dir, idx_bnd)
    !! get grid index bounds
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(out) :: idx_bnd(:,:)
      !! output lower/upper bound for each index (2, idx_dim)

    integer :: i, j0, j1, rdir

    m4_assert(this%idx_allowed(idx_type, idx_dir))
    m4_assert(size(idx_bnd,1) == 2)
    m4_assert(size(idx_bnd,2) == this%idx_dim)

    j1 = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      if ((idx_type == IDX_EDGE) .or. (idx_type == IDX_FACE)) then
        rdir = idx_dir - j0 + 1 ! relative direction
        if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
          if (idx_type == IDX_EDGE) then
            call this%g(i)%p%get_idx_bnd(IDX_VERTEX, 0, idx_bnd(:,j0:j1))
          else ! idx_type == IDX_FACE
            call this%g(i)%p%get_idx_bnd(IDX_CELL, 0, idx_bnd(:,j0:j1))
          end if
        else
          call this%g(i)%p%get_idx_bnd(idx_type, rdir, idx_bnd(:,j0:j1))
        end if
      else
        call this%g(i)%p%get_idx_bnd(idx_type, 0, idx_bnd(:,j0:j1))
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

    m4_assert(this%idx_allowed(IDX_VERTEX, 0, idx=idx))
    m4_assert(size(p) == this%dim)

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

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, 2]))

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
      !! output face coordinates (dim x face_nvert(idx_dir))

    integer :: i, j, k, l, m, n, c, i0, i1, j0, j1, rdir
    real    :: tmp(this%dim,this%face_nvert(idx_dir))

    m4_assert(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%face_nvert(idx_dir)]))

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
        n = this%g(i)%p%cell_nvert
        call this%g(i)%p%get_cell(idx(j0:j1), tmp(i0:i1,1:n))
      else
        n = this%g(i)%p%face_nvert(rdir)
        call this%g(i)%p%get_face(idx(j0:j1), rdir, tmp(i0:i1,1:n))
      end if

      ! tensor product of points (combine points from this grid with all points from other grids)
      m = 0
      do j = 1, this%face_nvert(idx_dir)/(n*c) ! repeat until result points are filled
        do k = 1, n                            ! loop over all n points from this grid
          do l = 1, c                          ! replicate point c times
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
      !! output cell coordinates (dim x cell_nvert)

    integer :: i, j, k, l, m, n, c, i0, i1, j0, j1
    real    :: tmp(this%dim,this%cell_nvert)

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))
    m4_assert(all(shape(p) == [this%dim, this%cell_nvert]))

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
      n = this%g(i)%p%cell_nvert
      call this%g(i)%p%get_cell(idx(j0:j1), tmp(i0:i1,1:n))

      ! tensor product of points (combine points from this grid with all points from other grids)
      m = 0
      do j = 1, this%cell_nvert/(n*c) ! repeat until result points are filled
        do k = 1, n                   ! loop over all n points from this grid
          do l = 1, c                 ! replicate point c times
            m = m + 1
            p(i0:i1,m) = tmp(i0:i1,k)
          end do
        end do
      end do

      ! update counter
      c = c * n
    end do
  end subroutine

  function tensor_grid_get_len(this, idx, idx_dir) result(len)
    !! get edge length
    class(tensor_grid), intent(in) :: this
    integer,            intent(in) :: idx(:)
      !! edge indices (idx_dim)
    integer,            intent(in) :: idx_dir
      !! edge direction
    real                           :: len
      !! return edge length

    integer :: i, j0, j1, rdir

    m4_assert(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))

    j1 = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! relative direction for i-th grid
      rdir = idx_dir - j0 + 1

      ! get edge length
      if ((rdir >= 1) .and. (rdir <= this%g(i)%p%idx_dim)) then
        len = this%g(i)%p%get_len(idx(j0:j1), rdir)
        return
      end if
    end do
  end function

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

    m4_assert(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))

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

    m4_assert(this%idx_allowed(IDX_CELL, 0, idx=idx))

    vol = 1.0
    j1  = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim
      vol = vol * this%g(i)%p%get_vol(idx(j0:j1))
    end do
  end function

  function tensor_grid_get_max_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir) result(max_neighb)
    !! get maximal number of nearest neighbours
    class(tensor_grid), intent(in) :: this
    integer,            intent(in) :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in) :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in) :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in) :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer                        :: max_neighb
      !! return maximal number of nearest neighbours

    m4_assert(this%idx_allowed(idx1_type, idx1_dir))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))

    max_neighb = this%max_neighb(idx1_type,idx1_dir,idx2_type,idx2_dir)
  end function

  subroutine tensor_grid_get_neighb(this, idx1_type, idx1_dir, idx2_type, idx2_dir, idx1, j, idx2, status)
    !! get j-th neighbour.
    !!
    !! j: we count neighbors from 1,2,...,N.
    !! N: depends on idx1 (e.g. boundary nodes might have fewer neighbors).
    !! status: indicates if j-th neighbor exists
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx1_type
      !! first index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx1_dir
      !! first index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in)  :: idx2_type
      !! neighbour index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,            intent(in)  :: idx2_dir
      !! neighbour index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,            intent(in)  :: idx1(:)
      !! first indices. size: (idx_dim)
    integer,            intent(in)  :: j
      !! j-th neighbor
    integer,            intent(out) :: idx2(:)
      !! output neighbour indices. size: (idx_dim)
    logical,            intent(out) :: status
      !! does j-th neighbor exist?

    integer :: i0, i1, i

    m4_assert(this%idx_allowed(idx1_type, idx1_dir, idx=idx1))
    m4_assert(this%idx_allowed(idx2_type, idx2_dir))
    m4_assert(size(idx2) == this%idx_dim)

    ! lookup neighbour index
    i0 = this%neighb_i0(idx1_type,idx1_dir,idx2_type,idx2_dir)%get(idx1)
    i1 = this%neighb_i1(idx1_type,idx1_dir,idx2_type,idx2_dir)%get(idx1)
    i  = i0 + j - 1

    ! check if valid neighbour
    if ((i < 1) .or. (i > i1)) then
      idx2 = 0
      status = .false.
      return
    end if

    ! extract j-th neighbour indices from precomputed table
    idx2   = idx1 + this%neighb(idx1_type,idx1_dir,idx2_type,idx2_dir)%d(:,i)
    status = .true.
  end subroutine

  subroutine tensor_grid_get_adjoint(this, idx, len, surf, vol)
    !! get adjoint grid information per cell
    class(tensor_grid), intent(in)  :: this
    integer,            intent(in)  :: idx(:)
      !! cell indices. size: (idx_dim)
    real, optional,     intent(out) :: len(:,:)
      !! edge lengths (max_cell_nedge, idx_dim)
    real, optional,     intent(out) :: surf(:,:)
      !! adjoint surface parts per edge (max_cell_nedge, idx_dim)
    real, optional,     intent(out) :: vol(:)
      !! adjoint volume parts per vertex (cell_nvert)

    integer :: i, i0, i1, isurf, ivol, j
    integer :: idx_dir, idx0(this%idx_dim), idx2(this%idx_dim)
    logical :: status
    real    :: gi_surf(this%max_cell_nedge,this%idx_dim), gi_vol(this%cell_nvert)

    if (present(len)) then
      m4_assert(size(len,1) == this%max_cell_nedge)
      m4_assert(size(len,2) == this%idx_dim)

      do idx_dir = 1, this%idx_dim
        do j = 1, this%cell_nedge(idx_dir)
          call this%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          len(j,idx_dir)  = this%get_len(idx2, idx_dir)
        end do
      end do
    end if

    if (present(surf)) then
      m4_assert(size(surf,1) == this%max_cell_nedge)
      m4_assert(size(surf,2) == this%idx_dim)

      surf = 1.0

      ! loop over subgrids, idx_dir range for i-th grid given by i0:i1
      i1 = 0
      do i = 1, size(this%g)
        associate (gi => this%g(i)%p)
          i0 = i1 + 1
          i1 = i1 + gi%idx_dim

          ! get adjoint cell information for i-th subgrid
          call gi%get_adjoint(idx(i0:i1), surf = gi_surf(1:gi%max_cell_nedge,1:gi%idx_dim), vol = gi_vol(1:gi%cell_nvert))

          ! loop over all edge directions
          do idx_dir = 1, this%idx_dim
            isurf = 0
            ivol  = 0
            idx2  = huge(idx2)

            ! loop over all edges in this direction
            do j = 1, this%cell_nedge(idx_dir)
              ! get edge
              idx0 = idx2
              call this%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)

              if ((idx_dir >= i0) .and. (idx_dir <= i1)) then
                ! update surf index if necessary
                if (any(idx2(i0:i1) /= idx0(i0:i1))) isurf = mod(isurf, gi%cell_nedge(idx_dir - i0 + 1)) + 1

                ! update surf
                surf(j,idx_dir) = surf(j,idx_dir) * gi_surf(isurf,idx_dir - i0 + 1)
              else
                ! update vol index if necessary
                if (any(idx2(i0:i1) /= idx0(i0:i1))) ivol  = mod(ivol, gi%cell_nvert) + 1

                ! update surf
                surf(j,idx_dir) = surf(j,idx_dir) * gi_vol(ivol)
              end if
            end do
          end do
        end associate
      end do
    end if

    if (present(vol)) then
      m4_assert(size(idx) == this%idx_dim)
      m4_assert(size(vol) == this%cell_nvert)

      vol = 1.0

      i1 = 0
      do i = 1, size(this%g)
        associate (gi => this%g(i)%p)
          i0 = i1 + 1
          i1 = i1 + gi%idx_dim

          ! get adjoint cell information for i-th subgrid
          call gi%get_adjoint(idx(i0:i1), vol = gi_vol(1:gi%cell_nvert))

          ivol = 0
          idx2 = huge(idx2)

          ! loop over all vertices
          do j = 1, this%cell_nvert
            idx0 = idx2
            call this%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx, j, idx2, status)

            if (any(idx2(i0:i1) /= idx0(i0:i1))) ivol = mod(ivol, gi%cell_nvert) + 1
            vol(j) = vol(j) * gi_vol(ivol)
          end do
        end associate
      end do
    end if
  end subroutine

  subroutine tensor_grid_output(this, of, unit)
    !! output tensor grid
    class(tensor_grid),     intent(in)    :: this
    type(output_file),      intent(inout) :: of
      !! output file handle
    character(*), optional, intent(in)    :: unit
      !! physical unit of coordinates (ignored)

    integer                    :: i
    type(json_array),  pointer :: subgrids
    type(json_object), pointer :: obj

    m4_ignore(unit)

    obj => of%new_object("Grids")
    call obj%add_string("Name", this%name)
    call obj%add_string("Type", "Tensor")
    call obj%add_array("Subgrids", subgrids)
    do i = 1, size(this%g)
      call subgrids%add_string(this%g(i)%p%name)
    end do
  end subroutine

  subroutine tensor_grid_init_neighb(this)
    !! init neighbour data
    class(tensor_grid), intent(inout) :: this

    integer          :: i, j, idx1_type, idx2_type, idx1_dir, idx2_dir, idir0(4), idir1(4), todo(4,4*4*(this%idx_dim+1)*(this%idx_dim+1)), n

    ! allocate memory
    allocate (this%max_neighb(4,0:this%idx_dim,4,0:this%idx_dim), source = 0)
    call allocate_grid_data(this%neighb_i0, this%idx_dim, [1,0,1,0], [4,this%idx_dim,4,this%idx_dim])
    call allocate_grid_data(this%neighb_i1, this%idx_dim, [1,0,1,0], [4,this%idx_dim,4,this%idx_dim])
    allocate (this%neighb(4,0:this%idx_dim,4,0:this%idx_dim))

    ! direction bounds
    idir0 = [0,            1,            1, 0]
    idir1 = [0, this%idx_dim, this%idx_dim, 0]

    ! use todo list (better granularity for openmp)
    n = 0
    do idx2_type = 1, 4
      do idx1_type = 1, 4
        do idx2_dir = idir0(idx2_type), idir1(idx2_type)
          do idx1_dir = idir0(idx1_type), idir1(idx1_type)
            n = n + 1
            todo(1,n) = idx1_type
            todo(2,n) = idx1_dir
            todo(3,n) = idx2_type
            todo(4,n) = idx2_dir
          end do
        end do
      end do
    end do

    ! init neighbours
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i, j, idx1_type, idx2_type, idx1_dir, idx2_dir) &
    !$omp shared(this, idir0, idir1, todo, n)
    do i = 1, n
      idx1_type = todo(1,i)
      idx1_dir  = todo(2,i)
      idx2_type = todo(3,i)
      idx2_dir  = todo(4,i)
      call this%init_neighb_loop(idx1_type, idx1_dir, idx2_type, idx2_dir)
    end do
    !$omp end parallel do
  end subroutine

  subroutine tensor_grid_init_neighb_loop(this, idx1_type, idx1_dir, idx2_type, idx2_dir)
    !! init neighbour data for one combination of idx1_type, idx1_dir, idx2_type idx2_dir
    class(tensor_grid), intent(inout) :: this
    integer,            intent(in)    :: idx1_type
    integer,            intent(in)    :: idx2_type
    integer,            intent(in)    :: idx1_dir
    integer,            intent(in)    :: idx2_dir

    integer, parameter :: ADDOP = 1, MULOP = 2, SELOP = 3, MULSELOP = 4
    integer, parameter :: OP_TABLE(4,4) = reshape([ADDOP, SELOP, MULOP, MULOP, &
      &                                            SELOP, ADDOP, MULOP, MULOP, &
      &                                            MULOP, MULOP, ADDOP, SELOP, &
      &                                            MULOP, MULOP, SELOP, ADDOP], [4, 4])

    integer           :: i, i0, i1, j, j0, j1, k, op, bnd(2,this%idx_dim), max_neighb, isearch
    integer           :: idx1(this%idx_dim), idx2(this%idx_dim), rtype1, rtype2, rdir1, rdir2
    integer           :: imul(size(this%g)), nmul(size(this%g))
    logical           :: select1, select2, status
    type(vector_int)  :: vmul(size(this%g),this%idx_dim), key, vneighb
    type(hashmap_int) :: hmap

    ! allocate memory
    do i = 1, this%idx_dim
      do j = 1, size(this%g)
        call vmul(j,i)%init(0, c = 16)
      end do
    end do
    call hmap%init()
    call key%init(0, c = 16)
    call vneighb%init(0, c = 16)

    ! allocate memory
    call this%neighb_i0(idx1_type,idx1_dir,idx2_type,idx2_dir)%init(this, idx1_type, idx1_dir)
    call this%neighb_i1(idx1_type,idx1_dir,idx2_type,idx2_dir)%init(this, idx1_type, idx1_dir)

    ! get operation (how to combine neighbours from sub-grids)
    op = OP_TABLE(idx1_type, idx2_type)
    if ((idx1_type == idx2_type) .and. (idx1_dir /= idx2_dir)) op = MULSELOP

    ! loop over all grid indices
    call this%get_idx_bnd(idx1_type, idx1_dir, bnd)
    idx1 = bnd(1,:)
    do while (idx1(this%idx_dim) <= bnd(2,this%idx_dim))
      ! clear temporary neighbour table (used as key in hmap)
      call key%reset()

      ! loop over sub-grids
      j1 = 0
      do i = 1, size(this%g)
        j0 = j1 + 1
        j1 = j1 + this%g(i)%p%idx_dim

        ! set relative types, directions and selection flags for result and dependency (uses i and j0 implicitly)
        call set_relative(idx1_type, idx1_dir, rtype1, rdir1, select1)
        call set_relative(idx2_type, idx2_dir, rtype2, rdir2, select2)

        ! selection operation: select neighbours only from valid sub-grid
        if ((op == SELOP) .and. .not. (select1 .and. select2)) cycle

        ! maximum number of neighbours for this sub-grid
        max_neighb = this%g(i)%p%get_max_neighb(rtype1, rdir1, rtype2, rdir2)

        ! perform operation
        if ((op == ADDOP) .or. (op == SELOP)) then ! add neighbours from sub-grids
          ! loop over neighbours for this sub-grid
          idx2 = idx1
          do j = 1, max_neighb
            call this%g(i)%p%get_neighb(rtype1, rdir1, rtype2, rdir2, idx1(j0:j1), j, idx2(j0:j1), status)
            if (.not. status) exit

            ! add neighbour to temporary table
            do k = 1, this%idx_dim
              call key%push(idx2(k) - idx1(k))
            end do
          end do

          ! selection operation: done, remaining sub-grids will not be selected and can be skipped
          if (op == SELOP) exit
        else ! MULOP, MULSELOP: tensor product of neighbours from all sub-grids
          do k = j0, j1
            call vmul(i,k)%reset()
          end do
          if ((op == MULSELOP) .and. (.not. select1 .and. .not. select2)) then
            do k = j0, j1
              call vmul(i,k)%push(idx1(k))
            end do
          else
            do j = 1, max_neighb
              call this%g(i)%p%get_neighb(rtype1, rdir1, rtype2, rdir2, idx1(j0:j1), j, idx2(j0:j1), status)
              if (.not. status) exit

              ! save neighbours temporarily
              do k = j0, j1
                call vmul(i,k)%push(idx2(k))
              end do
            end do
          end if

          ! number of neighbours for i-th sub-grid
          nmul(i) = vmul(i,j0)%n
        end if
      end do

      if ((op == MULOP) .or. (op == MULSELOP)) then
        ! perform tensor product
        imul = 1 ! select first combination
        do while (imul(size(this%g)) <= nmul(size(this%g)))
          ! set idx2 by combining indices from sub-grids
          j1 = 0
          do i = 1, size(this%g)
            j0 = j1 + 1
            j1 = j1 + this%g(i)%p%idx_dim
            do k = j0, j1
              ! select indices of imul(i)-th neighbour from i-th grid
              idx2(k) = vmul(i,k)%d(imul(i))
            end do
          end do

          ! add neighbour to temporary table
          do k = 1, this%idx_dim
            call key%push(idx2(k)-idx1(k))
          end do

          ! select next combination
          imul(1) = imul(1) + 1
          do i = 1, size(this%g) - 1
            if (imul(i) <= nmul(i)) exit
            imul(i  ) = 1
            imul(i+1) = imul(i+1) + 1
          end do
        end do
      end if

      if (key%n > 0) then
        ! check if same data already exists (use hashmap for constant time search)
        call hmap%get(key%d(1:key%n), isearch, status = status)

        ! save new neighbour data
        if (.not. status) then
          isearch = vneighb%n / this%idx_dim
          call vneighb%push(key%d(1:key%n))
          call hmap%set(key%d(1:key%n), isearch)
        end if
      else
        isearch = vneighb%n / this%idx_dim
      end if

      i0 = isearch + 1
      i1 = isearch + key%n / this%idx_dim

      ! save start/end indices of chunk
      call this%neighb_i0(idx1_type,idx1_dir,idx2_type,idx2_dir)%set(idx1, i0)
      call this%neighb_i1(idx1_type,idx1_dir,idx2_type,idx2_dir)%set(idx1, i1)

      ! update this%max_neighb
      associate (n => this%max_neighb(idx1_type,idx1_dir,idx2_type,idx2_dir))
        if (i1 - i0 + 1 > n) n = i1 - i0 + 1
      end associate

      ! next idx1
      idx1(1) = idx1(1) + 1
      do i = 1, this%idx_dim - 1
        if (idx1(i) <= bnd(2,i)) exit
        idx1(i  ) = bnd(1,i)
        idx1(i+1) = idx1(i+1) + 1
      end do
    end do

    ! save neighb
    this%neighb(idx1_type,idx1_dir,idx2_type,idx2_dir)%d = reshape(vneighb%d(1:vneighb%n), [this%idx_dim, vneighb%n/this%idx_dim])

  contains

    subroutine set_relative(idx_type, idx_dir, rtype, rdir, select)
      !! set relative type, direction and selection flag for i-th grid
      integer, intent(in)  :: idx_type
      integer, intent(in)  :: idx_dir
      integer, intent(out) :: rtype
      integer, intent(out) :: rdir
      logical, intent(out) :: select

      if ((idx_type == IDX_EDGE) .or. (idx_type == IDX_FACE)) then
        rdir = idx_dir - j0 + 1
        if ((rdir < 1) .or. (rdir > this%g(i)%p%idx_dim)) then
          if (idx_type == IDX_EDGE) then
            rtype = IDX_VERTEX
          else
            rtype = IDX_CELL
          end if
          rdir   = 0
          select = .false.
        else
          rtype  = idx_type
          select = .true.
        end if
      else
        rtype  = idx_type
        rdir   = 0
        select = .true.
      end if
    end subroutine

  end subroutine

end module
