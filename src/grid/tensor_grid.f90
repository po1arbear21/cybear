#include "../util/macro.f90.inc"

module tensor_grid_m

  use error_m
  use grid_m,   only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, grid, grid_ptr, grid_data_int, &
    &                 grid_data1_int, grid_data2_int, grid_data3_int, grid_data4_int, &
    &                 grid_data5_int, grid_data6_int, grid_data7_int, grid_data8_int
  use vector_m, only: vector_int

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
    integer,              allocatable :: neighb(:,:)
      !! neighbour data (idx2, i)
  contains
    procedure :: init           => tensor_grid_init
    procedure :: get_idx_bnd    => tensor_grid_get_idx_bnd
    procedure :: get_vertex     => tensor_grid_get_vertex
    procedure :: get_edge       => tensor_grid_get_edge
    procedure :: get_face       => tensor_grid_get_face
    procedure :: get_cell       => tensor_grid_get_cell
    procedure :: get_len        => tensor_grid_get_len
    procedure :: get_surf       => tensor_grid_get_surf
    procedure :: get_vol        => tensor_grid_get_vol
    procedure :: get_max_neighb => tensor_grid_get_max_neighb
    procedure :: get_neighb     => tensor_grid_get_neighb

    procedure, private :: init_neighb => tensor_grid_init_neighb
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

    ! init neighbour data
    call this%init_neighb()

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

    ASSERT(this%idx_allowed(IDX_EDGE, idx_dir, idx=idx))

    j1 = 0
    do i = 1, size(this%g)
      j0 = j1 + 1
      j1 = j1 + this%g(i)%p%idx_dim

      ! relative direction for i-th grid
      rdir = idx_dir - j0 + 1

      ! get edge length
      if ((rdir >= 1) .and. (rdir <= this%g(i)%p%idx_dim)) then
        len = this%g(i)%p%get_len(idx, rdir)
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

    ASSERT(this%idx_allowed(IDX_FACE, idx_dir, idx=idx))

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

    ASSERT(this%idx_allowed(IDX_CELL, 0, idx=idx))

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
    idx2   = this%neighb(:,i)
    status = .true.
  end subroutine

  subroutine tensor_grid_init_neighb(this)
    !! init neighbour data
    class(tensor_grid), intent(inout) :: this

    integer, parameter :: ADDOP = 1, MULOP = 2, SELOP = 3
    integer, parameter :: OP_TABLE(4,4) = reshape([ADDOP, SELOP, MULOP, MULOP, &
      &                                            SELOP, ADDOP, MULOP, MULOP, &
      &                                            MULOP, MULOP, ADDOP, SELOP, &
      &                                            MULOP, MULOP, SELOP, ADDOP], [4, 4])

    integer          :: idx1_type, idx2_type, idx1_dir, idx2_dir, idir0(4), idir1(4), op, rtype1, rtype2, rdir1, rdir2
    integer          :: i, j, k, j0, j1, idx1(this%idx_dim), idx2(this%idx_dim), bnd(this%idx_dim), max_neighb, i0, i1
    integer          :: idx2_mul_n(size(this%g)), imul(size(this%g))
    logical          :: select1, select2, status
    type(vector_int) :: idx2_vec(this%idx_dim), idx2_mul(size(this%g),this%idx_dim)

    ! allocate memory
    allocate (this%max_neighb(4,0:this%idx_dim,4,0:this%idx_dim), source = 0)
    select case (this%idx_dim)
      case (1)
        allocate (grid_data1_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (2)
        allocate (grid_data2_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (3)
        allocate (grid_data3_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (4)
        allocate (grid_data4_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (5)
        allocate (grid_data5_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (6)
        allocate (grid_data6_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (7)
        allocate (grid_data7_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case (8)
        allocate (grid_data8_int :: this%neighb_i0(4,0:this%idx_dim,4,0:this%idx_dim), this%neighb_i1(4,0:this%idx_dim,4,0:this%idx_dim))
      case default
        call program_error("idx_dim must be in range 1:8")
    end select
    do i = 1, this%idx_dim
      call idx2_vec(i)%init(0, c = 1024)
      do j = 1, size(this%g)
        call idx2_mul(j,i)%init(0, c = 16)
      end do
    end do

    ! direction bounds
    idir0 = [0,            1,            1, 0]
    idir1 = [0, this%idx_dim, this%idx_dim, 0]

    ! init neighbours
    do idx2_type = 1, 4
      do idx2_dir = idir0(idx2_type), idir1(idx2_type)
        do idx1_type = 1, 4
          do idx1_dir = idir0(idx1_type), idir1(idx1_type)
            call this%neighb_i0(idx1_type,idx1_dir,idx2_type,idx2_dir)%init(this, idx1_type, idx1_dir)
            call this%neighb_i1(idx1_type,idx1_dir,idx2_type,idx2_dir)%init(this, idx1_type, idx1_dir)

            ! get operation (how to combine neighbours from sub-grids)
            op = OP_TABLE(idx1_type, idx2_type)

            ! loop over all grid indices
            call this%get_idx_bnd(idx1_type, idx1_dir, bnd)
            idx1 = 1
            do while (idx1(this%idx_dim) <= bnd(this%idx_dim))
              ! start idx2 table index
              i0 = idx2_vec(1)%n + 1

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

                    ! add neighbours to global vector
                    do k = 1, this%idx_dim
                      call idx2_vec(k)%push(idx2(k))
                    end do
                  end do

                  ! selection operation: done, remaining sub-grids will not be selected and can be skipped
                  if (op == SELOP) exit
                else ! MULOP: tensor product of neighbours from all sub-grids
                  do k = j0, j1
                    call idx2_mul(i,k)%reset()
                  end do
                  do j = 1, max_neighb
                    call this%g(i)%p%get_neighb(rtype1, rdir1, rtype2, rdir2, idx1(j0:j1), j, idx2(j0:j1), status)
                    if (.not. status) exit

                    ! save neighbours temporarily
                    do k = j0, j1
                      call idx2_mul(i,k)%push(idx2(k))
                    end do
                  end do

                  ! number of neighbours for i-th sub-grid
                  idx2_mul_n(i) = idx2_mul(i,j0)%n
                end if
              end do

              if (op == MULOP) then
                ! perform tensor product
                imul = 1 ! select first combination
                do while (imul(size(this%g)) <= idx2_mul_n(size(this%g)))
                  ! set idx2 by combining indices from sub-grids
                  j1 = 0
                  do i = 1, size(this%g)
                    j0 = j1 + 1
                    j1 = j1 + this%g(i)%p%idx_dim
                    do k = j0, j1
                      ! select indices of imul(i)-th neighbour from i-th grid
                      idx2(k) = idx2_mul(i,k)%d(imul(i))
                    end do
                  end do

                  ! save idx2
                  do k = 1, this%idx_dim
                    call idx2_vec(k)%push(idx2(k))
                  end do

                  ! select next combination
                  imul(1) = imul(1) + 1
                  do i = 1, size(this%g) - 1
                    if (imul(1) <= idx2_mul_n(i)) exit
                    imul(i  ) = 1
                    imul(i+1) = imul(i+1) + 1
                  end do
                end do
              end if

              ! end idx2 table index
              i1 = idx2_vec(1)%n

              ! save start + end idx2 table indices
              call this%neighb_i0(idx1_type,idx1_dir,idx2_type,idx2_dir)%set(idx1, i0)
              call this%neighb_i1(idx1_type,idx1_dir,idx2_type,idx2_dir)%set(idx1, i1)

              ! update this%max_neighb
              associate (n => this%max_neighb(idx1_type,idx1_dir,idx2_type,idx2_dir))
                if (i1 - i0 + 1 > n) n = i1 - i0 + 1
              end associate

              ! next idx1
              idx1(1) = idx1(1) + 1
              do i = 1, this%idx_dim - 1
                if (idx1(i) <= bnd(i)) exit
                idx1(i  ) = 1
                idx1(i+1) = idx1(i+1) + 1
              end do
            end do
          end do
        end do
      end do
    end do

    ! convert and save idx2_vec
    allocate (this%neighb(this%idx_dim,idx2_vec(1)%n))
    do i = 1, this%idx_dim
      this%neighb(i,1:idx2_vec(i)%n) = idx2_vec(i)%d(1:idx2_vec(i)%n)
    end do

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
