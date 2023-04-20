m4_include(../util/macro.f90.inc)

module stencil_m

  use error_m, only: assert_failed
  use grid_m,  only: grid
  use math_m,  only: eye_int

  implicit none

  private
  public stencil, stencil_ptr
  public sparse_stencil
  public dynamic_stencil
  public full_stencil
  public empty_stencil
  public dirichlet_stencil
  public near_neighb_stencil


  type, abstract :: stencil
    !! abstract stencil
  contains
    procedure :: get_ptr => stencil_get_ptr
  end type

  type stencil_ptr
    class(stencil), pointer :: p => null()
  end type

  type, abstract, extends(stencil) :: sparse_stencil
    !! abstract sparse stencil
    !! derive from it for specific equations (or use dirichlet/nearest neighb stencils).
    !! use when exact stencil is known at initialization time.
    integer :: nmax
      !! maximum number of dependency points
  contains
    procedure                               :: sparse_stencil_init
    procedure(sparse_stencil_get), deferred :: get
  end type

  abstract interface
    subroutine sparse_stencil_get(this, idx1, j, idx2, status)
      !! get j-th dependency wrt idx1
      !! note that j-th dependency might not be used which is why there is a status:
      !!    e.g. for boundary boxes actual number of dependencies might be smaller than standard number of dependenc == nmax
      import sparse_stencil
      class(sparse_stencil), intent(in)  :: this
      integer,               intent(in)  :: idx1(:)
        !! result indices. size: idx_dim1
      integer,               intent(in)  :: j
        !! j-th dependency
      integer,               intent(out) :: idx2(:)
        !! output j-th dependency indices. size: idx_dim2
      logical, optional,     intent(out) :: status
        !! is j-th dependency used?
    end subroutine
  end interface

  type, extends(stencil) :: dynamic_stencil
    !! placeholder to distinguish from sparse and full stencil.
    !! use when exact stencil is only known at evaluation time
    !! (in contrast to sparse/full stencil which is known at initialization time.)
  end type

  type, extends(stencil) :: full_stencil
    !! placeholder to distinguish from sparse and dynamic stencil
    !! use when stencil is given by the complete grid => dense jacobian
  end type

  type, extends(sparse_stencil) :: empty_stencil
    !! empty stencil with no neighbors.
    !! indicates no dependence on this variable
  contains
    procedure :: init => empty_stencil_init
    procedure :: get  => empty_stencil_get
  end type

  type, extends(sparse_stencil) :: dirichlet_stencil
    !! dirichlet stencil
    !!
    !! computation by: idx2 = perm(idx1) + off1:off2
    !!    idx1: result index
    !!    idx2: dependency index
    !!
    !! example
    !!    connect two grids in the following way
    !!      idx1 = (ix, iy,  iz     )
    !!      idx2 = (ix, 1:N, iz+2, 3)
    !!    you need following arguments for that
    !!      perm = (1, -1, 3, -1)
    !!      off1 = (0,  1, 2,  3)
    !!      off2 = (0,  N, 2,  3)

    integer, allocatable :: perm(:)
      !! permutation array.
      !! -1 for empty index element.
    integer, allocatable :: off1(:)
      !! start offset vector
    integer, allocatable :: off2(:)
      !! end offset vector
  contains
    procedure :: init => dirichlet_stencil_init
    procedure :: get  => dirichlet_stencil_get
  end type

  type, extends(sparse_stencil) :: near_neighb_stencil
    !! nearest neighbours stencil
    class(grid), pointer :: g => null()
      !! grid this stencil is defined on

    integer :: idx1_type
      !! result index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer :: idx1_dir
      !! result index direction for edges and faces
    integer :: idx2_type
      !! dependency index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer :: idx2_dir
      !! dependency index direction for edges and faces
  contains
    procedure :: init => near_neighb_stencil_init
    procedure :: get  => near_neighb_stencil_get
  end type

contains

  function stencil_get_ptr(this) result(ptr)
    !! returns pointer type to this stencil
    class(stencil), target, intent(in) :: this
    type(stencil_ptr)                  :: ptr

    ptr%p => this
  end function

  subroutine sparse_stencil_init(this, nmax)
    !! initialize sparse stencil
    class(sparse_stencil), intent(out) :: this
    integer,               intent(in)  :: nmax
      !! maximum number of dependency points

    this%nmax = nmax
  end subroutine

  subroutine empty_stencil_init(this)
    !! initialize empty stencil
    class(empty_stencil), intent(out) :: this

    call this%sparse_stencil_init(0)
  end subroutine

  subroutine empty_stencil_get(this, idx1, j, idx2, status)
    !! get j-th dependency wrt idx1
    class(empty_stencil), intent(in)  :: this
    integer,              intent(in)  :: idx1(:)
      !! result indices. size: idx_dim1
    integer,              intent(in)  :: j
      !! j-th dependency
    integer,              intent(out) :: idx2(:)
      !! output j-th dependency indices. size: idx_dim2
    logical, optional,    intent(out) :: status
      !! is j-th dependency used?

    m4_ignore(this)
    m4_ignore(idx1)
    m4_ignore(j)

    idx2(:) = 0
    if (present(status)) status = .false.
  end subroutine

  subroutine dirichlet_stencil_init(this, g1, g2, perm, off1, off2)
    !! initialize dirichlet stencil
    !!
    !! computation by: idx2 = perm(idx1) + off1:off2
    !!    idx1: result index
    !!    idx2: dependency index
    class(dirichlet_stencil),      intent(out) :: this
    class(grid),           target, intent(in)  :: g1
      !! result grid
    class(grid), optional, target, intent(in)  :: g2
      !! dependency grid (default: g1)
    integer,     optional,         intent(in)  :: perm(:)
      !! permutation array (default: identity).
      !! -1 for empty index element.
    integer,     optional,         intent(in)  :: off1(:)
      !! start offset vector (default: 0)
    integer,     optional,         intent(in)  :: off2(:)
      !! end offset vector (default: off1)

    class(grid), pointer :: g2_
    integer              :: i
    integer, allocatable :: off1_(:), off2_(:)

    ! optional g2
    g2_ => g1
    if (present(g2)) g2_ => g2

    ! optional off1
    allocate (off1_(g2_%idx_dim), source=0)
    if (present(off1)) off1_ = off1

    ! optional off2
    off2_ = off1_
    if (present(off2)) then
      m4_assert(present(off1))
      m4_assert(size(off2) == g2_%idx_dim)
      off2_ = off2
    end if

    ! init base
    call this%sparse_stencil_init(product(off2_-off1_+1))

    ! optional perm
    if (associated(g2_, target=g1)) then
      allocate (this%perm(g1%idx_dim), source=[(i, i=1,g1%idx_dim)])
    else
      m4_assert(present(perm))
    end if
    if (present(perm)) then
      m4_assert(size(perm) == g2_%idx_dim)
      if (g2_%idx_dim > 0) then
        m4_assert(minval(perm) >= -1)
        m4_assert(maxval(perm) <= g1%idx_dim)
      end if
      this%perm = perm
    end if

    ! save off1_, off2_
    allocate (this%off1(size(off1_)), source=off1_)
    allocate (this%off2(size(off2_)), source=off2_)
  end subroutine

  subroutine dirichlet_stencil_get(this, idx1, j, idx2, status)
    !! get j-th dependency wrt idx1
    class(dirichlet_stencil), intent(in)  :: this
    integer,                  intent(in)  :: idx1(:)
      !! result indices. size: idx_dim1
    integer,                  intent(in)  :: j
      !! j-th dependency
    integer,                  intent(out) :: idx2(:)
      !! output j-th dependency indices. size: idx_dim2
    logical, optional,        intent(out) :: status
      !! is j-th dependency used?

    integer :: i, k, div, rem

    m4_assert(j > 0)

    ! status
    if (present(status)) then
      status = ((j >= 1) .and. (j <= this%nmax))
      if (.not. status) return
    end if

    ! idx2 <- perm(idx1)
    where (this%perm > 0)
      idx2 = idx1(this%perm)
    elsewhere
      idx2 = 0
    end where

    ! idx2 += off1
    idx2 = idx2 + this%off1

    ! idx2 += 0:(off2-off1)
    if (all(this%off1 == this%off2)) return
    k = j
    do i = 1, size(this%off1)
      div     = k / (  this%off2(i)-this%off1(i)+1)
      rem     = mod(k, this%off2(i)-this%off1(i)+1)
      idx2(i) = idx2(i) + rem
      k       = div
    end do
  end subroutine

  subroutine near_neighb_stencil_init(this, g, idx1_type, idx1_dir, idx2_type, idx2_dir)
    !! initialize nearest neighbours stencil
    class(near_neighb_stencil), intent(out) :: this
    class(grid), target,        intent(in)  :: g
      !! grid this stencil is defined on
    integer,                    intent(in)  :: idx1_type
      !! result index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,                    intent(in)  :: idx1_dir
      !! result index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    integer,                    intent(in)  :: idx2_type
      !! dependency index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,                    intent(in)  :: idx2_dir
      !! dependency index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)

    integer :: shift

    shift = 0
    if ((this%idx1_type == this%idx2_type) .and. (this%idx1_dir == this%idx2_dir)) shift = 1

    call this%sparse_stencil_init(g%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir) + shift)

    ! set members
    this%g => g
    this%idx1_type = idx1_type
    this%idx1_dir  = idx1_dir
    this%idx2_type = idx2_type
    this%idx2_dir  = idx2_dir
  end subroutine

  subroutine near_neighb_stencil_get(this, idx1, j, idx2, status)
    !! get j-th dependency wrt idx1
    !! note that j-th dependency might not be used which is why there is a status:
    !!    e.g. for boundary boxes actual number of dependencies might be smaller than standard number of dependenc == nmax
    class(near_neighb_stencil), intent(in)  :: this
    integer,                    intent(in)  :: idx1(:)
      !! result indices. size: idx_dim1
    integer,                    intent(in)  :: j
      !! j-th dependency
    integer,                    intent(out) :: idx2(:)
      !! output j-th dependency indices. size: idx_dim2
    logical, optional,          intent(out) :: status
      !! is j-th dependency used?

    integer :: shift
    logical :: status_

    m4_assert(j > 0)

    shift = 0
    if ((this%idx1_type == this%idx2_type) .and. (this%idx1_dir == this%idx2_dir)) shift = 1

    if (j-shift == 0) then
      ! couple to self
      idx2 = idx1
      status_ = .true.
    else
      ! grid neighbours do not contain self => subtract 1 from j
      call this%g%get_neighb(this%idx1_type, this%idx1_dir, this%idx2_type, this%idx2_dir, idx1, j-shift, idx2, status_)
    end if

    if (present(status)) status = status_
  end subroutine

end module
