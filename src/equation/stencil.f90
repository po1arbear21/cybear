#include "../util/macro.f90.inc"

module stencil_m

  use error_m, only: assert_failed
  use grid_m,  only: grid
  use math_m,  only: eye_int

  implicit none

  private
  public stencil, stencil_ptr
  public dynamic_stencil
  public base_static_stencil
  public empty_stencil
  public static_stencil
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

  type, abstract, extends(stencil) :: base_static_stencil
    !! abstract static stencil.
    !! stencil w/ constant numbers of neighbors
    integer :: nmax
      !! maximum number of dependency points
  contains
    procedure, private :: base_static_stencil_init
  end type

  type, extends(base_static_stencil) :: empty_stencil
    !! empty stencil with no neighbors.
    !! indicates no dependence on this variable
  contains
    procedure :: init => empty_stencil_init
  end type

  type, abstract, extends(base_static_stencil) :: static_stencil
    !! abstract stationary stencil
    !! derive from it for specific equations (or use dirichlet/nearest neighb stencils).
    !! useful when exact stencil is known and actual neighbors exist (in contrast to empty_stencil).
  contains
    procedure(static_stencil_get), deferred :: get
  end type

  abstract interface
    subroutine static_stencil_get(this, idx1, j, idx2, status)
      !! get j-th dependency wrt idx1
      !! note that j-th dependency might not be used which is why there is a status:
      !!    e.g. for boundary boxes actual number of dependencies might be smaller than standard number of dependenc == nmax
      import static_stencil
      class(static_stencil), intent(in)  :: this
      integer,               intent(in)  :: idx1(:)
        !! result indices. size: idx_dim1
      integer,               intent(in)  :: j
        !! j-th dependency
      integer,               intent(out) :: idx2(:)
        !! output j-th dependency indices. size: idx_dim2
      logical,               intent(out) :: status
        !! is j-th dependency used?
    end subroutine
  end interface

  type, extends(stencil) :: dynamic_stencil
    !! placeholder to distinguish from static stencils.
  end type

  type, extends(static_stencil) :: dirichlet_stencil
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
    !!      off1 = (0, 0,   2, 2)
    !!      off2 = (0, N-1, 2, 2)

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

  type, extends(static_stencil) :: near_neighb_stencil
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

  subroutine base_static_stencil_init(this, nmax)
    !! initialize base stencil
    class(base_static_stencil), intent(out) :: this
    integer,                    intent(in)  :: nmax
      !! maximum number of dependency points

    this%nmax = nmax
  end subroutine

  subroutine empty_stencil_init(this)
    !! initialize empty stencil
    class(empty_stencil), intent(out) :: this

    call this%base_static_stencil_init(0)
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
    if (present(off1)) then
      off1_ = off1
    end if

    ! optional off2
    off2_ = off1_
    if (present(off2)) then
      ASSERT(present(off1))
      ASSERT(size(off2) == g2_%idx_dim)
      off2_ = off2
    end if

    ! init base
    call this%base_static_stencil_init(product(off2_-off1_+1))

    ! optional perm
    if (associated(g2_, target=g1)) then
      allocate (this%perm(g1%idx_dim), source=[(i, i=1,g1%idx_dim)])
    else
      ASSERT(present(perm))
    end if
    if (present(perm)) then
      ASSERT(size(perm) == g2_%idx_dim)
      ASSERT(minval(perm) >= -1)
      ASSERT(maxval(perm) <= g1%idx_dim)
      ASSERT(.not. any(perm == 0))
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
    logical,                  intent(out) :: status
      !! is j-th dependency used?

    integer :: i, k, div, rem

    ASSERT(j > 0)

    ! status
    status = ((j >= 1) .and. (j <= this%nmax))
    if (.not. status) return

    ! idx2 <- perm(idx1)
    do i = 1, size(this%perm)
      if (this%perm(i) == -1) then
        idx2(i) = 1
      else
        idx2(i) = idx1(this%perm(i))
      end if
    end do

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

    call this%base_static_stencil_init(g%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir))

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
    logical,                    intent(out) :: status
      !! is j-th dependency used?

    call this%g%get_neighb(this%idx1_type, this%idx1_dir, this%idx2_type, this%idx2_dir, idx1, j, idx2, status)
  end subroutine

end module
