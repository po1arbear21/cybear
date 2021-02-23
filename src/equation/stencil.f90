#include "../util/macro.f90.inc"

module stencil_m

  use grid_m, only: grid

  implicit none

  private
  public stencil, stencil_ptr
  public dynamic_stencil
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

  type, abstract, extends(stencil) :: static_stencil
    !! abstract stationary stencil
    !! derive from it for specific equations (or use dirichlet/nearest neighb stencils).
    !! useful when exact stencil is known.
    integer :: nmax
      !! maximum number of dependency points
  contains
    procedure                               :: static_stencil_init
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
    !! dirichlet stencil, selects only one grid point
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

  subroutine static_stencil_init(this, nmax)
    !! initialize base stencil
    class(static_stencil), intent(out) :: this
    integer,               intent(in)  :: nmax
      !! maximum number of dependency points

    this%nmax = nmax
  end subroutine

  subroutine dirichlet_stencil_init(this)
    !! initialize dirichlet stencil
    class(dirichlet_stencil), intent(out) :: this

    ! init base
    call this%static_stencil_init(1)
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

    IGNORE(this)

    ! only one dependency point
    idx2   = idx1
    status = (j == 1)
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

    call this%static_stencil_init(g%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir))

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
