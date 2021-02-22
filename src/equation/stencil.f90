#include "../util/macro.f90.inc"

module stencil_m

  use grid_m, only: grid

  implicit none

  private
  public stencil, stencil_ptr
  public dirichlet_stencil

  type, abstract :: stencil
    !! abstract stencil
    integer :: nmax
      !! maximum number of dependency points
  contains
    procedure                        :: stencil_init
    procedure                        :: get_ptr => stencil_get_ptr
    procedure(stencil_get), deferred :: get
  end type

  abstract interface
    subroutine stencil_get(this, idx1, idx2, ndep)
      !! get list of dependency points
      import stencil
      class(stencil), intent(in)  :: this
      integer,        intent(in)  :: idx1(:)
        !! result indices (idx_dim1)
      integer,        intent(out) :: idx2(:,:)
        !! output dependency indices (idx_dim2 x nmax)
      integer,        intent(out) :: ndep
        !! output actual number of dependency points
    end subroutine
  end interface

  type stencil_ptr
    class(stencil), pointer :: p => null()
  end type

  type, extends(stencil) :: dirichlet_stencil
    !! dirichlet stencil, selects only one grid point
  contains
    procedure :: init => dirichlet_stencil_init
    procedure :: get  => dirichlet_stencil_get
  end type

  type, extends(stencil) :: near_neighb_stencil
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

  subroutine stencil_init(this, nmax)
    !! initialize base stencil
    class(stencil), intent(out) :: this
    integer,        intent(in)  :: nmax
      !! maximum number of dependency points

    this%nmax = nmax
  end subroutine

  function stencil_get_ptr(this) result(ptr)
    !! returns pointer type to this stencil
    class(stencil), target, intent(in) :: this
    type(stencil_ptr)                  :: ptr

    ptr%p => this
  end function

  subroutine dirichlet_stencil_init(this)
    !! initialize dirichlet stencil
    class(dirichlet_stencil), intent(out) :: this

    ! init base
    call this%stencil_init(1)
  end subroutine

  subroutine dirichlet_stencil_get(this, idx1, idx2, ndep)
    !! get list of dependency points
    class(dirichlet_stencil), intent(in)  :: this
    integer,                  intent(in)  :: idx1(:)
      !! result indices (idx_dim1)
    integer,                  intent(out) :: idx2(:,:)
      !! output dependency indices (idx_dim2 x nmax)
    integer,                  intent(out) :: ndep
      !! output actual number of dependency points

    IGNORE(this)

    ! simply copy indices
    idx2(:,1) = idx1
    ndep = 1
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

    call this%stencil_init(g%get_max_neighb(idx1_type, idx1_dir, idx2_type, idx2_dir))

    ! set members
    this%g => g
    this%idx1_type = idx1_type
    this%idx1_dir  = idx1_dir
    this%idx2_type = idx2_type
    this%idx2_dir  = idx2_dir
  end subroutine

  subroutine near_neighb_stencil_get(this, idx1, idx2, ndep)
    !! get list of dependency points
    class(near_neighb_stencil), intent(in)  :: this
    integer,                    intent(in)  :: idx1(:)
      !! result indices (idx_dim1)
    integer,                    intent(out) :: idx2(:,:)
      !! output dependency indices (idx_dim2 x nmax)
    integer,                    intent(out) :: ndep
      !! output actual number of dependency points

    call this%g%get_neighb(this%idx1_type, this%idx1_dir, this%idx2_type, this%idx2_dir, idx1, idx2, ndep)
  end subroutine

end module
