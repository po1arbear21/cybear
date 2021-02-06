#include "../util/macro.f90.inc"

module stencil_m
  implicit none

  type, abstract :: stencil
    !! abstract stencil
    integer :: max_ndep
      !! maximum number of dependency points
  contains
    procedure                            :: stencil_init
    procedure(stencil_get_dep), deferred :: get_dep
  end type

  abstract interface
    subroutine stencil_get_dep(this, idx1, idx2, ndep)
      !! get list of dependency points
      import stencil
      class(stencil), intent(in)  :: this
      integer,        intent(in)  :: idx1(:)
        !! result indices (idx_dim1)
      integer,        intent(out) :: idx2(:,:)
        !! output dependency indices (idx_dim2 x max_ndep)
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
    procedure :: init    => dirichlet_stencil_init
    procedure :: get_dep => dirichlet_stencil_get_dep
  end type

  ! TODO: nearest_neighbour_stencil

contains

  subroutine stencil_init(this, max_ndep)
    !! initialize base stencil
    class(stencil), intent(out) :: this
    integer,        intent(in)  :: max_ndep
      !! maximum number of dependency points

    this%max_ndep = max_ndep
  end subroutine

  subroutine dirichlet_stencil_init(this)
    !! initialize dirichlet stencil
    class(dirichlet_stencil), intent(out) :: this

    ! init base
    call this%stencil_init(1)
  end subroutine

  subroutine dirichlet_stencil_get_dep(this, idx1, idx2, ndep)
    !! get list of dependency points
    class(dirichlet_stencil), intent(in)  :: this
    integer,                  intent(in)  :: idx1(:)
      !! result indices (idx_dim1)
    integer,                  intent(out) :: idx2(:,:)
      !! output dependency indices (idx_dim2 x max_ndep)
    integer,                  intent(out) :: ndep
      !! output actual number of dependency points

    IGNORE(this)

    ! simply copy indices
    idx2(:,1) = idx1
    ndep = 1
  end subroutine

end module