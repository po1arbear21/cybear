module example_continuity_m

  use matrix_m,                  only: matrix_convert, sparse_real
  use example_contact_m,         only: contacts, uncontacted, grd_contacts
  use example_current_density_m, only: current_dens
  use example_density_m,         only: dens
  use example_device_m,          only: adj_v, grd, n_intrin
  use example_imref_m,           only: iref
  use example_voltage_m,         only: volt
  use grid_m,                    only: IDX_VERTEX, IDX_EDGE
  use jacobian_m,                only: jacobian, jacobian_ptr
  use res_equation_m,            only: res_equation
  use stencil_m,                 only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,               only: vselector

  implicit none

  private
  public contin

  type, extends(res_equation) :: continuity
    !! dn/dt + dj/dx = 0 <=> dx * dn/dt + dj = 0
    type(vselector)           :: dens
      !! main variable
    type(vselector)           :: current_dens
      !! dependency

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn
    type(empty_stencil)       :: st_em

    type(jacobian), pointer   :: jaco_current_dens => null()
    type(jacobian), pointer   :: jaco_dens         => null()
    type(jacobian), pointer   :: jaco_dens_t       => null()
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

  type(continuity) :: contin

contains

  subroutine continuity_init(this)
    class(continuity), intent(out) :: this

    integer :: i, idx0(1), idx1(1)

    ! init res_equation
    call this%equation_init("continuity")

    ! vselect variables differing between contacts and uncontacted
    call this%current_dens%init(current_dens)
    call this%dens%init(dens,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])

    ! set main variable
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_nn%init( grd, IDX_VERTEX, 0, IDX_EDGE, 1)
    call this%st_em%init()

    ! init jacos
    this%jaco_current_dens => this%init_jaco_f(this%depend(this%current_dens),  [this%st_nn%get_ptr(),  (this%st_em%get_ptr(),  i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_dens         => this%init_jaco_f(this%depend(this%dens),          [this%st_em%get_ptr(),  (this%st_dir%get_ptr(), i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_dens_t       => this%init_jaco_f(this%depend(this%dens),          [this%st_dir%get_ptr(), (this%st_em%get_ptr(),  i=1, size(contacts))], const = .true.,  dtime = .true.)

    ! loop over vertices
    do i = 1, size(grd%x)
      idx0 = [i-1]
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        ! setting contacted const jacobian entries
        call this%jaco_current_dens%add( idx1, idx0, -1.0)
        call this%jaco_current_dens%add( idx1, idx1,  1.0)
        call this%jaco_dens_t%set(       idx1, idx1, adj_v%get(idx1))
      else
        ! setting uncontacted const jacobian entries
        call this%jaco_dens%set(idx1, idx1, 1.0)
      end if
    end do
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    class(continuity), intent(inout) :: this

    real, allocatable :: tmp(:)
    integer           :: i, j, idx(1)

    allocate(tmp(this%dens%n))

    ! calculating the density
    call this%jaco_current_dens%matr%mul_vec(this%current_dens%get(), tmp)
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp)

    do i = 1, size(contacts)
      do j = 1, contacts(i)%conts%n
        idx = contacts(i)%conts%get_idx(j)
        call this%f%update(idx, [-n_intrin * exp(contacts(i)%phi_ms)])
      end do
    end do
  end subroutine

end module
