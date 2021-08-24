module example_continuity_m

  use matrix_m,                  only: matrix_convert, sparse_real
  use example_contact_m,         only: contacts, uncontacted, grd_cont
  use example_current_density_m, only: curr
  use example_density_m,         only: dens
  use example_device_m,          only: grd, n_intrin, adj_v
  use example_imref_m,           only: iref
  use example_voltage_m,         only: cont_v
  use grid_m,                    only: IDX_VERTEX, IDX_EDGE
  use jacobian_m,                only: jacobian, jacobian_ptr
  use res_equation_m,            only: res_equation
  use stencil_m,                 only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,               only: vselector

  implicit none

  private
  public cont

  type, extends(res_equation) :: continuity
    !! dn/dt + dj/dx = 0 <=> dx * dn/dt + dj = 0
    type(vselector)           :: dens
      !! main variable
    type(vselector)           :: curr
      !! dependency

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn
    type(empty_stencil)       :: st_em

    type(jacobian), pointer   :: jaco_curr   => null()
    type(jacobian), pointer   :: jaco_dens   => null()
    type(jacobian), pointer   :: jaco_dens_t => null()
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

  type(continuity) :: cont

contains

  subroutine continuity_init(this)
    class(continuity), intent(out) :: this

    integer           :: i, idx0(1), idx1(1)
    type(sparse_real) :: sp

    ! init res_equation
    call this%equation_init("continuity")

    ! vselect variables differing between contacts and uncontacted
    call this%curr%init(curr)
    call this%dens%init(dens,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])

    ! set main variable
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_nn%init( grd, IDX_VERTEX, 0, IDX_EDGE, 1)
    call this%st_em%init()

    ! init jaco
    this%jaco_curr    => this%init_jaco_f(this%depend(this%curr),  [this%st_nn%get_ptr(),  (this%st_em%get_ptr(),  i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_dens    => this%init_jaco_f(this%depend(this%dens),  [this%st_em%get_ptr(), (this%st_dir%get_ptr(), i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_dens_t  => this%init_jaco_f(this%depend(this%dens),  [this%st_dir%get_ptr(), (this%st_em%get_ptr(), i=1, size(contacts))], const = .true.,  dtime = .true.)

    ! loop over vertices
    do i = 1, size(grd%x)
      idx0 = [i-1]
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        ! setting uncontacted const jacobs
        call this%jaco_dens%set(idx1, idx1, 1.0)
      else
        ! setting contacted const jacobs
        call this%jaco_curr%add( idx1, idx0, -1.0)
        call this%jaco_curr%add( idx1, idx1,  1.0)
        call this%jaco_dens_t%set( idx1, idx1, adj_v%get(idx1))
      end if
    end do
    call this%init_final()

    call matrix_convert(this%jaco_curr%matr, sp)
    call sp%output("jaco_curr.csv")
    call matrix_convert(this%jaco_dens%matr, sp)
    call sp%output("jaco_dens.csv")
  end subroutine

  subroutine continuity_eval(this)
    class(continuity), intent(inout) :: this

    real, allocatable :: tmp(:)
    integer :: i, j, idx(1)

    allocate(tmp(this%dens%n))
    call this%jaco_curr%matr%mul_vec(this%curr%get(), tmp)
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
