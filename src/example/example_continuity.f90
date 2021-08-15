module example_continuity_m

  use example_contact_m,        only: contacts, uncontacted
  use example_current_m,        only: curr
  use example_device_m,         only: grd, n_intrin, adj_v
  use example_imref_m,          only: iref
  use example_voltage_m,        only: cont_v
  use grid_m,                   only: IDX_VERTEX, IDX_EDGE
  use jacobian_m,               only: jacobian, jacobian_ptr
  use res_equation_m,           only: res_equation
  use stencil_m,                only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,              only: vselector

  implicit none

  private
  public cont

  type, extends(res_equation) :: continuity
    !! dn/dt + dj/dx = 0 <=> dx * dn/dt + dj = 0
    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn
    type(empty_stencil)       :: st_em
    type(jacobian), pointer   :: jaco_adj_v => null()
    type(jacobian), pointer   :: jaco_curr  => null()
    type(jacobian), pointer   :: jaco_dens  => null()
    type(jacobian), pointer   :: jaco_iref  => null()
    type(jacobian), pointer   :: jaco_volt  => null()
    type(vselector)           :: adj_v
    type(vselector)           :: curr
    type(vselector)           :: dens
    type(vselector)           :: iref
    type(vselector)           :: volt
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

  type(continuity) :: cont

contains

  subroutine continuity_init(this)
    class(continuity), intent(out) :: this

    integer :: i, idx1, idx2

    ! init res_equation
    call this%equation_init("continuity")

    ! vselect variables differing between contacts and uncontacted
    call this%dens%init(adj_v,  [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%curr%init(curr,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%dens%init(dens,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%curr%init(iref,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%dens%init(cont_v, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])

    ! set main variable
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_nn%init( grd, IDX_VERTEX, 0, IDX_EDGE, 1)
    call this%st_em%init()

    ! init jaco
    this%jaco_adj_v => this%init_jaco_f(this%depend(this%adj_v), [this%st_dir%get_ptr(), (this%st_em%get_ptr(),  i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_curr  => this%init_jaco_f(this%depend(this%curr),  [this%st_nn%get_ptr(),  (this%st_em%get_ptr(),  i=1, size(contacts))], const = .true.,  dtime = .false.)
    this%jaco_dens  => this%init_jaco_f(this%depend(this%dens),  [this%st_dir%get_ptr(), (this%st_dir%get_ptr(), i=1, size(contacts))], const = .true.,  dtime = .true.)
    this%jaco_iref  => this%init_jaco_f(this%depend(this%iref),  [this%st_em%get_ptr(),  (this%st_dir%get_ptr(), i=1, size(contacts))], const = .false., dtime = .false.)
    this%jaco_volt  => this%init_jaco_f(this%depend(this%volt),  [this%st_em%get_ptr(),  (this%st_dir%get_ptr(), i=1, size(contacts))], const = .false., dtime = .false.)

    ! loop over cells
    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      if (contacted%flags%get(idx1)) then
        ! setting contacted const jacobs
        call this%jaco_adj_v%set(idx1, idx1, this%dens%get([idx1]))
        call this%jaco_curr%add( idx1, idx1, -1)
        call this%jaco_curr%add( idx1, idx2,  1)
        call this%jaco_dens%set( idx1, idx1, this%adj_v%get([idx1]))
      else
        ! setting uncontacted const jacobs
        call this%jaco_dens%set(idx1, idx1, 1)
      end if

      if (contacted%flags%get(idx2)) then
         ! setting contacted const jacobs
        call this%jaco_adj_v%set(idx2, idx2, this%dens%get([idx1]))
        call this%jaco_curr%add( idx2, idx1, -1)
        call this%jaco_curr%add( idx2, idx2,  1)
        call this%jaco_dens%set( idx2, idx2, this%adj_v%get([idx1]))
      else
        ! setting uncontacted const jacobs
        call this%jaco_dens%set(idx2, idx2, 1)
      end if
    end do

    ! loop over vertices
    do i = 1, size(grd%x)
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        call this%jaco_rho%set(idx1, idx1, -adj_v%get(idx1))
      else
        call this%jaco_pot%set(idx1, idx1, 1.0)
      end if
    end do
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    class(continuity), intent(out) :: this

    real, allocatable :: tmp(:)
    integer :: i, idx1

    allocate(tmp(this%dense%n))

    do i = 1, size(grd%x)-1
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        ! setting uncontacted jacobs
        call this%jaco_iref%set(idx1, idx1, -n_intrin*exp(iref%get(idx1)+cont_v(c_index)))
        call this%jaco_volt%set(idx1, idx1, -n_intrin*exp(iref%get(idx1)+cont_v(c_index)))

        ! setting for uncontacted dens
        call dens
      end if
    end do
  end subroutine

end module
