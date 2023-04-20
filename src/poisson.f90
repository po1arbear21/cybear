module poisson_m

  use charge_density_m, only: charge_density
  use device_params_m,  only: device_params
  use grid_m,           only: IDX_VERTEX, IDX_CELL, IDX_EDGE
  use grid0D_m,         only: get_dummy_grid
  use jacobian_m,       only: jacobian
  use math_m,           only: eye_real
  use potential_m,      only: potential
  use res_equation_m,   only: res_equation
  use stencil_m,        only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use voltage_m,        only: voltage
  use vselector_m,      only: vselector

  implicit none

  private
  public poisson

  type, extends(res_equation) :: poisson
    !! Quasi-Stationary Poisson Equation

    type(device_params), pointer :: par => null()
      !! pointer to device parameters

    type(vselector) :: pot
      !! main variable
    type(vselector) :: rho
      !! dependency
    type(vselector) :: volt
      !! dependency

    type(dirichlet_stencil)   :: st_dir
    type(dirichlet_stencil)   :: st_dir_volt
    type(empty_stencil)       :: st_em
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_rho  => null()
    type(jacobian), pointer :: jaco_volt => null()
  contains
    procedure :: init => poisson_init
    procedure :: eval => poisson_eval
  end type

contains

  subroutine poisson_init(this, par, pot, rho, volt)
    !! initialize poisson equation
    class(poisson),              intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(potential),             intent(in)  :: pot
      !! potential variable
    type(charge_density),        intent(in)  :: rho
      !! charge density variable
    type(voltage),               intent(in)  :: volt(:)
      !! voltages for all contacts

    integer               :: i,  idx_dir, ict, dum(0)
    real                  :: cap
    real, allocatable     :: d_volt(:,:), eye(:,:)
    integer, allocatable  :: idx1(:), idx2(:), idx(:)
    logical               :: status

    ! init base
    call this%equation_init("poisson")
    this%par => par
    eye = eye_real(size(par%contacts))

    ! create variable selectors
    call this%pot%init(pot, [par%oxide_vct(     0)%get_ptr(), &                              ! uncontacted oxide vertices
      &                      par%transport_vct( 0)%get_ptr(), &                              ! uncontacted transport vertices
      &                      (par%poisson_vct(ict)%get_ptr(), ict = 1, size(par%contacts))]) ! contacted vertices
    call this%rho%init(rho, par%transport(IDX_VERTEX,0))
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, size(volt))], "voltages")

    ! init residuals using this%pot as main variable
    call this%init_f(this%pot)

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_dir_volt%init(par%g, g2 = get_dummy_grid(), perm = dum)
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()

    ! init jacos
    this%jaco_pot  => this%init_jaco_f(this%depend(this%pot ), st = [this%st_nn%get_ptr(), this%st_nn%get_ptr(),  &
      & (this%st_dir%get_ptr(),      ict = 1, size(par%contacts))], const = .true.)
    this%jaco_rho  => this%init_jaco_f(this%depend(this%rho ), st = [this%st_em%get_ptr(), this%st_dir%get_ptr(), &
      & (this%st_em%get_ptr(),       ict = 1, size(par%contacts))], const = .true.)
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), st = [this%st_em%get_ptr(), this%st_em%get_ptr(),  &
      & (this%st_dir_volt%get_ptr(), ict = 1, size(par%contacts))], const = .true.)

    allocate (idx1(par%g%idx_dim), idx2(par%g%idx_dim))

    ! loop over poisson edges
    do idx_dir = 1, par%g%idx_dim
      do i = 1, par%poisson(IDX_EDGE,idx_dir)%n
        idx = par%poisson(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! capacitance
        cap = par%eps(IDX_EDGE,idx_dir)%get(idx) * par%surf(idx_dir)%get(idx) / par%g%get_len(idx, idx_dir)

        if (par%ict%get(idx1) == 0) then
          call this%jaco_pot%add(idx1, idx1,  cap)
          call this%jaco_pot%add(idx1, idx2, -cap)
        end if
        if (par%ict%get(idx2) == 0) then
          call this%jaco_pot%add(idx2, idx1, -cap)
          call this%jaco_pot%add(idx2, idx2,  cap)
        end if
       end do
    end do

    ! loop over uncontacted transport vertices
    do i = 1, par%transport_vct(0)%n
      idx1 = par%transport_vct(0)%get_idx(i)

      ! set jaco_rho entries
      call this%jaco_rho%set(idx1, idx1, - par%tr_vol%get(idx1))
    end do

    ! loop over contacted vertices
    do ict = 1, size(par%contacts)
      d_volt = - eye(ict:ict,:)
      do i = 1, par%poisson_vct(ict)%n
        idx1 = par%poisson_vct(ict)%get_idx(i)

        ! set jaco_volt entry
        call this%jaco_volt%set(idx1, dum, d_volt)

        ! set jaco_pot entry
        call this%jaco_pot%set(idx1, idx1, 1.0)
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine poisson_eval(this)
    !! evaluate poisson equation
    class(poisson), intent(inout) :: this

    real,    allocatable  :: tmp(:)
    integer               :: ict, i
    integer, allocatable  :: idx(:)

    allocate (tmp(this%pot%n))
    ! calculate residuals (without phims)
    call this%jaco_pot%matr%mul_vec( this%pot%get(),  tmp              )
    call this%jaco_rho%matr%mul_vec( this%rho%get(),  tmp, fact_y = 1.0)
    call this%jaco_volt%matr%mul_vec(this%volt%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp)

    ! phims
    do ict = 1, size(this%par%contacts)
      do i = 1, this%par%poisson_vct(ict)%n
        idx = this%par%poisson_vct(ict)%get_idx(i)
        call this%f%update(idx, [-this%par%contacts(ict)%phims])
      end do
    end do
  end subroutine

end module
