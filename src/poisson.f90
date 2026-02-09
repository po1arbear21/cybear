module poisson_m

  use charge_density_m,  only: charge_density
  use device_params_m,   only: device_params
  use grid_m,            only: IDX_VERTEX, IDX_CELL, IDX_EDGE
  use grid0D_m,          only: get_dummy_grid
  use jacobian_m,        only: jacobian
  use math_m,            only: eye_real
  use potential_m,       only: potential
  use res_equation_m,    only: res_equation
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use surface_charge_m,  only: surface_charge
  use voltage_m,         only: voltage
  use vselector_m,       only: vselector

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
    type(vselector) :: scharge
      !! surface charge dependency (optional, for charged surface model)

    type(dirichlet_stencil)   :: st_dir
    type(dirichlet_stencil)   :: st_dir_volt
    type(empty_stencil)       :: st_em
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot     => null()
    type(jacobian), pointer :: jaco_rho     => null()
    type(jacobian), pointer :: jaco_volt    => null()
    type(jacobian), pointer :: jaco_scharge => null()
      !! Jacobian for surface charge contribution (optional)
  contains
    procedure :: init => poisson_init
    procedure :: eval => poisson_eval
  end type

contains

  subroutine poisson_init(this, par, pot, rho, volt, scharge)
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
    type(surface_charge), optional, intent(in) :: scharge
      !! surface charge variable (for charged surface model)

    integer               :: i, idx_dir, ict, dum(0), n_surf_vert
    real                  :: cap, A_surf
    real, allocatable     :: d_volt(:,:), eye(:,:)
    integer, allocatable  :: idx1(:), idx2(:), idx(:)
    logical               :: status, has_scharge

    print "(A)", "poisson_init"
    has_scharge = present(scharge)

    ! init base
    call this%equation_init("poisson")
    this%par => par

    allocate (idx(par%g%idx_dim), idx1(par%g%idx_dim), idx2(par%g%idx_dim), eye(par%nct,par%nct))
    eye = eye_real(par%nct)

    ! create variable selectors
    call this%pot%init(pot, [par%oxide_vct(     0)%get_ptr(), &                   ! uncontacted oxide vertices
      &                      par%transport_vct( 0)%get_ptr(), &                   ! uncontacted transport vertices
      &                      (par%poisson_vct(ict)%get_ptr(), ict = 1, par%nct)]) ! contacted vertices
    call this%rho%init(rho, par%transport_vct(0))
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, par%nct)], "voltages")

    ! init residuals using this%pot as main variable
    call this%init_f(this%pot)

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_dir_volt%init(par%g, g2 = get_dummy_grid(), perm = dum)
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()

    ! init jacobians
    this%jaco_pot  => this%init_jaco_f(this%depend(this%pot ), st = [this%st_nn%get_ptr(), this%st_nn%get_ptr(),  &
      & (this%st_dir%get_ptr(),      ict = 1, par%nct)], const = .true.)
    this%jaco_rho  => this%init_jaco_f(this%depend(this%rho ), st = [this%st_em%get_ptr(), this%st_dir%get_ptr(), &
      & (this%st_em%get_ptr(),       ict = 1, par%nct)], const = .true.)
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), st = [this%st_em%get_ptr(), this%st_em%get_ptr(),  &
      & (this%st_dir_volt%get_ptr(), ict = 1, par%nct)], const = .true.)

    ! surface charge jacobian (optional, for charged surface model)
    if (has_scharge) then
      call this%scharge%init(scharge, par%transport_vct(0))
      this%jaco_scharge => this%init_jaco_f(this%depend(this%scharge), &
        & st = [this%st_em%get_ptr(), this%st_dir%get_ptr(), &
        &       (this%st_em%get_ptr(), ict = 1, par%nct)], const = .true.)
    end if

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
      call this%jaco_rho%add(idx1, idx1, - par%tr_vol%get(idx1))
    end do

    ! set surface charge jacobian entries at surface vertices
    ! Surface charge contributes: F += -rho_surf * A_surf
    ! (negative because positive surface charge increases potential, reducing residual)
    if (has_scharge) then
      n_surf_vert = 0
      do i = 1, par%transport(IDX_VERTEX, 0)%n
        idx1 = par%transport(IDX_VERTEX, 0)%get_idx(i)
        ! Surface vertices: x=1 or x=nx, not on contact
        if ((idx1(1) == 1 .or. idx1(1) == par%g1D(1)%n) .and. par%ict%get(idx1) == 0) then
          n_surf_vert = n_surf_vert + 1

          ! Calculate adjoint surface area (y-extent for x-boundary)
          if (idx1(2) == 1) then
            A_surf = 0.5 * (par%g1D(2)%x(2) - par%g1D(2)%x(1))
          elseif (idx1(2) == par%g1D(2)%n) then
            A_surf = 0.5 * (par%g1D(2)%x(par%g1D(2)%n) - par%g1D(2)%x(par%g1D(2)%n - 1))
          else
            A_surf = 0.5 * (par%g1D(2)%x(idx1(2) + 1) - par%g1D(2)%x(idx1(2) - 1))
          end if

          ! dF/d(rho_surf) = -A_surf
          call this%jaco_scharge%add(idx1, idx1, -A_surf)
        end if
      end do
      print "(A,I0,A)", "  Surface charge: ", n_surf_vert, " surface vertices in Poisson"
    end if

    ! loop over contacted vertices
    do ict = 1, par%nct
      d_volt = - eye(ict:ict,:)
      do i = 1, par%poisson_vct(ict)%n
        idx1 = par%poisson_vct(ict)%get_idx(i)

        ! set jaco_volt entry
        call this%jaco_volt%add(idx1, dum, d_volt)

        ! set jaco_pot entry
        call this%jaco_pot%add(idx1, idx1, 1.0)
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

    allocate (idx(this%par%g%idx_dim), tmp(this%pot%n))

    ! calculate residuals (without phims)
    call this%jaco_pot%matr%mul_vec( this%pot%get(),  tmp              )
    call this%jaco_rho%matr%mul_vec( this%rho%get(),  tmp, fact_y = 1.0)
    call this%jaco_volt%matr%mul_vec(this%volt%get(), tmp, fact_y = 1.0)
    if (associated(this%jaco_scharge)) then
      call this%jaco_scharge%matr%mul_vec(this%scharge%get(), tmp, fact_y = 1.0)
    end if
    call this%f%set(tmp)

    ! phims
    do ict = 1, this%par%nct
      do i = 1, this%par%poisson_vct(ict)%n
        idx = this%par%poisson_vct(ict)%get_idx(i)
        call this%f%update(idx, [-this%par%contacts(ict)%phims])
      end do
    end do
  end subroutine

end module
