m4_include(util/macro.f90.inc)

module approx_m

  use charge_density_m, only: charge_density
  use contact_m,        only: CT_OHMIC
  use device_params_m,  only: device_params
  use error_m,          only: assert_failed, program_error
  use esystem_m,        only: esystem
  use grid_m,           only: IDX_VERTEX, IDX_CELL, IDX_EDGE
  use grid0D_m,         only: get_dummy_grid
  use ieee_arithmetic,  only: ieee_is_finite
  use imref_m,          only: imref
  use jacobian_m,       only: jacobian
  use math_m,           only: eye_real
  use matrix_m,         only: block_real
  use normalization_m,  only: norm
  use poisson_m,        only: poisson
  use potential_m,      only: potential
  use res_equation_m,   only: res_equation
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON
  use stencil_m,        only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use voltage_m,        only: voltage
  use vselector_m,      only: vselector

  implicit none

  private
  public :: approx_imref, approx_potential

  type, extends(res_equation) :: imref_approx_eq
    !! equation to approximate initial imref using density ~ doping concentration and mob ~ mob_0
    !! div (mob_0 * max(n_intrin,dcon) * grad(iref)) == 0

    type(device_params), pointer :: par => null()
      !! pointer to device parameters

    type(vselector) :: iref
      !! main variable
    type(vselector) :: volt
      !! dependency

    type(dirichlet_stencil)   :: st_dir
    type(dirichlet_stencil)   :: st_dir_volt
    type(empty_stencil)       :: st_em
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer   :: jaco_iref => null()
    type(jacobian), pointer   :: jaco_volt => null()
  contains
    procedure :: init => imref_approx_eq_init
    procedure :: eval => imref_approx_eq_eval
  end type

contains

  subroutine approx_imref(par, iref, volt)
    !! approximate imref based on doping and low-field mobility
    type(device_params), intent(in)    :: par
      !! device parameters
    type(imref),         intent(inout) :: iref
      !! imref variable
    type(voltage),       intent(in)    :: volt(:)
      !! voltages for all contacts

    integer                   :: ict
    real, allocatable         :: dx(:), f(:)
    type(imref_approx_eq)     :: eq
    type(esystem)             :: sys
    type(block_real), pointer :: df

    ! initialize approximation equation
    call eq%init(par, iref, volt)

    ! init system
    call sys%init("imref approximation")
    call sys%add_equation(eq)
    do ict = 1, size(par%contacts)
      call sys%provide(volt(ict))
    end do
    call sys%init_final()

    ! solve for imref
    allocate (dx(sys%n), f(sys%n))
    call sys%eval(f = f, df = df)
    call df%factorize()
    call df%solve_vec(f, dx)
    call sys%set_x(sys%get_x() - dx)

    ! free memory
    call eq%destruct()
    call sys%destruct()
  end subroutine

  subroutine approx_potential(par, pot, iref)
    !! approximate potential in transport region based on doping and previously approximated imrefs
    type(device_params),  intent(in)    :: par
      !! device parameters
    type(potential),      intent(inout) :: pot
      !! potential variable
    type(imref),          intent(in)    :: iref(:)


    integer               :: i, ci
    real                  :: ireff(2), dop(2), dop_eff, L, p
    integer, allocatable  :: idx(:)

    ! loop over uncontacted transport vertices
    allocate (idx(par%g%idx_dim))
    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx = par%transport(IDX_VERTEX,0)%get_idx(i)

      ! get doping and imrefs
      dop(DOP_DCON) = par%dop(IDX_VERTEX,0,DOP_DCON)%get(idx)
      dop(DOP_ACON) = par%dop(IDX_VERTEX,0,DOP_ACON)%get(idx)
      do ci = par%ci0, par%ci1
        ireff(ci) = iref(ci)%get(idx)
      end do
      dop_eff = dop(DOP_DCON) - dop(DOP_ACON)

      ! estimate potential
      if ((par%ci0 == CR_ELEC) .and. (par%ci1 == CR_HOLE)) then
        p = 0
        if (dop_eff /= 0) then
          L = 0.5 * par%smc%band_gap + log(abs(dop_eff) / sqrt(par%smc%edos(1) * par%smc%edos(2)))
          if (L < 9) then
            p = asinh(0.5 * exp(L))
          else
            p = L + exp(-2 * L)
          end if
        end if
        call pot%set(idx, 0.5 * (ireff(CR_ELEC) + ireff(CR_HOLE)) + sign(p, dop_eff))
      elseif (par%ci0 == CR_ELEC) then
        if (dop(DOP_DCON) > 0) then
          call pot%set(idx, ireff(CR_ELEC) + 0.5 * par%smc%band_gap + log(dop(DOP_DCON) / sqrt(par%smc%edos(CR_ELEC) * par%smc%edos(CR_HOLE))))
        else
          call pot%set(idx, ireff(CR_ELEC))
        end if
      elseif (par%ci0 == CR_HOLE) then
        if (dop(DOP_ACON) > 0) then
          call pot%set(idx, ireff(CR_HOLE) + 0.5 * par%smc%band_gap - log(dop(DOP_ACON) / sqrt(par%smc%edos(CR_ELEC) * par%smc%edos(CR_HOLE))))
        else
          call pot%set(idx, ireff(CR_HOLE))
        end if
      end if
    end do

    ! safety check
    m4_assert(all(ieee_is_finite(pot%get())))
  end subroutine

  subroutine imref_approx_eq_init(this, par, iref, volt)
    !! initialize imref approximation equation
    class(imref_approx_eq),      intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(imref),                 intent(in)  :: iref
      !! imref variable
    type(voltage),               intent(in)  :: volt(:)
      !! voltages for all contacts

    integer              :: i, idx_dir, ict, dum(0)
    real                 :: eps, cap, ndum
    real,    allocatable :: d_volt(:,:), eye(:,:)
    integer, allocatable :: idx1(:), idx2(:), idx(:)
    logical              :: status

    ndum = norm(1e-10, "1/cm^3")

    ! init base
    call this%equation_init("imref_approx_eq")
    this%par => par
    allocate (eye(par%nct,par%nct))
    eye = eye_real(par%nct)

    ! create variable selectors
    call this%iref%init(iref, [par%transport_vct( 0)%get_ptr(), &                   ! uncontacted transport vertices
      &                      (par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)]) ! contacted vertices
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, par%nct)], "voltages")

    ! init residuals using this%iref as main variable
    call this%init_f(this%iref)

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_dir_volt%init(par%g, g2 = get_dummy_grid(), perm = dum)
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()

    ! init jacobians
    this%jaco_iref  => this%init_jaco_f(this%depend(this%iref ), st = [this%st_nn%get_ptr(), &
    & (this%st_dir%get_ptr(),      ict = 1, par%nct)], const = .true.)
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), st = [this%st_em%get_ptr(),   &
    & (this%st_dir_volt%get_ptr(), ict = 1, par%nct)], const = .true.)

    ! loop over transport cells
    allocate (idx(par%g%idx_dim), idx1(par%g%idx_dim), idx2(par%g%idx_dim))
    do idx_dir = 1, par%g%idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get "permittivity" (= mob * dens)
        eps = par%mob0(IDX_EDGE,idx_dir,iref%ci)%get(idx) * max(ndum, par%dop(IDX_EDGE,idx_dir,iref%ci)%get(idx))

        ! "capacitance" (= surf * mob * dens / len)
        cap = par%tr_surf(idx_dir)%get(idx) * eps / par%g%get_len(idx, idx_dir)

        if (par%ict%get(idx1) == 0) then
          call this%jaco_iref%add(idx1, idx1,  cap)
          call this%jaco_iref%add(idx1, idx2, -cap)
        end if
        if (par%ict%get(idx2) == 0) then
          call this%jaco_iref%add(idx2, idx1, -cap)
          call this%jaco_iref%add(idx2, idx2,  cap)
        end if
      end do
    end do

    ! loop over contacted vertices
    do ict = 1, par%nct
      d_volt = - eye(ict:ict,:)
      do i = 1, par%transport_vct(ict)%n
        idx1 = par%transport_vct(ict)%get_idx(i)

        ! set jaco_volt entry
        call this%jaco_volt%set(idx1, dum, d_volt)

        ! set jaco_iref entry
        call this%jaco_iref%set(idx1, idx1, 1.0)
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine imref_approx_eq_eval(this)
    !! evaluate imref approximation equation
    class(imref_approx_eq), intent(inout) :: this

    real, allocatable :: tmp(:)

    allocate (tmp(this%iref%n))

    ! calculate residuals
    call this%jaco_iref%matr%mul_vec(this%iref%get(), tmp              )
    call this%jaco_volt%matr%mul_vec(this%volt%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp)
  end subroutine

end module
