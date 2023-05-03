module approx_m

  use device_params_m, only: device_params, CR_ELEC, CR_HOLE
  use esystem_m,       only: esystem
  use grid_m,          only: IDX_VERTEX, IDX_CELL, IDX_EDGE
  use grid0D_m,        only: get_dummy_grid
  use imref_m,         only: imref
  use jacobian_m,      only: jacobian
  use math_m,          only: eye_real
  use matrix_m,        only: sparse_real
  use potential_m,     only: potential
  use res_equation_m,  only: res_equation
  use stencil_m,       only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use voltage_m,       only: voltage
  use vselector_m,     only: vselector

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

    integer               :: ict
    real, allocatable     :: x(:), f(:)
    type(imref_approx_eq) :: eq
    type(esystem)         :: sys
    type(sparse_real)     :: df

    ! initialize approximation equation
    call eq%init(par, iref, volt)

    ! init system
    call sys%init("imref approximation")
    call sys%add_equation(eq)
    do ict = 1, size(par%contacts)
      call sys%provide(volt(ict))
    end do
    call sys%init_final()

    ! memory
    allocate (x(sys%n), f(sys%n))

    ! evaluate system
    call sys%eval(f = f, df = df)

    ! solve
    call df%factorize()
    call df%solve_vec(-f, x)
    call df%destruct()

    ! save values
    call sys%set_x(x)
  end subroutine

  subroutine approx_potential(par, pot, iref)
    !! approximate potential in transport region based on doping and previously approximated imrefs
    type(device_params), intent(in)    :: par
      !! device parameters
    type(potential),     intent(inout) :: pot
      !! potential variable
    type(imref),         intent(in)    :: iref(:)
      !! electron/hole quasi-fermi potential

    integer               :: i, ci
    real                  :: ireff(2), dop(2)
    integer, allocatable  :: idx(:)

    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx = par%transport(IDX_VERTEX,0)%get_idx(i)

      ! get doping and imrefs
      do ci = par%ci0, par%ci1
        dop(ci) = par%dop(IDX_VERTEX,0,ci)%get(idx)
        ireff(ci) = iref(ci)%get(idx)
      end do

      ! estimate potential
      if ((par%ci0 == CR_ELEC) .and. (par%ci1 == CR_HOLE)) then
        call pot%set(idx, 0.5 * (ireff(CR_ELEC) + ireff(CR_HOLE)) + asinh(0.5 * (dop(CR_ELEC) - dop(CR_HOLE)) / par%smc%n_intrin))
      elseif (par%ci0 == CR_ELEC) then
        call pot%set(idx, ireff(CR_ELEC) + log(max(dop(CR_ELEC)/par%smc%n_intrin, 1.0)))
      elseif (par%ci0 == CR_HOLE) then
        call pot%set(idx, ireff(CR_HOLE) - log(max(dop(CR_HOLE)/par%smc%n_intrin, 1.0)))
      end if
    end do
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
    real                 :: eps, cap
    real,    allocatable :: d_volt(:,:), eye(:,:)
    integer, allocatable :: idx1(:), idx2(:), idx(:)
    logical              :: status

    ! init base
    call this%equation_init("imref_approx_eq")
    this%par => par
    eye = eye_real(size(par%contacts))

    ! create variable selectors
    call this%iref%init(iref, [par%transport_vct( 0)%get_ptr(), &                              ! uncontacted transport vertices
      &                      (par%transport_vct(ict)%get_ptr(), ict = 1, size(par%contacts))]) ! contacted vertices
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, size(volt))], "voltages")

    ! init residuals using this%iref as main variable
    call this%init_f(this%iref)

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_dir_volt%init(par%g, g2 = get_dummy_grid(), perm = dum)
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()

    ! init jacos
    this%jaco_iref  => this%init_jaco_f(this%depend(this%iref ), st = [this%st_nn%get_ptr(), &
    & (this%st_dir%get_ptr(),      ict = 1, size(par%contacts))], const = .true.)
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), st = [this%st_em%get_ptr(),   &
    & (this%st_dir_volt%get_ptr(), ict = 1, size(par%contacts))], const = .true.)

    ! loop over transport cells
    allocate (idx1(par%g%idx_dim), idx2(par%g%idx_dim))
    do idx_dir = 1, par%g%idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get "permittivity" (= mob * dens)
        eps = par%mob0(IDX_EDGE,idx_dir,iref%ci)%get(idx) * max(par%smc%n_intrin, par%dop(IDX_EDGE,idx_dir,iref%ci)%get(idx))

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
    do ict = 1, size(par%contacts)
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
