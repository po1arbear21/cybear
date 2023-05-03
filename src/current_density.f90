m4_include(util/macro.f90.inc)

module current_density_m

  use density_m,       only: density
  use device_params_m, only: device_params, CR_NAME, CR_CHARGE, DIR_NAME
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use jacobian_m,      only: jacobian
  use math_m,          only: ber, dberdx
  use mobility_m,      only: mobility
  use potential_m,     only: potential
  use stencil_m,       only: dirichlet_stencil, near_neighb_stencil
  use variable_m,      only: variable_real
  use error_m,         only: assert_failed, program_error
  use radau5_m,        only: ode_options, ode_result, radau5
  use newton_m,        only: newton1D, newton1D_opt
  use distributions_m, only: fermi_dirac_integral_1h, fermi_dirac_integral_m1h, inv_fermi_dirac_integral_1h
  use dual_m,          only: dual_1

  implicit none

  private
  public current_density, calc_current_density

  type, extends(variable_real) :: current_density
    !! electron/hole current density
    integer                    :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
    real, pointer :: x2(:,:)   => null()
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_current_density
    !! calculate current density by drift-diffusion model

    type(device_params),   pointer :: par   => null()
    type(potential),       pointer :: pot   => null()
    type(density),         pointer :: dens  => null()
    type(current_density), pointer :: cdens => null()
    type(mobility),        pointer :: mob   => null()

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_dens => null()
    type(jacobian), pointer :: jaco_mob  => null()
  contains
    procedure :: init => calc_current_density_init
    procedure :: eval => calc_current_density_eval

    procedure, private :: eval_sg    => calc_current_density_eval_sg
    procedure, private :: eval_degen => calc_current_density_eval_degen
  end type

contains

  subroutine current_density_init(this, par, ci, idx_dir)
    !! initialize current density
    class(current_density), intent(out) :: this
    type(device_params),    intent(in)  :: par
      !! device parameters
    integer,                intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    integer,                intent(in)  :: idx_dir
      !! edge direction

    integer                        :: idx_dim
    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
    this%ci = ci

    idx_dim = par%g%idx_dim
    select case (idx_dim)
    case (1)
      p1 => this%data%get_ptr1()
      this%x1 => p1%data
    case (2)
      p2 => this%data%get_ptr2()
      this%x2 => p2%data
    case (3)
      p3 => this%data%get_ptr3()
      this%x3 => p3%data
    case default
      call program_error("Maximal 3 dimensions allowed")
    end select
  end subroutine

  subroutine calc_current_density_init(this, par, pot, dens, cdens, mob)
    !! initialize drift-diffusion equation
    class(calc_current_density),   intent(out) :: this
    type(device_params),   target, intent(in)  :: par
      !! device parameters
    type(potential),       target, intent(in)  :: pot
      !! potential variable
    type(density),         target, intent(in)  :: dens
      !! density variable
    type(current_density), target, intent(in)  :: cdens
      !! current density variable
    type(mobility),        target, intent(in)  :: mob
      !! mobility variable

    integer :: idx_dir, ci, iprov

    idx_dir = cdens%idx_dir
    ci      = cdens%ci

    ! init base
    call this%equation_init("calc_"//cdens%name)
    this%par   => par
    this%pot   => pot
    this%dens  => dens
    this%cdens => cdens
    this%mob   => mob

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_nn%init(par%g, IDX_EDGE, idx_dir, IDX_VERTEX, 0)

    ! provide current density
    iprov = this%provide(cdens, par%transport(IDX_EDGE, idx_dir))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on density
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on mobility
    this%jaco_mob => this%init_jaco(iprov, this%depend(mob, par%transport(IDX_EDGE, idx_dir)), st = [this%st_dir%get_ptr()])

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    !! evaluate drift-diffusion equation
    class(calc_current_density), intent(inout) :: this

    integer               :: i, idx_dir, idx_dim
    real                  :: pot(2), dens(2), j, djdpot(2), djddens(2), djdmob, len, mob
    integer, allocatable  :: idx(:), idx1(:), idx2(:)
    logical               :: status

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    ! ci      = this%cdens%ci
    ! ch      = CR_CHARGE(ci)

    allocate (idx1(idx_dim), idx2(idx_dim))

    ! loop over transport edges
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(i,pot,dens,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    !$omp shared(this,idx_dir,idx_dim)
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

      ! parameters
      len   = this%par%g%get_len(idx1, idx_dir)
      pot(1) = this%pot%get(idx1)
      pot(2) = this%pot%get(idx2)
      dens(1) = this%dens%get(idx1)
      dens(2) = this%dens%get(idx2)
      mob   = this%mob%get(idx1)

      if (this%par%smc%degen) then
        ! generalized Scharfetter-Gummel
        call this%eval_degen(len, pot, dens, mob, j, djdpot, djddens, djdmob)
      else
        ! Scharfetter-Gummel
        call this%eval_sg(len, pot, dens, mob, j, djdpot, djddens, djdmob)
      end if

      ! set current density + derivatives
      call this%cdens%set(idx, j)
      call this%jaco_pot%set( idx, idx1, djdpot(1))
      call this%jaco_pot%set( idx, idx2, djdpot(2))
      call this%jaco_dens%set(idx, idx1, djddens(1))
      call this%jaco_dens%set(idx, idx2, djddens(2))
      call this%jaco_mob%set( idx, idx,  djdmob)
    end do
     !$omp end parallel do
  end subroutine

  subroutine calc_current_density_eval_sg(this, len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! Scharfetter-Gummel stabilization
    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    real :: ch, ber1, ber2, dber1, dber2

    ch = CR_CHARGE(this%cdens%ci)

    ber1  = ber(ch * (pot(1) - pot(2)))
    ber2  = ber(ch * (pot(2) - pot(1)))
    dber1 = ch * dberdx(ch * (pot(1) - pot(2)))
    dber2 = ch * dberdx(ch * (pot(2) - pot(1)))

    j = - mob * (ber1 * dens(2) - ber2 * dens(1)) / len

    djdpot(1) = - mob * (dber1 * dens(2) + dber2 * dens(1)) / len
    djdpot(2) =   mob * (dber1 * dens(2) + dber2 * dens(1)) / len

    djddens(1)  =   mob * ber2 / len
    djddens(2)  = - mob * ber1 / len

    djdmob = -(ber1 * dens(2) - ber2 * dens(1)) / len
  end subroutine

  subroutine calc_current_density_eval_degen(this, len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! Generalized Scharfetter-Gummel stabilization for degenerate case
    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    integer            :: ci
    real               :: ch, edos, jsg, djsgdpot(2), djsgddens(2), djsgdmob
    real               :: jj0, jj, djjdp(3), nn(2)
    type(newton1D_opt) :: newt_opt
    type(ode_options)  :: ode_opt
    type(ode_result)   :: ode_res1, ode_res2

    ci   = this%cdens%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)

    ! initial guess
    call this%eval_sg(len, pot, dens, mob, jsg, djsgdpot, djsgddens, djsgdmob)
    jj0 = jsg * len / (mob * edos)
    nn  = dens / edos

    ! solve with newton iteration
    call ode_opt%init(1, atol = [minval(nn*1e-14)], rtol = [1e-10], max_rejected = 20)
    call newt_opt%init()
    call newton1D(newton_fun, [pot(2) - pot(1), nn(1), nn(2)], newt_opt, jj0, jj, djjdp)

    ! extract solution + derivatives
    j       = jj * mob * edos / len
    djdpot  = [-1.0, 1.0] * djjdp(1) * mob * edos / len
    djddens = djjdp(2:3) * mob / len
    djdmob  = jj * edos / len

  contains

    subroutine newton_fun(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (jj)
      real,              intent(in)  :: p(:)
        !! parameters (pot(2) - pot(1), nn(1), nn(2))
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: nnr(2), djj(2), dpot(2), dnn0(2)

      if (p(2) < p(3) * 1e-3) then
        ! solve ode from left to right
        call radau5(ode_fun, 0.0, 1.0, [1.0], [p(2)], [p(1), x], ode_opt, ode_res1)
        nnr( 1) = ode_res1%Usmp(    1,  1)
        djj( 1) = ode_res1%dUsmpdP( 1,2,1)
        dpot(1) = ode_res1%dUsmpdP( 1,1,1)
        dnn0(1) = ode_res1%dUsmpdU0(1,1,1)
        nnr( 2) = p(3)
        djj( 2) = 0
        dpot(2) = 0
        dnn0(2) = 1
      elseif (p(3) < p(2) * 1e-3) then
        ! solve ode from right to left
        call radau5(ode_fun, 1.0, 0.0, [0.0], [p(3)], [p(1), x], ode_opt, ode_res2)
        nnr( 1) = p(2)
        djj( 1) = 0
        dpot(1) = 0
        dnn0(1) = 1
        nnr( 2) = ode_res2%Usmp(    1,  1)
        djj( 2) = ode_res2%dUsmpdP( 1,2,1)
        dpot(2) = ode_res2%dUsmpdP( 1,1,1)
        dnn0(2) = ode_res2%dUsmpdU0(1,1,1)
      else
        ! solve ode from left to center and from right to center
        call radau5(ode_fun, 0.0, 0.5, [0.5], [p(2)], [p(1), x], ode_opt, ode_res1)
        call radau5(ode_fun, 1.0, 0.5, [0.5], [p(3)], [p(1), x], ode_opt, ode_res2)
        nnr( 1) = ode_res1%Usmp(    1,  1)
        djj( 1) = ode_res1%dUsmpdP( 1,2,1)
        dpot(1) = ode_res1%dUsmpdP( 1,1,1)
        dnn0(1) = ode_res1%dUsmpdU0(1,1,1)
        nnr( 2) = ode_res2%Usmp(    1,  1)
        djj( 2) = ode_res2%dUsmpdP( 1,2,1)
        dpot(2) = ode_res2%dUsmpdP( 1,1,1)
        dnn0(2) = ode_res2%dUsmpdU0(1,1,1)
      end if

      ! residual
      f = nnr(1) - nnr(2)
      if (present(dfdx)) dfdx = djj(1) - djj(2)
      if (present(dfdp)) dfdp = [dpot(1) - dpot(2), dnn0(1), -dnn0(2)]
    end subroutine

    subroutine ode_fun(x, U, P, f, dfdU, dfdP)
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (nn)
      real,           intent(in)  :: P(:)
        !! parameters (pot(2) - pot(1), jj)
      real, optional, intent(out) :: f(:)
        !! output dnn/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt nn
      real, optional, intent(out) :: dfdP(:,:)
        !! output derivative of f wrt P

      real :: eta, deta, Fm1h, dFm12, alpha, dalpha, fsg, dfsg

      m4_ignore(x)

      ! get "degeneracy-factor" (~1 for small densities; < 1 for large densities)
      call inv_fermi_dirac_integral_1h(U(1), eta, deta)
      call fermi_dirac_integral_m1h(eta, Fm1h, dFm12)
      alpha = Fm1h / U(1)
      dalpha = (dFm12 * deta - alpha) / U(1)

      ! Scharfetter-Gummel
      fsg  = - ch * P(1) * U(1) - P(2)
      dfsg = - ch * P(1)

      ! scale Scharfetter-Gummel by degeneracy factor
      if (present(f)) then
        f(1) = alpha * fsg
      end if
      if (present(dfdU)) then
        dfdU(1,1) = dalpha * fsg + alpha * dfsg
      end if
      if (present(dfdp)) then
        dfdp(1,1) = - alpha * ch * U(1)
        dfdp(1,2) = - alpha
      end if
    end subroutine

  end subroutine

end module
