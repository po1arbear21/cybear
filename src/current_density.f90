m4_include(util/macro.f90.inc)

module current_density_m

  use density_m,        only: density
  use device_params_m,  only: device_params
  use distributions_m,  only: fermi_dirac_integral_1h, fermi_dirac_integral_m1h, inv_fermi_dirac_integral_1h
  use dual_m
  use equation_m,       only: equation
  use error_m,          only: assert_failed, program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m, only: DIR_NAME
  use ieee_arithmetic,  only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_value
  use jacobian_m,       only: jacobian
  use math_m,           only: ber, dberdx, log1p
  use mobility_m,       only: mobility
  use newton_m,         only: newton1D, newton1D_opt
  use potential_m,      only: potential
  use radau5_m,         only: ode_options, ode_result, radau5
  use semiconductor_m,  only: CR_NAME, CR_CHARGE
  use stencil_m,        only: dirichlet_stencil, near_neighb_stencil
  use variable_m,       only: variable_real
use ode_m
  implicit none

  private
  public current_density, calc_current_density

  type, extends(variable_real) :: current_density
    !! electron/hole current density

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
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

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
    this%ci = ci

    ! get pointer to data
    select case (par%g%idx_dim)
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

! block
!   len  = 2.0177215314875774E-01
!   pot  = [3.0775684460780085E+03, 2.7697214714903826E+03]
!   dens = [3.1291100007260222E+00, 1.0000000000000001E-30]
!   mob  = 3.5473379003530772E-02
!   call this%eval_degen(len, pot, dens, mob, j, djdpot, djddens, djdmob)
!   print *, j, djdpot, djddens, djdmob
!   stop
! end block

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    allocate (idx1(idx_dim), idx2(idx_dim))

    ! loop over transport edges in parallel
    ! $omp parallel do default(none) schedule(dynamic) &
    ! $omp private(i,pot,dens,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    ! $omp shared(this,idx_dir,idx_dim)
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
        ! generalized Scharfetter-Gummel (slow)
        call this%eval_degen(len, pot, dens, mob, j, djdpot, djddens, djdmob)
      else
        ! Scharfetter-Gummel (fast)
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
    ! $omp end parallel do
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

! integer, save :: irun = 0

    ! 1/FD1H(x) = exp(-x) + sqrt(1/8) for x < -16
    real, parameter :: alpha = sqrt(0.125)
    real, parameter :: etaF  = -16

    integer            :: ci, dir
    logical            :: l2r
    real               :: ch, djdp(3), dpot, edos, j0, jmin, jmax, jsgn, n(2), eta(2), deta, detadn(2), t, tol, xsmp
    type(newton1D_opt) :: newt_opt
    type(ode_options)  :: ode_opt
    type(ode_result)   :: ode_res1, ode_res2

    ! abbreviations
    ci   = this%cdens%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)

    ! normalize
    dpot = - ch * (pot(2) - pot(1))
    n    = dens / edos

    ! eta
    call inv_fermi_dirac_integral_1h(n(1), eta(1), detadn(1))
    call inv_fermi_dirac_integral_1h(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

block
  use math_m
  integer, parameter :: NN = 201
  integer            :: ii
  real, allocatable  :: ee(:)

  allocate (ee(NN))
  ee = linspace(-100.0, 5.0, NN)

  eta(1) = 0
  dpot   = -100
  call fermi_dirac_integral_1h(eta(1), n(1), t)
  do ii = 1, NN
    eta(2) = ee(ii)
    deta = eta(2) - eta(1)
    call fermi_dirac_integral_1h(eta(2), n(2), t)



! irun = irun + 1
! print "(A,I6)", "irun = ", irun

! if (irun == 720) then
!   print "(A,ES24.16)", "eta1 = ", eta(1)
!   print "(A,ES24.16)", "eta2 = ", eta(2)
!   print "(A,ES24.16)", "dpot = ", dpot
! end if

    if (all(eta <= etaF)) then
      ! special case: 1/F12(eta) ~ exp(-eta) + sqrt(1/8) -> Bernoulli iteration
      call this%eval_sg(len, pot, dens, mob, j0, djdpot, djddens, djdmob)
      j0 = len * j0 / (mob * edos)
      call newt_opt%init(atol = 1e-16 * abs(j0), rtol = 1e-14, log = .true.)
      call newton1D(small_eta_residual, [dpot, eta(1), eta(2)], newt_opt, j0, j, djdp)
      djdp(2:3) = djdp(2:3) * detadn
    else
      ! jmin, jmax by slope (detadx must be equal to deta for some x in [0, 1])
      jmin = min(n(1), n(2)) * abs(dpot - deta)
      jmax = max(n(1), n(2)) * abs(dpot - deta)
      jsgn = sign(1.0, dpot - deta)
      if (jsgn < 0) then
        t    = jmin
        jmin = - jmax
        jmax = - t
      end if

      ! reduce j range further if possible
      if (dpot < 0) then
        if (eta(2) > eta(1)) then
          dir  = 1
          jmax = min(jmax, dpot * n(2))
        elseif (eta(2) >= eta(1) + dpot) then
          dir  = 0
          jmin = max(jmin, dpot * n(2))
          jmax = min(jmax, 0.0)
        else
          dir  = -1
          jmin = max(jmin, 0.0)
        end if
      else
        if (eta(2) < eta(1)) then
          dir  = -1
          jmin = max(jmin, dpot * n(1))
        elseif (eta(2) <= eta(1) + dpot) then
          dir  = 0
          jmin = max(jmin, 0.0)
          jmax = min(jmax, dpot * n(1))
        else
          dir  = 1
          jmax = min(jmax, 0.0)
        end if
      end if

      ! print "(3ES24.16)", eta(2), jmin, jmax

      tol = 1e-10 * (max(abs(eta(1)), abs(eta(2)), abs(dpot)) + 1)
      if ((abs(deta - dpot) < tol) .or. (abs(deta) < tol)) then
        ! eta changes approximately linear
        j       = 0.5 * (n(1) + n(2)) * (dpot - deta)
        djdp(1) = 0.5 * (n(1) + n(2))
        djdp(2) = 0.5 * (dpot - deta) + 0.5 * (n(1) + n(2)) * detadn(1)
        djdp(3) = 0.5 * (dpot - deta) - 0.5 * (n(1) + n(2)) * detadn(2)
      else
        if (dir == 0) then
          if (dpot < 0) then
            j0   = jmin + abs(jmin) * 1e-14
            jmin = jmin - abs(jmin) * 1e-14
            ! jmax = jmax - abs(jmax) * 1e-14
          else
            j0   = jmax - abs(jmax) * 1e-14
            ! jmin = jmin + abs(jmin) * 1e-14
            jmax = jmax + abs(jmax) * 1e-14
          end if
        else
          call this%eval_sg(len, pot, dens, mob, j0, djdpot, djddens, djdmob)
          j0 = len * j0 / (mob * edos)
          if ((j0 < jmin) .or. (j0 > jmax)) then
            j0 = 0.5 * (jmin + jmax)
          end if
        end if
        ! if ((dir == 0) .and. (abs(deta) > ?)) then
        !   if (dpot < 0) then
        !     j = jmin
        !     ....
        !   else
        !     j = jmax
        !     ....
        !   end if
        ! else
          ! ode-solver + newton iteration
          call ode_opt%init(1, atol = [1e-16], rtol = [1e-14], max_rejected = 200)
          call newt_opt%init(atol = 1e-16 * max(abs(jmin), abs(jmax)), rtol = 1e-14, xmin = jmin, xmax = jmax, log = .false.)

! block
!   integer, parameter :: NNN = 201
!   integer            :: iii
!   real, allocatable :: jjj(:), fff(:)

!   allocate (jjj(NNN), fff(NNN))
!   jjj = linspace(jmin - abs(jmin)*1e-15, jmax + abs(jmax)*1e-15, NNN)
!   do iii = 1, NNN
!     call shooting_residual(jjj(iii), [dpot, eta(1), eta(2)], fff(iii))

!     print "(2ES24.16)", jjj(iii), fff(iii)
!   end do

! end block

          call newton1D(shooting_residual, [dpot, eta(1), eta(2)], newt_opt, j0, j, djdp)
          djdp(2:3) = djdp(2:3) * detadn


        ! end if
      end if
    end if

    print "(4ES24.16)", eta(2), j, jmin, jmax
  end do

  stop
end block
    ! extract solution + derivatives
    j       = j * mob * edos / len
    djdpot  = [-1.0, 1.0] * ch * djdp(1) * mob * edos / len
    djddens = djdp(2:3) * mob / len
    djdmob  = j / mob

! if (irun == 720) then
!   stop
! end if

  contains

    subroutine small_eta_residual(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (j)
      real,              intent(in)  :: p(:)
        !! parameters (pot(2) - pot(1), eta(1), eta(2))
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: arg, ber1, ber2, dber1, dber2

      ! shifted bernoulli argument
      arg = p(1) - alpha * x

      ber1  = ber(   -arg)
      ber2  = ber(    arg)
      dber1 = dberdx(-arg)
      dber2 = dberdx( arg)

      f = x - ber1 * exp(p(2)) + ber2 * exp(p(3))
      if (present(dfdx)) then
        dfdx = 1 + alpha * (dber1 * exp(p(2)) + dber2 * exp(p(3)))
      end if
      if (present(dfdp)) then
        dfdp(1) = - dber1 * exp(p(2)) + dber2 * exp(p(3))
        dfdp(2) = - ber1 * exp(p(2))
        dfdp(3) =   ber2 * exp(p(3))
      end if
    end subroutine

    subroutine shooting_residual(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (j)
      real,              intent(in)  :: p(:)
        !! parameters (pot(2) - pot(1), eta(1), eta(2))
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: eta1, deta1deta0, deta1ddpot, deta1dj, eta2, deta2deta0, deta2ddpot, deta2dj, xsmp

      if (dir > 0) then ! left to right
        l2r = .true.
        call solve_ode(0.0, 1.0, p(2), p(1), x, eta1, deta1deta0, deta1ddpot, deta1dj)
        f = eta1 - p(3)
        if (present(dfdx)) then
          dfdx = deta1dj
        end if
        if (present(dfdp)) then
          dfdp(1) = deta1ddpot
          dfdp(2) = deta1deta0
          dfdp(3) = -1.0
        end if
      elseif (dir < 0) then ! right to left
        l2r = .false.
        call solve_ode(1.0, 0.0, p(3), p(1), x, eta1, deta1deta0, deta1ddpot, deta1dj)
        f = eta1 - p(2)
        if (present(dfdx)) then
          dfdx = deta1dj
        end if
        if (present(dfdp)) then
          dfdp(1) = deta1ddpot
          dfdp(2) = -1.0
          dfdp(3) = deta1deta0
        end if
      else ! left to center, right to center
        if (dpot < 0) then
          xsmp = (eta(2) - eta(1)) / dpot
        elseif (dpot > 0) then
          xsmp = 1 - (eta(2) - eta(1)) / dpot
        else
          xsmp = 0.5
        end if
        if ((xsmp < 0) .or. (xsmp > 1)) then
          print *, xsmp
          error stop "xsmp not in [0, 1]"
        end if

        l2r = .true.
        call solve_ode(0.0, xsmp, p(2), p(1), x, eta1, deta1deta0, deta1ddpot, deta1dj)
        l2r = .false.
        call solve_ode(1.0, xsmp, p(3), p(1), x, eta2, deta2deta0, deta2ddpot, deta2dj)
        f = eta1 - eta2
        if (present(dfdx)) then
          dfdx = deta1dj - deta2dj
        end if
        if (present(dfdp)) then
          dfdp(1) = deta1ddpot - deta2ddpot
          dfdp(2) = deta1deta0
          dfdp(3) =            - deta2deta0
        end if
      end if
    end subroutine

    recursive subroutine solve_ode(x0, x1, eta0, dpot, j, eta1, deta1deta0, deta1ddpot, deta1dj)
      real, intent(in)  :: x0
      real, intent(in)  :: x1
      real, intent(in)  :: eta0
      real, intent(in)  :: dpot
      real, intent(in)  :: j
      real, intent(out) :: eta1
      real, intent(out) :: deta1deta0
      real, intent(out) :: deta1ddpot
      real, intent(out) :: deta1dj

      real             :: deta1dx(1)
      type(dual_3)     :: dpot_, j_, eta0_, eta1_, t, xF, B, t1, t2
      type(ode_result) :: result

! if (irun == 720) then
!   print *
!   print "(A,ES24.16)", "x0   = ", x0
!   print "(A,ES24.16)", "x1   = ", x1
!   print "(A,ES24.16)", "eta0 = ", eta0
!   print "(A,ES24.16)", "dpot = ", dpot
!   print "(A,ES24.16)", "j    = ", j
! end if

      if (eta0 < etaF) then
        ! use dual numbers to avoid most of the manual derivative calculations
        call eta0_%init(eta0, i = 1)
        call dpot_%init(dpot, i = 2)
        call j_%init(      j, i = 3)
        t    = dpot_ - alpha * j_

        ! find x where eta crosses etaF
        t1 = t * exp(etaF ) - j_
        t2 = t * exp(eta0_) - j_
        if (((t1%x > 0) .and. (t2%x > 0)) .or. ((t1%x < 0) .and. (t2%x < 0))) then
          if (abs(t%x) < 1e-6) then
            xF = x0 + (exp(eta0_) - exp(etaF)) / j_ + t * (exp(2 * eta0_) - exp(2 * etaF)) / (2 * j_**2) + t**2 * (exp(3 * eta0_) - exp(3 * etaF)) / (3 * j_**3)
          else
            xF = x0 + log(t1 / t2) / t
          end if

          if (((x1 > x0) .and. (xF%x >= x0) .and. (xF%x < x1)) .or. ((x1 < x0) .and. (xF%x <= x0) .and. (xF%x > x1))) then
            ! go from (xF,etaF) to (x1,eta1) using ode solver
            call solve_ode(xF%x, x1, etaF, dpot, j, eta1, deta1deta0, deta1ddpot, deta1dj)
            call detadx(x1, [eta1], [dpot, j], f = deta1dx)
            deta1deta0 =            - deta1dx(1) * xF%dx(1)
            deta1ddpot = deta1ddpot - deta1dx(1) * xF%dx(2)
            deta1dj    = deta1dj    - deta1dx(1) * xF%dx(3)
          else
            ! stay below -16
            B%x  = ber(   t%x * (x1 - x0))
            B%dx = dberdx(t%x * (x1 - x0)) * t%dx
            eta1_ = eta0_ + log1p((t - j_ * exp(-eta0_)) * (x1 - x0) / B)

            eta1       = eta1_%x
            deta1deta0 = eta1_%dx(1)
            deta1ddpot = eta1_%dx(2)
            deta1dj    = eta1_%dx(3)
          end if
        else
          ! stay constant (avoid numerical issues)
          eta1       = eta0
          deta1deta0 = 1
          deta1ddpot = 0
          deta1dj    = 0
        end if
      else
        ! go from (x0,eta0) to (x1,eta1) using ode solver
!FIXME        ! call radau5(detadx, x0, x1, [x1], [eta0], [dpot, j], ode_opt, result)
        eta1       = result%Usmp(    1,  1)
        deta1deta0 = result%dUsmpdU0(1,1,1)
        deta1ddpot = result%dUsmpdP( 1,1,1)
        deta1dj    = result%dUsmpdP( 1,2,1)
      end if

! if (irun == 720) then
!   print "(A,ES24.16)", "eta1 = ", eta1
! end if
    end subroutine

    subroutine detadx(x, U, p, f, dfdU, dfdp)
      !! ode right-hand side
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (eta)
      real,           intent(in)  :: p(:)
        !! parameters (dpot, j)
      real, optional, intent(out) :: f(:)
        !! output deta/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt eta
      real, optional, intent(out) :: dfdp(:,:)
        !! output derivative of f wrt P

      real :: F12, dF12

      ! if (dir == 0) then
      !   if ((dpot < 0) .and. (.not. l2r)) then
      !     if ((U(1) < eta(2)) .or. (U(1) > eta(2) - dpot * (1.0 - x))) then
      !       if (present(f)) f(1) = (eta(2) - U(1)) / (1.0 - x)
      !       if (present(dfdU)) dfdU(1,1) = - 1.0 / (1.0 - x)
      !       if (present(dfdp)) dfdp = 0 ! FIXME: put eta(2) in p ?
      !       return
      !     end if
      !   elseif ((dpot > 0) .and. (l2r)) then
      !     if ((U(1) < eta(1)) .or. (U(1) > eta(1) + dpot * x)) then
      !       if (present(f)) f(1) = (U(1) - eta(1)) / x
      !       if (present(dfdU)) dfdU(1,1) = 1.0 / x
      !       if (present(dfdp)) dfdp = 0 ! FIXME: put eta(1) in p ?
      !       return
      !     end if
      !   end if
      ! end if

      call fermi_dirac_integral_1h(U(1), F12, dF12)

      if (present(f)) then
        f(1) = p(1) - p(2) / F12

        if (.not. ieee_is_finite(f(1))) then
          print "(A,ES24.16)", "x    = ", x
          print "(A,ES24.16)", "eta  = ", U(1)
          print "(A,ES24.16)", "dpot = ", p(1)
          print "(A,ES24.16)", "j    = ", p(2)
          print "(A,ES24.16)", "F12  = ", F12
          print "(A,ES24.16)", "dF12 = ", dF12
          print "(A,ES24.16)", "f    = ", f(1)
          error stop "f not finite"
        end if
      end if
      if (present(dfdU)) then
        dfdU(1,1) = p(2) / F12**2 * dF12

        if (.not. ieee_is_finite(dfdU(1,1))) then
          print "(A,ES24.16)", "x    = ", x
          print "(A,ES24.16)", "eta  = ", U(1)
          print "(A,ES24.16)", "dpot = ", p(1)
          print "(A,ES24.16)", "j    = ", p(2)
          print "(A,ES24.16)", "F12  = ", F12
          print "(A,ES24.16)", "dF12 = ", dF12
          print "(A,ES24.16)", "f    = ", f(1)
          print "(A,ES24.16)", "dfdU = ", dfdU(1,1)
          error stop "dfdU not finite"
        end if
      end if
      if (present(dfdp)) then
        dfdp(1,1) = 1.0
        dfdp(1,2) = - 1.0 / F12
        if (any(.not. ieee_is_finite(dfdp))) then
          print "(A,ES24.16)", "x    = ", x
          print "(A,ES24.16)", "eta  = ", U(1)
          print "(A,ES24.16)", "dpot = ", p(1)
          print "(A,ES24.16)", "j    = ", p(2)
          print "(A,ES24.16)", "F12  = ", F12
          print "(A,ES24.16)", "dF12 = ", dF12
          print "(A,ES24.16)", "f    = ", f(1)
          print "(A,2ES24.16)", "dfdp = ", dfdp(1,1), dfdp(1,2)
          error stop "dfdp not finite"
        end if
      end if
    end subroutine

  end subroutine

end module
