m4_include(util/macro.f90.inc)

module current_density_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_value

  use current_integral_m, only: current_integral_get, CURRENT_INTEGRAL_DEBUG
  use density_m,          only: density
  use device_params_m,    only: device_params
  use equation_m,         only: equation
  use error_m,            only: assert_failed, program_error
  use grid_m,             only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,        only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m,   only: DIR_NAME
  use jacobian_m,         only: jacobian
  use math_m,             only: ber, dberdx, dberdx2, log1p
  use mobility_m,         only: mobility
  use potential_m,        only: potential
  use radau5_m,           only: ode_options, ode_result, radau5
  use semiconductor_m,    only: CR_NAME, CR_CHARGE, DOS_PARABOLIC, DIST_MAXWELL, &
    &                           STAB_ED, STAB_MED, STAB_EXACT, STAB_SG, semiconductor
  use stencil_m,          only: dirichlet_stencil, near_neighb_stencil
  use variable_m,         only: variable

  implicit none

  private
  public current_density, calc_current_density, time

  type, extends(variable) :: current_density
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

    logical :: stat

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

    procedure, private :: get_curr     => calc_current_density_get_curr
    procedure, private :: get_curr_sg  => calc_current_density_get_curr_sg
    procedure, private :: get_curr_ed  => calc_current_density_get_curr_ed
    procedure, private :: get_curr_ed2 => calc_current_density_get_curr_ed2
  end type

  real :: time = 0.0

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
    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", &
      &                     g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
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

    print "(A)", "calc_current_density_init"

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
    if (par%smc%mob) this%jaco_mob => this%init_jaco(iprov, this%depend(mob, par%transport(IDX_EDGE, idx_dir)), st = [this%st_dir%get_ptr()])

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    !! evaluate drift-diffusion equation
    class(calc_current_density), intent(inout) :: this

    integer               :: i, idx_dir, idx_dim, start, end, rate
    real                  :: pot(2), dens(2), j, djdpot(2), djddens(2), djdmob, len, mob
    integer, allocatable  :: idx(:), idx1(:), idx2(:)
    logical               :: status

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim))
    call system_clock(start, rate)

    ! loop over transport edges in parallel
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(i,pot,dens,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    !$omp shared(this,idx_dir,idx_dim)
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

      ! parameters
      len     = this%par%g%get_len(idx, idx_dir)
      pot( 1) = this%pot%get( idx1)
      pot( 2) = this%pot%get( idx2)
      dens(1) = this%dens%get(idx1)
      dens(2) = this%dens%get(idx2)
      if (this%par%smc%mob) then
        mob = this%mob%get(idx)
      else
        mob = this%par%mob0(IDX_EDGE, idx_dir, this%cdens%ci)%get(idx)
      end if

      ! get current along edge
      call this%get_curr(this%par%smc, len, pot, dens, mob, j, djdpot, djddens, djdmob)

      ! set current density + derivatives
      call this%cdens%set(idx, j)
      call this%jaco_pot%add( idx, idx1, djdpot(1))
      call this%jaco_pot%add( idx, idx2, djdpot(2))
      call this%jaco_dens%add(idx, idx1, djddens(1))
      call this%jaco_dens%add(idx, idx2, djddens(2))
      if (this%par%smc%mob) call this%jaco_mob%add(idx, idx, djdmob)
    end do
    !$omp end parallel do

    call system_clock(end)
    time = time + real(end-start)/real(rate)
  end subroutine

  subroutine calc_current_density_get_curr(this, smc, len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! get current along edge
    class(calc_current_density), intent(in)  :: this
    type(semiconductor),         intent(in)  :: smc
      !! semiconductor material
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

    integer :: ci
    real    :: ch, dpot, edos, n(2), djdn(2), djddpot

    ! abbreviations
    ci   = this%cdens%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)

    ! normalize potential drop and density
    dpot = - ch * (pot(2) - pot(1))
    n    = dens / edos

    if (smc%dist == DIST_MAXWELL) then
      ! non-degenerate case: always use Scharfetter-Gummel
      call this%get_curr_sg(n, dpot, j, djdn, djddpot)
    else
      ! get current based on different stabilization method
      select case (smc%stab)
      case (STAB_SG)
        call this%get_curr_sg(n, dpot, j, djdn, djddpot)
      case (STAB_ED)
        call this%get_curr_ed(smc, n, dpot, j, djdn, djddpot)
      case (STAB_MED)
        call this%get_curr_ed2(smc, n, dpot, j, djdn, djddpot)
      case (STAB_EXACT)
        call current_integral_get(dist, inv_dist, int_dist, n, dpot, j, djdn, djddpot)
      end select

      if (.not. ieee_is_finite(j)) then
        print "(A,ES25.16E3)", "n1   = ", n(1)
        print "(A,ES25.16E3)", "n2   = ", n(2)
        print "(A,ES25.16E3)", "dpot = ", dpot
        error stop "j is not finite"
      end if
    end if

    ! denormalize
    j       = j * mob * edos / len
    djdpot  = ch * [1.0, -1.0] * djddpot * mob * edos / len
    djddens = djdn * mob / len
    djdmob  = j / mob
  contains

    subroutine dist(eta, k, F, dFdeta)
      real,    intent(in)  :: eta
      integer, intent(in)  :: k
      real,    intent(out) :: F
      real,    intent(out) :: dFdeta

      call smc%get_dist(eta, k, F, dFdeta)
    end subroutine

    subroutine inv_dist(F, eta, detadF)
      real, intent(in)  :: F
      real, intent(out) :: eta
      real, intent(out) :: detadF

      call smc%get_inv_dist(F, eta, detadF)
    end subroutine

    subroutine int_dist(eta, k, I, dIdeta)
      real,    intent(in)  :: eta(2)
      integer, intent(in)  :: k
      real,    intent(out) :: I
      real,    intent(out) :: dIdeta(2)

      call smc%get_int_dist(eta, k, I, dIdeta)
    end subroutine

  end subroutine

  subroutine calc_current_density_get_curr_sg(this, n, dpot, j, djdn, djddpot)
    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real :: B2, B1, dB2, dB1

    m4_ignore(this)

    B1  = ber(   -dpot)
    B2  = ber(    dpot)
    dB1 = dberdx(-dpot)
    dB2 = dberdx( dpot)

    j       = B1 * n(1) - B2 * n(2)
    djdn    = [B1, -B2]
    djddpot = - dB1 * n(1) - dB2 * n(2)
  end subroutine

  subroutine calc_current_density_get_curr_ed(this, smc, n, dpot, j, djdn, djddpot)
    class(calc_current_density), intent(in)  :: this
    type(semiconductor),         intent(in)  :: smc
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real :: eta(2), detadn(2), deta, h, dhdn(2), B1, B2, dB1, dB2, djdh
    real :: etam, F0, dF0, F1, dF1, F2, dF2, F3, dF3, h0, dh0dn(2), h2, dh2dn(2)

    m4_ignore(this)

    call smc%get_inv_dist(n(1), eta(1), detadn(1))
    call smc%get_inv_dist(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

    if (abs(deta) > 1e-5) then
      h    = log(n(2) / n(1)) / deta
      dhdn = [1.0, -1.0] * (h * detadn - 1 / n) / deta
    else
      etam = 0.5 * (eta(1) + eta(2))
      call smc%get_dist(etam, 0, F0, dF0)
      call smc%get_dist(etam, 1, F1, dF1)
      call smc%get_dist(etam, 2, F2, dF2)
      call smc%get_dist(etam, 3, F3, dF3)

      ! scale F1, F2, F3 by 1/F0 to simplify further calculations and avoid over-/underflow
      F1  = F1 / F0
      dF1 = (dF1 - F1 * dF0) / F0
      F2  = F2 / F0
      dF2 = (dF2 - F2 * dF0) / F0
      F3  = F3 / F0
      dF3 = (dF3 - F3 * dF0) / F0

      ! taylor series up to second order
      h0    = F1
      dh0dn = dF1 * 0.5 * detadn
      h2    = 0.25  * ((0.5 * F3 + F1**3) / 3 - 0.5 * F1 * F2)
      dh2dn = 0.125 * ((0.5 * dF3 + 3*F1**2*dF1) / 3 - 0.5 * (dF1 * F2 + F1 * dF2)) * detadn
      h     = h0 + deta**2 * h2
      dhdn  = dh0dn + 2*deta*[-1.0,1.0]*detadn * h2 + deta**2 * dh2dn
    end if

    B1  = ber(- dpot * h)
    B2  = ber(  dpot * h)
    dB1 = dberdx(-dpot * h)
    dB2 = dberdx( dpot * h)

    j       = (B1 * n(1) - B2 * n(2)) / h
    djdh    = - (j + dpot * (dB1 * n(1) + dB2 * n(2))) / h
    djdn    = [B1, -B2] / h + djdh * dhdn
    djddpot = - (dB1 * n(1) + dB2 * n(2))
  end subroutine

  recursive subroutine calc_current_density_get_curr_ed2(this, smc, n, dpot, j, djdn, djddpot)
    class(calc_current_density), intent(in)  :: this
    type(semiconductor),         intent(in)  :: smc
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real, parameter :: DETA_TOL = 1e-5, H_TOL = 1e-6

    real :: eta(2), detadn(2), dndeta(2), d2ndeta2(2), d2etadn2(2), deta, etam
    real :: I1, I2, J1, dI1deta(2), dI2deta(2), dJ1deta(2), dI1dn(2), dI2dn(2), dJ1dn(2)
    real :: F0, F1, F2, F3, F4, dF0, dF1, dF2, dF3, dF4
    real :: h00, h02, d01, d03, dh00dn(2), dh02dn(2), dd01dn(2), dd03dn(2)
    real :: h10, h12, d11, d13, dh10dn(2), dh12dn(2), dd11dn(2), dd13dn(2)
    real :: h0, d0, h1, d1, dh0dn(2), dd0dn(2), dh1dn(2), dd1dn(2)
    real :: t1, dt1dn(2), B1, B2, dB1, dB2, dB1dn(2), dB2dn(2)
    real :: c1, c2, dc1dn(2), dc2dn(2), d, dddn(2), p, dpdn(2), dpddpot, q, dqdn(2), dqddpot
    real :: hL, dhLdn(2), hR, dhRdn(2)
    real :: h, dhdn(2), dhddpot, dhdhL, dhdh0, dhdhR, dhdh1, dhdd0, dhdd1
    real :: djdh

    if (n(1) > n(2)) then
      call this%get_curr_ed2(smc, [n(2), n(1)], -dpot, j, djdn, djddpot)
      j    = -j
      djdn = - [djdn(2), djdn(1)]
      return
    end if

    call smc%get_inv_dist(n(1), eta(1), detadn(1))
    call smc%get_inv_dist(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

    if (abs(deta) <= DETA_TOL) then
      etam = 0.5 * (eta(1) + eta(2))
      call smc%get_dist(etam, 0, F0, dF0)
      call smc%get_dist(etam, 1, F1, dF1)
      call smc%get_dist(etam, 2, F2, dF2)
      call smc%get_dist(etam, 3, F3, dF3)
      call smc%get_dist(etam, 4, F4, dF4)
    end if

    if (dpot <= deta) then
      if (abs(deta) <= DETA_TOL) then
        h00    = F1 / F0
        dh00dn = (dF1 - h00 * dF0) / F0 * 0.5 * detadn
        h02    = (F0*F3 - F1*F2) / (12 * F0**2)
        dh02dn = (dF0 * F3 - F0 * dF3 - dF1 * F2 - F1 * dF2 - h02 * (24 * F0 * dF0)) / (12 * F0**2) * 0.5 * detadn
        d01    = F1 * (F1**2 - F0*F2) / (12 * F0**3)
        dd01dn = (3*F1**2*dF1 - dF0*F1*F2 - F0*dF1*F2 - F0*F1*dF2 - d01*36*F0**2*dF0) / (12*F0**3) * 0.5*detadn
        d03    = (-15*F1**3*F2 + 12*F0*F1*F2**2 + 11*F0*F1**2*F3 - 3*F0**2*F1*F4 - 5*F0**2*F2*F3) / (240*F0**4)
        dd03dn = (-45*F1**2*dF1*F2 - 15*F1**3*dF2 + 12*dF0*F1*F2**2 + 12*F0*dF1*F2**2 + 24*F0*F1*F2*dF2 &
          &    + 11*dF0*F1**2*F3 + 22*F0*F1*dF1*F3 + 11*F0*F1**2*dF3 - 6*F0*dF0*F1*F4 - 3*F0**2*dF1*F4 &
          &    - 3*F0**2*F1*dF4 - 10*F0*dF0*F2*F3 - 5*F0**2*dF2*F3 - 5*F0**2*F2*dF3 - d03*960*F0**3*dF0) &
          &    / (240*F0**4) * 0.5 * detadn

        h0    = h00 + deta**2 * h02
        dh0dn = dh00dn + 2*deta*[-1.0,1.0]*detadn * h02 + deta**2 * dh02dn

        d0    = deta * (d01 + deta**2 * d03)
        dd0dn = [-1.0,1.0]*detadn * (d01 + deta**2*d03) + deta*(dd01dn + 2*deta*[-1.0,1.0]*detadn*d03 + deta**2*dd03dn)
      else
        call smc%get_int_dist(eta, 1, I1, dI1deta)
        call smc%get_int_dist(eta, 2, I2, dI2deta)
        dI1dn = dI1deta * detadn
        dI2dn = dI2deta * detadn

        h0    = (n(2) - n(1)) / I1
        dh0dn = ([-1, 1] - h0 * dI1dn) / I1

        d0    = h0 / I1 * (I2 / I1 - 0.5 * (n(1) + n(2)))
        dd0dn = (dh0dn*I1 - h0*dI1dn)/I1**2 * (I2/I1 - 0.5*(n(1)+n(2))) + h0/I1 * ((dI2dn*I1 - I2*dI1dn) / I1**2 - 0.5)
      end if
    end if

    if (dpot > 0) then
      if (abs(deta) <= DETA_TOL) then
        h10    = F1 / F0
        dh10dn = (dF1 - h10 * dF0) / F0 * 0.5 * detadn
        h12    = (F0**2*F3 + 2*F1**3 - 3*F0*F1*F2) / (12*F0**3)
        dh12dn = (2*F0*dF0*F3 + F0**2*dF3 + 6*F1**2*dF1 - 3*(dF0*F1*F2 + F0*dF1*F2 + F0*F1*dF2) - h12*36*F0**2*dF0) &
          &    / (12*F0**3) * 0.5*detadn
        d11    = F1 * (F1**2 - F0*F2) / (12 * F0**3)
        dd11dn = (3*F1**2*dF1 - dF0*F1*F2 - F0*dF1*F2 - F0*F1*dF2 - d11*36*F0**2*dF0) / (12*F0**3) * 0.5*detadn
        d13    = (34*F1**5 - 73*F0*F1**3*F2 - 3*F0**3*F1*F4 - 5*F0**3*F2*F3 + 28*F0**2*F1*F2**2 &
          &    + 19*F0**2*F1**2*F3)/(240*F0**5)
        dd13dn = (170*F1**4*dF1 - 73*dF0*F1**3*F2 - 219*F0*F1**2*dF1*F2 - 73*F0*F1**3*dF2 - 9*F0**2*dF0*F1*F4 &
          &    - 3*F0**3*dF1*F4 - 3*F0**3*F1*dF4 - 15*F0**2*dF0*F2*F3 - 5*F0**3*dF2*F3 - 5*F0**3*F2*dF3 &
          &    + 56*F0*dF0*F1*F2**2 + 28*F0**2*dF1*F2**2 + 56*F0**2*F1*F2*dF2 + 38*F0*dF0*F1**2*F3 &
          &    + 38*F0**2*F1*dF1*F3 + 19*F0**2*F1**2*dF3 - d13*1200*F0**4*dF0) / (240*F0**5)

        h1    = h10 + deta**2 * h12
        dh1dn = dh10dn + 2*deta*[-1.0,1.0]*detadn * h12 + deta**2 * dh12dn

        d1    = deta * (d11 + deta**2 * d13)
        dd1dn = [-1.0,1.0]*detadn * (d11 + deta**2*d13) + deta*(dd11dn + 2*deta*[-1.0,1.0]*detadn*d13 + deta**2*dd13dn)
      else
        call smc%get_int_dist(eta, -1, J1, dJ1deta)
        dJ1dn = dJ1deta * detadn

        h1    = log(n(2) / n(1)) / deta
        dh1dn = [-1, 1] * (1 / n - h1 * detadn) / deta

        B1  = dberdx( -deta * h1)
        dB1 = dberdx2(-deta * h1)
        B2  = dberdx(  deta * h1)
        dB2 = dberdx2( deta * h1)

        dB1dn = dB1 * ([ 1, -1] * h1 * detadn - deta * dh1dn)
        dB2dn = dB2 * ([-1,  1] * h1 * detadn + deta * dh1dn)

        t1    = - B1 * n(1) - B2 * n(2)
        dt1dn = - dB1dn * n(1) - dB2dn * n(2) - [B1, B2]

        d1    = h1 * (1/(J1*t1) - 1/deta)
        dd1dn = dh1dn * (1/(J1*t1) - 1/deta) + h1*(-(dJ1dn * t1 + J1 * dt1dn)/(J1 * t1)**2 + ([-1, 1]*detadn/deta**2))
      end if
    end if

    if (dpot <= 0.0) then
      hL = detadn(2) / n(2)
      call smc%get_dist(eta(2), 1, dndeta(2), d2ndeta2(2))
      d2etadn2(2) = - d2ndeta2(2) / dndeta(2)**3
      dhLdn(1) = 0
      dhLdn(2) = d2etadn2(2) / n(2) - detadn(2) / n(2)**2

      h = (h0 - d0 * dpot * hL / (h0 - hL)) / (1 - d0 * dpot / (h0 - hL))
      dhddpot = d0 / (1 - d0 * dpot / (h0 - hL))**2
      dhdhL   = (d0 * dpot)*2 / (hL - h0 + d0 * dpot)**2
      dhdh0   = (h0 - hL) * (h0 - 2 * d0 * dpot - hL) / (hL - h0 + d0 * dpot)**2
      dhdd0   = dpot / (1 - d0 * dpot / (h0 - hL))**2
      dhdn    = dhdhL * dhLdn + dhdh0 * dh0dn + dhdd0 * dd0dn

    elseif (dpot >= deta) then
      hR = detadn(1) / n(1)
      call smc%get_dist(eta(1), 1, dndeta(1), d2ndeta2(1))
      d2etadn2(1) = - d2ndeta2(1) / dndeta(1)**3
      dhRdn(1) = d2etadn2(1) / n(1) - detadn(1) / n(1)**2
      dhRdn(2) = 0

      h       = (h1 + d1 * (dpot - deta) * hR / (hR - h1)) / (1 + d1 * (dpot - deta) / (hR - h1))
      dhddpot = d1 / (1 + d1 * (dpot - deta) / (hR - h1))**2
      dhdhR   = (d1 * (dpot - deta))**2 / (hR - h1 + d1 * (dpot - deta))**2
      dhdh1   = (h1 - hR) * (h1 - 2 * d1 * (dpot - deta) - hR) / (hR - h1 + d1 * (dpot - deta))**2
      dhdd1   = (dpot - deta) / (1 + d1 * (dpot - deta) / (hR - h1))**2
      dhdn    = dhdhR * dhRdn + dhdh1 * dh1dn + dhdd1 * dd1dn - dhddpot * [-1, 1] * detadn
    elseif (h0 == h1) then
      h       = h0
      dhdn    = dh0dn
      dhddpot = 0
    else
      c1    = (h1 - h0) / deta
      dc1dn = (dh1dn - dh0dn - (h1 - h0) * [-1, 1] * detadn / deta) / deta
      c2    = (h1 * d0 + h0 * d1) / (deta * (h1 - h0)) - (h0 + h1) / deta**2
      dc2dn = (dh1dn * d0 + h1 * dd0dn + dh0dn * d1 + h0 * dd1dn) / (deta * (h1 - h0)) &
        &   - (h1 * d0 + h0 * d1) / (deta * (h1 - h0))**2 * ([-1, 1] * detadn * (h1 - h0) + deta * (dh1dn - dh0dn)) &
        &   + (- (dh0dn + dh1dn) + 2 * (h0 + h1) * [-1, 1] * detadn / deta) / deta**2
      d     = (d0 + d1) / (deta * (h1 - h0)) - 2 / deta**2
      dddn  = (dd0dn + dd1dn) / (deta * (h1 - h0)) &
        &   - (d0 + d1) / (deta * (h1 - h0))**2 * ([-1,1]*detadn * (h1 - h0) + deta * (dh1dn - dh0dn)) &
        &   + 4 / deta**3 * [-1, 1] * detadn

      p       = h0 + c1 * dpot + c2 * dpot * (deta - dpot)
      dpdn    = dh0dn + dpot * (dc1dn + dc2dn * (deta - dpot) + c2 * [-1, 1] * detadn)
      dpddpot = c1 + c2 * (deta - 2 * dpot)

      q       = 1 + d * dpot * (deta - dpot)
      dqdn    = dddn * dpot * (deta - dpot) + d * dpot * [-1, 1] * detadn
      dqddpot = d * (deta - 2 * dpot)

      h       = p / q
      dhdn    = (dpdn - p * dqdn / q) / q
      dhddpot = (dpddpot - p * dqddpot / q) / q
    end if

    B1  = ber(-dpot * h)
    B2  = ber( dpot * h)
    dB1 = dberdx(-dpot * h)
    dB2 = dberdx( dpot * h)

    j       = (B1 * n(1) - B2 * n(2)) / h
    djdh    = - (j + dpot * (dB1 * n(1) + dB2 * n(2))) / h
    djdn    = [B1, - B2] / h + djdh * dhdn
    djddpot = - (dB1 * n(1) + dB2 * n(2)) + djdh * dhddpot
  end subroutine

end module
