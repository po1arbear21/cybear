m4_include(util/macro.f90.inc)

module current_density_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_value

  use current_integral_m, only: current_integral_get
  use density_m,          only: density
  use device_params_m,    only: device_params
  use equation_m,         only: equation
  use error_m,            only: assert_failed, program_error
  use grid_m,             only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,        only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m,   only: DIR_NAME
  use jacobian_m,         only: jacobian
  use math_m,             only: ber, dberdx
  use mobility_m,         only: mobility
  use potential_m,        only: potential
  use semiconductor_m,    only: CR_NAME, CR_CHARGE, DIST_MAXWELL, STAB_ED, STAB_EXACT, STAB_SG, semiconductor
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

    procedure, private :: get_curr    => calc_current_density_get_curr
    procedure, private :: get_curr_sg => calc_current_density_get_curr_sg
    procedure, private :: get_curr_ed => calc_current_density_get_curr_ed
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
    if (par%smc(par%smc_default)%mob_sat) this%jaco_mob => this%init_jaco(iprov, this%depend(mob, par%transport(IDX_EDGE, idx_dir)), st = [this%st_dir%get_ptr()])

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    !! evaluate drift-diffusion equation
    class(calc_current_density), intent(inout) :: this

    integer               :: i, idx_dir, idx_dim, start, end, rate, ci
    real                  :: pot(2), dens(2), be(2), edos_e(2), j, djdpot(2), djddens(2), djdmob, len, mob
    integer, allocatable  :: idx(:), idx1(:), idx2(:)
    logical               :: status

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    ci      = this%cdens%ci
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim))
    call system_clock(start, rate)

    ! loop over transport edges in parallel
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(i,pot,dens,be,edos_e,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    !$omp shared(this,idx_dir,idx_dim,ci)
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

      ! parameters
      len       = this%par%g%get_len(idx, idx_dir)
      pot(   1) = this%pot%get(idx1)
      pot(   2) = this%pot%get(idx2)
      dens(  1) = this%dens%get(idx1)
      dens(  2) = this%dens%get(idx2)
      be(    1) = this%par%band_edge_v(ci)%get(idx1)
      be(    2) = this%par%band_edge_v(ci)%get(idx2)
      edos_e(1) = this%par%edos_v(     ci)%get(idx1)
      edos_e(2) = this%par%edos_v(     ci)%get(idx2)
      if (this%par%smc(this%par%smc_default)%mob_sat) then
        mob = this%mob%get(idx)
      else
        mob = this%par%mob0(IDX_EDGE, idx_dir, this%cdens%ci)%get(idx)
      end if

      ! get current along edge
      call this%get_curr(this%par%smc(this%par%smc_default), len, pot, dens, be, edos_e, mob, j, djdpot, djddens, djdmob)

      ! set current density + derivatives
      call this%cdens%set(idx, j)
      call this%jaco_pot%add( idx, idx1, djdpot(1))
      call this%jaco_pot%add( idx, idx2, djdpot(2))
      call this%jaco_dens%add(idx, idx1, djddens(1))
      call this%jaco_dens%add(idx, idx2, djddens(2))
      if (this%par%smc(this%par%smc_default)%mob_sat) call this%jaco_mob%add(idx, idx, djdmob)
    end do
    !$omp end parallel do

    call system_clock(end)
    time = time + real(end-start)/real(rate)
  end subroutine

  subroutine calc_current_density_get_curr(this, smc, len, pot, dens, be, edos_e, mob, j, djdpot, djddens, djdmob)
    !! get current along edge. The Scharfetter-Gummel argument uses the per-vertex
    !! band edge values be(:) so that the conduction/valence-band step at a
    !! heterojunction edge contributes to the drift term. The density is
    !! normalized per-vertex by edos_e(:) so n = dens/edos collapses to the
    !! single-material expression when both endpoints share the same material.
    !! Reference: Selberherr §4.3 / Yu & Dutton.
    class(calc_current_density), intent(in)  :: this
    type(semiconductor),         intent(in)  :: smc
      !! semiconductor material (used for dist function, stabilization method; global per v1)
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: be(2)
      !! band edge at edge endpoints (Ec for electrons, Ev for holes)
    real,                        intent(in)  :: edos_e(2)
      !! effective DOS at edge endpoints
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
    real    :: ch, dpot, edos_avg, n(2), djdn(2), djddpot

    ! abbreviations
    ci       = this%cdens%ci
    ch       = CR_CHARGE(ci)

    ! normalize potential drop and density. The band-edge term must carry the
    ! per-carrier sign `ch`: the SG zero condition is Δln(n) = dpot, and from
    ! the simulator's imref convention (iref = pot - be + ch·eta, with
    ! eta = ln(dens/edos)) at flat quasi-Fermi level Δln(n) = -ch·(Δpot - Δbe).
    ! For electrons this gives Δpot - Δbe; for holes, Δbe - Δpot. The (be(2) -
    ! be(1)) term collapses to zero when the two endpoints share the same
    ! material, recovering the single-material expression bit-identically.
    dpot = - ch * ((pot(2) - pot(1)) - (be(2) - be(1)))
    n(1) = dens(1) / edos_e(1)
    n(2) = dens(2) / edos_e(2)
    ! geometric mean of edos for denormalization — better paired with the
    ! per-vertex n normalization than arithmetic mean (both reduce to edos
    ! in the single-material limit).
    edos_avg = sqrt(edos_e(1) * edos_e(2))

    if (smc%dist == DIST_MAXWELL .and. (.not. smc%reg)) then
      ! non-degenerate case: always use Scharfetter-Gummel
      call this%get_curr_sg(n, dpot, j, djdn, djddpot)
    else
      ! get current based on different stabilization method
      select case (smc%stab)
      case (STAB_SG)
        call this%get_curr_sg(n, dpot, j, djdn, djddpot)
      case (STAB_ED)
        call this%get_curr_ed(smc, n, dpot, j, djdn, djddpot)
      case (STAB_EXACT)
        call current_integral_get(dist, inv_dist, int_dist, n, dpot, j, djdn, djddpot)
      end select

      if (.not. ieee_is_finite(j)) then
        print "(A,ES25.16E3)", "n1   = ", n(1)
        print "(A,ES25.16E3)", "n2   = ", n(2)
        print "(A,ES25.16E3)", "dpot = ", dpot
        print "(A,I6)",        "stab = ", smc%stab
        error stop "j is not finite"
      end if
    end if

    ! denormalize. edos_avg is the natural choice for the homogeneous limit:
    ! edos_e(1)=edos_e(2)=edos => edos_avg=edos and djddens stays edos-free.
    j         = j * mob * edos_avg / len
    djdpot    = ch * [1.0, -1.0] * djddpot * mob * edos_avg / len
    djddens(1) = djdn(1) * mob * edos_avg / (edos_e(1) * len)
    djddens(2) = djdn(2) * mob * edos_avg / (edos_e(2) * len)
    djdmob    = j / mob
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

      real, parameter :: ITOL = 1e-13

      call smc%get_int_dist(eta, k, ITOL, I, dIdeta)
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

    real, parameter :: DETA_TOL = 1e-5

    real :: eta(2), detadn(2), deta, etam, F0, dF0, F1, dF1, F2, dF2, F3, dF3
    real :: a0, da0dn(2), a2, da2dn(2), a, dadn(2)

    call smc%get_inv_dist(n(1), eta(1), detadn(1))
    call smc%get_inv_dist(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

    if (abs(deta) > DETA_TOL) then
      a    = log(n(2) / n(1)) / deta
      dadn = [1.0, -1.0] * (a * detadn - 1 / n) / deta
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
      a0    = F1
      da0dn = dF1 * 0.5 * detadn
      a2    = 0.25  * ((0.5 * F3 + F1**3) / 3 - 0.5 * F1 * F2)
      da2dn = 0.125 * ((0.5 * dF3 + 3*F1**2*dF1) / 3 - 0.5 * (dF1 * F2 + F1 * dF2)) * detadn
      a     = a0 + deta**2 * a2
      dadn  = da0dn + 2*deta*[-1.0,1.0]*detadn * a2 + deta**2 * da2dn
    end if

    call this%get_curr_sg(n, dpot * a, j, djdn, djddpot)
    j    = j / a
    djdn = (djdn + (djddpot * dpot - j)  * dadn) / a
  end subroutine

end module
