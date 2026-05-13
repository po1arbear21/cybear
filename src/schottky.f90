module schottky_m

  use contact_m,        only: CT_SCHOTTKY
  use density_m,        only: density
  use device_params_m,  only: device_params
  use electric_field_m, only: electric_field
  use equation_m,       only: equation
  use error_m,          only: program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE
  use jacobian_m,       only: jacobian
  use math_m,           only: PI, expm1, log1p
  use normalization_m,  only: norm, denorm
  use quad_m,           only: quad
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, CR_NAME
  use stencil_m,        only: dirichlet_stencil, empty_stencil, stencil_ptr
  use variable_m,       only: variable

  implicit none
  private

  public :: schottky_n0b
  public :: schottky_tunneling
  public :: schottky_bc
  public :: calc_schottky_bc

  type, extends(variable) :: schottky_bc
    !! Schottky boundary carrier flux at one Schottky contact.
    !!   j_bc(idx) = v_th * exp(dphi) * dens - v_th * n0b - J_t
    !! Units: 1/cm^2/s. Continuity folds in -A_ct * j_bc per Schottky vertex.

    integer :: ci    !! carrier index (CR_ELEC / CR_HOLE)
    integer :: ict   !! contact index (must be CT_SCHOTTKY)
  contains
    procedure :: init => schottky_bc_init
  end type

  type, extends(equation) :: calc_schottky_bc
    !! Compute schottky_bc + Jacobian w.r.t. dens at one Schottky contact.
    !! Lagged dF/dE_normal: empty stencil on efield (value read at eval, no Jacobian).

    type(device_params),  pointer :: par           => null()
    type(density),        pointer :: dens          => null()
    type(electric_field), pointer :: efield_normal => null()
    type(schottky_bc),    pointer :: bc            => null()

    integer :: normal_dir         !! cached at init via get_normal_dir
    real    :: v_th               !! cached at init via schottky_velocity
    real, allocatable :: eps_sc(:)
      !! per-contact-vertex semiconductor permittivity (from the adjacent edge in
      !! the normal direction); sized transport_vct(ict)%n. Cached at init so
      !! IFBL sees the local permittivity at heterointerfaces along the contact.

    type(dirichlet_stencil) :: st_dir
    type(empty_stencil)     :: st_em

    type(jacobian), pointer :: jaco_dens   => null()
    type(jacobian), pointer :: jaco_efield => null()
  contains
    procedure :: init => calc_schottky_bc_init
    procedure :: eval => calc_schottky_bc_eval
  end type

contains


  function get_normal_dir(par, ict) result(dir)
    !! Determine which coordinate direction is normal to a contact surface.
    !! The normal direction is the one with zero coordinate range across the
    !! contact vertices (i.e., the flat direction). Hard-fails if no such
    !! direction exists, since an oblique contact would silently produce
    !! a meaningless IFBL geometry.

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict
    integer :: dir

    integer :: i, d, n
    integer, allocatable :: idx(:)
    real :: cmin(3), cmax(3), range(3), pos(3)
    real, parameter :: REL_TOL = 1.0e-6

    dir = 1
    cmin =  huge(0.0)
    cmax = -huge(0.0)

    allocate(idx(par%g%idx_dim))
    n = par%transport_vct(ict)%n
    if (n == 0) return

    do i = 1, n
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, pos(1:par%g%dim))
      do d = 1, par%g%dim
        cmin(d) = min(cmin(d), pos(d))
        cmax(d) = max(cmax(d), pos(d))
      end do
    end do

    range = cmax - cmin
    do d = 2, par%g%dim
      if (range(d) < range(dir)) dir = d
    end do

    if (range(dir) > REL_TOL * maxval(range(1:par%g%dim))) &
      call program_error("get_normal_dir: contact "//par%contacts(ict)%name// &
        & " is not axis-aligned (no flat direction found)")
  end function get_normal_dir


  function schottky_velocity(par, ci, ict) result(v_th)
    !! Calculate thermionic injection velocity for Robin BC.
    !!
    !! v_th = A* T^2 / (q N_c)
    !!
    !! This is the surface recombination velocity coefficient in the
    !! Robin BC: J_TE = v_th * (n - n0b)

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci   !! carrier index (CR_ELEC/CR_HOLE)
    integer, intent(in) :: ict  !! contact index
    real :: v_th

    real :: A_richardson

    if (ci == CR_ELEC) then
      A_richardson = par%contacts(ict)%A_richardson_n
    else
      A_richardson = par%contacts(ict)%A_richardson_p
    end if

    v_th = A_richardson * norm(par%T, "K")**2 / par%smc(par%smc_default)%edos(ci)
  end function schottky_velocity


  pure subroutine ifbl_delta_phi(ifbl, E_normal, eps_sc, dphi)
    !! Image-force barrier lowering at a Schottky contact.
    !!   dphi = sqrt(|E_normal| / (4 pi eps_sc))
    !! Returns 0 when IFBL is disabled or |E_normal| is below a numerical floor
    !! (avoids a singular sqrt at zero field). All quantities normalized.

    logical, intent(in)  :: ifbl
    real,    intent(in)  :: E_normal, eps_sc
    real,    intent(out) :: dphi

    real, parameter :: E_TOL = 1.0e-10

    if (ifbl .and. abs(E_normal) > E_TOL) then
      dphi = sqrt(abs(E_normal) / (4.0 * PI * eps_sc))
    else
      dphi = 0.0
    end if
  end subroutine ifbl_delta_phi


  subroutine schottky_n0b(par, ci, ict, E_normal, n0b, eps_sc, delta_phi)
    !! Calculate equilibrium boundary density at Schottky barrier.
    !!
    !! n0b = N_c exp(-phi_b_eff / kT),  phi_b_eff = max(phi_b - dphi, 0)
    !! where dphi is the IFBL barrier lowering (zero when IFBL is disabled).

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci, ict
    real, intent(in) :: E_normal
    real, intent(out) :: n0b
    real, intent(in) :: eps_sc
      !! normalized semiconductor permittivity (used only when IFBL is enabled)
    real, optional, intent(out) :: delta_phi
      !! effective barrier lowering, clamped to phi_b (normalized to kT).
      !! Clamping keeps n0b and the v_th*exp(delta_phi) prefactor in eval consistent,
      !! so detailed balance holds even when |E| is large enough that the raw
      !! image-force formula would exceed phi_b.

    real :: phi_b, phi_b_eff, dphi

    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
    else
      phi_b = par%smc(par%smc_default)%band_gap - par%contacts(ict)%phi_b
    end if

    call ifbl_delta_phi(par%contacts(ict)%ifbl, E_normal, eps_sc, dphi)
    dphi = min(dphi, phi_b)
    phi_b_eff = phi_b - dphi

    n0b = par%smc(par%smc_default)%edos(ci) * exp(-phi_b_eff)

    if (present(delta_phi)) delta_phi = dphi

  end subroutine schottky_n0b


  subroutine schottky_tunneling(par, ci, ict, E_normal, dens, J_t, dJ_t_ddens, eps_sc)
    !! Calculate sub-barrier tunneling current via Tsu-Esaki integral.
    !!
    !! J_t = (m* / 2 pi^2) * integral_0^{phi_b} T(E) * supply(E) dE
    !!
    !! The above-barrier thermionic component is handled separately
    !! by the Robin BC velocity term.
    !!
    !! NOTE: validated for electrons (CR_ELEC). The hole branch reuses the
    !! same integrand with m_tunnel_p and phi_b_hole = band_gap - phi_b;
    !! this symmetric form has not yet been benchmarked against an
    !! independent hole-tunneling reference.

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict, ci
    real, intent(in) :: E_normal, dens
    real, intent(out) :: J_t
    real, intent(out) :: dJ_t_ddens
    real, intent(in) :: eps_sc
      !! normalized semiconductor permittivity (used only when IFBL is enabled)

    real :: phi_b, m_tunnel, eta, eta_m, detadF, dphi
    real :: params(4), integral, dIda, dIdb, dIdp(4), err
    real :: prefactor
    integer :: ncalls

    if (.not. par%contacts(ict)%tunneling .or. dens <= 0.0) then
      J_t = 0.0
      dJ_t_ddens = 0.0
      return
    end if

    ! barrier height and tunneling mass
    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_n
    else
      phi_b = par%smc(par%smc_default)%band_gap - par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_p
    end if

    prefactor = m_tunnel / (2.0 * PI**2)

    ! compute eta from density via inverse distribution function
    call par%smc(par%smc_default)%get_inv_dist(dens / par%smc(par%smc_default)%edos(ci), eta, detadF)

    ! quasi-Fermi level relative to metal Fermi level (unmodified phi_b:
    ! eta_m = V_n/kT represents the Fermi level splitting, independent of IFBL)
    eta_m = eta + phi_b

    ! apply IFBL to phi_b for WKB barrier shape and integration bounds
    call ifbl_delta_phi(par%contacts(ict)%ifbl, E_normal, eps_sc, dphi)
    phi_b = max(phi_b - dphi, 0.0)

    ! integration parameters: [phi_b, |E|, m*, eta_m]
    params = [phi_b, E_normal, m_tunnel, eta_m]

    ! integrate from 0 to phi_b (tunneling regime only)
    call quad(tsu_esaki_integrand, 0.0, phi_b, params, integral, &
              dIda, dIdb, dIdp, rtol=1.0e-6, err=err, ncalls=ncalls)

    ! current density (negative = injection into semiconductor)
    J_t = -prefactor * integral

    ! derivative dJ_t/d(dens) via chain rule:
    ! dJ/d(dens) = dJ/d(eta_m) * d(eta_m)/d(eta) * d(eta)/d(dens)
    dJ_t_ddens = -prefactor * dIdp(4) * detadF / par%smc(par%smc_default)%edos(ci)
  end subroutine schottky_tunneling


  subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
    !! Integrand: f(E) = T(E) * supply(E)
    !! where supply = ln(1+exp(eta_m-E)) - ln(1+exp(-E))
    !! and T = WKB transmission coefficient

    real, intent(in)  :: E, p(:)
    real, intent(out) :: f, dfdE, dfdp(:)

    real :: phi_b, efield, m_star, eta_m
    real :: T, dT_dE, dT_dphi, dT_dF
    real :: N_supply, dN_dE, f_metal, f_semi

    phi_b  = p(1)
    efield = p(2)
    m_star = p(3)
    eta_m  = p(4)

    call wkb_transmission(E, phi_b, efield, m_star, T, dT_dE, dT_dphi, dT_dF)

    f_metal = 1.0 / (1.0 + exp(E))
    f_semi  = 1.0 / (1.0 + exp(E - eta_m))
    dN_dE   = f_metal - f_semi

    ! supply function: rewritten for numerical stability
    ! log1p(expm1(eta_m) * f_metal) avoids catastrophic cancellation
    N_supply = log1p(expm1(eta_m) * f_metal)

    f = T * N_supply
    dfdE = dT_dE * N_supply + T * dN_dE

    dfdp(1) = dT_dphi * N_supply
    dfdp(2) = dT_dF * N_supply
    dfdp(3) = 0.0
    dfdp(4) = T * f_semi
  end subroutine tsu_esaki_integrand


  pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
    !! Standard WKB transmission coefficient for triangular barrier.
    !!
    !! T = exp(-gamma)  for E < phi_b (tunneling)
    !! T = 1            for E >= phi_b (thermionic, above barrier)
    !!
    !! gamma = (4/3) * sqrt(2 m*) * (phi_b - E)^(3/2) / F

    real, intent(in)  :: E, phi_b, F, m_star
    real, intent(out) :: T, dT_dE, dT_dphi, dT_dF

    real :: gamma, E_diff, F_safe, coeff, exp_gamma
    real, parameter :: F_reg = 1e-10

    E_diff = phi_b - E
    F_safe = sqrt(F**2 + F_reg**2)
    coeff = (4.0/3.0) * sqrt(2.0 * m_star)

    if (E_diff < 0.0) then
      T = 1.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
    else
      gamma = coeff * E_diff**1.5 / F_safe
      exp_gamma = exp(-gamma)
      T = exp_gamma
      dT_dE   =  exp_gamma * coeff * 1.5 * sqrt(E_diff) / F_safe
      dT_dphi = -dT_dE
      dT_dF   =  exp_gamma * coeff * E_diff**1.5 * F / (F_safe**3)
    end if
  end subroutine wkb_transmission


  subroutine schottky_bc_init(this, par, ict, ci)
    !! initialize Schottky boundary current variable at contact ict for carrier ci
    class(schottky_bc),  intent(out) :: this
    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ict, ci

    if (par%contacts(ict)%type /= CT_SCHOTTKY) &
      call program_error("schottky_bc_init: contact "//par%contacts(ict)%name//" is not Schottky")

    call this%variable_init("schottky_bc_"//CR_NAME(ci)//"_"//par%contacts(ict)%name, "1/cm^2/s", &
      &                     g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci  = ci
    this%ict = ict
  end subroutine


  subroutine calc_schottky_bc_init(this, par, dens, efield, bc)
    !! initialize the calc equation for one Schottky contact
    class(calc_schottky_bc),       intent(out) :: this
    type(device_params),   target, intent(in)  :: par
    type(density),         target, intent(in)  :: dens
    type(electric_field),  target, intent(in)  :: efield(:)
      !! all electric field components; the contact-normal one is selected internally
    type(schottky_bc),     target, intent(in)  :: bc

    integer              :: iprov, idep_dens, idep_efield, ict, ci, i
    integer, allocatable :: idx(:), idx1(:)
    logical              :: status

    ict = bc%ict
    ci  = bc%ci

    print "(A)", "calc_schottky_bc_init for "//bc%name

    call this%equation_init("calc_"//bc%name)
    this%par  => par
    this%dens => dens
    this%bc   => bc

    ! resolve the contact-normal direction and bind the matching efield component
    this%normal_dir    =  get_normal_dir(par, ict)
    this%efield_normal => efield(this%normal_dir)
    this%v_th          =  schottky_velocity(par, ci, ict)

    ! per-vertex semiconductor permittivity from the adjacent edge in the normal
    ! direction; populated once at init so IFBL in eval just indexes by i
    allocate(idx(par%g%idx_dim), idx1(par%g%idx_dim))
    allocate(this%eps_sc(par%transport_vct(ict)%n))
    do i = 1, par%transport_vct(ict)%n
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, this%normal_dir, idx, 1, idx1, status)
      if (.not. status) call program_error("calc_schottky_bc_init: contact "//par%contacts(ict)%name// &
        &                                  " has a vertex with no neighboring edge in the normal direction")
      this%eps_sc(i) = par%eps(IDX_EDGE, this%normal_dir)%get(idx1)
    end do

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_em%init()

    ! provide schottky_bc on this contact's vertices
    iprov = this%provide(bc, par%transport_vct(ict))

    ! depend on density at this contact's vertices (non-const diagonal Jacobian)
    idep_dens = this%depend(dens, par%transport_vct(ict))
    this%jaco_dens => this%init_jaco(iprov, idep_dens, &
      & st = [this%st_dir%get_ptr()], const = .false.)

    ! depend on E-field for evaluation order; lagged (empty stencil = no Jacobian entries)
    idep_efield = this%depend(efield(this%normal_dir), par%transport_vct(ict))
    this%jaco_efield => this%init_jaco(iprov, idep_efield, &
      & st = [this%st_em%get_ptr()], const = .true.)

    call this%init_final()
  end subroutine


  subroutine calc_schottky_bc_eval(this)
    !! evaluate Schottky boundary current and its derivative w.r.t. density
    class(calc_schottky_bc), intent(inout) :: this

    integer :: i, ict, ci
    integer :: idx(this%par%g%idx_dim)
    real    :: E_normal, dens_val, n0b, J_t, dJ_t_ddens, delta_phi, j_bc

    ict = this%bc%ict
    ci  = this%bc%ci

    do i = 1, this%par%transport_vct(ict)%n
      idx = this%par%transport_vct(ict)%get_idx(i)

      E_normal = this%efield_normal%get(idx)
      dens_val = this%dens%get(idx)

      ! equilibrium boundary density (with optional IFBL)
      call schottky_n0b(this%par, ci, ict, E_normal, n0b, eps_sc = this%eps_sc(i), delta_phi = delta_phi)

      ! tunneling current and its derivative w.r.t. density
      call schottky_tunneling(this%par, ci, ict, E_normal, dens_val, J_t, dJ_t_ddens, eps_sc = this%eps_sc(i))

      ! boundary current density: j_bc = v_th * exp(dphi) * dens - v_th * n0b - J_t
      j_bc = this%v_th * exp(delta_phi) * dens_val - this%v_th * n0b - J_t
      call this%bc%set(idx, j_bc)

      ! diagonal Jacobian entry: dj_bc/d(dens) = v_th * exp(dphi) - dJ_t/d(dens)
      call this%jaco_dens%add(idx, idx, this%v_th * exp(delta_phi) - dJ_t_ddens)
    end do
  end subroutine

end module schottky_m
