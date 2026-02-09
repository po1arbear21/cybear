module surface_charge_m
  !! Surface charge density from Fermi-level pinning at surface defects
  !!
  !! Following Haney et al. 2016, Sec. IV:
  !!   rho_surf = q * N_surf/2 * (1 - 2*f_surf)
  !!   f_surf = (n_s + p_trap) / (n_s + p_s + n_trap + p_trap)
  !!
  !! where:
  !!   n_trap = N_c * exp((E_surf - E_c) / kT)
  !!   p_trap = N_v * exp((E_v - E_surf) / kT)
  !!
  !! The surface charge enters the Poisson equation as a boundary source term.
  !!
  !! NOTE: Currently only 2D simulations are supported. The surface is defined
  !! as the x-boundaries (x=0 and x=Lx), with area integration along y.

  use density_m,       only: density
  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_table_m,    only: grid_table
  use jacobian_m,      only: jacobian
  use semiconductor_m, only: CR_ELEC, CR_HOLE
  use variable_m,      only: variable
  use normalization_m, only: denorm

  implicit none

  private
  public surface_charge, calc_surface_charge

  type, extends(variable) :: surface_charge
    !! Surface charge density variable
    !!
    !! Unit: C/cm^2 (charge per unit area, normalized)
    !! Positive = net positive charge at surface

    real, pointer :: x1(:)     => null()
    real, pointer :: x2(:,:)   => null()
    real, pointer :: x3(:,:,:) => null()
  contains
    procedure :: init => surface_charge_init
  end type

  type, extends(equation) :: calc_surface_charge
    !! Calculate surface charge from carrier densities and Fermi-level pinning

    type(device_params), pointer :: par => null()
    type(surface_charge), pointer :: scharge => null()
    type(density), pointer :: dens_n => null()
    type(density), pointer :: dens_p => null()

    ! Physics parameters (normalized)
    real :: E_surf      ! Surface defect level relative to midgap (kT units)
    real :: Ns_dens     ! Surface state density (cm^-2, normalized)
    real :: n_trap      ! n_trap = Nc * exp((E_surf - Ec) / kT)
    real :: p_trap      ! p_trap = Nv * exp((Ev - E_surf) / kT)

    ! Surface vertex information
    integer :: n_surf_vert
      !! number of surface vertices
    integer, allocatable :: surf_idx(:,:)
      !! surface vertex indices (idx_dim x n_surf_vert)
    type(grid_table) :: surf_tab
      !! grid table for surface vertices only (for efficient Jacobian)

    ! Jacobians
    type(jacobian), pointer :: jaco_dens_n => null()
      !! d(rho_surf)/dn
    type(jacobian), pointer :: jaco_dens_p => null()
      !! d(rho_surf)/dp
  contains
    procedure :: init => calc_surface_charge_init
    procedure :: eval => calc_surface_charge_eval
  end type

contains

  subroutine surface_charge_init(this, par)
    !! Initialize surface charge variable
    class(surface_charge), intent(out) :: this
    type(device_params),   intent(in)  :: par

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    call this%variable_init("rho_surf", "C/cm^2", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)

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
    end select
  end subroutine

  subroutine calc_surface_charge_init(this, par, dens_n, dens_p, scharge)
    !! Initialize surface charge calculation equation
    class(calc_surface_charge),   intent(out)   :: this
    type(device_params),  target, intent(in)    :: par
    type(density),        target, intent(in)    :: dens_n
    type(density),        target, intent(in)    :: dens_p
    type(surface_charge), target, intent(in)    :: scharge

    integer              :: i, iprov, idep_n, idep_p, n_alloc
    integer, allocatable :: idx(:)
    real                 :: Eg

    print "(A)", "calc_surface_charge_init"

    ! Dimension check: charged surface currently only supports 2D
    if (par%g%idx_dim /= 2) then
      print "(A)", "ERROR: charged_surf requires exactly 2D simulation"
      print "(A,I0)", "       Current dimension: ", par%g%idx_dim
      error stop "calc_surface_charge_init: charged surface only implemented for 2D"
    end if

    ! init equation
    call this%equation_init("calc_rho_surf")
    this%par     => par
    this%scharge => scharge
    this%dens_n  => dens_n
    this%dens_p  => dens_p

    ! Store physics parameters
    this%E_surf = par%smc%E_surf
    this%Ns_dens = par%smc%N_surf

    ! Calculate n_trap and p_trap from E_surf
    ! E_surf is relative to midgap, so:
    !   n_trap = Nc * exp((E_surf - (Ec - E_mid)) / kT) = Nc * exp(E_surf - Eg/2)
    !   p_trap = Nv * exp(((E_mid - Ev) - E_surf) / kT) = Nv * exp(-E_surf - Eg/2)
    ! Note: band_gap is already in kT units, E_surf is also in kT units
    Eg = par%smc%band_gap
    this%n_trap = par%smc%edos(CR_ELEC) * exp(this%E_surf - 0.5 * Eg)
    this%p_trap = par%smc%edos(CR_HOLE) * exp(-this%E_surf - 0.5 * Eg)

    print "(A,ES12.4,A)", "  E_surf = ", denorm(this%E_surf, 'eV'), " eV (relative to midgap)"
    print "(A,ES12.4,A)", "  N_surf = ", denorm(this%Ns_dens, 'cm^-2'), " cm^-2"
    print "(A,ES12.4,A)", "  n_trap = ", denorm(this%n_trap, 'cm^-3'), " cm^-3"
    print "(A,ES12.4,A)", "  p_trap = ", denorm(this%p_trap, 'cm^-3'), " cm^-3"

    ! Count surface vertices: x=1 or x=nx, AND not on contact
    allocate(idx(par%g%idx_dim))
    this%n_surf_vert = 0
    do i = 1, par%transport(IDX_VERTEX, 0)%n
      idx = par%transport(IDX_VERTEX, 0)%get_idx(i)
      if ((idx(1) == 1 .or. idx(1) == par%g1D(1)%n) .and. par%ict%get(idx) == 0) then
        this%n_surf_vert = this%n_surf_vert + 1
      end if
    end do

    print "(A,I0,A)", "  Found ", this%n_surf_vert, " surface vertices for charged surface"

    ! Allocate and fill surface vertex arrays + create grid_table for surface vertices
    n_alloc = max(this%n_surf_vert, 1)
    allocate(this%surf_idx(par%g%idx_dim, n_alloc))

    ! Initialize grid_table for surface vertices (all flags initially false)
    call this%surf_tab%init("surf_charge", par%g, IDX_VERTEX, 0)

    this%n_surf_vert = 0
    do i = 1, par%transport(IDX_VERTEX, 0)%n
      idx = par%transport(IDX_VERTEX, 0)%get_idx(i)
      if ((idx(1) == 1 .or. idx(1) == par%g1D(1)%n) .and. par%ict%get(idx) == 0) then
        this%n_surf_vert = this%n_surf_vert + 1
        this%surf_idx(:, this%n_surf_vert) = idx
        ! Set flag in grid_table for this surface vertex
        call this%surf_tab%flags%set(idx, .true.)
      end if
    end do

    ! Finalize grid_table (builds internal lookup tables)
    call this%surf_tab%init_final()

    ! provide surface charge on transport vertices (Poisson depends on this)
    iprov = this%provide(scharge, par%transport_vct(0))

    ! dependencies on carrier densities (must match provide selector for Jacobian dimensions)
    idep_n = this%depend(dens_n, par%transport_vct(0))
    idep_p = this%depend(dens_p, par%transport_vct(0))

    ! initialize Jacobians (non-constant since they depend on carrier densities)
    this%jaco_dens_n => this%init_jaco(iprov, idep_n, const = .false.)
    this%jaco_dens_p => this%init_jaco(iprov, idep_p, const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_surface_charge_eval(this)
    !! Evaluate surface charge density
    !!
    !! rho_surf = q * Ns_dens/2 * (1 - 2*f_surf)
    !! f_surf = (n_s + p_trap) / (n_s + p_s + n_trap + p_trap)
    !!
    !! Jacobians:
    !! df/dn = (p_s + n_trap) / D^2
    !! df/dp = -(n_s + p_trap) / D^2
    !! d(rho_surf)/dn = -q * Ns_dens * df/dn
    !! d(rho_surf)/dp = -q * Ns_dens * df/dp
    class(calc_surface_charge), intent(inout) :: this

    integer              :: j
    integer, allocatable :: idx(:)
    real                 :: n_s, p_s, D, f, rho, df_dn, df_dp, drho_dn, drho_dp
    real                 :: qN  ! q * Ns_dens (in normalized units, q=1)

    allocate(idx(this%par%g%idx_dim))

    ! q = 1 in normalized units, so qN = Ns_dens
    qN = this%Ns_dens

    ! Note: scharge is initialized to zero, and we only set non-zero values at surface vertices.
    ! Jacobian entries are only set for surface vertices (sparse).
    ! Non-surface Jacobian entries remain at their initialized values (zero).

    ! Evaluate at surface vertices only
    do j = 1, this%n_surf_vert
      idx = this%surf_idx(:, j)

      n_s = this%dens_n%get(idx)
      p_s = this%dens_p%get(idx)

      ! Denominator: D = n_s + p_s + n_trap + p_trap
      D = n_s + p_s + this%n_trap + this%p_trap

      if (D > 1e-30) then
        ! Surface occupancy: f = (n_s + p_trap) / D
        f = (n_s + this%p_trap) / D

        ! Surface charge: rho = q * Ns_dens/2 * (1 - 2*f)
        rho = 0.5 * qN * (1.0 - 2.0 * f)

        ! Jacobians
        ! df/dn = (D - (n_s + p_trap)) / D^2 = (p_s + n_trap) / D^2
        df_dn = (p_s + this%n_trap) / (D * D)
        ! df/dp = -(n_s + p_trap) / D^2
        df_dp = -(n_s + this%p_trap) / (D * D)

        ! d(rho)/dn = -q * Ns_dens * df/dn
        drho_dn = -qN * df_dn
        ! d(rho)/dp = -q * Ns_dens * df/dp
        drho_dp = -qN * df_dp
      else
        rho     = 0.0
        drho_dn = 0.0
        drho_dp = 0.0
      end if

      call this%scharge%set(idx, rho)
      call this%jaco_dens_n%add(idx, idx, drho_dn)
      call this%jaco_dens_p%add(idx, idx, drho_dp)

      ! Debug output for first few surface vertices (only on first call)
      ! if (j <= 3) then
      !   print "(A,I3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)", &
      !     "  surf_vert ", j, ": n_s=", denorm(n_s, 'cm^-3'), &
      !     ", p_s=", denorm(p_s, 'cm^-3'), ", f=", f, ", rho=", denorm(rho, '1/cm^2')
      ! end if
    end do

  end subroutine

end module
