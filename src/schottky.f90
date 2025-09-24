module schottky_m

  use device_params_m,  only: device_params
  use normalization_m,  only: norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, CR_NAME
  use grid_m,           only: IDX_VERTEX
  use variable_m,       only: variable_real
  use equation_m,       only: equation
  use jacobian_m,       only: jacobian
  use electric_field_m, only: electric_field
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use vselector_m,      only: vselector
  use stencil_m,        only: dirichlet_stencil, empty_stencil
  use error_m,          only: program_error
  use contact_m,        only: CT_SCHOTTKY

  implicit none

  private
  public :: schottky_injection_mb, schottky_velocity
  public :: schottky_injection_mb_bias, schottky_barrier_lowering
  public :: get_schottky_contact_normal_dir
  public :: schottky_injection, calc_schottky_injection

  type, extends(variable_real) :: schottky_injection
    !! Bias-dependent injection density n0B at Schottky contacts
    !! Only exists at Schottky contact vertices

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => schottky_injection_init
  end type

  type, extends(equation) :: calc_schottky_injection
    !! Calculate bias-dependent injection density n0B from electric field
    !! n0B = n0B_base * exp(delta_phi_b(E)) for electrons
    !! n0B = n0B_base * exp(-delta_phi_b(E)) for holes

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    integer :: ci
      !! carrier index

    type(vselector) :: efield_normal
      !! electric field component normal to contacts (combined selector)
    type(vselector) :: n0b
      !! injection density variable

    integer, allocatable :: contact_normal(:)
      !! normal direction for each contact (0 if not Schottky)

    type(dirichlet_stencil) :: st_dir
      !! stencil for contact vertices
    type(empty_stencil) :: st_em
      !! empty stencil

    type(jacobian), pointer :: jaco_efield => null()
      !! jacobian for normal E-field dependency
  contains
    procedure :: init => calc_schottky_injection_init
    procedure :: eval => calc_schottky_injection_eval
  end type

contains

  function get_schottky_contact_normal_dir(par, ict) result(normal_dir)
    !! Determine the normal direction of a Schottky contact for barrier lowering
    !! Analyzes vertex coordinates to find the constant dimension
    !! Returns: 1 for x-normal, 2 for y-normal, 3 for z-normal
    !! Returns: 0 if contact is not axis-aligned

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ict     ! Contact index
    integer                         :: normal_dir

    integer :: i, dir, idx(par%g%idx_dim)
    integer :: n_vertices
    real    :: coord_min(3), coord_max(3), coord_range(3)
    real    :: vertex_pos(3)
    real    :: tolerance, min_range, char_length

    ! Initialize
    normal_dir = 0
    coord_min = 1e30
    coord_max = -1e30
    tolerance = 1e-6  ! Relative tolerance for planarity

    ! Get number of transport vertices for this contact
    n_vertices = par%transport_vct(ict)%n

    if (n_vertices == 0) then
      print *, "WARNING: Schottky contact ", ict, " has no transport vertices"
      return
    end if

    ! Analyze coordinate ranges
    do i = 1, n_vertices
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, vertex_pos(1:par%g%dim))

      do dir = 1, par%g%dim
        coord_min(dir) = min(coord_min(dir), vertex_pos(dir))
        coord_max(dir) = max(coord_max(dir), vertex_pos(dir))
      end do
    end do

    ! Calculate range in each direction
    do dir = 1, par%g%dim
      coord_range(dir) = coord_max(dir) - coord_min(dir)
    end do

    ! Find direction with minimum range (normal to contact)
    min_range = huge(1.0)
    do dir = 1, par%g%dim
      if (coord_range(dir) < min_range) then
        min_range = coord_range(dir)
        normal_dir = dir
      end if
    end do

    ! Validate planarity
    if (par%g%dim > 1) then
      char_length = maxval(coord_range(1:par%g%dim))
      if (min_range > tolerance * char_length .and. char_length > 0.0) then
        print *, "WARNING: Schottky contact ", ict, " not perfectly axis-aligned for barrier lowering"
        print *, "  Range in normal direction: ", denorm(min_range, "nm"), " nm"
      end if
    end if

  end function

  subroutine schottky_injection_mb(par, ci, ict, ninj)
    !! Calculate equilibrium density n0 for Schottky contact Robin BC
    !! This is the thermionic emission density at zero bias
    !! Electrons: n0 = Nc * exp(-phi_Bn)
    !! Holes:     p0 = Nv * exp(-phi_Bp) where phi_Bp = E_g - phi_Bn

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(out) :: ninj    ! Equilibrium density n0 (normalized)

    real :: phi_Bn, phi_Bp

    ! Get normalized barrier height (already converted in device_params)
    phi_Bn = par%contacts(ict)%phi_b

    if (ci == CR_ELEC) then
      ! Electrons: n0 = Nc * exp(-phi_Bn)
      ninj = par%smc%edos(CR_ELEC) * exp(-phi_Bn)
      print "(A,ES12.5)", "DEBUG: ninj = ", denorm(ninj,"cm^-3")
    else  ! CR_HOLE
      ! Holes: barrier from valence band
      phi_Bp = par%smc%band_gap - phi_Bn
      ninj = par%smc%edos(CR_HOLE) * exp(-phi_Bp)
    end if
  end subroutine

  function schottky_velocity(par, ci, ict) result(s)
    !! Calculate thermionic emission velocity at Schottky contact
    !! Using Richardson constant: v_surf = A*T^2/(q*Nc)

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! Carrier index
    integer,             intent(in) :: ict  ! Contact index
    real                            :: s

    ! Check if Richardson constant is provided and > 0
    if (par%contacts(ict)%A_richardson > 0.0) then
      ! Calculate normalized surface velocity
      ! v_surf = A*T^2/(q*Nc) where q is handled by normalization
      ! T must be normalized, Nc is already normalized
      s = par%contacts(ict)%A_richardson * norm(par%T, "K") * norm(par%T, "K") / par%smc%edos(ci)
      print "(A,ES12.5)", "DEBUG: s = ", denorm(s,"cm/s")
    else
      ! Default thermal velocity estimate (v_th/4)
      s = 0.25  ! v_th/4 in normalized units
    end if
  end function

  subroutine schottky_injection_mb_bias(par, ci, ict, E_field, ninj, dninj_dE)
    !! Calculate bias-dependent equilibrium density and derivative for Schottky contact
    !! Includes image force barrier lowering effect

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(in)  :: E_field ! Electric field at interface (normalized)
    real,                intent(out) :: ninj    ! Injection density (normalized)
    real,                intent(out) :: dninj_dE ! Derivative d(ninj)/dE (normalized)

    real :: ninj_base, delta_phi_b, d_delta_phi_dE
    real :: eps_r = 11.7  ! Silicon relative permittivity (TODO: get from par)

    ! Get base injection without field
    call schottky_injection_mb(par, ci, ict, ninj_base)

    ! Calculate barrier lowering and its derivative
    call schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)

    ! Apply barrier lowering to injection density
    if (ci == CR_ELEC) then
      ! Electrons: lower barrier increases injection
      ninj = ninj_base * exp(delta_phi_b)
      dninj_dE = ninj * d_delta_phi_dE
    else  ! CR_HOLE
      ! Holes: opposite effect
      ninj = ninj_base * exp(-delta_phi_b)
      dninj_dE = -ninj * d_delta_phi_dE
    end if
  end subroutine

  subroutine schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)
    !! Calculate image force barrier lowering (Schottky effect)
    !! Implements Δφ_b = sqrt(q*|E|/(4π*ε))

    use math_m, only: PI

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ict           ! Contact index
    real,                intent(in)  :: E_field       ! Electric field (normalized)
    real,                intent(in)  :: eps_r         ! Relative permittivity
    real,                intent(out) :: delta_phi_b   ! Barrier lowering (normalized)
    real,                intent(out) :: d_delta_phi_dE ! Derivative d(Δφ_b)/dE

    real :: E_smooth, eps_smooth
    real :: gamma, q_norm
    integer, save :: debug_count = 0

    ! Check if IFBL is disabled for this contact
    if (.not. par%contacts(ict)%ifbl) then
      delta_phi_b = 0.0
      d_delta_phi_dE = 0.0
      return
    end if

    ! Smoothing parameter to avoid kink at E=0
    eps_smooth = 1e-10

    ! Smooth |E| to avoid numerical issues: |E| ≈ sqrt(E² + ε²)
    E_smooth = sqrt(E_field**2 + eps_smooth**2)

    ! Normalized charge (already in normalized units)
    q_norm = 1.0

    ! Image force coefficient: γ = sqrt(q/(4π*ε₀*εᵣ))
    ! In normalized units
    gamma = sqrt(q_norm / (4.0 * PI * eps_r))

    ! Barrier lowering: Δφ_b = γ*sqrt(|E|)
    delta_phi_b = gamma * sqrt(E_smooth)

    ! Derivative: d(Δφ_b)/dE = γ * 0.5 * E / (|E| * sqrt(|E|))
    ! Using smoothed version to avoid division by zero
    if (E_smooth > eps_smooth) then
      d_delta_phi_dE = gamma * 0.5 * E_field / (E_smooth * sqrt(E_smooth))
    else
      d_delta_phi_dE = 0.0
    end if

    ! Debug output for first few calls
    if (par%contacts(1)%type == 3) then  ! Only debug if first contact is Schottky
      if (debug_count < 3) then
        debug_count = debug_count + 1
        print *, "DEBUG schottky_barrier_lowering:"
        print *, "  E_field (normalized) = ", E_field
        print *, "  E_field (V/cm) = ", denorm(abs(E_field), "V/cm")
        print *, "  eps_r = ", eps_r
        print *, "  gamma = ", gamma
        print *, "  delta_phi_b (normalized) = ", delta_phi_b
        print *, "  delta_phi_b (eV) = ", denorm(delta_phi_b, "eV")
        print *, "  d_delta_phi_dE = ", d_delta_phi_dE
      end if
    end if
  end subroutine

  subroutine schottky_injection_init(this, par, ci)
    !! Initialize schottky_injection variable for storing n0B
    class(schottky_injection), intent(out) :: this
    type(device_params),       intent(in)  :: par
      !! device parameters
    integer,                   intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base - only exists at vertices
    call this%variable_init(CR_NAME(ci)//"n0b", "cm^-3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
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

  subroutine calc_schottky_injection_init(this, par, ci, efield, n0b)
    !! Initialize equation for calculating n0B from electric field
    class(calc_schottky_injection), intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    integer,                        intent(in)  :: ci
      !! carrier index
    type(electric_field), target,   intent(in)  :: efield(:)
      !! electric field components
    type(schottky_injection), target, intent(in) :: n0b
      !! injection density variable

    integer :: ict, i, iprov, idx(par%g%idx_dim), dir, idep
    integer :: normal_dir
    logical :: has_schottky
    integer, allocatable :: iprov_array(:), idep_array(:,:)

    print "(A)", "calc_schottky_injection_init for "//CR_NAME(ci)

    ! Check if we have any Schottky contacts
    has_schottky = .false.
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        has_schottky = .true.
        exit
      end if
    end do

    if (.not. has_schottky) then
      print *, "  No Schottky contacts found, skipping initialization"
      return
    end if

    ! init base equation
    call this%equation_init("calc_"//CR_NAME(ci)//"n0b")
    this%par => par
    this%ci = ci

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_em%init()

    ! Allocate and determine normal directions for each contact
    allocate(this%contact_normal(par%nct))
    this%contact_normal = 0  ! Initialize all to 0 (not Schottky or invalid)

    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        this%contact_normal(ict) = get_schottky_contact_normal_dir(par, ict)
        if (this%contact_normal(ict) > 0 .and. this%contact_normal(ict) <= par%g%dim) then
          print *, "  Schottky contact ", ict, " normal direction: ", this%contact_normal(ict)
        else
          print *, "  WARNING: Could not determine normal for Schottky contact ", ict
        end if
      end if
    end do

    ! Create n0b selector for all Schottky contact vertices
    ! We'll include all contacts but only Schottky ones will have non-zero values
    call this%n0b%init(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide n0B at Schottky contact vertices (provide for all, but only Schottky will be set)
    iprov = this%provide(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Create a combined selector for normal E-field components
    ! For each contact, we need the E-field component in its normal direction
    ! This is tricky - we need to combine different E-field components based on contact normal

    ! For simplicity, let's start with just the first valid normal direction
    ! TODO: This needs refinement for multiple contacts with different normals
    do dir = 1, par%g%dim
      if (any(this%contact_normal == dir)) then
        ! Use this direction's E-field for all Schottky contacts
        ! Simplified: use all contact vertices, framework will handle zero contributions
        call this%efield_normal%init(efield(dir), [(par%transport_vct(ict)%get_ptr(), &
                                                    ict = 1, par%nct)])

        ! Create dependency and Jacobian
        idep = this%depend(efield(dir), [(par%transport_vct(ict)%get_ptr(), &
                                         ict = 1, par%nct)])

        this%jaco_efield => this%init_jaco(iprov, idep, &
          st = [(this%st_dir%get_ptr(), ict = 1, par%nct)], &
          const = .false.)

        exit  ! For now, handle only one direction
      end if
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_schottky_injection_eval(this)
    !! Evaluate n0B = n0B_base * exp(±delta_phi_b(E))
    class(calc_schottky_injection), intent(inout) :: this

    integer :: ict, i, j, idx(this%par%g%idx_dim)
    real :: E_field, n0b_base, n0b_val, delta_phi_b, d_delta_phi_dE
    real :: eps_r
    real, allocatable :: tmp(:), E_field_arr(:)

    ! Check if properly initialized
    if (.not. allocated(this%contact_normal)) return
    if (.not. associated(this%jaco_efield)) return

    allocate(tmp(this%n0b%n))
    tmp = 0.0

    ! Counter for position in n0b array
    j = 0

    ! Loop over Schottky contacts
    do ict = 1, this%par%nct
      if (this%par%contacts(ict)%type /= CT_SCHOTTKY) cycle

      ! Skip if no valid normal direction
      if (this%contact_normal(ict) <= 0) cycle

      ! Get permittivity (TODO: make this material-dependent)
      eps_r = 11.7  ! Silicon

      ! Process each vertex of this contact
      do i = 1, this%par%transport_vct(ict)%n
        idx = this%par%transport_vct(ict)%get_idx(i)
        j = j + 1  ! Increment position in n0b array

        ! Get normal E-field at this vertex
        E_field_arr = this%efield_normal%get(idx)
        if (size(E_field_arr) > 0) then
          E_field = E_field_arr(1)  ! Extract scalar value
        else
          E_field = 0.0  ! No field if not properly set
        end if

        ! Get base injection density (without field)
        call schottky_injection_mb(this%par, this%ci, ict, n0b_base)

        ! Calculate barrier lowering
        call schottky_barrier_lowering(this%par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)

        ! Apply barrier lowering to get n0B
        if (this%ci == CR_ELEC) then
          n0b_val = n0b_base * exp(delta_phi_b)
          ! Set Jacobian entry
          call this%jaco_efield%set(idx, idx, n0b_val * d_delta_phi_dE)
        else  ! CR_HOLE
          n0b_val = n0b_base * exp(-delta_phi_b)
          ! Set Jacobian entry
          call this%jaco_efield%set(idx, idx, -n0b_val * d_delta_phi_dE)
        end if

        ! Store in temporary array at position j
        tmp(j) = n0b_val
      end do
    end do

    ! Materialize Jacobian
    call this%jaco_efield%set_matr(const = .false., nonconst = .true.)

    ! Set the n0B values
    call this%n0b%set(tmp)

    deallocate(tmp)
    if (allocated(E_field_arr)) deallocate(E_field_arr)
  end subroutine

end module
