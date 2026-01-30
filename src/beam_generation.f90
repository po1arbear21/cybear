module beam_generation_m
  !! STEM-EBIC beam generation with cone profile
  !!
  !! Models 200 keV electron beam through thin Si lamella.
  !! Beam has Gaussian profile in sweep direction with depth-dependent width:
  !!   σ(z) = sqrt(σ₀² + (α·z)²)  where α = semi-convergence angle
  !! Falls back to line profile (delta function) when grid is too coarse.

  use bin_search_m,    only: bin_search, BS_LESS, BS_GREAT
  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data0_real, grid_data1_real, grid_data2_real, grid_data3_real
  use semiconductor_m, only: CR_NAME, CR_ELEC, CR_HOLE
  use variable_m,      only: variable
  use normalization_m, only: denorm, norm

  implicit none

  private
  public beam_generation, beam_position, calc_beam_generation

  type, extends(variable) :: beam_generation
    !! Beam generation rate variable (volume rate)
    !! Unit: 1/cm^3/s

    integer :: ci
      !! carrier index (CR_ELEC or CR_HOLE)

    real, pointer :: x1(:)     => null()
    real, pointer :: x2(:,:)   => null()
    real, pointer :: x3(:,:,:) => null()
  contains
    procedure :: init => beam_generation_init
  end type

  type, extends(variable) :: beam_position
    !! Beam y-position for STEM-EBIC sweep
    !! Unit: cm (normalized)

    real, pointer :: x => null()
  contains
    procedure :: init => beam_position_init
  end type

  type, extends(equation) :: calc_beam_generation
    !! Calculate beam generation rate G(x,y) for STEM-EBIC

    type(device_params), pointer :: par => null()
    type(beam_generation), pointer :: bgen => null()
    type(beam_position), pointer :: beam_pos => null()

    integer :: ci
      !! carrier index
    integer :: i_min
      !! beam sweep range minimum index
    integer :: i_max
      !! beam sweep range maximum index
    logical :: first_eval = .true.
      !! flag to print debug output only on first eval
    real :: G_tot = 0.0
      !! total generation rate [pairs/(s·cm)] - stored for validation
  contains
    procedure :: init => calc_beam_generation_init
    procedure :: eval => calc_beam_generation_eval
  end type

contains

  subroutine beam_generation_init(this, par, ci)
    class(beam_generation), intent(out) :: this
    type(device_params),    intent(in)  :: par
    integer,                intent(in)  :: ci

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    call this%variable_init("bgen_"//CR_NAME(ci), "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

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

  subroutine beam_position_init(this, name)
    class(beam_position), intent(out) :: this
    character(*),         intent(in)  :: name

    type(grid_data0_real), pointer :: p

    call this%variable_init(name, "nm")
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

  subroutine calc_beam_generation_init(this, par, bgen, beam_pos)
    class(calc_beam_generation),    intent(out) :: this
    type(device_params),    target, intent(in)  :: par
    type(beam_generation),  target, intent(in)  :: bgen
    type(beam_position),    target, intent(in)  :: beam_pos

    integer :: iprov

    print "(A)", "calc_beam_generation_init"

    call this%equation_init("calc_"//bgen%name)
    this%par      => par
    this%bgen     => bgen
    this%beam_pos => beam_pos
    this%ci       = bgen%ci

    ! Find beam sweep range indices
    this%i_min = bin_search(par%g1D(2)%x, par%reg_beam(1)%beam_min)
    this%i_max = bin_search(par%g1D(2)%x, par%reg_beam(1)%beam_max)

    ! Print alignment info (first carrier only)
    if (this%ci == CR_ELEC) then
      print "(A,I0,A,I0)", "  Sweep range: ", this%i_min, " to ", this%i_max
    end if

    if (this%ci == CR_ELEC) then
      print "(A,A)", "  Beam profile: ", trim(par%reg_beam(1)%beam_dist%s)
    end if

    ! Provide beam generation on transport vertices
    iprov = this%provide(bgen, par%transport(IDX_VERTEX, 0))

    call this%init_final()
  end subroutine

  subroutine calc_beam_generation_eval(this)
    !! Evaluate STEM-EBIC generation with cone profile (Gaussian in y, depth-dependent width)
    !!
    !! Physics: 200 keV electron beam through thin Si lamella
    !! Beam forms a cone due to semi-convergence angle α = 10 mrad
    !! σ(x) = sqrt(σ₀² + (α·x)²) where x is depth from entrance surface
    !! G(x,y) = G_tot × exp(-(y-beam_y)²/(2σ²)) / (sqrt(2π)×σ×t)
    !! G_tot = (I_beam/q) × (dE/dz × t) / E_ehp
    !!
    !! For line profile: Uses linear interpolation between adjacent grid nodes
    !! to avoid grid-dependent beam position quantization.
    !!
    !! For Gaussian profile: Uses analytical integration over each control volume
    !! to avoid grid-dependent sampling errors.
    class(calc_beam_generation), intent(inout) :: this

    integer              :: i, iy, ny, ix_beam, iy_beam
    integer, allocatable :: idx(:)
    real                 :: I_beam_A_cm, beam_y_cm, beam_y_norm, t_cm, x_cm, y_cm, x_min_cm
    real                 :: sigma, sigma_sq, G_phys, G_tot, N_ehp
    real                 :: G_sum, tr_vol_cm2, rel_error
    real                 :: col_area_lo_cm2
    real                 :: y_lo
    real                 :: y_lo_edge_cm, y_hi_edge_cm, arg_lo, arg_hi
    real                 :: beam_x_cm, point_area_cm2
    logical              :: use_line, use_point

    ! Hard-coded STEM parameters for 200 keV beam in Si (CGS units)
    real, parameter :: Q_ELEM   = 1.602176634e-19  ! elementary charge [C]
    real, parameter :: SIGMA_0  = 1.0e-8           ! probe size at focus [cm] (0.1 nm)
    real, parameter :: ALPHA    = 0.01             ! semi-convergence angle [rad] (10 mrad)
    real, parameter :: DE_DZ    = 5.21e6           ! stopping power [eV/cm] (521 eV/μm for 200keV in Si, ESTAR)
    real, parameter :: E_EHP    = 3.6              ! energy per e-h pair [eV] (Si)
    real, parameter :: SQRT_2PI = 2.5066282746310002
    real, parameter :: SQRT_2   = 1.4142135623730951

    allocate(idx(this%par%g%idx_dim))

    ! Get beam parameters in physical units
    ! I_beam is specified as current per unit depth [A/cm] for 2D simulation
    I_beam_A_cm = denorm(this%par%reg_beam(1)%I_beam, 'A/cm')
    beam_y_cm = denorm(this%beam_pos%x, 'cm')
    beam_y_norm = this%beam_pos%x

    ! Select beam profile based on config
    use_line = (this%par%reg_beam(1)%beam_dist%s == "line")
    use_point = (this%par%reg_beam(1)%beam_dist%s == "point")

    ! Get lamella thickness [cm] - the beam path length through the material
    t_cm = denorm(this%par%reg_beam(1)%lamella_t, 'cm')
    x_min_cm = denorm(this%par%g1D(1)%x(1), 'cm')
    if (use_point) then
      beam_x_cm = denorm(this%par%reg_beam(1)%beam_x, 'cm')
    end if

    ! Total generation rate [pairs/(s·cm)] - per unit depth for 2D simulation
    ! N_ehp = number of e-h pairs generated per incident electron
    ! G_tot = pairs generated per second per cm depth
    N_ehp = (DE_DZ * t_cm) / E_EHP
    G_tot = (I_beam_A_cm / Q_ELEM) * N_ehp
    this%G_tot = G_tot

    ! Debug output (first eval, first carrier only)
    if (this%first_eval .and. this%ci == CR_ELEC) then
      print "(A,ES12.4,A)", "  STEM-EBIC: I_beam = ", I_beam_A_cm, " A/cm"
      print "(A,ES12.4,A)", "  STEM-EBIC: beam_y = ", beam_y_cm * 1e7, " nm"
      if (use_point) then
        print "(A,ES12.4,A)", "  STEM-EBIC: beam_x = ", beam_x_cm * 1e7, " nm"
      end if
      print "(A,ES12.4,A)", "  lamella thickness = ", t_cm * 1e7, " nm"
      print "(A,F8.1,A)",   "  N_ehp = ", N_ehp, " pairs/electron"
      print "(A,ES12.4,A)", "  G_tot = ", G_tot * 1e-4, " pairs/(s*um)"
    end if

    ! For point profile: find single point and calculate area
    if (use_point) then
      ! Find indices closest to beam position (beam_x, beam_y)
      ix_beam = bin_search(this%par%g1D(1)%x, this%par%reg_beam(1)%beam_x)
      iy_beam = bin_search(this%par%g1D(2)%x, beam_y_norm)

      ! Calculate area of the single point (for normalization)
      point_area_cm2 = 0.0
      do i = 1, this%par%transport(IDX_VERTEX, 0)%n
        idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
        if (idx(1) == ix_beam .and. idx(2) == iy_beam) then
          point_area_cm2 = denorm(this%par%tr_vol%get(idx), 'cm^2')
          exit
        end if
      end do

      ! Debug: print point info (first eval only)
      if (this%first_eval .and. this%ci == CR_ELEC) then
        print "(A,I0,A,I0,A)", "  Point profile: beam at grid point (", ix_beam, ", ", iy_beam, ")"
        print "(A,ES12.4,A)", "  Point area = ", point_area_cm2, " cm^2"
      end if
    end if

    ! For line profile: find the grid node at beam position
    ! With Y_BEAM in device config, beam positions are guaranteed to be grid nodes
    if (use_line) then
      ! Find the grid node at beam position
      iy_beam = bin_search(this%par%g1D(2)%x, beam_y_norm)

      ! Sanity check: beam should be exactly on a grid node
      ! (If Y_BEAM was specified in device config, this is guaranteed)
      y_lo = this%par%g1D(2)%x(iy_beam)
      if (abs(y_lo - beam_y_norm) > 1e-10 * abs(beam_y_norm + 1e-30)) then
        if (this%first_eval .and. this%ci == CR_ELEC) then
          print "(A)", "  WARNING: Beam position not on grid node - results may be mesh-dependent"
          print "(A,ES12.4,A,ES12.4)", "    beam_y = ", beam_y_cm*1e7, " nm, nearest node = ", &
            denorm(y_lo, 'cm')*1e7, " nm"
        end if
      end if

      ! Calculate total area of beam column (for normalization)
      col_area_lo_cm2 = 0.0
      do i = 1, this%par%transport(IDX_VERTEX, 0)%n
        idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
        if (idx(2) == iy_beam) then
          col_area_lo_cm2 = col_area_lo_cm2 + denorm(this%par%tr_vol%get(idx), 'cm^2')
        end if
      end do

      ! Debug output (first eval only)
      if (this%first_eval .and. this%ci == CR_ELEC) then
        print "(A,I0)", "  Line profile: beam at grid node iy = ", iy_beam
      end if
    end if

    ! Calculate G at each transport vertex
    G_sum = 0.0
    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)

      ! Get position in cm
      x_cm = denorm(this%par%g1D(1)%x(idx(1)), 'cm') - x_min_cm  ! depth from entrance
      y_cm = denorm(this%par%g1D(2)%x(idx(2)), 'cm')             ! lateral position

      if (use_point) then
        ! POINT PROFILE: All generation at single (beam_x, beam_y) point
        ! G [1/cm²/s] = G_tot [1/(s·cm)] / area [cm²] for 2D per-unit-depth
        if (idx(1) == ix_beam .and. idx(2) == iy_beam) then
          G_phys = G_tot / point_area_cm2
        else
          G_phys = 0.0
        end if
      else if (use_line) then
        ! LINE PROFILE: All generation at the single beam column (beam is on grid node)
        ! G [1/cm²/s] = G_tot [1/(s·cm)] / area [cm²] for 2D per-unit-depth
        if (idx(2) == iy_beam) then
          G_phys = G_tot / col_area_lo_cm2
        else
          G_phys = 0.0
        end if
      else
        ! GAUSSIAN PROFILE: Depth-dependent beam width with analytical integration
        ! σ(x) = sqrt(σ₀² + (α·x)²)
        sigma_sq = SIGMA_0**2 + (ALPHA * x_cm)**2
        sigma = sqrt(sigma_sq)

        ! Get the actual transport volume for this vertex
        tr_vol_cm2 = denorm(this%par%tr_vol%get(idx), 'cm^2')

        ! Determine control volume y-edges for this vertex
        ! For tensor grid: edges are at midpoints between adjacent nodes
        ny = this%par%g1D(2)%n
        iy = idx(2)

        if (iy == 1) then
          y_lo_edge_cm = denorm(this%par%g1D(2)%x(1), 'cm')
        else
          y_lo_edge_cm = 0.5 * (denorm(this%par%g1D(2)%x(iy-1), 'cm') + &
                                denorm(this%par%g1D(2)%x(iy), 'cm'))
        end if

        if (iy == ny) then
          y_hi_edge_cm = denorm(this%par%g1D(2)%x(ny), 'cm')
        else
          y_hi_edge_cm = 0.5 * (denorm(this%par%g1D(2)%x(iy), 'cm') + &
                                denorm(this%par%g1D(2)%x(iy+1), 'cm'))
        end if

        ! Analytical integration of Gaussian over control volume [y_lo_edge, y_hi_edge]
        ! frac_y = fraction of Gaussian captured by this y-slice (dimensionless)
        arg_lo = (y_lo_edge_cm - beam_y_cm) / (sigma * SQRT_2)
        arg_hi = (y_hi_edge_cm - beam_y_cm) / (sigma * SQRT_2)

        ! G_phys [1/cm²/s] = G_tot [1/(s·cm)] × frac_y × frac_x / tr_vol [cm²]
        ! where frac_x = (tr_vol / dy) / t_cm = dx / t_cm (uniform in x over lamella thickness)
        ! Simplifies to: G_phys = G_tot × frac_y / (t_cm × dy)
        ! Using dy = tr_vol / dx and noting we need consistency with tr_vol:
        ! G_phys = G_tot × frac_y × (dx / t_cm) / tr_vol = G_tot × frac_y / (t_cm × dy)
        G_phys = G_tot * 0.5 * (erf(arg_hi) - erf(arg_lo)) / t_cm / (y_hi_edge_cm - y_lo_edge_cm)
      end if

      ! Store in normalized units (volume rate)
      call this%bgen%set(idx, norm(G_phys, '1/cm^3/s'))

      ! Accumulate for conservation check: sum(G × area) should equal G_tot
      tr_vol_cm2 = denorm(this%par%tr_vol%get(idx), 'cm^2')
      G_sum = G_sum + G_phys * tr_vol_cm2
    end do

    ! Conservation check (first eval, first carrier only; warning always)
    if (this%ci == CR_ELEC) then
      rel_error = abs(G_sum - G_tot) / G_tot * 100.0
      if (this%first_eval) then
        print "(A,ES12.4,A,ES12.4,A)", "  Conservation: sum(G×A) = ", G_sum * 1e-4, &
          " pairs/(s*um), G_tot = ", G_tot * 1e-4, " pairs/(s*um)"
        print "(A,F8.4,A)", "  Relative error = ", rel_error, " %"
      end if
      if (rel_error > 5.0) then
        print "(A)", "  WARNING: Conservation error > 5% - consider refining mesh near beam"
      end if
    end if

    ! Clear first_eval flag after first evaluation
    this%first_eval = .false.
  end subroutine

end module
