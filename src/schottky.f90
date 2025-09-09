module schottky_m

  use device_params_m,  only: device_params
  use normalization_m,  only: norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE
  use grid_m,           only: IDX_VERTEX

  implicit none

  private
  public :: schottky_injection_mb, schottky_velocity
  public :: schottky_injection_mb_bias, schottky_barrier_lowering
  public :: get_schottky_contact_normal_dir

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

      ! Debug output
      print *, "DEBUG schottky_velocity for contact ", ict, ":"
      print *, "  A_richardson = ", par%contacts(ict)%A_richardson, " A/cm^2/K^2"
      print *, "  Temperature = ", par%T, " K"
      print *, "  Nc (normalized) = ", par%smc%edos(ci)
      print *, "  v_surf (normalized) = ", s
      print *, "  v_surf (physical) = ", denorm(s, "cm/s"), " cm/s"
    else
      ! Default thermal velocity estimate (v_th/4)
      s = 0.25  ! v_th/4 in normalized units
    end if
  end function

  subroutine schottky_injection_mb_bias(par, ci, ict, E_field, ninj, dninj_dE)
    !! Calculate bias-dependent equilibrium density and derivative for Schottky contact
    !! PLACEHOLDER: Will include image force barrier lowering effect
    
    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(in)  :: E_field ! Electric field at interface (normalized)
    real,                intent(out) :: ninj    ! Injection density (normalized)
    real,                intent(out) :: dninj_dE ! Derivative d(ninj)/dE (normalized)
    
    ! PLACEHOLDER: For now, just call the original function and set derivative to zero
    call schottky_injection_mb(par, ci, ict, ninj)
    dninj_dE = 0.0
    
    ! TODO: Implement barrier lowering physics
    ! - Calculate barrier lowering from E_field
    ! - Modify ninj based on lowered barrier
    ! - Calculate derivative for Jacobian
  end subroutine

  subroutine schottky_barrier_lowering(par, E_field, delta_phi_b, d_delta_phi_dE)
    !! Calculate image force barrier lowering (Schottky effect)
    !! PLACEHOLDER: Will implement Δφ_b = sqrt(q*|E|/(4π*ε))
    
    type(device_params), intent(in)  :: par
    real,                intent(in)  :: E_field       ! Electric field (normalized)
    real,                intent(out) :: delta_phi_b   ! Barrier lowering (normalized)
    real,                intent(out) :: d_delta_phi_dE ! Derivative d(Δφ_b)/dE
    
    ! PLACEHOLDER: No barrier lowering for now
    delta_phi_b = 0.0
    d_delta_phi_dE = 0.0
    
    ! TODO: Implement proper physics
    ! - Calculate barrier lowering: Δφ_b = sqrt(q*|E|/(4π*ε))
    ! - Calculate derivative for Jacobian
  end subroutine

end module
