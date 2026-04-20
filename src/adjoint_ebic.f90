module adjoint_ebic_m
  !! Adjoint-method EBIC kernel.
  !!
  !! Solves A^T phi = c once per bias point, then evaluates
  !! I_EBIC(r_beam) = phi^T s(r_beam) per beam position. Reuses existing
  !! Ramo-Shockley current DOFs (c as unit indicator) and the continuity
  !! beam-generation coupling (s via sys_full residual at u* with G switched on).
  !!
  !! Usage:
  !!   call ss%run(...)                          ! forward DC solve at dark state u*
  !!   call refactorize_at_converged(dev, ss%solver)
  !!   call build_collection_functional(dev, k, c)
  !!   call solve_adjoint(ss%solver, c, phi)
  !!   do i_beam = 1, n_beam
  !!     call build_adjoint_source(dev, r_beam(i_beam), s)
  !!     I_EBIC(i_beam) = evaluate_ebic(phi, s)
  !!   end do

  use bin_search_m,    only: bin_search
  use block_m,         only: block_real
  use device_m,        only: device
  use error_m,         only: program_error
  use normalization_m, only: norm, denorm
  use semiconductor_m, only: CR_ELEC, CR_HOLE
  use solver_base_m,   only: solver_real

  implicit none

  private
  public build_collection_functional
  public build_adjoint_source
  public refactorize_at_converged
  public solve_adjoint
  public evaluate_ebic
  public extract_phi_for_var
  public adjoint_fast_path
  public init_adjoint_fast_path
  public eval_I_EBIC_fast

  type adjoint_fast_path
    !! Precomputed state for the profile-aware fast evaluation path.
    !! Built once per (device, bias) at setup, reused across every beam position.
    logical :: enabled = .false.
    logical :: is_3D   = .false.
    logical :: has_elec, has_hole
    integer :: i_lo_elec = -1   !! sys_full DOF offset for ndens interior tab (0-indexed: +k-1)
    integer :: i_lo_hole = -1   !! sys_full DOF offset for pdens interior tab
    integer, allocatable :: rev_map(:,:,:)   !! (ix,iy,iz) -> k in transport_vct(0); 0 on contact
    real,    allocatable :: tr_vol_norm(:)   !! cached V_k (normalized), length = n_interior
  end type

contains

  subroutine build_collection_functional(dev, contact_id, c)
    !! c_k = unit indicator vector at the sys_full DOF index of dev%curr(k)%x.
    !!
    !! Uses the Ramo-Shockley terminal-current DOF already present in sys_full
    !! (the "currents" vselector main variable, one block with nct values in
    !! contact-id order). No hand-assembly of dI/du needed.
    class(device), target, intent(inout) :: dev
    integer,               intent(in)    :: contact_id
    real, allocatable,     intent(out)   :: c(:)

    integer :: iimvar, ibl, dof

    if (contact_id < 1 .or. contact_id > dev%par%nct) then
      call program_error("build_collection_functional: contact_id out of range")
    end if

    iimvar = dev%sys_full%search_main_var("currents")
    ibl    = dev%sys_full%res2block(iimvar)%d(1)

    if (dev%sys_full%i1(ibl) - dev%sys_full%i0(ibl) + 1 /= dev%par%nct) then
      call program_error("build_collection_functional: currents block size mismatch")
    end if

    allocate (c(dev%sys_full%n), source = 0.0)
    dof = dev%sys_full%i0(ibl) + contact_id - 1
    c(dof) = 1.0
  end subroutine

  subroutine refactorize_at_converged(dev, solver)
    !! Re-evaluate Jacobian at the current state and re-factorize the solver.
    !!
    !! The Newton loop in steady_state stores the LU at x_{final-1} (the state
    !! before the converged update). For strict adjoint correctness we want
    !! LU(J(u*)), so one extra eval + factorize pins the LU to the converged
    !! state. Negligible cost vs. the full Newton loop.
    class(device),       target, intent(inout) :: dev
    class(solver_real),          intent(inout) :: solver

    type(block_real), pointer :: dfdx, dfdx_precon

    call dev%sys_full%eval(df = dfdx, dfp = dfdx_precon)
    call solver%factorize(dfdx, P = dfdx_precon)
  end subroutine

  subroutine build_adjoint_source(dev, r_beam, s, x_beam, z_beam)
    !! Linearized source s = V * G(r; r_beam), scattered into the continuity
    !! rows of sys_full, zero elsewhere.
    !!
    !! Fast path: only the pieces that depend on beam_pos get re-evaluated.
    !!   1. Set dev%beam_pos%x = r_beam (y-coordinate). If x_beam is present,
    !!      also set dev%par%reg_beam(1)%beam_x for point-beam 2D sweeps.
    !!   2. Per active carrier: refresh G via calc_bgen%eval (local update
    !!      of dev%bgen(ci) — does NOT touch Poisson, Ramo-Shockley, etc.).
    !!   3. Per active carrier: tmp = jaco_bgen * bgen = -V * G on interior
    !!      vertex rows, 0 on contact-vertex rows.
    !!   4. Scatter -tmp (= +V*G) into s at the interior-tab DOF slice of
    !!      ndens/pdens in sys_full.
    !!
    !! Equivalent in output to the earlier implementation that called
    !! dev%sys_full%eval(f = f); s = -f  --- but 10x faster, because the
    !! Poisson residual, Ramo-Shockley residual, mobility/imref/current-
    !! density re-evaluations that sys_full%eval would trigger all yield
    !! zero at u* and contribute nothing to s. We skip them entirely.
    !!
    !! Assumes the esystem state is at u* (i.e. called right after Newton
    !! convergence, no set_x in between).
    class(device), target, intent(inout) :: dev
    real,                  intent(in)    :: r_beam
    real, allocatable,     intent(out)   :: s(:)
    real,        optional, intent(in)    :: x_beam
    real,        optional, intent(in)    :: z_beam

    integer              :: ci, iimvar, ibl_int, i_lo, i_hi, n_int
    integer              :: contin_f_n
    real,    allocatable :: tmp(:)
    character(len=5)     :: dens_name

    allocate (s(dev%sys_full%n), source = 0.0)

    dev%beam_pos%x = r_beam
    if (present(x_beam)) dev%par%reg_beam(1)%beam_x = x_beam
    if (present(z_beam)) dev%par%reg_beam(1)%beam_z = z_beam

    do ci = dev%par%ci0, dev%par%ci1
      ! Refresh G for this carrier at the new beam position (touches only bgen).
      call dev%calc_bgen(ci)%eval()

      ! tmp = jaco_bgen * bgen = -V * G  (jaco_bgen entries are -V at interior
      ! vertices, 0 at contact vertices). fact_y=0 zeros out any prior contents.
      contin_f_n = dev%contin(ci)%f%n
      if (allocated(tmp)) deallocate (tmp)
      allocate (tmp(contin_f_n), source = 0.0)
      call dev%contin(ci)%jaco_bgen%matr%mul_vec(            &
        &      dev%contin(ci)%bgen%get(), tmp, fact_y = 0.0)

      ! Locate the interior-tab DOF slice for this carrier's density in sys_full.
      if (ci == CR_ELEC) then
        dens_name = "ndens"
      else
        dens_name = "pdens"
      end if
      iimvar = dev%sys_full%search_main_var(trim(dens_name))
      ibl_int = dev%sys_full%res2block(iimvar)%d(1)
      i_lo    = dev%sys_full%i0(ibl_int)
      i_hi    = dev%sys_full%i1(ibl_int)
      n_int   = i_hi - i_lo + 1

      ! Scatter -tmp[1:n_int] (the V*G values) into s. Entries past n_int in tmp
      ! are zero (contact-vertex rows of jaco_bgen have no non-zero cols).
      s(i_lo:i_hi) = -tmp(1:n_int)
    end do

    if (allocated(tmp)) deallocate (tmp)
  end subroutine

  subroutine solve_adjoint(solver, c, phi)
    !! Transpose-solve A^T phi = c using the existing factorization.
    !! PARDISO/MUMPS both support this via a flag; no re-factorization.
    class(solver_real), intent(inout) :: solver
    real,               intent(in)    :: c(:)
    real, allocatable,  intent(out)   :: phi(:)

    allocate (phi(size(c)))
    call solver%solve(c, phi, trans = 'T')
  end subroutine

  function evaluate_ebic(phi, s) result(I_EBIC)
    !! I_EBIC = phi * s (scalar dot product in sys_full flat indexing).
    !! phi's Poisson rows are naturally skipped because s vanishes there.
    real, intent(in) :: phi(:), s(:)
    real             :: I_EBIC

    if (size(phi) /= size(s)) then
      call program_error("evaluate_ebic: phi and s size mismatch")
    end if
    I_EBIC = dot_product(phi, s)
  end function

  subroutine init_adjoint_fast_path(dev, fp)
    !! Precompute state for eval_I_EBIC_fast_line. Runs once per bias point.
    !!
    !! Builds:
    !!   - rev_map(ix,iy,iz) -> k in transport_vct(0), or 0 if that vertex is on
    !!     a contact (no continuity equation there, so no contribution to I_EBIC).
    !!   - tr_vol_norm(k): cached normalized control volume for interior vertex k.
    !!   - i_lo_elec / i_lo_hole: starting DOF offsets in sys_full for each
    !!     carrier's interior ndens/pdens tab. phi(i_lo + k - 1) is the adjoint
    !!     value at interior vertex k.
    class(device), target, intent(in)  :: dev
    type(adjoint_fast_path), intent(out) :: fp

    integer              :: k, iimvar, ibl, idx_dim
    integer, allocatable :: idx(:)
    integer              :: Nx, Ny, Nz

    idx_dim = dev%par%g%idx_dim
    fp%is_3D = (idx_dim == 3)
    fp%has_elec = (dev%par%ci0 <= CR_ELEC) .and. (dev%par%ci1 >= CR_ELEC)
    fp%has_hole = (dev%par%ci0 <= CR_HOLE) .and. (dev%par%ci1 >= CR_HOLE)

    Nx = dev%par%g1D(1)%n
    Ny = dev%par%g1D(2)%n
    if (fp%is_3D) then
      Nz = dev%par%g1D(3)%n
    else
      Nz = 1
    end if
    allocate (fp%rev_map(Nx, Ny, Nz), source = 0)
    allocate (fp%tr_vol_norm(dev%par%transport_vct(0)%n))
    allocate (idx(idx_dim))

    do k = 1, dev%par%transport_vct(0)%n
      idx = dev%par%transport_vct(0)%get_idx(k)
      if (fp%is_3D) then
        fp%rev_map(idx(1), idx(2), idx(3)) = k
      else
        fp%rev_map(idx(1), idx(2), 1) = k
      end if
      fp%tr_vol_norm(k) = dev%par%tr_vol%get(idx)
    end do

    if (fp%has_elec) then
      iimvar = dev%sys_full%search_main_var("ndens")
      ibl    = dev%sys_full%res2block(iimvar)%d(1)
      fp%i_lo_elec = dev%sys_full%i0(ibl)
    end if
    if (fp%has_hole) then
      iimvar = dev%sys_full%search_main_var("pdens")
      ibl    = dev%sys_full%res2block(iimvar)%d(1)
      fp%i_lo_hole = dev%sys_full%i0(ibl)
    end if

    fp%enabled = .true.
  end subroutine

  subroutine eval_I_EBIC_fast(dev, fp, phi, r_beam, I_EBIC, x_beam, z_beam)
    !! Profile-aware fast evaluation of I_EBIC. Handles both "line" and "point"
    !! beam profiles. Falls through to program_error for anything else.
    !!
    !! Math: given the set S of interior vertices where the beam deposits
    !! nonzero generation,
    !!     I_EBIC = sum_{k in S} V_k * G_k * phi(dof_k)
    !!           = G_tot / V_total * sum_{k in S} V_k * phi(dof_k)
    !! where V_total = sum_{k in S} V_k_phys (col_vol for line, point_vol for
    !! point). Skips calc_bgen, jaco_bgen mat-vec, scatter, and the full
    !! sys_full dot product entirely — O(|S|) operations per position instead
    !! of O(N_dof).
    !!
    !! S is:
    !!   line 2D:   {(ix, iy_beam)        : ix = 1..Nx}              (~Nx terms)
    !!   line 3D:   {(ix, iy_beam, iz_beam): ix = 1..Nx}             (~Nx terms)
    !!   point 2D:  {(ix_beam, iy_beam)}                             (1 term)
    !!   point 3D:  {(ix_beam, iy_beam, iz_beam)}                    (1 term)
    class(device),           target, intent(inout) :: dev
    type(adjoint_fast_path),         intent(in)    :: fp
    real,                            intent(in)    :: phi(:)
    real,                            intent(in)    :: r_beam     !! y beam position (normalized)
    real,                            intent(out)   :: I_EBIC     !! result (normalized units)
    real,                  optional, intent(in)    :: x_beam     !! x beam position (normalized), point profile
    real,                  optional, intent(in)    :: z_beam     !! z beam position (normalized), 3D only

    ! Physics constants (must match beam_generation.f90)
    real, parameter :: Q_ELEM = 1.602176634e-19
    real, parameter :: DE_DZ  = 5.21e6
    real, parameter :: E_EHP  = 3.6

    integer :: ix, iy, iz, k, dof, ix_lo, ix_hi
    real    :: t_cm, I_beam_phys, N_ehp, G_tot, vol_total_phys, V_k_phys
    real    :: G_phys, G_norm, wsum
    logical :: is_line, is_point

    if (.not. fp%enabled) &
      & call program_error("eval_I_EBIC_fast: fast path not initialized")

    is_line  = (dev%par%reg_beam(1)%beam_dist%s == "line")
    is_point = (dev%par%reg_beam(1)%beam_dist%s == "point")
    if (.not. (is_line .or. is_point)) &
      & call program_error("eval_I_EBIC_fast: only 'line' and 'point' profiles supported")

    ! Update the device-level beam coords so downstream printouts stay consistent
    ! and bin_search on the current values works.
    dev%beam_pos%x = r_beam
    if (present(x_beam)) dev%par%reg_beam(1)%beam_x = x_beam
    if (present(z_beam)) dev%par%reg_beam(1)%beam_z = z_beam

    ! Snap y to grid (always needed; mirrors beam_generation.f90)
    iy = bin_search(dev%par%g1D(2)%x, r_beam)

    ! Snap z for 3D; irrelevant for 2D
    if (fp%is_3D) then
      if (dev%par%reg_beam(1)%beam_z >= 0.0) then
        iz = bin_search(dev%par%g1D(3)%x, dev%par%reg_beam(1)%beam_z)
      else
        iz = (dev%par%g1D(3)%n + 1) / 2
      end if
    else
      iz = 1
    end if

    ! x iteration range: all x for line, single x for point
    if (is_line) then
      ix_lo = 1
      ix_hi = dev%par%g1D(1)%n
    else
      ix_lo = bin_search(dev%par%g1D(1)%x, dev%par%reg_beam(1)%beam_x)
      ix_hi = ix_lo
    end if

    ! Pass 1: vol_total_phys = sum of control volumes of beam-active vertices
    vol_total_phys = 0.0
    do ix = ix_lo, ix_hi
      k = fp%rev_map(ix, iy, iz)
      if (k == 0) cycle   ! vertex on a contact — skip
      if (fp%is_3D) then
        V_k_phys = denorm(fp%tr_vol_norm(k), 'cm^3')
      else
        V_k_phys = denorm(fp%tr_vol_norm(k), 'cm^2')
      end if
      vol_total_phys = vol_total_phys + V_k_phys
    end do

    if (vol_total_phys <= 0.0) then
      I_EBIC = 0.0
      return
    end if

    ! G_tot from beam params (matches beam_generation.f90:209-210)
    t_cm  = denorm(dev%par%reg_beam(1)%lamella_t, 'cm')
    if (fp%is_3D) then
      I_beam_phys = denorm(dev%par%reg_beam(1)%I_beam, 'A')
    else
      I_beam_phys = denorm(dev%par%reg_beam(1)%I_beam, 'A/cm')
    end if
    N_ehp = DE_DZ * t_cm / E_EHP
    G_tot = (I_beam_phys / Q_ELEM) * N_ehp

    ! Mirror the slow path's side effect: calc_bgen%eval stores G_tot on each
    ! calc_bgen for the driver's eta_absolute reporting. Without this the value
    ! left over from the dark solve (I_beam=0 -> G_tot=0) would divide-by-zero.
    if (fp%has_elec) dev%calc_bgen(CR_ELEC)%G_tot = G_tot
    if (fp%has_hole) dev%calc_bgen(CR_HOLE)%G_tot = G_tot

    ! Uniform G per active vertex (normalized units, matching bgen storage)
    G_phys = G_tot / vol_total_phys
    G_norm = norm(G_phys, '1/cm^3/s')

    ! Pass 2: I_EBIC = sum_{ci} sum_{k active} V_k_norm * G_norm * phi(dof_k)
    wsum = 0.0
    do ix = ix_lo, ix_hi
      k = fp%rev_map(ix, iy, iz)
      if (k == 0) cycle
      if (fp%has_elec) then
        dof = fp%i_lo_elec + k - 1
        wsum = wsum + phi(dof) * fp%tr_vol_norm(k)
      end if
      if (fp%has_hole) then
        dof = fp%i_lo_hole + k - 1
        wsum = wsum + phi(dof) * fp%tr_vol_norm(k)
      end if
    end do

    I_EBIC = wsum * G_norm
  end subroutine

  subroutine extract_phi_for_var(dev, phi, var_name, tab_idx, phi_slice)
    !! Extract the portion of phi corresponding to main variable var_name,
    !! tab tab_idx. For "ndens"/"pdens" the tabs are:
    !!   tab 1       = interior vertices (transport_vct(0))
    !!   tab 2..nct+1 = contact vertices (transport_vct(1..nct))
    !! The interior slice is what's needed for the collection-probability map.
    class(device), target, intent(inout) :: dev
    real,                  intent(in)    :: phi(:)
    character(*),          intent(in)    :: var_name
    integer,               intent(in)    :: tab_idx
    real, allocatable,     intent(out)   :: phi_slice(:)

    integer :: iimvar, ibl

    iimvar = dev%sys_full%search_main_var(var_name)
    if (tab_idx < 1 .or. tab_idx > size(dev%sys_full%res2block(iimvar)%d)) then
      call program_error("extract_phi_for_var: tab_idx out of range")
    end if

    ibl      = dev%sys_full%res2block(iimvar)%d(tab_idx)
    phi_slice = phi(dev%sys_full%i0(ibl):dev%sys_full%i1(ibl))
  end subroutine

end module
