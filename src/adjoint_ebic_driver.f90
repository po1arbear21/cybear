program adjoint_ebic_driver
  !! Adjoint-method EBIC driver.
  !!
  !! Workflow per [adjoint ebic] section:
  !!   1. Dark DC solve (V=0 on all contacts, I_beam=0, G=0) → u*.
  !!   2. Re-evaluate J at u* and re-factorize (pin LU to the converged state).
  !!   3. Build c = unit indicator at the specified collecting contact.
  !!   4. Transpose-solve J^T phi = c.
  !!   5. Restore physical I_beam; for each r_beam in the list:
  !!        form s(r_beam) via sys_full%eval, evaluate I_EBIC = phi . s.
  !!   6. Dump phi^n + phi^p interior slices for the collection-probability map.
  !!
  !! Use via: fargo run adjoint_ebic --device D --run R

  use adjoint_ebic_m,  only: build_collection_functional, build_adjoint_source, &
                           & refactorize_at_converged, solve_adjoint,            &
                           & evaluate_ebic, extract_phi_for_var,                 &
                           & adjoint_fast_path, init_adjoint_fast_path,          &
                           & eval_I_EBIC_fast
  use approx_m,        only: approx_imref, approx_potential
  use cl_options_m,    only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,        only: dev
  use error_m,         only: program_error
  use input_m,         only: input_file
  use input_src_m,     only: polygon_src
  use math_m,          only: linspace
  use normalization_m, only: init_normconst, denorm
  use semiconductor_m, only: CR_NAME, CR_ELEC, CR_HOLE
  use steady_state_m,  only: steady_state
  use storage_m,       only: storage, STORAGE_WRITE
  use string_m,        only: string
  use util_m,          only: get_hostname
  use vselector_m,     only: vselector

  implicit none

  logical                   :: gummel_restart, gummel_once, gummel_enabled
  real                      :: temperature
  type(input_file)          :: runfile
  character(:), allocatable :: device_file, run_file
  integer                   :: k
  integer,      allocatable :: sids(:)

  print "(A)", "=============================================="
  print "(A)", " adjoint_ebic — collection-map via reciprocity"
  print "(A)", "=============================================="
  print "(A)", "Host: " // get_hostname()

  call command_line()

  if (.not. allocated(device_file)) call program_error("--device is required")
  if (.not. allocated(run_file))    call program_error("--run is required")

  call dev%init(device_file, temperature)
  call runfile%init(run_file)

  call runfile%get_sections("adjoint ebic", sids)
  if (size(sids) == 0) call program_error("no [adjoint ebic] sections found in run file")

  do k = 1, size(sids)
    call run_adjoint_ebic(sids(k))
  end do

contains

  subroutine command_line()
    !! Parse --device, --run, --temperature.
    integer                      :: idesc, i, j
    integer,         allocatable :: iclopt(:), jclopt(:)
    type(cl_option), allocatable :: clopt(:)
    type(cl_option_descriptor)   :: desc(3) = [ &
      cl_option_descriptor('T', "temperature", .true.,  .false., .true., .true.), &
      cl_option_descriptor('d', "device",      .true.,  .false., .true., .true.), &
      cl_option_descriptor('r', "run",         .true.,  .false., .true., .true.)  &
    ]

    call get_cl_options(desc, clopt, iclopt, jclopt)
    do idesc = 1, size(desc)
      do i = iclopt(idesc), iclopt(idesc + 1) - 1
        j = jclopt(i)
        select case (clopt(j)%short)
        case ('T')
          read (clopt(j)%arg, *) temperature
          call init_normconst(temperature)
          print "(A,ES25.16E3,A)", "T = ", temperature, " K"
        case ('d')
          device_file = clopt(j)%arg
        case ('r')
          run_file = clopt(j)%arg
        end select
      end do
    end do
  end subroutine

  subroutine run_adjoint_ebic(si)
    !! Execute one [adjoint ebic] section. Expected keys:
    !!   name               — run label
    !!   collecting_contact — name of the contact whose current we solve for
    !!   r_beam             — list of beam y-positions (in length units of the device grid)
    integer, intent(in) :: si

    type(steady_state)         :: ss
    type(string)               :: suite_name, contact_name
    integer                    :: sj, ict_collect, nk, k, n_n_int, n_p_int, i, ict
    real                       :: I_beam_saved
    real,          allocatable :: r_beam_list(:), x_beam_list(:), z_beam_list(:), I_EBIC(:)
    real,          allocatable :: c(:), phi(:), s(:)
    real,          allocatable :: phi_n_int(:), phi_p_int(:)
    real,          allocatable :: V_contact(:), I_dark(:)
    real,          allocatable :: u_dark_state(:), delta_u_primal(:), s_primal_save(:)
    logical                    :: have_holes, sweep_has_x, sweep_has_z
    logical                    :: status, has_V, any_biased

    ! Absolute collection efficiency accounting
    real, parameter            :: Q_ELEM_C = 1.602176634e-19   ! Coulombs
    real                       :: G_tot_per_cm, curr_fact_cm, total_pairs_per_sec
    real                       :: I_EBIC_peak_A, I_EBIC_A, eta_abs, eta_rel

    ! Timing
    integer(8)                 :: t_setup_start, t_setup_end, t_sweep_start, t_sweep_end
    integer(8)                 :: clock_rate
    real                       :: setup_sec, sweep_sec

    call runfile%get(si, "name", suite_name)
    print "(A)", ""
    print "(A)", "=== Adjoint EBIC case: " // suite_name%s // " ==="

    call runfile%get(si, "collecting_contact", contact_name)

    ! Beam-position specification. Supports four modes by key presence:
    !   (a) r_beam / r_beam_min/_max/_n       -> 1D y-sweep (legacy)
    !   (b) x_beam_*, y_beam_*                -> 2D (x, y) Cartesian for 2D point
    !   (c) y_beam_*, z_beam_*                -> 2D (y, z) Cartesian for 3D line
    !   (d) x_beam_*, y_beam_*, z_beam_*      -> 3D (x, y, z) Cartesian for 3D point
    ! Flat ordering: outermost loop over y, then z, then x inner-most.
    block
      real    :: r_beam_min, r_beam_max
      real    :: x_beam_min, x_beam_max, y_beam_min, y_beam_max
      real    :: z_beam_min, z_beam_max
      integer :: r_beam_n, x_beam_n, y_beam_n, z_beam_n, i_x, i_y, i_z, k_idx
      integer :: n_axis_x, n_axis_y, n_axis_z
      logical :: has_min, has_max, has_n
      logical :: hx_min, hx_max, hx_n, hy_min, hy_max, hy_n
      logical :: hz_min, hz_max, hz_n
      real, allocatable :: x_axis(:), y_axis(:), z_axis(:)

      call runfile%get(si, "x_beam_min", x_beam_min, status = hx_min)
      call runfile%get(si, "x_beam_max", x_beam_max, status = hx_max)
      call runfile%get(si, "x_beam_n",   x_beam_n,   status = hx_n)
      call runfile%get(si, "y_beam_min", y_beam_min, status = hy_min)
      call runfile%get(si, "y_beam_max", y_beam_max, status = hy_max)
      call runfile%get(si, "y_beam_n",   y_beam_n,   status = hy_n)
      call runfile%get(si, "z_beam_min", z_beam_min, status = hz_min)
      call runfile%get(si, "z_beam_max", z_beam_max, status = hz_max)
      call runfile%get(si, "z_beam_n",   z_beam_n,   status = hz_n)
      sweep_has_x = hx_min .and. hx_max .and. hx_n
      sweep_has_z = hz_min .and. hz_max .and. hz_n

      if (sweep_has_x .or. sweep_has_z) then
        if (.not. (hy_min .and. hy_max .and. hy_n)) &
          & call program_error("sweep requires y_beam_min/_max/_n when x_beam_* or z_beam_* is present")

        if (sweep_has_x) then
          x_axis = linspace(x_beam_min, x_beam_max, x_beam_n)
          n_axis_x = x_beam_n
        else
          allocate (x_axis(1))
          x_axis = 0.0
          n_axis_x = 1
        end if
        y_axis = linspace(y_beam_min, y_beam_max, y_beam_n)
        n_axis_y = y_beam_n
        if (sweep_has_z) then
          z_axis = linspace(z_beam_min, z_beam_max, z_beam_n)
          n_axis_z = z_beam_n
        else
          allocate (z_axis(1))
          z_axis = 0.0
          n_axis_z = 1
        end if

        allocate (r_beam_list(n_axis_x * n_axis_y * n_axis_z))
        if (sweep_has_x) allocate (x_beam_list(n_axis_x * n_axis_y * n_axis_z))
        if (sweep_has_z) allocate (z_beam_list(n_axis_x * n_axis_y * n_axis_z))
        k_idx = 0
        do i_y = 1, n_axis_y
          do i_z = 1, n_axis_z
            do i_x = 1, n_axis_x
              k_idx = k_idx + 1
              r_beam_list(k_idx) = y_axis(i_y)
              if (sweep_has_x) x_beam_list(k_idx) = x_axis(i_x)
              if (sweep_has_z) z_beam_list(k_idx) = z_axis(i_z)
            end do
          end do
        end do
      else
        call runfile%get(si, "r_beam_min", r_beam_min, status = has_min)
        call runfile%get(si, "r_beam_max", r_beam_max, status = has_max)
        call runfile%get(si, "r_beam_n",   r_beam_n,   status = has_n)
        if (has_min .and. has_max .and. has_n) then
          r_beam_list = linspace(r_beam_min, r_beam_max, r_beam_n)
        else
          call runfile%get(si, "r_beam", r_beam_list)
        end if
      end if
      if (.not. allocated(x_beam_list)) allocate (x_beam_list(0))
      if (.not. allocated(z_beam_list)) allocate (z_beam_list(0))
    end block

    call dev%par%contact_map%get(contact_name, ict_collect, status)
    if (.not. status) call program_error("contact not found: " // contact_name%s)
    print "(A,A,A,I0,A)", "  Collecting contact: ", trim(contact_name%s), &
                       &  " (index ", ict_collect, ")"
    if (sweep_has_x .and. sweep_has_z) then
      print "(A)",    "  Sweep mode: 3D Cartesian (x, y, z) product"
    else if (sweep_has_x) then
      print "(A)",    "  Sweep mode: 2D Cartesian (x, y) product"
    else if (sweep_has_z) then
      print "(A)",    "  Sweep mode: 2D Cartesian (y, z) product"
    else
      print "(A)",    "  Sweep mode: 1D (y only)"
    end if
    print "(A,I0)", "  Number of beam positions: ", size(r_beam_list)

    ! Per-contact bias: optional keys V_<contact_name> (unit ": V"). Missing keys
    ! default to 0 V (short-circuit). Non-zero entries enable "operando" EBIC: the
    ! adjoint linearizes around the biased DC operating point, and the forward
    ! validation subtracts the dark contact current before comparing.
    allocate (V_contact(dev%par%nct))
    allocate (I_dark(dev%par%nct))
    any_biased = .false.
    do ict = 1, dev%par%nct
      call runfile%get(si, "V_" // trim(dev%par%contacts(ict)%name), V_contact(ict), &
        & status = has_V)
      if (.not. has_V) V_contact(ict) = 0.0
      if (abs(V_contact(ict)) > 0.0) any_biased = .true.
    end do
    if (any_biased) then
      print "(A)", "  Contact bias (operando — adjoint linearizes around the biased DC state):"
      do ict = 1, dev%par%nct
        print "(A,A,A,ES12.4,A)", "    V_", trim(dev%par%contacts(ict)%name), " = ", &
          & denorm(V_contact(ict), 'V'), " V"
      end do
    end if

    call runfile%get_section("full newton params", sj)
    call ss%init(dev%sys_full)
    call ss%set_params(runfile%sections%d(sj))
    ss%msg = "Newton: "

    call system_clock(count_rate = clock_rate)
    call system_clock(t_setup_start)

    I_beam_saved = dev%par%reg_beam(1)%I_beam
    ! Pass r_beam_list(1) as the dark-solve beam_y so that u_dark and u_fwd
    ! agree on the Y_BEAM input-DOF value. I_beam = 0 during dark, so the
    ! physics is unaffected (no generation regardless of where we park the
    ! beam), but norm(u_fwd - u_dark) no longer picks up the 7000 nm input
    ! difference that would otherwise dominate the diagnostic norm.
    call run_dark_solve(ss, V_contact, I_dark, beam_y_park = r_beam_list(1))

    ! Capture flat u_dark for the Jacobian-linearization residual diagnostic
    ! (compared later against u_fwd - u_dark from the forward linescan).
    u_dark_state = dev%sys_full%get_x()

    print "(A)", ""
    print "(A)", "--- Refactorize J at converged u* ---"
    call refactorize_at_converged(dev, ss%solver)

    print "(A)", ""
    print "(A)", "--- Build c, solve adjoint J^T phi = c ---"
    call build_collection_functional(dev, ict_collect, c)
    call solve_adjoint(ss%solver, c, phi)
    print "(A,I0)", "  adjoint DOFs: ", size(phi)

    ! Restore physical I_beam; G enters s linearly so sign of the beam matters
    ! but the magnitude cleanly multiplies the final I_EBIC.
    dev%par%reg_beam(1)%I_beam = I_beam_saved

    nk = size(r_beam_list)
    allocate (I_EBIC(nk))

    print "(A)", ""
    print "(A)", "--- Beam sweep via adjoint dot-product ---"
    ! Profile-aware fast path for "line" and "point" beams: O(N_active_vertices)
    ! per position instead of O(N_dof). Bypasses calc_bgen, jaco_bgen mat-vec,
    ! scatter, and dense dot_product. Falls back to the general path for
    ! "gaussian" or any future profile.
    block
      type(adjoint_fast_path) :: fp
      logical                 :: use_fast, use_fast_requested, has_fast_key
      character(:), allocatable :: dist_s

      ! Read optional override from run INI. Default: auto (use fast when profile
      ! is line/point, slow otherwise). Set `use_fast_path = false` to force the
      ! general mat-vec path — useful as a correctness fallback or for debugging.
      call runfile%get(si, "use_fast_path", use_fast_requested, status = has_fast_key)
      if (.not. has_fast_key) use_fast_requested = .true.

      dist_s = dev%par%reg_beam(1)%beam_dist%s
      use_fast = use_fast_requested .and. (dist_s == "line" .or. dist_s == "point")

      if (use_fast) then
        print "(A,A,A)", "  Sweep kernel: profile-aware fast path (",                 &
                       & trim(dist_s), " profile, O(N_active_vertices) per position)"
        call init_adjoint_fast_path(dev, fp)
      else if (.not. use_fast_requested) then
        print "(A,A,A)", "  Sweep kernel: general slow path (use_fast_path=false in ",&
                       & "run INI, O(N_dof) per position)"
      else
        print "(A,A,A)",  "  Sweep kernel: general slow path (profile '",             &
                       & trim(dist_s), "' not handled by fast path, O(N_dof) per position)"
      end if

      ! Close setup timing AFTER fast-path init so setup_time includes all
      ! one-time amortized work. Everything between here and t_sweep_start is
      ! per-run negligible (print statements).
      call system_clock(t_setup_end)
      call system_clock(t_sweep_start)
      do k = 1, nk
        if (use_fast) then
          if (sweep_has_x .and. sweep_has_z) then
            call eval_I_EBIC_fast(dev, fp, phi, r_beam_list(k), I_EBIC(k), &
              & x_beam = x_beam_list(k), z_beam = z_beam_list(k))
          else if (sweep_has_x) then
            call eval_I_EBIC_fast(dev, fp, phi, r_beam_list(k), I_EBIC(k), &
              & x_beam = x_beam_list(k))
          else if (sweep_has_z) then
            call eval_I_EBIC_fast(dev, fp, phi, r_beam_list(k), I_EBIC(k), &
              & z_beam = z_beam_list(k))
          else
            call eval_I_EBIC_fast(dev, fp, phi, r_beam_list(k), I_EBIC(k))
          end if
        else
          if (sweep_has_x .and. sweep_has_z) then
            call build_adjoint_source(dev, r_beam_list(k), s, &
              & x_beam = x_beam_list(k), z_beam = z_beam_list(k))
          else if (sweep_has_x) then
            call build_adjoint_source(dev, r_beam_list(k), s, x_beam = x_beam_list(k))
          else if (sweep_has_z) then
            call build_adjoint_source(dev, r_beam_list(k), s, z_beam = z_beam_list(k))
          else
            call build_adjoint_source(dev, r_beam_list(k), s)
          end if
          I_EBIC(k) = evaluate_ebic(phi, s)
        end if
      end do
      call system_clock(t_sweep_end)
    end block

    call run_primal_diagnostic(si, ss, c, phi, r_beam_list, x_beam_list, z_beam_list, &
                             & sweep_has_x, sweep_has_z, u_dark_state,                &
                             & delta_u_primal, s_primal_save)

    ! Total pair generation rate — needed for absolute collection efficiency.
    ! G_tot (stored on calc_bgen during eval) is in pairs/(s·cm) for 2D; multiply
    ! by the implicit simulation depth curr_fact to get total pairs/s.
    G_tot_per_cm       = dev%calc_bgen(CR_ELEC)%G_tot
    curr_fact_cm       = denorm(dev%par%curr_fact, 'cm')
    total_pairs_per_sec = G_tot_per_cm * curr_fact_cm
    I_EBIC_peak_A       = maxval(abs([ (denorm(I_EBIC(k), 'A'), k = 1, nk) ]))

    print "(A,ES12.4,A)", "  Total pair generation rate: ", total_pairs_per_sec, " pairs/s"
    print "(A,ES12.4,A)", "  q * G_total              : ", Q_ELEM_C * total_pairs_per_sec, &
                        &                                  " A  (= I_EBIC if eta=1)"
    print "(A)", ""
    if (sweep_has_x .or. sweep_has_z) then
      print "(A,ES16.6,A)", "  Multi-D sweep peak |I_EBIC|: ", I_EBIC_peak_A, " A"
      print "(A,ES16.6)",   "  Multi-D sweep peak eta_abs : ", &
                          & I_EBIC_peak_A / (Q_ELEM_C * total_pairs_per_sec)
      print "(A)",          "  (Per-position values written to fbs; table skipped for multi-D sweep.)"
    else
      print "(A)", "  ┌──────────┬──────────────────┬────────────────┬────────────────┐"
      print "(A)", "  │ y [nm]   │   I_EBIC [A]     │  eta_absolute  │  eta / peak    │"
      print "(A)", "  ├──────────┼──────────────────┼────────────────┼────────────────┤"
      do k = 1, nk
        I_EBIC_A = denorm(I_EBIC(k), 'A')
        eta_abs  = I_EBIC_A / (Q_ELEM_C * total_pairs_per_sec)
        eta_rel  = I_EBIC_A / I_EBIC_peak_A
        print "(A,F8.1,A,ES16.6,A,F14.6,A,F14.6,A)",                                  &
          "  │ ", denorm(r_beam_list(k), 'nm'),                                       &
          " │ ", I_EBIC_A,                                                            &
          " │ ", eta_abs,                                                             &
          " │ ", eta_rel, " │"
      end do
      print "(A)", "  └──────────┴──────────────────┴────────────────┴────────────────┘"
    end if

    setup_sec = real(t_setup_end - t_setup_start) / real(clock_rate)
    sweep_sec = real(t_sweep_end - t_sweep_start) / real(clock_rate)
    print "(A)", ""
    print "(A)", "--- Timing (demonstrates amortization) ---"
    print "(A,F10.4,A)",    "  Setup (Newton + refactorize + adjoint solve):  ", setup_sec, " s"
    print "(A,I0,A,F10.4,A)","  Sweep across ", nk, " positions (eval + dot): ", sweep_sec, " s"
    print "(A,ES12.4,A)",   "  Per-position cost (eval + dot):                ", sweep_sec / real(nk), " s"
    print "(A,F8.1,A)",     "  Per-position / per-setup ratio:                ", sweep_sec / real(nk) / setup_sec * 100.0, " %"

    ! Collection-probability map: interior vertex slice of phi for ndens and pdens.
    have_holes = (dev%par%ci1 >= CR_HOLE)
    print "(A)", ""
    print "(A)", "--- Collection probability phi^n (+ phi^p) on interior vertices ---"

    call extract_phi_for_var(dev, phi, "ndens", 1, phi_n_int)
    n_n_int = size(phi_n_int)
    print "(A,I0,A,ES14.6,A,ES14.6)", "  ndens interior DOFs: ", n_n_int, &
                                 &     "   min(phi^n) = ", minval(phi_n_int), &
                                 &     "   max(phi^n) = ", maxval(phi_n_int)
    print "(A,ES14.6)",            "  |phi^n| = ", norm2(phi_n_int)
    print "(A)", ""
    print "(A)", "  phi^n stride-through (first 30 entries):"
    do i = 1, min(30, n_n_int)
      print "(A,I3,A,ES14.6)", "    phi^n[", i, "] = ", phi_n_int(i)
    end do

    if (have_holes) then
      call extract_phi_for_var(dev, phi, "pdens", 1, phi_p_int)
      n_p_int = size(phi_p_int)
      print "(A,I0,A,ES14.6,A,ES14.6)", "  pdens interior DOFs: ", n_p_int, &
                                   &     "   min(phi^p) = ", minval(phi_p_int), &
                                   &     "   max(phi^p) = ", maxval(phi_p_int)
      print "(A,ES14.6)",              "  |phi^p| = ", norm2(phi_p_int)
    end if

    print "(A)", ""
    print "(A,ES14.6)", "  |phi|   = ", norm2(phi)
    print "(A,ES14.6)", "  max|c|  = ", maxval(abs(c))
    if (sweep_has_x .and. sweep_has_z) then
      call build_adjoint_source(dev, r_beam_list((nk + 1) / 2), s, &
        & x_beam = x_beam_list((nk + 1) / 2), z_beam = z_beam_list((nk + 1) / 2))
    else if (sweep_has_x) then
      call build_adjoint_source(dev, r_beam_list((nk + 1) / 2), s, &
        & x_beam = x_beam_list((nk + 1) / 2))
    else if (sweep_has_z) then
      call build_adjoint_source(dev, r_beam_list((nk + 1) / 2), s, &
        & z_beam = z_beam_list((nk + 1) / 2))
    else
      call build_adjoint_source(dev, r_beam_list((nk + 1) / 2), s)
    end if
    print "(A,ES14.6)", "  |s|     = ", norm2(s)
    print "(A,ES14.6)", "  max|s|  = ", maxval(abs(s))

    ! Persist the scalar sweep, collection-probability interior slices, and a
    ! few scalar summaries to an fbs file for downstream plotting. The grid
    ! (coordinates for the phi^n/phi^p interior entries) lives in device.fbs,
    ! which dev%init already emits alongside this file.
    block
      type(storage)             :: st
      character(:), allocatable :: fbs_path
      real,         allocatable :: eta_abs_arr(:)
      integer                   :: kk

      fbs_path = suite_name%s // ".fbs"
      allocate (eta_abs_arr(nk))
      do kk = 1, nk
        eta_abs_arr(kk) = denorm(I_EBIC(kk), 'A') / (Q_ELEM_C * total_pairs_per_sec)
      end do

      call st%open(fbs_path, flag = STORAGE_WRITE)
      call st%write("adjoint/r_beam",             r_beam_list,   unit = "nm")
      if (sweep_has_x) &
        call st%write("adjoint/x_beam",           x_beam_list,   unit = "nm")
      if (sweep_has_z) &
        call st%write("adjoint/z_beam",           z_beam_list,   unit = "nm")
      call st%write("adjoint/I_EBIC",             I_EBIC,        unit = "A")
      call st%write("adjoint/eta_absolute",       eta_abs_arr)
      call st%write("adjoint/phi_n_interior",     phi_n_int)
      if (have_holes) &
        call st%write("adjoint/phi_p_interior",   phi_p_int)
      call st%write("adjoint/G_total",            [total_pairs_per_sec])
      call st%write("adjoint/q_times_G_total",    [Q_ELEM_C * total_pairs_per_sec])
      call st%write("adjoint/setup_time",         [setup_sec])
      call st%write("adjoint/sweep_time",         [sweep_sec])
      call st%write("adjoint/V_contact",          V_contact, unit = "V")
      call st%write("adjoint/I_dark",             I_dark,    unit = "A")
      call st%close()
      print "(A)", ""
      print "(A)", "  Results written to " // fbs_path
    end block

    print "(A)", ""
    print "(A)", "=== Case " // suite_name%s // " complete ==="

    block
      logical :: run_fwd, has_fwd_key
      call runfile%get(si, "run_forward_validation", run_fwd, status = has_fwd_key)
      if (.not. has_fwd_key) run_fwd = .true.
      if (.not. run_fwd) then
        print "(A)", ""
        print "(A)", "--- Forward validation skipped (run_forward_validation=false) ---"
      else
        print "(A)", ""
        print "(A)", "--- Forward EBIC validation at middle r_beam ---"
        if (sweep_has_x .and. sweep_has_z) then
          call run_forward_validation(r_beam_list((nk + 1) / 2), I_beam_saved, ict_collect, &
            & V_contact, I_dark, &
            & x_beam = x_beam_list((nk + 1) / 2), z_beam = z_beam_list((nk + 1) / 2))
        else if (sweep_has_x) then
          call run_forward_validation(r_beam_list((nk + 1) / 2), I_beam_saved, ict_collect, &
            & V_contact, I_dark, x_beam = x_beam_list((nk + 1) / 2))
        else if (sweep_has_z) then
          call run_forward_validation(r_beam_list((nk + 1) / 2), I_beam_saved, ict_collect, &
            & V_contact, I_dark, z_beam = z_beam_list((nk + 1) / 2))
        else
          call run_forward_validation(r_beam_list((nk + 1) / 2), I_beam_saved, ict_collect, &
            & V_contact, I_dark)
        end if
      end if
    end block

    ! Optional forward-linescan (direct reference for adjoint validation). Runs
    ! N_fwd forward Newton solves with the beam swept along y at fixed x and
    ! (for 3D) fixed z. Results saved to <name>_fwd_linescan.fbs. Controlled by
    ! forward_linescan_n (+ optional forward_linescan_x, _y_min, _y_max, _z).
    block
      integer              :: n_fwd, kk
      real, allocatable    :: fwd_y(:), fwd_I(:)
      real                 :: fx, fz, fy_min, fy_max, y_lo, y_hi
      logical              :: has_n, has_x, has_z, has_ymin, has_ymax
      type(storage)        :: st_fwd
      integer(8)           :: t_ls_start, t_ls_end, clock_rate_ls

      call runfile%get(si, "forward_linescan_n",     n_fwd,  status = has_n)
      if (has_n .and. n_fwd > 0) then
        call runfile%get(si, "forward_linescan_x",     fx,     status = has_x)
        call runfile%get(si, "forward_linescan_z",     fz,     status = has_z)
        call runfile%get(si, "forward_linescan_y_min", fy_min, status = has_ymin)
        call runfile%get(si, "forward_linescan_y_max", fy_max, status = has_ymax)
        y_lo = fy_min
        y_hi = fy_max
        if (.not. has_ymin) y_lo = minval(r_beam_list)
        if (.not. has_ymax) y_hi = maxval(r_beam_list)
        if (.not. has_x)    fx   = dev%par%reg_beam(1)%beam_x

        fwd_y = linspace(y_lo, y_hi, n_fwd)
        allocate (fwd_I(n_fwd))

        print "(A)", ""
        print "(A,I0,A)", "--- Forward linescan: ", n_fwd, " points ---"

        call system_clock(count_rate = clock_rate_ls)
        call system_clock(t_ls_start)
        do kk = 1, n_fwd
          if (has_z) then
            call run_forward_validation(fwd_y(kk), I_beam_saved, ict_collect,       &
              & V_contact, I_dark, x_beam = fx, z_beam = fz,                        &
              & I_ebic_out = fwd_I(kk), quiet = .true.)
          else
            call run_forward_validation(fwd_y(kk), I_beam_saved, ict_collect,       &
              & V_contact, I_dark, x_beam = fx,                                     &
              & I_ebic_out = fwd_I(kk), quiet = .true.)
          end if
          print "(A,I4,A,F8.1,A,ES16.6,A)", "    k=", kk, "  y = ",                 &
            & denorm(fwd_y(kk), 'nm'), " nm  I_EBIC = ",                            &
            & denorm(fwd_I(kk), 'A'), " A"
        end do
        call system_clock(t_ls_end)
        print "(A,F10.4,A)", "  Forward linescan total time: ",                      &
          & real(t_ls_end - t_ls_start) / real(clock_rate_ls), " s"

        call st_fwd%open(suite_name%s // "_fwd_linescan.fbs", flag = STORAGE_WRITE)
        call st_fwd%write("fwd/y",      fwd_y, unit = "nm")
        call st_fwd%write("fwd/x",      [fx],  unit = "nm")
        call st_fwd%write("fwd/I_EBIC", fwd_I, unit = "A")
        if (has_z) call st_fwd%write("fwd/z", [fz], unit = "nm")
        call st_fwd%close()
        print "(A)", "  Results written to " // suite_name%s // "_fwd_linescan.fbs"
      end if
    end block

    if (allocated(delta_u_primal)) &
      & call compare_du_fwd_vs_primal(c, ict_collect, u_dark_state, delta_u_primal)
  end subroutine

  subroutine run_primal_diagnostic(si, ss, c, phi, r_beam_list, x_beam_list,           &
                                 & z_beam_list, sweep_has_x, sweep_has_z,              &
                                 & u_dark_state, delta_u_primal, s_primal_save)
    !! Adjoint self-consistency + Jacobian-correctness diagnostic. Gated by the
    !! debug_primal_solve key in the run INI; if the key is absent or false,
    !! returns silently with delta_u_primal / s_primal_save unallocated, which
    !! the caller uses as a "did the diagnostic run?" flag for the downstream
    !! compare_du_fwd_vs_primal call.
    !!
    !! Three independent tests, in order:
    !!   1. Transpose identity:   c · (A^-1 s) ?= phi · s    (must be ~1e-12).
    !!   2. Source equivalence:   adjoint s = -jaco_bgen.bgen vs forward
    !!                            implicit -F(u_dark, beam on) at continuity rows.
    !!   3. Taylor probe:         |F(u_dark + A^-1 s, beam on)| / |s|. If A is
    !!                            the true linearization, residual is O(|s|^2);
    !!                            an O(|s|) residual localizes the bug to A.
    integer,            intent(in)    :: si
    type(steady_state), intent(inout) :: ss
    real,               intent(in)    :: c(:), phi(:)
    real,               intent(in)    :: r_beam_list(:), x_beam_list(:), z_beam_list(:)
    logical,            intent(in)    :: sweep_has_x, sweep_has_z
    real,               intent(in)    :: u_dark_state(:)
    real,  allocatable, intent(out)   :: delta_u_primal(:), s_primal_save(:)

    logical              :: do_primal, has_key
    real,    allocatable :: s_primal(:), delta_u(:)
    real,    allocatable :: f_beam_on(:), s_implicit(:)
    real,    allocatable :: u_test(:), f_test(:)
    real                 :: I_primal, I_adjoint, rel_diff
    real                 :: ratio, diff_norm, s_adj_norm, s_imp_norm, f_test_norm
    integer              :: iimvar, ibl, ilo, ihi

    call runfile%get(si, "debug_primal_solve", do_primal, status = has_key)
    if (.not. (has_key .and. do_primal)) return

    print "(A)", ""
    print "(A)", "--- DEBUG: primal-vs-adjoint self-consistency (r_beam_list(1)) ---"

    if (sweep_has_x .and. sweep_has_z) then
      call build_adjoint_source(dev, r_beam_list(1), s_primal,                          &
        & x_beam = x_beam_list(1), z_beam = z_beam_list(1))
    else if (sweep_has_x) then
      call build_adjoint_source(dev, r_beam_list(1), s_primal, x_beam = x_beam_list(1))
    else if (sweep_has_z) then
      call build_adjoint_source(dev, r_beam_list(1), s_primal, z_beam = z_beam_list(1))
    else
      call build_adjoint_source(dev, r_beam_list(1), s_primal)
    end if

    allocate (delta_u(size(s_primal)))
    call ss%solver%solve(s_primal, delta_u, trans = 'N')

    allocate (delta_u_primal, source = delta_u)
    allocate (s_primal_save,  source = s_primal)

    I_primal  = dot_product(c, delta_u)
    I_adjoint = dot_product(phi, s_primal)
    if (abs(I_adjoint) > 0.0) then
      rel_diff = (I_primal - I_adjoint) / I_adjoint
    else
      rel_diff = 0.0
    end if

    print "(A,ES20.12,A)", "  primal  c*(A^-1 s):  ", denorm(I_primal,  'A'), " A"
    print "(A,ES20.12,A)", "  adjoint phi*s     :  ", denorm(I_adjoint, 'A'), " A"
    print "(A,ES12.4)",    "  (primal - adjoint) / adjoint (should be ~1e-12): ", rel_diff

    ! Source equivalence: adjoint s vs forward implicit -F(u_dark, beam on).
    ! At u = u_dark, F_dark(u_dark) = 0; turning the beam on perturbs only the
    ! continuity rows (-V*G), so -F is exactly what the first Newton step solves.
    dev%beam_pos%x = r_beam_list(1)

    allocate (f_beam_on(dev%sys_full%n))
    call dev%sys_full%eval(f = f_beam_on)
    s_implicit = -f_beam_on

    print "(A)", ""
    print "(A)", "--- DEBUG: adjoint's s vs forward's implicit s (-F at u_dark + beam) ---"
    print "(A)", "  main_var / tab      |  |s_adj|       |s_imp|       |diff|        diff/|s_adj|"

    do iimvar = 1, dev%sys_full%g%imvar%n
      associate (mv => dev%sys_full%get_main_var(iimvar))
        do ibl = 1, mv%ntab
          ilo = dev%sys_full%i0(dev%sys_full%res2block(iimvar)%d(ibl))
          ihi = dev%sys_full%i1(dev%sys_full%res2block(iimvar)%d(ibl))
          s_adj_norm = norm2(s_primal_save(ilo:ihi))
          s_imp_norm = norm2(s_implicit(ilo:ihi))
          diff_norm  = norm2(s_primal_save(ilo:ihi) - s_implicit(ilo:ihi))
          if (s_adj_norm > 0.0 .or. s_imp_norm > 0.0) then
            print "(A,A16,A,I1,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)",                   &
              & "  ", adjustl(mv%name), "/t", ibl,                                       &
              & "   ", s_adj_norm, "  ", s_imp_norm, "  ", diff_norm,                    &
              & "   ", diff_norm / max(s_adj_norm, tiny(1.0))
          end if
        end do
      end associate
    end do

    if (norm2(s_primal_save) > 0.0) then
      ratio = dot_product(s_primal_save, s_implicit) / (norm2(s_primal_save)**2)
      print "(A,ES20.12)", "  <s_adj, s_imp> / |s_adj|^2 (= s_imp / s_adj if parallel): ", ratio
    end if

    ! Taylor probe at u_dark + A^-1 s. If A = J(u_dark), residual is O(|s|^2).
    u_test = u_dark_state + delta_u_primal
    call dev%sys_full%set_x(u_test)
    allocate (f_test(dev%sys_full%n))
    call dev%sys_full%eval(f = f_test)
    f_test_norm = norm2(f_test)
    print "(A,ES14.6)", "  |F(u_dark + A^-1 s, beam on)| = ", f_test_norm
    print "(A,ES14.6)", "  |s_adj|                       = ", norm2(s_primal_save)
    print "(A,ES14.6)", "  ratio (should be ~|s|^2-scale): ",                            &
      & f_test_norm / max(norm2(s_primal_save), tiny(1.0))
    ! Restore u_dark so the downstream forward Newton starts from a clean state.
    call dev%sys_full%set_x(u_dark_state)
  end subroutine

  subroutine compare_du_fwd_vs_primal(c, ict_collect, u_dark_state, delta_u_primal)
    !! Forward-vs-primal linearization comparison. Must be called AFTER the
    !! forward Newton solve so dev%sys_full%get_x() returns u_fwd. At V=0 and
    !! small I_beam, u_fwd - u_dark must equal A^-1 s to machine precision
    !! (up to O(|s|^2)). Discrepancy here is the smoking gun for A != J(u_dark)
    !! as seen by the forward -- the open ~5.24e-4 bug.
    real,    intent(in) :: c(:)
    integer, intent(in) :: ict_collect
    real,    intent(in) :: u_dark_state(:), delta_u_primal(:)

    real,    allocatable      :: u_fwd_state(:), delta_u_fwd(:)
    real                      :: du_fwd_norm, du_pri_norm, du_diff_norm
    real                      :: I_fwd_dot, I_adj_dot, I_dark_at_ict
    real                      :: d_fwd_norm, d_pri_norm, d_dif_norm
    integer                   :: ict_dof, iimvar, itab, ibl, ilo, ihi, n_imvar
    class(vselector), pointer :: mv

    u_fwd_state = dev%sys_full%get_x()
    delta_u_fwd = u_fwd_state - u_dark_state

    du_fwd_norm  = norm2(delta_u_fwd)
    du_pri_norm  = norm2(delta_u_primal)
    du_diff_norm = norm2(delta_u_fwd - delta_u_primal)

    iimvar  = dev%sys_full%search_main_var("currents")
    ibl     = dev%sys_full%res2block(iimvar)%d(1)
    ict_dof = dev%sys_full%i0(ibl) + ict_collect - 1
    I_fwd_dot     = dot_product(c, delta_u_fwd)
    I_adj_dot     = dot_product(c, delta_u_primal)
    I_dark_at_ict = u_dark_state(ict_dof)

    print "(A)", ""
    print "(A)", "--- DEBUG: forward vs adjoint linearization of u ---"
    print "(A,ES14.6)",    "  |delta_u_forward|           = ", du_fwd_norm
    print "(A,ES14.6)",    "  |delta_u_primal (A^-1 s)|   = ", du_pri_norm
    print "(A,ES14.6)",    "  |fwd - primal|              = ", du_diff_norm
    print "(A,ES14.6)",    "  relative error              = ",                          &
      & du_diff_norm / max(du_pri_norm, tiny(1.0))
    print "(A)",           ""
    print "(A,ES20.12,A)", "  c . delta_u_fwd  (should = I_fwd - I_dark):  ",            &
      & denorm(I_fwd_dot, 'A'), " A"
    print "(A,ES20.12,A)", "  c . delta_u_pri  (= phi . s, adjoint linear):",            &
      & denorm(I_adj_dot, 'A'), " A"
    print "(A,ES20.12,A)", "  u_dark (I_ct_dof) [check ~= 0 at V=0]:       ",            &
      & denorm(I_dark_at_ict, 'A'), " A"
    if (abs(I_adj_dot) > 0.0) &
      print "(A,ES12.4)",  "  (c.dU_fwd - c.dU_pri) / c.dU_pri:             ",           &
        & (I_fwd_dot - I_adj_dot) / I_adj_dot

    ! Per-main-variable breakdown of |dU_fwd - dU_pri|. Identifies which main
    ! variable (pot / ndens / pdens / iref_* / cdens_* / currents) carries the
    ! gauge direction.
    print "(A)", ""
    print "(A)", "--- Per-main-variable breakdown ---"
    print "(A)", "  main_var / tab      | blk   i0..i1       |dU_fwd|      |dU_pri|      |diff|        diff/|pri|"

    n_imvar = dev%sys_full%g%imvar%n
    do iimvar = 1, n_imvar
      mv => dev%sys_full%get_main_var(iimvar)
      do itab = 1, mv%ntab
        ibl = dev%sys_full%res2block(iimvar)%d(itab)
        ilo = dev%sys_full%i0(ibl)
        ihi = dev%sys_full%i1(ibl)
        d_fwd_norm = norm2(delta_u_fwd(ilo:ihi))
        d_pri_norm = norm2(delta_u_primal(ilo:ihi))
        d_dif_norm = norm2(delta_u_fwd(ilo:ihi) - delta_u_primal(ilo:ihi))
        print "(A,A16,A,I1,A,I3,A,I6,A,I6,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)",         &
          & "  ", adjustl(mv%name), "/t", itab,                                           &
          & " | ", ibl, "  ", ilo, "..", ihi,                                             &
          & "  ", d_fwd_norm, "  ", d_pri_norm, "  ", d_dif_norm,                         &
          & "  ", d_dif_norm / max(d_pri_norm, tiny(1.0))
      end do
    end do
  end subroutine

  subroutine run_forward_validation(r_beam, I_beam_phys, ict_collect, V_contact, I_dark, &
    & x_beam, z_beam, I_ebic_out, quiet)
    !! Run a fresh forward DC solve with the beam turned on at r_beam, at the same
    !! per-contact bias as the adjoint's DC operating point. The beam-induced
    !! signal is I(beam on) − I_dark — that is what the adjoint's φ·s matches.
    !! For 2D point-beam sweeps pass x_beam; for 3D line sweeps pass z_beam; for 3D
    !! point pass both. Pass I_ebic_out to capture the dark-subtracted current for
    !! a calling linescan. Pass quiet=true to suppress prints when looping.
    real,    intent(in) :: r_beam, I_beam_phys
    integer, intent(in) :: ict_collect
    real,    intent(in) :: V_contact(:), I_dark(:)
    real, optional, intent(in)  :: x_beam, z_beam
    real, optional, intent(out) :: I_ebic_out
    logical, optional, intent(in) :: quiet

    logical :: be_quiet

    type(steady_state)   :: ss_fwd
    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:), frac(:)
    integer              :: sj, ict, n_ramp, kk
    real                 :: I_collect_fwd, I_ebic_fwd, max_bias
    integer(8)           :: t_fwd_start, t_fwd_end, clock_rate
    real                 :: t_fwd_sec

    call runfile%get_section("full newton params", sj)
    call ss_fwd%init(dev%sys_full)
    call ss_fwd%set_params(runfile%sections%d(sj))
    ss_fwd%msg = "Fwd-Newton: "

    ! Same ramping discipline as run_dark_solve — beam is on from the start so
    ! G also scales with the ramp, but that just means the final leg lands on
    ! the target DC state. (For small beam currents the correction is linear
    ! anyway; ramping robustness dominates.)
    max_bias = maxval(abs(V_contact))
    if (denorm(max_bias, 'V') < 0.026) then
      n_ramp = 0
    else
      n_ramp = ceiling(denorm(max_bias, 'V') / 0.1)
    end if

    allocate (V_all(dev%par%nct + 1, n_ramp + 1))
    if (n_ramp == 0) then
      frac = [ 1.0 ]
    else
      frac = [ (real(kk) / real(n_ramp), kk = 0, n_ramp) ]
    end if
    do kk = 1, n_ramp + 1
      do ict = 1, dev%par%nct
        V_all(ict, kk) = frac(kk) * V_contact(ict)
      end do
      V_all(dev%par%nct + 1, kk) = r_beam
    end do

    t_inp = frac
    t_sim = frac
    call input%init(t_inp, V_all)

    dev%par%reg_beam(1)%I_beam = I_beam_phys
    if (present(x_beam)) dev%par%reg_beam(1)%beam_x = x_beam
    if (present(z_beam)) dev%par%reg_beam(1)%beam_z = z_beam

    gummel_restart = .true.
    gummel_once    = .true.
    gummel_enabled = .true.

    call system_clock(count_rate = clock_rate)
    call system_clock(t_fwd_start)
    call ss_fwd%run(input=input, t_input=t_sim, gummel=gummel)
    call system_clock(t_fwd_end)
    t_fwd_sec = real(t_fwd_end - t_fwd_start) / real(clock_rate)

    I_collect_fwd = dev%curr(ict_collect)%x
    I_ebic_fwd    = I_collect_fwd - I_dark(ict_collect)

    be_quiet = .false.
    if (present(quiet)) be_quiet = quiet

    if (.not. be_quiet) then
      if (present(x_beam)) then
        print "(A,F9.2,A,F9.2,A,ES16.6,A)", "  Forward I_collect at (x, y) = (", &
          & denorm(x_beam, 'nm'), " nm, ", denorm(r_beam, 'nm'), " nm) is ", &
          & denorm(I_collect_fwd, 'A'), " A"
      else
        print "(A,F9.2,A,ES16.6,A)", "  Forward I_collect at r_beam = ", denorm(r_beam, 'nm'), &
                             &      " nm  is  ", denorm(I_collect_fwd, 'A'), " A"
      end if
      print "(A,ES16.6,A)",         "  Dark I_collect (beam off, bias on)        = ", &
                             &      denorm(I_dark(ict_collect), 'A'), " A"
      print "(A,ES16.6,A)",         "  Forward I_EBIC = I_collect − I_dark       = ", &
                             &      denorm(I_ebic_fwd, 'A'),                       " A"
      print "(A,ES16.6)",           "  Forward eta = I_EBIC / I_beam             = ", &
                             &      I_ebic_fwd / I_beam_phys
      print "(A,F10.4,A)",          "  Forward Newton time (= brute-force cost per beam position):", &
                             &      t_fwd_sec, " s"
    end if

    if (present(I_ebic_out)) I_ebic_out = I_ebic_fwd
  end subroutine

  subroutine run_dark_solve(ss, V_contact, I_dark, beam_y_park)
    !! Reference DC operating point u* with beam off. Contacts are held at the
    !! (possibly non-zero) V_contact biases. Captured per-contact currents
    !! I_dark are the baseline that the forward validation subtracts to isolate
    !! the beam-induced signal.
    !!
    !! Bias ramping: if any |V_contact| exceeds V_T, solve at intermediate bias
    !! fractions to keep Newton inside its basin of attraction. Without this,
    !! Newton oscillates for V_P ≳ 0.2V forward on a pn junction — the V=0
    !! nlpe guess is exponentially far from the strongly-injected solution.
    !!
    !! Optional beam_y_park: where to park the beam-y input during dark solve.
    !! I_beam = 0 so physics is unaffected, but setting this to match the
    !! forward's target r_beam keeps the Y_BEAM input-DOF value equal between
    !! u_dark and u_fwd (removes a spurious ~input-scale term from the
    !! |u_fwd - u_dark| diagnostic).
    type(steady_state), intent(inout) :: ss
    real,               intent(in)    :: V_contact(:)
    real,               intent(out)   :: I_dark(:)
    real,  optional,    intent(in)    :: beam_y_park

    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:), frac(:)
    integer              :: ict, n_ramp, k
    real                 :: max_bias, y_park_val

    print "(A)", ""
    print "(A)", "--- Operating-point solve (user bias, I_beam=0) ---"

    max_bias = maxval(abs(V_contact))
    ! At ~0 bias, a single-point solve is enough; no ramp. Otherwise use
    ! ~0.1 V per Newton leg (≈ 4 V_T at 300K), producing n_ramp+1 time points
    ! from 0 to V_contact. Threshold is 1 V_T — below that the dark nlpe/gummel
    ! initial guess is already in Newton's basin of attraction.
    if (denorm(max_bias, 'V') < 0.026) then
      n_ramp = 0
    else
      n_ramp = ceiling(denorm(max_bias, 'V') / 0.1)
      print "(A,I0,A,ES10.3,A)", "  Ramping bias in ", n_ramp, " legs (max |V| = ", &
        & denorm(max_bias, 'V'), " V)"
    end if

    allocate (V_all(dev%par%nct + 1, n_ramp + 1))
    if (n_ramp == 0) then
      frac = [ 1.0 ]   ! single-step solve directly at target (zero bias is a
                       ! trivial 1.0 * 0 = 0; small bias is one-shot without ramp)
    else
      frac = [ (real(k) / real(n_ramp), k = 0, n_ramp) ]
    end if
    y_park_val = 0.0
    if (present(beam_y_park)) y_park_val = beam_y_park
    do k = 1, n_ramp + 1
      do ict = 1, dev%par%nct
        V_all(ict, k) = frac(k) * V_contact(ict)
      end do
      V_all(dev%par%nct + 1, k) = y_park_val   ! beam y (beam is off anyway)
    end do

    t_inp = frac
    t_sim = frac
    call input%init(t_inp, V_all)

    dev%par%reg_beam(1)%I_beam = 0.0

    gummel_restart = .true.
    gummel_once    = .true.
    gummel_enabled = .true.

    call ss%run(input=input, t_input=t_sim, gummel=gummel)

    do ict = 1, dev%par%nct
      I_dark(ict) = dev%curr(ict)%x
      print "(A,A,A,ES12.4,A)", "  I_", trim(dev%par%contacts(ict)%name), &
        " = ", denorm(I_dark(ict), 'A'), " A"
    end do
  end subroutine

  subroutine gummel()
    !! Gummel iteration — adapted from linearity_test.f90.
    !! Runs once at the start of the Newton loop (gummel_once = .true.):
    !! approx imref/potential → nlpe pre-solve → (optional SRH beam-excess hack).
    integer            :: ci, si_nlpe
    type(steady_state) :: ss_nlpe

    if (.not. gummel_enabled) return
    if (gummel_once) gummel_enabled = .false.

    if (gummel_restart) then
      gummel_restart = .false.

      do ci = dev%par%ci0, dev%par%ci1
        call approx_imref(dev%par, dev%iref(ci), dev%volt)
      end do
      call approx_potential(dev%par, dev%pot, dev%iref)
      call dev%sys_nlpe%eval()

      ! Nonlinear Poisson pre-solve: essential for Schottky contacts, harmless for ohmic.
      call runfile%get_section("nlpe params", si_nlpe)
      call ss_nlpe%init(dev%sys_nlpe)
      call ss_nlpe%set_params(runfile%sections%d(si_nlpe))
      ss_nlpe%msg = "nlpe: "
      call ss_nlpe%run()
    end if
  end subroutine

end program
