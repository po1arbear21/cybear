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
                           & evaluate_ebic, extract_phi_for_var
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
  use string_m,        only: string
  use util_m,          only: get_hostname

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
    integer                    :: sj, ict_collect, nk, k, n_n_int, n_p_int, i
    real                       :: I_beam_saved
    real,          allocatable :: r_beam_list(:), I_EBIC(:)
    real,          allocatable :: c(:), phi(:), s(:)
    real,          allocatable :: phi_n_int(:), phi_p_int(:)
    logical                    :: have_holes
    logical                    :: status

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

    ! r_beam can be specified two ways:
    !   Explicit list:   r_beam = v1, v2, ... : nm
    !   Linspace:        r_beam_min = X : nm, r_beam_max = Y : nm, r_beam_n = N
    ! If all three linspace keys are present, they take precedence.
    block
      real    :: r_beam_min, r_beam_max
      integer :: r_beam_n
      logical :: has_min, has_max, has_n

      call runfile%get(si, "r_beam_min", r_beam_min, status = has_min)
      call runfile%get(si, "r_beam_max", r_beam_max, status = has_max)
      call runfile%get(si, "r_beam_n",   r_beam_n,   status = has_n)
      if (has_min .and. has_max .and. has_n) then
        r_beam_list = linspace(r_beam_min, r_beam_max, r_beam_n)
      else
        call runfile%get(si, "r_beam", r_beam_list)
      end if
    end block

    call dev%par%contact_map%get(contact_name, ict_collect, status)
    if (.not. status) call program_error("contact not found: " // contact_name%s)
    print "(A,A,A,I0,A)", "  Collecting contact: ", trim(contact_name%s), &
                       &  " (index ", ict_collect, ")"
    print "(A,I0)", "  Number of beam positions: ", size(r_beam_list)

    call runfile%get_section("full newton params", sj)
    call ss%init(dev%sys_full)
    call ss%set_params(runfile%sections%d(sj))
    ss%msg = "Newton: "

    call system_clock(count_rate = clock_rate)
    call system_clock(t_setup_start)

    I_beam_saved = dev%par%reg_beam(1)%I_beam
    call run_dark_solve(ss)

    print "(A)", ""
    print "(A)", "--- Refactorize J at converged u* ---"
    call refactorize_at_converged(dev, ss%solver)

    print "(A)", ""
    print "(A)", "--- Build c, solve adjoint J^T phi = c ---"
    call build_collection_functional(dev, ict_collect, c)
    call solve_adjoint(ss%solver, c, phi)
    print "(A,I0)", "  adjoint DOFs: ", size(phi)

    call system_clock(t_setup_end)

    ! Restore physical I_beam; G enters s linearly so sign of the beam matters
    ! but the magnitude cleanly multiplies the final I_EBIC.
    dev%par%reg_beam(1)%I_beam = I_beam_saved

    nk = size(r_beam_list)
    allocate (I_EBIC(nk))

    print "(A)", ""
    print "(A)", "--- Beam sweep via adjoint dot-product ---"
    call system_clock(t_sweep_start)
    do k = 1, nk
      call build_adjoint_source(dev, r_beam_list(k), s)
      I_EBIC(k) = evaluate_ebic(phi, s)
    end do
    call system_clock(t_sweep_end)

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
    call build_adjoint_source(dev, r_beam_list((nk + 1) / 2), s)
    print "(A,ES14.6)", "  |s|     = ", norm2(s)
    print "(A,ES14.6)", "  max|s|  = ", maxval(abs(s))

    print "(A)", ""
    print "(A)", "=== Case " // suite_name%s // " complete ==="
    print "(A)", ""
    print "(A)", "--- Forward EBIC validation at middle r_beam ---"
    call run_forward_validation(r_beam_list((nk + 1) / 2), I_beam_saved, ict_collect)
  end subroutine

  subroutine run_forward_validation(r_beam, I_beam_phys, ict_collect)
    !! Run a fresh forward DC solve with the beam turned on at r_beam.
    !! Print the resulting I_<collect> for direct comparison against the adjoint value.
    real,    intent(in) :: r_beam, I_beam_phys
    integer, intent(in) :: ict_collect

    type(steady_state)   :: ss_fwd
    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:)
    integer              :: sj, ict

    call runfile%get_section("full newton params", sj)
    call ss_fwd%init(dev%sys_full)
    call ss_fwd%set_params(runfile%sections%d(sj))
    ss_fwd%msg = "Fwd-Newton: "

    allocate (V_all(dev%par%nct + 1, 1))
    do ict = 1, dev%par%nct
      V_all(ict, 1) = 0.0
    end do
    V_all(dev%par%nct + 1, 1) = r_beam

    t_inp = [0.0]
    t_sim = [0.0]
    call input%init(t_inp, V_all)

    dev%par%reg_beam(1)%I_beam = I_beam_phys

    gummel_restart = .true.
    gummel_once    = .true.
    gummel_enabled = .true.

    call ss_fwd%run(input=input, t_input=t_sim, gummel=gummel)

    print "(A,F9.2,A,ES16.6,A)", "  Forward I_collect at r_beam = ", denorm(r_beam, 'nm'), &
                           &      " nm  is  ", denorm(dev%curr(ict_collect)%x, 'A'), " A"
    print "(A,ES16.6)",           "  Forward eta = I_collect / I_beam = ", &
                           &      dev%curr(ict_collect)%x / I_beam_phys
  end subroutine

  subroutine run_dark_solve(ss)
    !! Dark reference u_0: V=0 on all contacts, I_beam=0, beam parked arbitrarily.
    type(steady_state), intent(inout) :: ss

    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:)
    integer              :: ict

    print "(A)", ""
    print "(A)", "--- Dark reference solve (V=0, I_beam=0) ---"

    allocate (V_all(dev%par%nct + 1, 1))
    do ict = 1, dev%par%nct
      V_all(ict, 1) = 0.0
    end do
    V_all(dev%par%nct + 1, 1) = 0.0   ! beam y (arbitrary; beam is off)

    t_inp = [0.0]
    t_sim = [0.0]
    call input%init(t_inp, V_all)

    dev%par%reg_beam(1)%I_beam = 0.0

    gummel_restart = .true.
    gummel_once    = .true.
    gummel_enabled = .true.

    call ss%run(input=input, t_input=t_sim, gummel=gummel)

    do ict = 1, dev%par%nct
      print "(A,A,A,ES12.4,A)", "  I_", trim(dev%par%contacts(ict)%name), &
        " = ", denorm(dev%curr(ict)%x, 'A'), " A"
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
