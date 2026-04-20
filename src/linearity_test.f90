program linearity_test
  !! EBIC Linearity Test Driver
  !!
  !! Runs the three-test suite from docs/adjoint_method/linearity_test.md:
  !!   Test 1 — Homogeneity (I_beam scaling at fixed r_A)
  !!   Test 2 — Superposition (G_A, G_B, G_A+G_B)
  !!   Test 3 — Small-signal ratios (post-processing of Test 1 solutions)
  !!
  !! Does NOT use the adjoint. Forward-only; prerequisite for any adjoint work.
  !!
  !! Two operating modes:
  !!   --device D --run R       : single test case
  !!   --suite M                : read manifest M, spawn subprocess per [case]
  !!                              (re-launches self because dev is a module singleton
  !!                              and full re-init is not supported cleanly)

  use approx_m,           only: approx_imref, approx_potential
  use cl_options_m,       only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,           only: dev
  use error_m,            only: program_error
  use input_m,            only: input_file
  use input_src_m,        only: polygon_src
  use normalization_m,    only: init_normconst, denorm
  use semiconductor_m,    only: CR_NAME, CR_ELEC, CR_HOLE
  use steady_state_m,     only: steady_state
  use string_m,           only: string
  use util_m,             only: get_hostname

  implicit none

  logical                   :: gummel_restart, gummel_once, gummel_enabled
  real                      :: temperature
  type(input_file)          :: runfile
  character(:), allocatable :: suite_file, device_file, run_file

  print "(A)", "=============================================="
  print "(A)", " linearity_test — EBIC adjoint prerequisite"
  print "(A)", "=============================================="
  print "(A)", "Host: " // get_hostname()

  call command_line()

  if (allocated(suite_file)) then
    ! Suite mode: the main process does NOT init dev; it forks children.
    call run_suite_mode()
  else
    ! Single-case mode: init dev + runfile, run one suite.
    if (.not. allocated(device_file)) call program_error("--device is required (or use --suite)")
    if (.not. allocated(run_file))    call program_error("--run is required (or use --suite)")
    call dev%init(device_file, temperature)
    call runfile%init(run_file)
    call run_linearity_suite()
  end if

contains

  subroutine command_line()
    !! Parse command line. Capture device/run/suite into strings without side effects,
    !! then the main block dispatches based on what was provided.
    integer                      :: idesc, i, j
    integer,         allocatable :: iclopt(:), jclopt(:)
    type(cl_option), allocatable :: clopt(:)
    type(cl_option_descriptor)   :: desc(4) = [ &
      cl_option_descriptor('T', "temperature", .true.,  .false., .true., .true.), &
      cl_option_descriptor('d', "device",      .false., .false., .true., .true.), &
      cl_option_descriptor('r', "run",         .false., .false., .true., .true.), &
      cl_option_descriptor('s', "suite",       .false., .false., .true., .true.)  &
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
        case ('s')
          suite_file = clopt(j)%arg
        end select
      end do
    end do
  end subroutine

  subroutine run_suite_mode()
    !! Read manifest and launch self as a subprocess per [case] section.
    !! On first subprocess failure, abort (user-requested strict mode).
    type(input_file)          :: manifest
    integer,      allocatable :: sids(:)
    integer                   :: k, status, exitstat, cmdstat
    type(string)              :: case_name, case_device, case_run
    character(:), allocatable :: self_path, temp_str, cmd
    character(8)              :: temp_buf

    print "(A)", ""
    print "(A)", "--- Suite mode: reading manifest ---"
    print "(A)", "  manifest = " // suite_file

    call manifest%init(suite_file)
    call manifest%get_sections("case", sids)
    if (size(sids) == 0) call program_error("suite manifest contains no [case] sections")
    print "(A,I0,A)", "  Found ", size(sids), " case(s)"

    ! Locate our own executable path for the subprocess invocation.
    ! get_command_argument(0) gives the invocation path; assume it's usable as-is.
    call get_self_path(self_path)

    ! Format temperature as a plain string for the subprocess CLI.
    write (temp_buf, "(F8.3)") temperature
    temp_str = trim(adjustl(temp_buf))

    do k = 1, size(sids)
      call manifest%get(sids(k), "name",   case_name)
      call manifest%get(sids(k), "device", case_device)
      call manifest%get(sids(k), "run",    case_run)

      print "(A)", ""
      print "(A)", "=============================================="
      print "(A,I0,A,I0,A)", " Case ", k, "/", size(sids), ": " // case_name%s
      print "(A)", "=============================================="

      ! Clean up cybear output artifacts from previous cases. device_params_init
      ! opens device.fbs in STORAGE_WRITE mode which EXTENDS an existing file
      ! rather than truncating — so the second case crashes with "variable
      ! sys/grids/grid already exists" unless we remove the file first.
      call execute_command_line("rm -f device.fbs device.fbs.log", wait=.true.)

      cmd = self_path // &
            " --temperature " // temp_str // &
            " --device "      // case_device%s // &
            " --run "         // case_run%s

      print "(A)", "  $ " // cmd

      call execute_command_line(cmd, wait=.true., exitstat=exitstat, cmdstat=cmdstat)

      if (cmdstat /= 0) then
        print "(A,I0)", "  ERROR: failed to spawn subprocess. cmdstat = ", cmdstat
        call program_error("suite aborted: could not launch case " // case_name%s)
      end if
      if (exitstat /= 0) then
        print "(A,I0)", "  ERROR: case exited with status ", exitstat
        call program_error("suite aborted: case " // case_name%s // " failed")
      end if

      print "(A)", "  Case " // case_name%s // " PASSED"
    end do

    print "(A)", ""
    print "(A)", "=============================================="
    print "(A,I0,A)", " Suite complete: all ", size(sids), " case(s) passed"
    print "(A)", "=============================================="
  end subroutine

  subroutine get_self_path(path)
    !! Return the path used to invoke this binary (argv[0] equivalent).
    character(:), allocatable, intent(out) :: path
    character(4096) :: buf
    integer         :: arglen
    call get_command_argument(0, buf, length=arglen)
    path = buf(1:arglen)
  end subroutine

  subroutine run_linearity_suite()
    !! Top-level dispatch for a single test case (one device + one run INI).
    integer              :: sid, sj
    integer, allocatable :: sids(:)
    type(steady_state)   :: ss
    type(string)         :: suite_name
    real                 :: r_A, r_B, I_exp
    real, allocatable    :: decade_multipliers(:), superposition_levels(:)
    logical              :: run_t1, run_t2, run_t3

    call runfile%get_sections("linearity test", sids)
    if (size(sids) /= 1) call program_error("expected exactly one [linearity test] section in run INI")
    sid = sids(1)
    call runfile%get_section("full newton params", sj)

    call runfile%get(sid, "name", suite_name)
    call runfile%get(sid, "r_A",   r_A)
    call runfile%get(sid, "r_B",   r_B)
    call runfile%get(sid, "I_exp", I_exp)
    call runfile%get(sid, "decade_multipliers",   decade_multipliers)
    call runfile%get(sid, "superposition_levels", superposition_levels)
    call runfile%get(sid, "run_test_1", run_t1)
    call runfile%get(sid, "run_test_2", run_t2)
    call runfile%get(sid, "run_test_3", run_t3)

    print "(A)", ""
    print "(A)",       "Linearity case: "  // suite_name%s
    print "(A,F8.3,A)", "  r_A   = ", denorm(r_A, 'um'), " um"
    print "(A,F8.3,A)", "  r_B   = ", denorm(r_B, 'um'), " um"
    print "(A,ES12.4,A)", "  I_exp = ", denorm(I_exp, 'A/cm'), " A/cm"

    call ss%init(dev%sys_full)
    call ss%set_params(runfile%sections%d(sj))
    ss%msg = "Newton: "

    call run_dark_solve(ss, r_A)

    if (run_t1) call run_test_1_scaling(ss, r_A, I_exp, decade_multipliers)
    if (run_t2) print "(A)", "[Test 2 not yet implemented — stub]"
    if (run_t3) print "(A)", "[Test 3 not yet implemented — stub]"

    print "(A)", ""
    print "(A)", "=== Case " // suite_name%s // " complete ==="
  end subroutine

  subroutine run_dark_solve(ss, r_A)
    !! Dark reference u_0: V=0 on all contacts, I_beam = 0, beam parked at r_A.
    type(steady_state), intent(inout) :: ss
    real,               intent(in)    :: r_A

    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:)
    real                 :: I_beam_saved
    integer              :: ict

    print "(A)", ""
    print "(A)", "--- Dark reference solve (V=0, I_beam=0) ---"

    allocate (V_all(dev%par%nct + 1, 1))
    do ict = 1, dev%par%nct
      V_all(ict, 1) = 0.0
    end do
    V_all(dev%par%nct + 1, 1) = r_A

    t_inp = [0.0]
    t_sim = [0.0]
    call input%init(t_inp, V_all)

    I_beam_saved = dev%par%reg_beam(1)%I_beam
    dev%par%reg_beam(1)%I_beam = 0.0

    gummel_restart = .true.
    gummel_once    = .true.
    gummel_enabled = .true.

    call ss%run(input=input, t_input=t_sim, gummel=gummel)

    dev%par%reg_beam(1)%I_beam = I_beam_saved

    do ict = 1, dev%par%nct
      print "(A,A,A,ES12.4,A)", "  I_", trim(dev%par%contacts(ict)%name), &
        " = ", denorm(dev%curr(ict)%x, 'A'), " A"
    end do
  end subroutine

  subroutine run_test_1_scaling(ss, r_A, I_exp, decades)
    !! Test 1 — Homogeneity. Sweep I_beam = 10^decades × I_exp at fixed r_A.
    type(steady_state), intent(inout) :: ss
    real,               intent(in)    :: r_A, I_exp
    real,               intent(in)    :: decades(:)

    type(polygon_src)    :: input
    real,    allocatable :: V_all(:,:), t_inp(:), t_sim(:)
    real,    allocatable :: I_beam_sweep(:), I_ebic_sweep(:,:)
    integer              :: ict, k, nk

    nk = size(decades)
    allocate (I_beam_sweep(nk), I_ebic_sweep(dev%par%nct, nk))

    print "(A)", ""
    print "(A)", "--- Test 1: Homogeneity scaling ---"

    do k = 1, nk
      I_beam_sweep(k) = I_exp * 10.0**decades(k)
    end do

    allocate (V_all(dev%par%nct + 1, 1))
    do ict = 1, dev%par%nct
      V_all(ict, 1) = 0.0
    end do
    V_all(dev%par%nct + 1, 1) = r_A

    t_inp = [0.0]
    t_sim = [0.0]
    call input%init(t_inp, V_all)

    do k = 1, nk
      dev%par%reg_beam(1)%I_beam = I_beam_sweep(k)

      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.

      print "(A,I0,A,I0,A,ES12.4,A)", "  Level ", k, "/", nk, &
        ": I_beam = ", denorm(I_beam_sweep(k), 'A/cm'), " A/cm"
      call ss%run(input=input, t_input=t_sim, gummel=gummel)

      do ict = 1, dev%par%nct
        I_ebic_sweep(ict, k) = dev%curr(ict)%x
        print "(A,A,A,ES12.4,A)", "    I_", trim(dev%par%contacts(ict)%name), &
          " = ", denorm(dev%curr(ict)%x, 'A'), " A"
      end do
    end do

    ! Summary table — prints collection efficiency (I_EBIC / I_beam) for each contact.
    print "(A)", ""
    print "(A)", "  Test 1 summary (collection efficiency I_EBIC / I_beam):"
    print "(A)", "  ────────────────────────────────────────────────────────"
    write (*, "(A,A)", advance="no") "   decade    I_beam [A/cm]     "
    do ict = 1, dev%par%nct
      write (*, "(A,A,A)", advance="no") "I_", trim(dev%par%contacts(ict)%name), "/I_beam      "
    end do
    print *
    do k = 1, nk
      write (*, "(F8.2,3X,ES13.4)", advance="no") decades(k), denorm(I_beam_sweep(k), 'A/cm')
      do ict = 1, dev%par%nct
        write (*, "(3X,ES13.4)", advance="no") I_ebic_sweep(ict, k) / I_beam_sweep(k)
      end do
      print *
    end do
    print "(A)", "  ────────────────────────────────────────────────────────"
    print "(A)", "  Perfect linearity → ratio constant across all decades."
  end subroutine

  subroutine gummel()
    !! Gummel iteration adapted from stem_ebic.f90
    integer            :: ci, si, si_nlpe
    type(steady_state) :: ss_dd(2), ss_nlpe

    if (.not. gummel_enabled) return
    if (gummel_once) gummel_enabled = .false.

    if (gummel_restart) then
      gummel_restart = .false.

      do ci = dev%par%ci0, dev%par%ci1
        call approx_imref(dev%par, dev%iref(ci), dev%volt)
      end do
      call approx_potential(dev%par, dev%pot, dev%iref)
      call dev%sys_nlpe%eval()

      ! Nonlinear Poisson pre-solve. Establishes a potential consistent with
      ! the Schottky boundary conditions at fixed imref, so that the subsequent
      ! full-Newton DD solve doesn't have to walk uphill through ~φ_b of
      ! mismatched potential against a dx_lim step cap. Critical for Schottky;
      ! harmless for ohmic (converges in 0–1 iterations).
      call runfile%get_section("nlpe params", si_nlpe)
      call ss_nlpe%init(dev%sys_nlpe)
      call ss_nlpe%set_params(runfile%sections%d(si_nlpe))
      ss_nlpe%msg = "nlpe: "
      call ss_nlpe%run()

      if (dev%par%has_beam_gen .and. dev%par%smc%srh) then
        do ci = dev%par%ci0, dev%par%ci1
          call dev%calc_bgen(ci)%eval()
        end do
        block
          integer              :: i_bv
          integer, allocatable :: idx_bv(:)
          real                 :: G_bv, tau_bv, excess_bv

          allocate (idx_bv(dev%par%g%idx_dim))
          do i_bv = 1, dev%par%transport_vct(0)%n
            idx_bv = dev%par%transport_vct(0)%get_idx(i_bv)
            G_bv = dev%bgen(CR_ELEC)%get(idx_bv)
            if (G_bv > 0.0) then
              tau_bv = dev%par%srh_tau(CR_ELEC)%get(idx_bv)
              excess_bv = G_bv * tau_bv
              call dev%dens(CR_ELEC)%set(idx_bv, dev%dens(CR_ELEC)%get(idx_bv) + excess_bv)
              tau_bv = dev%par%srh_tau(CR_HOLE)%get(idx_bv)
              excess_bv = G_bv * tau_bv
              call dev%dens(CR_HOLE)%set(idx_bv, dev%dens(CR_HOLE)%get(idx_bv) + excess_bv)
            end if
          end do
          deallocate (idx_bv)
        end block
        do ci = dev%par%ci0, dev%par%ci1
          call dev%calc_iref(ci)%eval()
        end do
      end if
    end if

    call runfile%get_section("dd params", si)
    do ci = dev%par%ci0, dev%par%ci1
      call ss_dd(ci)%init(dev%sys_dd(ci))
      call ss_dd(ci)%set_params(runfile%sections%d(si))
      ss_dd(ci)%msg = CR_NAME(ci) // "DD: "
      call ss_dd(ci)%run()
    end do
  end subroutine

end program
