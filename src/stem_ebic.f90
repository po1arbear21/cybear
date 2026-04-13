program stem_ebic

  use approx_m,           only: approx_imref, approx_potential
  use bin_search_m,       only: bin_search
  use block_m,            only: block_real
  use cl_options_m,       only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,           only: dev
  use error_m,            only: program_error
  use input_m,            only: input_file, input_section
  use input_src_m,        only: polygon_src, harmonic_src
  use math_m,             only: linspace, logspace, PI
  use normalization_m,    only: init_normconst, norm, denorm
  use semiconductor_m,    only: CR_NAME, CR_ELEC, CR_HOLE
  use small_signal_m,     only: small_signal
  use solver_base_m,      only: solver_real
  use solver_m,           only: default_solver_params, init_solver_real
  use steady_state_m,     only: steady_state
  use storage_m,          only: storage, STORAGE_WRITE, DYNAMIC_EXT, DYNAMIC_APP
  use string_m,           only: string
  use transient_m,        only: transient
  use util_m,             only: get_hostname

  use current_density_m, only: time

  implicit none

  logical          :: gummel_restart, gummel_once, gummel_enabled
  real             :: temperature
  type(input_file) :: runfile

  print "(A)", "Start simulation on " // get_hostname()

  ! parse command line arguments
  call command_line()

  ! ! solve
  call solve_steady_state()
  ! call solve_small_signal()
  ! call solve_transient()
  ! call solve_harmonic_balance()
  ! call solve_responsivity()

contains

  subroutine command_line()
    !! parse command line arguments, init normalization, device and run file

    integer                      :: idesc, i, j
    integer,         allocatable :: iclopt(:), jclopt(:)
    type(cl_option), allocatable :: clopt(:)
    type(cl_option_descriptor)   :: desc(3) = [ &
      cl_option_descriptor('T', "temperature", .true., .false., .true., .true.), &
      cl_option_descriptor('d', "device",      .true., .false., .true., .true.), &
      cl_option_descriptor('r', "run",         .true., .false., .true., .true.)  &
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
          call dev%init(clopt(j)%arg, temperature)

        case ('r')
          call runfile%init(clopt(j)%arg)

        end select
      end do
    end do
  end subroutine

  ! subroutine load_iteration_params(section, n, opt)
  !   character(*),     intent(in)  :: section
  !     !! section name in input file
  !   integer,          intent(in)  :: n
  !     !! system size
  !   type(newton_opt), intent(out) :: opt
  !     !! output newton options

  !   integer :: min_it, max_it
  !   logical :: log
  !   real    :: rtol, atol, lim

  !   call runfile%get(section, "min_it", min_it)
  !   call runfile%get(section, "max_it", max_it)
  !   call runfile%get(section, "rtol",   rtol  )
  !   call runfile%get(section, "atol",   atol  )
  !   call runfile%get(section, "lim",    lim   )
  !   call runfile%get(section, "log",    log   )

  !   call opt%init(n, atol = atol, rtol = rtol, dx_lim = lim, min_it = min_it, max_it = max_it, log = log)
  ! end subroutine

  subroutine voltage_input_ss(sid, t_inp, V, t_sim)
    !! read in voltage configuration of contacts for steady-state from runfile
    !! If beam generation is enabled, includes constant beam position as last row
    integer,           intent(in)  :: sid
      !! section id for runfile
    real, allocatable, intent(out) :: t_inp(:)
      !! pseudo time points for polygon input source
    real, allocatable, intent(out) :: V(:,:)
      !! input values at each pseudo time point
      !! (nct, size(t_inp)) if no beam, (nct+1, size(t_inp)) if beam enabled
    real, allocatable, intent(out) :: t_sim(:)
      !! pseudo time points at which steady-state simulations are performed

    integer              :: ict, ict_sweep, i, ninput
    integer, allocatable :: nsweep_ct(:), nsweep(:)
    logical              :: status
    real,    allocatable :: Vbounds(:), tmp(:), beam_bounds(:)

    ! count inputs: nct voltages + 1 beam (if enabled)
    ninput = dev%par%nct
    if (dev%par%has_beam_gen) ninput = ninput + 1

    ! find out which, if any, contact is swept and how many points are used
    ict_sweep = 0
    do ict = 1, dev%par%nct
      call runfile%get(sid, "N_"//dev%par%contacts(ict)%name, nsweep_ct, status = status)
      if (status) status = (size(nsweep_ct) > 1 .or. nsweep_ct(1) > 1)
      if (status) then
        if (ict_sweep /= 0) call program_error("only one voltage sweep per steady-state section allowed")
        ict_sweep = ict
        nsweep = nsweep_ct
      end if
    end do
    if (.not. allocated(nsweep)) allocate (nsweep(0))

    ! t_inp
    t_inp = linspace(0.0, real(size(nsweep)), size(nsweep)+1)

    ! allocate input array
    allocate(V(ninput, size(nsweep)+1))

    ! contact voltages
    do ict = 1, dev%par%nct
      call runfile%get(sid, "V_"//dev%par%contacts(ict)%name, Vbounds)
      if (ict == ict_sweep) then
        V(ict,:) = Vbounds
      else
        if (size(Vbounds) > 1) call program_error("Voltage bounds given for a contact without sweep")
        V(ict,:) = Vbounds(1)
      end if
    end do

    ! beam position (constant, if beam enabled)
    if (dev%par%has_beam_gen) then
      call runfile%get(sid, "Y_BEAM", beam_bounds, status = status)
      if (status) then
        ! use center of range if bounds given
        if (size(beam_bounds) > 1) then
          V(ninput,:) = (beam_bounds(1) + beam_bounds(2)) / 2.0
        else
          V(ninput,:) = beam_bounds(1)
        end if
      else
        ! default: use device file y-range center
        V(ninput,:) = (dev%par%reg_beam(1)%beam_min + dev%par%reg_beam(1)%beam_max) / 2.0
      end if
    end if

    ! t_sim
    if (ict_sweep == 0) then
      t_sim = [0.0]
    else
      t_sim = [0.0]
      do i = 1, size(nsweep)
        tmp = linspace(i-1.0, i+0.0, nsweep(i))
        t_sim = [t_sim, tmp(2:size(tmp))]
      end do
    end if
  end subroutine

  subroutine beam_input_ss(sid, t_inp, inputs, t_sim)
    !! read in beam sweep configuration for STEM-EBIC steady-state from runfile
    !! Returns full input array: constant voltages + swept beam position
    integer,           intent(in)  :: sid
      !! section id for runfile
    real, allocatable, intent(out) :: t_inp(:)
      !! pseudo time points for polygon input source
    real, allocatable, intent(out) :: inputs(:,:)
      !! all input values: (nct+1, size(t_inp))
      !! rows 1..nct = constant voltages, row nct+1 = beam position
    real, allocatable, intent(out) :: t_sim(:)
      !! pseudo time points at which steady-state simulations are performed

    integer              :: i, ict, ninput
    integer, allocatable :: nsweep(:)
    logical              :: status
    real,    allocatable :: beam_bounds(:), tmp(:), Vbounds(:)

    ninput = dev%par%nct + 1

    ! get beam sweep points
    call runfile%get(sid, "N_BEAM", nsweep, status = status)
    if (.not. status) nsweep = [1]

    ! t_inp (bounds for interpolation)
    t_inp = linspace(0.0, real(size(nsweep)), size(nsweep)+1)

    ! allocate full input array
    allocate(inputs(ninput, size(t_inp)))

    ! read constant voltages
    do ict = 1, dev%par%nct
      call runfile%get(sid, "V_"//dev%par%contacts(ict)%name, Vbounds)
      inputs(ict, :) = Vbounds(1)  ! constant voltage
    end do

    ! read beam position bounds (supports multi-segment: Y_BEAM = y0, y1, y2, ... : um)
    call runfile%get(sid, "Y_BEAM", beam_bounds, status = status)
    if (.not. status) then
      ! default: use device file y-range
      beam_bounds = [dev%par%reg_beam(1)%beam_min, dev%par%reg_beam(1)%beam_max]
    end if
    ! validate: Y_BEAM must have size(N_BEAM)+1 boundary values
    if (size(beam_bounds) /= size(nsweep) + 1) then
      error stop "Y_BEAM must have size(N_BEAM)+1 boundary values"
    end if
    inputs(ninput, :) = beam_bounds  ! swept beam position

    ! t_sim (actual simulation points)
    t_sim = [0.0]
    do i = 1, size(nsweep)
      tmp = linspace(i-1.0, i+0.0, nsweep(i))
      t_sim = [t_sim, tmp(2:size(tmp))]
    end do
  end subroutine

  subroutine solve_steady_state()
    integer              :: si, sj, ict_n, ixb, n_xbeam
    integer, allocatable :: sids(:), nsweep_beam(:)
    logical              :: use_beam_sweep, has_xbeam, status
    real,    allocatable :: V(:,:), t_inp(:), t(:), x_beam_arr(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    call runfile%get_sections("steady state", sids)
    call runfile%get_section("full newton params", sj)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call runfile%get(sids(si), "name", name)

      ! check if beam sweep is active
      use_beam_sweep = .false.
      if (dev%par%has_beam_gen) then
        call runfile%get(sids(si), "N_BEAM", nsweep_beam, status = status)
        if (status) use_beam_sweep = (nsweep_beam(1) > 1)
      end if

      ! input config: either voltage sweep or beam sweep
      if (use_beam_sweep) then
        call beam_input_ss(sids(si), t_inp, V, t)
      else
        call voltage_input_ss(sids(si), t_inp, V, t)
      end if
      call input%init(t_inp, V)

      ! 2D beam sweep: X_BEAM provides multiple x-positions for the beam
      ! If absent, single x-position from device file (standard 1D sweep)
      has_xbeam = .false.
      n_xbeam = 1
      if (use_beam_sweep) then
        call runfile%get(sids(si), "X_BEAM", x_beam_arr, status = has_xbeam)
        if (has_xbeam) n_xbeam = size(x_beam_arr)
      end if

      ! solve steady-state
      call ss%init(dev%sys_full)
      call ss%set_params(runfile%sections%d(sj))
      ss%msg = "Newton: "
      block
        type(string), allocatable :: out_vars(:)
        integer :: ict_out
        allocate(out_vars(6 + 2 * dev%par%nct))
        out_vars(1) = string("pot")
        out_vars(2) = string("ndens")
        out_vars(3) = string("pdens")
        out_vars(4) = string("Ex")
        out_vars(5) = string("Ey")
        out_vars(6) = string("bgen_n")
        do ict_out = 1, dev%par%nct
          out_vars(6 + ict_out) = string("V_" // dev%par%contacts(ict_out)%name)
          out_vars(6 + dev%par%nct + ict_out) = string("I_" // dev%par%contacts(ict_out)%name)
        end do
        call ss%init_output(out_vars, name%s // ".fbs")
      end block

      ! outer loop over beam x-positions (1 iteration if no X_BEAM)
      do ixb = 1, n_xbeam
        if (has_xbeam) then
          dev%par%reg_beam(1)%beam_x = x_beam_arr(ixb)
          print "(A,I0,A,I0,A,F8.2,A)", "beam_x ", ixb, "/", n_xbeam, &
            & " = ", denorm(x_beam_arr(ixb), 'um'), " um"
        end if

        gummel_restart = (ixb == 1)
        gummel_once    = .true.
        gummel_enabled = .true.

        ! Y-sweep for this x-position
        call ss%run(input = input, t_input = t, gummel = gummel, output_hook = eval_efield)

        ! Print current at end of each x-line
        do ict_n = 1, dev%par%nct
          print "(A,A,A,ES12.4,A)", "I_", trim(dev%par%contacts(ict_n)%name), &
            & " = ", denorm(dev%curr(ict_n)%x, 'A'), " A"
        end do
      end do
    end do
  end subroutine

  subroutine solve_small_signal()
    integer              :: si, sj, Nf, i
    integer, allocatable :: sids(:)
    logical              :: flog
    real                 :: f0, f1
    real,    allocatable :: f(:), V(:,:), t_inp(:), t(:)
    complex, allocatable :: s(:), result(:,:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(small_signal)   :: ac
    type(steady_state)   :: ss
    type(storage)        :: st
    type(input_section)  :: params

    call runfile%get_sections("small signal", sids)
    call runfile%get_section("full newton params", sj)
    params = runfile%sections%d(sj)
    call params%set("currents_atol", 1e300)
    call params%set("ndens_dx_lim_rel", 0.2)
    call params%set("pdens_dx_lim_rel", 0.2)

    do si = 1, size(sids)
      print "(A)", "small signal"

      call runfile%get(sids(si), "name", name)

      ! get frequencies
      call runfile%get(sids(si), "f0", f0)
      call runfile%get(sids(si), "f1", f1)
      call runfile%get(sids(si), "Nf", Nf)
      call runfile%get(sids(si), "flog", flog)
      if (flog) then
        f = logspace(f0, f1, Nf)
      else
        f = linspace(f0, f1, Nf)
      end if
      s = 2 * PI * (0.0, 1.0) * f

      ! steady-state config
      call voltage_input_ss(sids(si), t_inp, V, t)
      call input%init(t_inp, V)

      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      ! solve steady-state
      call ss%init(dev%sys_full)
      call ss%set_params(params)
      ss%msg = "Newton: "
      call ss%run(input = input, t_input = t, gummel = gummel)

      ! run small-signal analysis for a single working point
      call ac%init(dev%sys_full)
      call ac%init_output([string("I_GAT")], name%s // ".fbs")
      call ac%run(s)

      deallocate (f, s)
    end do

    call runfile%get_sections("small signal voltage sweep", sids)
    do si = 1, size(sids)
      print "(A)", "small signal voltage sweep"

      call runfile%get(sids(si), "name", name)

      ! get frequencies
      call runfile%get(sids(si), "f0", f0)
      call runfile%get(sids(si), "f1", f1)
      call runfile%get(sids(si), "Nf", Nf)
      call runfile%get(sids(si), "flog", flog)
      if (flog) then
        f = logspace(f0, f1, Nf)
      else
        f = linspace(f0, f1, Nf)
      end if
      s = 2 * PI * (0.0, 1.0) * f

      ! steady-state config
      call voltage_input_ss(sids(si), t_inp, V, t)
      call input%init(t_inp, V)

      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      ! solve steady-state
      call ss%init(dev%sys_full)
      call ss%set_params(params)
      ss%msg = "Newton: "

      ! run small-signal analysis at each working point
      allocate (result(size(t), Nf), source = (0.0,0.0))
      do i = 1, size(t)
        print *, "steady-state step: ", i
        call ss%run(input = input, t_input = [t(i)], gummel = gummel)
        call ac%init(dev%sys_full)
        call ac%run(s)
        result(i,:) = ac%get_scalar("I_GAT", "V_GAT")
        call st%open(name%s // ".fbs", flag = STORAGE_WRITE)
        call st%write("small-signal/s", [s(i)], unit = "Hz", dynamic = DYNAMIC_APP)
        call st%write("small-signal/dI_GAT_dV_GAT", result(i,:), unit = "A/V", dynamic = DYNAMIC_EXT)
        call st%close()
      end do

      deallocate (f, s)
    end do
  end subroutine

  subroutine solve_transient()
    integer              :: ict, Nt, si, si_full_newton, si_transient, start, end, rate
    integer, allocatable :: sids(:)
    real                 :: dt0
    real,    allocatable :: ti(:), Vtmp(:), Vi(:,:)
    type(string)         :: name
    type(steady_state)   :: ss
    type(transient)      :: trans
    type(polygon_src)    :: input

    call runfile%get_sections("transient", sids)
    call runfile%get_section("full newton params", si_full_newton)
    call runfile%get_section("transient params", si_transient)
    do si = 1, size(sids)
      print "(A)", "transient"

      call runfile%get(sids(si), "name", name)
      call runfile%get(sids(si), "dt0", dt0)
      call runfile%get(sids(si), "t", ti)
      Nt = size(ti)
      allocate (Vi(dev%par%nct, Nt))
      do ict = 1, dev%par%nct
        call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vtmp)
        if (size(Vtmp) == 1) then
          Vi(ict,:) = Vtmp(1)
        else
          Vi(ict,:) = Vtmp
        end if
      end do

      call input%init(ti, Vi)

      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.

      ! solve steady-state
      call ss%init(dev%sys_full)
      call ss%set_params(runfile%sections%d(si_full_newton))
      ss%msg = "Newton: "
      call ss%run(input = input, gummel = gummel)

      ! run transient simulation
      call trans%init(dev%sys_full)
      call trans%set_params(runfile%sections%d(si_transient))
      trans%msg = "Transient: "
      call trans%init_output([string("pot"), string("ndens"), string("pdens"), string("ionD"), string("ionA"), &
        &                     string("V_GAT"), string("I_DRN"), string("I_GAT"), string("I_SRC"), string("I_BLK"), &
        &                     string("ncdensx"), string("pcdensx"), string("ncdensy"), string("pcdensy") &
        &                    ], name%s // ".fbs")
      call system_clock(start, rate)
      call trans%run(ti, dt0=dt0, input=input, start_steady_state=.true.)
      call system_clock(end)
      print *, "Transient time: ", real(end-start)/real(rate)
      print *, "Time used for current calculation: ", time
      print *, "which is ", time/(real(end-start)/real(rate)) * 100, "% of the total runtime"
    end do
  end subroutine

  ! subroutine solve_harmonic_balance()
  !   integer                :: si, ict, Nf, NH, Nt
  !   integer, allocatable   :: sids(:)
  !   logical                :: flog, log
  !   real                   :: f0, f1
  !   real,    allocatable   :: volt(:), f(:), c(:,:), s(:,:)
  !   type(string)           :: name

  !   type(harmonic_src)     :: input
  !   type(harmonic_balance) :: hb
  !   type(steady_state)     :: ss
  !   type(newton_opt)       :: opt_hb

  !   allocate (c(dev%par%nct,0:1), s(dev%par%nct,1))

  !   call runfile%get_sections("harmonic balance", sids)
  !   do si = 1, size(sids)
  !     print "(A)", "harmonic balance"

  !     call runfile%get(sids(si), "name", name)

  !     ! get input source
  !     do ict = 1, dev%par%nct
  !       call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, volt)
  !       if (size(volt) /= 3) call program_error("3 values for voltage expected (c0, c1, s1 coefficients)")

  !       c(ict,0) = volt(1)
  !       c(ict,1) = volt(2)
  !       s(ict,1) = volt(3)
  !     end do
  !     call input%init(0.0, c, s)

  !     ! get frequency
  !     call runfile%get(sids(si), "f0", f0)
  !     call runfile%get(sids(si), "f1", f1)
  !     call runfile%get(sids(si), "Nf", Nf)
  !     call runfile%get(sids(si), "flog", flog)
  !     if (flog) then
  !       f = logspace(f0, f1, Nf)
  !     else
  !       f = linspace(f0, f1, Nf)
  !     end if

  !     ! number of harmonics
  !     call runfile%get(sids(si), "NH", NH)

  !     ! number of time evaluation points
  !     call runfile%get(sids(si), "Nt", Nt)

  !     ! solve steady-state
  !     gummel_restart = .true.
  !     gummel_once    = .false.
  !     gummel_enabled = .true.

  !     call ss%init(dev%sys_full)
  !     call ss%input_newton_params(runfile, "full newton params")
  !     call ss%input_var_params(runfile, "full newton params")
  !     call ss%run(input = input, gummel = gummel)

  !     ! run harmonic balance
  !     call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
  !     call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

  !     deallocate (f)
  !   end do
  ! end subroutine



  ! subroutine solve_responsivity()
  !   integer                :: i, si, ict, isrc, idrn, Nf, NH, Nt, ofunit
  !   integer,   allocatable :: sids(:)
  !   logical                :: log
  !   real                   :: VA, f0, f1, power, curr, resp
  !   real,      allocatable :: f(:), c(:,:), s(:,:)
  !   type(string)           :: source, drain, output

  !   type(harmonic_src)     :: input
  !   type(harmonic_balance) :: hb
  !   type(steady_state)     :: ss
  !   type(newton_opt)       :: opt_hb

  !   call runfile%get_sections("responsivity", sids)
  !   do si = 1, size(sids)
  !     ! get DC voltages
  !     do ict = 1, dev%par%nct
  !       call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
  !     end do

  !     ! get source and drain contact ids
  !     call runfile%get(sids(si), "source", source)
  !     call runfile%get(sids(si), "drain",  drain )
  !     isrc = dev%par%contact_map%get(source)
  !     idrn = dev%par%contact_map%get(drain)

  !     ! get amplitude
  !     call runfile%get(sids(si), "VA", VA)

  !     ! get frequency
  !     call runfile%get(sids(si), "f0", f0)
  !     call runfile%get(sids(si), "f1", f1)
  !     call runfile%get(sids(si), "Nf", Nf)
  !     f = linspace(f0, f1, Nf)

  !     ! number of harmonics
  !     call runfile%get(sids(si), "NH", NH)

  !     ! number of time evaluation points
  !     call runfile%get(sids(si), "Nt", Nt)

  !     ! output file
  !     call runfile%get(sids(si), "output", output)

  !     ! init input source
  !     allocate (c(dev%par%nct,0:1))
  !     allocate (s(dev%par%nct,1))
  !     c(   :,0) = [(dev%volt(ict)%x, ict = 1, dev%par%nct)]
  !     c(   :,1) = 0
  !     s(   :,1) = 0
  !     s(isrc,1) = VA
  !     call input%init(0.0, c, s)

  !     ! solve steady-state
  !     gummel_restart = .true.
  !     gummel_once    = .false.
  !     gummel_enabled = .true.

  !     call ss%init(dev%sys_full)
  !     call ss%input_newton_params(runfile, "full newton params")
  !     call ss%input_var_params(runfile, "full newton params")
  !     call ss%run(input = input, gummel = gummel)

  !     call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
  !     call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

  !     open (newunit = ofunit, file = output%s, status = "replace", action = "write")
  !     do i = 1, Nf
  !       call hb%select_harmonic(1, i, sine = .true.)
  !       power = 0.5 * VA * dev%curr(isrc)%x
  !       call hb%select_harmonic(0, i)
  !       curr  = dev%curr(idrn)%x
  !       resp = curr / power
  !       write (ofunit, "(4ES25.16E3)") denorm(f(i), "Hz"), denorm(power, "W/um"), denorm(dev%curr(idrn)%x, "A/um"), denorm(resp, "A/W")
  !     end do
  !     close (ofunit)

  !     deallocate (f, c, s)
  !   end do
  ! end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: it, ci, min_it, max_it, si, imax_dd, ibl_dd(1), iloc_dd
    integer, allocatable :: idx_dd(:)
    logical            :: log
    real               :: error, err_pot, err_iref(2), atol
    real, allocatable  :: pot0(:), iref0(:,:), x_before(:), x_after(:)
    type(steady_state) :: ss_dd(2)

    if (.not. gummel_enabled) return
    if (gummel_once) gummel_enabled = .false.

    allocate (pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n,2))

    if (gummel_restart) then
      gummel_restart = .false.
      ! approximate imrefs
      do ci = dev%par%ci0, dev%par%ci1
        call approx_imref(dev%par, dev%iref(ci), dev%volt)
      end do

      ! approximate potential
      call approx_potential(dev%par, dev%pot, dev%iref)

      ! initialize density from iref + potential (fixes zero-density at Schottky vertices)
      call dev%sys_nlpe%eval()

      ! add beam-generation excess: Δn ≈ G * τ_SRH at each vertex
      if (dev%par%has_beam_gen) then
        do ci = dev%par%ci0, dev%par%ci1
          call dev%calc_bgen(ci)%eval()
        end do
        block
          integer :: i_bv
          integer, allocatable :: idx_bv(:)
          real :: G_bv, tau_bv, excess_bv

          allocate(idx_bv(dev%par%g%idx_dim))
          do i_bv = 1, dev%par%transport_vct(0)%n
            idx_bv = dev%par%transport_vct(0)%get_idx(i_bv)
            G_bv = dev%bgen(CR_ELEC)%get(idx_bv)
            if (G_bv > 0.0) then
              ! add excess electrons: Δn = G * τ_n
              tau_bv = dev%par%srh_tau(CR_ELEC)%get(idx_bv)
              excess_bv = G_bv * tau_bv
              call dev%dens(CR_ELEC)%set(idx_bv, dev%dens(CR_ELEC)%get(idx_bv) + excess_bv)
              ! add excess holes: Δp = G * τ_p
              tau_bv = dev%par%srh_tau(CR_HOLE)%get(idx_bv)
              excess_bv = G_bv * tau_bv
              call dev%dens(CR_HOLE)%set(idx_bv, dev%dens(CR_HOLE)%get(idx_bv) + excess_bv)
            end if
          end do
          deallocate(idx_bv)
        end block
        ! update imrefs to be consistent with the new density
        do ci = dev%par%ci0, dev%par%ci1
          call dev%calc_iref(ci)%eval()
        end do

        ! diagnostic
        block
          integer :: ix_b, iy_b, idx_b(dev%par%g%idx_dim)
          ix_b = bin_search(dev%par%g1D(1)%x, dev%par%reg_beam(1)%beam_x)
          iy_b = bin_search(dev%par%g1D(2)%x, dev%beam_pos%x)
          idx_b = [ix_b, iy_b]
          print "(A,I4,A,I4,A,ES10.3,A,ES10.3,A,ES10.3,A)", &
            & "  init+beam beam(", ix_b, ",", iy_b, &
            & ") pot=", denorm(dev%pot%get(idx_b), 'V'), &
            & " V  n=", denorm(dev%dens(CR_ELEC)%get(idx_b), '1/cm^3'), &
            & "  p=", denorm(dev%dens(CR_HOLE)%get(idx_b), '1/cm^3'), " cm^-3"
        end block
      end if
    end if

    call runfile%get_section("dd params", si)
    do ci = dev%par%ci0, dev%par%ci1
      call ss_dd(ci)%init(dev%sys_dd(ci))
      call ss_dd(ci)%set_params(runfile%sections%d(si))
      ss_dd(ci)%msg = CR_NAME(ci) // "DD: "
    end do

    ! get gummel params
    call runfile%get("gummel params", "log", log)
    call runfile%get("gummel params", "atol", atol)
    call runfile%get("gummel params", "min_it", min_it)
    call runfile%get("gummel params", "max_it", max_it)

    ! gummel iteration
    it = 0
    error = huge(1.0)
    do while (((error > atol) .and. (it < max_it)) .or. (it < min_it))
      it = it + 1

      ! After first iteration, don't restart for subsequent beam positions
      ! (reuse converged solution as initial guess)
      gummel_restart = .false.

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! diagnostic: potential + density at beam point after NLPE
      allocate(idx_dd(dev%par%g%idx_dim))
      idx_dd = dev%par%transport_vct(0)%get_idx(1)  ! just to get idx_dim
      ! find beam grid point
      block
        integer :: ix_b, iy_b
        ix_b = bin_search(dev%par%g1D(1)%x, dev%par%reg_beam(1)%beam_x)
        iy_b = bin_search(dev%par%g1D(2)%x, dev%beam_pos%x)
        idx_dd = [ix_b, iy_b]
        print "(A,I3,A,I4,A,I4,A,ES10.3,A,ES10.3,A,ES10.3,A)", &
          & "  after NLPE gum", it, " beam(", ix_b, ",", iy_b, &
          & ") pot=", denorm(dev%pot%get(idx_dd), 'V'), &
          & " V  n=", denorm(dev%dens(1)%get(idx_dd), '1/cm^3'), &
          & "  p=", denorm(dev%dens(2)%get(idx_dd), '1/cm^3'), " cm^-3"
      end block

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()

        ! snapshot density before DD
        x_before = dev%sys_dd(ci)%get_x()

        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))

        ! DD oscillation diagnostic: find DOF with largest density change
        x_after = dev%sys_dd(ci)%get_x()
        imax_dd = maxloc(abs(x_after - x_before), dim=1)
        ibl_dd = findloc(dev%sys_dd(ci)%i1 >= imax_dd, .true.)
        iloc_dd = imax_dd - dev%sys_dd(ci)%i0(ibl_dd(1)) + 1
        if (ibl_dd(1) == 1) then
          idx_dd = dev%par%transport_vct(0)%get_idx(iloc_dd)
        else
          idx_dd = dev%par%transport_vct(ibl_dd(1) - 1)%get_idx(iloc_dd)
        end if
        print "(A,I0,A,I3,A,I4,A,I4,A,ES10.3,A,ES10.3,A,ES10.3,A)", &
          & "  DD", ci, " gum", it, " worst(", idx_dd(1), ",", idx_dd(2), &
          & ") before=", denorm(x_before(imax_dd), '1/cm^3'), &
          & " after=", denorm(x_after(imax_dd), '1/cm^3'), &
          & " delta=", denorm(x_after(imax_dd) - x_before(imax_dd), '1/cm^3'), " cm^-3"
        deallocate(x_before, x_after)
      end do
      deallocate(idx_dd)

      ! log
      if (log) then
        print "(A,I6,ES25.16E3)", "Gummel: ", it, denorm(error, "V")
      end if
    end do

    if ((it > max_it) .and. (error > atol)) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if
  end subroutine

  subroutine solve_nlpe()
    !! solve non-linear poisson equation

    integer                         :: it, nx, min_it, max_it, si
    logical                         :: status
    real                            :: err, res0, res1, damping, dx0, atol, dx_lim
    real,               allocatable :: x0(:), f(:), dx(:)
    type(block_real),   pointer     :: dfdx
    type(string)                    :: solver_name
    type(input_section)             :: solver_params, solver_params_tmp
    class(solver_real), allocatable :: solver

    ! memory
    nx = dev%sys_nlpe%n
    allocate (x0(nx), f(nx), dx(nx), source = 0.0)

    ! init iteration params
    it   = 0
    err  = huge(err)
    res0 = huge(res0)
    damping = 1.0

    ! initial approximation
    x0 = dev%sys_nlpe%get_x()

    ! get nlpe params
    call runfile%get_section("nlpe params", si)
    call runfile%get(si, "atol", atol)
    call runfile%get(si, "dx_lim", dx_lim)
    call runfile%get(si, "min_it", min_it)
    call runfile%get(si, "max_it", max_it)

    ! solver name
    call runfile%get(si, "solver", solver_name, status = status)
    if (.not. status) solver_name%s = "pardiso"

    ! solver parameters
    solver_params = default_solver_params(solver_name%s)
    call runfile%sections%d(si)%get_subsection(solver_name%s, solver_params_tmp)
    if (solver_params_tmp%params%n > 0) call solver_params%set_subsection(solver_params_tmp)

    ! init solver
    call init_solver_real(solver_name%s, solver_params, solver)

    ! newton iteration
    do while (((err > atol) .and. (it <= max_it)) .or. (it < min_it))
      it = it + 1

      ! evaluate system
      call dev%sys_nlpe%eval(f = f, df = dfdx)
      res1 = dot_product(f, f)
      write (*, "(A,I6,2ES25.16E3)", advance = "no") "NLPE: ", it, res1, damping

      ! repeat step with 0.5*dx if residual gets larger
      if (res1 > res0) then
        it = it - 1
        damping = damping * 0.5
        call dev%sys_nlpe%set_x(x0 + damping * dx)
        print *
        cycle
      end if

      ! accept step
      x0      = dev%sys_nlpe%get_x()
      res0    = res1
      damping = 1.0

      ! solve
      call solver%factorize(dfdx)
      call solver%solve(-f, dx)

      ! limit update
      dx0 = maxval(abs(dx) / dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! absolute error
      err = maxval(abs(dx))

      ! update variables
      call dev%sys_nlpe%set_x(x0 + dx)

      write (*, "(ES25.16E3)") denorm(err, "V")
    end do

    call solver%destruct()
  end subroutine

  subroutine eval_efield()
    !! Evaluate electric field components before output
    !! (calc_efield is not in the Newton eval order because nothing depends on it)
    integer :: dir

    do dir = 1, dev%par%g%dim
      call dev%calc_efield(dir)%eval()
    end do
  end subroutine

end program stem_ebic
