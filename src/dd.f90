program dd

  use approx_m,           only: approx_imref, approx_potential
  use cl_options_m,       only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,           only: dev
  use error_m,            only: program_error
  use harmonic_balance_m, only: harmonic_balance
  use input_m,            only: input_file
  use input_src_m,        only: polygon_src, harmonic_src
  use math_m,             only: linspace, logspace, PI
  use matrix_m,           only: block_real
  use newton_m,           only: newton_opt
  use normalization_m,    only: init_normconst, norm, denorm
  use semiconductor_m,    only: CR_NAME
  use small_signal_m,     only: small_signal
  use steady_state_m,     only: steady_state
  use storage_m,          only: storage, STORAGE_WRITE, DYNAMIC_EXT, DYNAMIC_APP
  use string_m,           only: string, new_string
  use transient_m,        only: transient, TRANS_TRBDF2
  use util_m,             only: get_hostname

  use current_density_m, only: time

  implicit none

  logical          :: gummel_restart, gummel_once, gummel_enabled
  real             :: temperature
  type(input_file) :: runfile

  print "(A)", "Start simulation on " // get_hostname()

  ! parse command line arguments
  call command_line()

  ! solve
  call solve_steady_state()
  call solve_small_signal()
  call solve_harmonic_balance()
  call solve_transient()
  call solve_responsivity()

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

  subroutine load_iteration_params(section, n, opt)
    character(*),     intent(in)  :: section
      !! section name in input file
    integer,          intent(in)  :: n
      !! system size
    type(newton_opt), intent(out) :: opt
      !! output newton options

    integer :: min_it, max_it
    logical :: log
    real    :: rtol, atol, lim

    call runfile%get(section, "min_it", min_it)
    call runfile%get(section, "max_it", max_it)
    call runfile%get(section, "rtol",   rtol  )
    call runfile%get(section, "atol",   atol  )
    call runfile%get(section, "lim",    lim   )
    call runfile%get(section, "log",    log   )

    call opt%init(n, atol = atol, rtol = rtol, dx_lim = lim, min_it = min_it, max_it = max_it, log = log)
  end subroutine

  subroutine voltage_input_ss(sid, t_inp, V, t_sim)
    !! read in voltage configuration of contacts for steady-state from runfile
    integer,           intent(in)  :: sid
      !! section id for runfile
    real, allocatable, intent(out) :: t_inp(:)
      !! pseudo time points for polygon input source
    real, allocatable, intent(out) :: V(:,:)
      !! voltages at each pseudo time point (dev%par%nct, size(t_inp))
    real, allocatable, intent(out) :: t_sim(:)
      !! pseudo time points at which steady-state simulations are performed

    integer              :: ict, ict_sweep, i
    integer, allocatable :: nsweep_ct(:), nsweep(:)
    logical              :: status
    real,    allocatable :: Vbounds(:), tmp(:)

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

    ! contact voltages
    allocate(V(dev%par%nct, size(nsweep)+1))
    do ict = 1, dev%par%nct
      call runfile%get(sid, "V_"//dev%par%contacts(ict)%name, Vbounds)
      if (ict == ict_sweep) then
        V(ict,:) = Vbounds
      else
        if (size(Vbounds) > 1) call program_error("Voltage bounds given for a contact without sweep")
        V(ict,:) = Vbounds(1)
      end if
    end do

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

  subroutine solve_steady_state()
    integer              :: si
    integer, allocatable :: sids(:)
    logical              :: log
    real,    allocatable :: V(:,:), t_inp(:), t(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    call runfile%get_sections("steady state", sids)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call runfile%get(sids(si), "name", name)

      ! input config
      call voltage_input_ss(sids(si), t_inp, V, t)
      call input%init(t_inp, V)

      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      ! solve steady-state
      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%init_output([new_string("pot"), new_string("ndens"), new_string("pdens"), new_string("ionD"), &
                         & new_string("ionA"), new_string("V_GAT"), new_string("I_DRN")], name%s // ".fbs")
      call ss%run(input = input, t_input = t, gummel = gummel)
    end do
  end subroutine

  subroutine solve_small_signal()
    integer              :: si, Nf, i
    integer, allocatable :: sids(:)
    logical              :: flog, log
    real                 :: f0, f1
    real,    allocatable :: f(:), V(:,:), t_inp(:), t(:)
    complex, allocatable :: s(:), result(:,:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(small_signal)   :: ac
    type(steady_state)   :: ss
    type(storage)        :: st

    call runfile%get_sections("small signal", sids)
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
      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%run(input = input, t_input = t, gummel = gummel)

      ! run small-signal analysis for a single working point
      call ac%init(dev%sys_full, log=.true.)
      call ac%init_output([new_string("I_GAT")], name%s // ".fbs")
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
      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%set_var_params("currents", atol=1e300)
      call ss%set_var_params("ndens", dx_lim_rel=0.2)
      call ss%set_var_params("pdens", dx_lim_rel=0.2)

      ! run small-signal analysis at each working point
      allocate (result(size(t), Nf), source = (0.0,0.0))
      do i = 1, size(t)
        print *, "steady-state step: ", i
        call ss%run(input = input, t_input = [t(i)], gummel = gummel)
        call ac%init(dev%sys_full, log = .true.)
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

  subroutine solve_harmonic_balance()
    integer                :: si, ict, Nf, NH, Nt
    integer, allocatable   :: sids(:)
    logical                :: flog, log
    real                   :: f0, f1
    real,    allocatable   :: volt(:), f(:), c(:,:), s(:,:)
    type(string)           :: name

    type(harmonic_src)     :: input
    type(harmonic_balance) :: hb
    type(steady_state)     :: ss
    type(newton_opt)       :: opt_hb

    allocate (c(dev%par%nct,0:1), s(dev%par%nct,1))

    call runfile%get_sections("harmonic balance", sids)
    do si = 1, size(sids)
      print "(A)", "harmonic balance"

      call runfile%get(sids(si), "name", name)

      ! get input source
      do ict = 1, dev%par%nct
        call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, volt)
        if (size(volt) /= 3) call program_error("3 values for voltage expected (c0, c1, s1 coefficients)")

        c(ict,0) = volt(1)
        c(ict,1) = volt(2)
        s(ict,1) = volt(3)
      end do
      call input%init(0.0, c, s)

      ! get frequency
      call runfile%get(sids(si), "f0", f0)
      call runfile%get(sids(si), "f1", f1)
      call runfile%get(sids(si), "Nf", Nf)
      call runfile%get(sids(si), "flog", flog)
      if (flog) then
        f = logspace(f0, f1, Nf)
      else
        f = linspace(f0, f1, Nf)
      end if

      ! number of harmonics
      call runfile%get(sids(si), "NH", NH)

      ! number of time evaluation points
      call runfile%get(sids(si), "Nt", Nt)

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%run(input = input, gummel = gummel)

      ! run harmonic balance
      call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
      call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

      deallocate (f)
    end do
  end subroutine

  subroutine solve_transient()
    integer              :: ict, Nt, si, start, end, rate
    integer, allocatable :: sids(:)
    logical              :: log
    real                 :: dt0
    real,    allocatable :: ti(:), Vtmp(:), Vi(:,:)
    type(string)         :: name
    type(steady_state)   :: ss
    type(transient)      :: trans
    type(polygon_src)    :: input

    call runfile%get_sections("transient", sids)
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
      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%run(input = input, gummel = gummel)

      ! run transient simulation
      call runfile%get("transient params", "log", log)
      call trans%init(dev%sys_full, log = log, msg = "Transient: ")
      call trans%set_ode_params(method = TRANS_TRBDF2, adaptive = .true., eabs = 0.01)
      call trans%input_newton_params(runfile, "transient params")
      call trans%input_var_params(runfile, "transient params")
      call trans%init_output([new_string("pot"), new_string("ndens"), new_string("pdens"), new_string("ionD"), &
                            & new_string("ionA"), new_string("V_GAT"), new_string("I_DRN"), new_string("I_GAT"), &
                            & new_string("I_SRC"), new_string("I_BLK"), new_string("ncdensx"), new_string("pcdensx"), &
                            & new_string("ncdensy"), new_string("pcdensy")], name%s // ".fbs")
      call system_clock(start, rate)
      call trans%run(ti, dt0=dt0, input=input, start_steady_state=.true.)
      call system_clock(end)
      print *, "Transient time: ", real(end-start)/real(rate)
      print *, "Time used for current calculation: ", time
      print *, "which is ", time/(real(end-start)/real(rate)) * 100, "% of the total runtime"
    end do
  end subroutine

  subroutine solve_responsivity()
    integer                :: i, si, ict, isrc, idrn, Nf, NH, Nt, ofunit
    integer,   allocatable :: sids(:)
    logical                :: log
    real                   :: VA, f0, f1, power, curr, resp
    real,      allocatable :: f(:), c(:,:), s(:,:)
    type(string)           :: source, drain, output

    type(harmonic_src)     :: input
    type(harmonic_balance) :: hb
    type(steady_state)     :: ss
    type(newton_opt)       :: opt_hb

    call runfile%get_sections("responsivity", sids)
    do si = 1, size(sids)
      ! get DC voltages
      do ict = 1, dev%par%nct
        call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do

      ! get source and drain contact ids
      call runfile%get(sids(si), "source", source)
      call runfile%get(sids(si), "drain",  drain )
      isrc = dev%par%contact_map%get(source)
      idrn = dev%par%contact_map%get(drain)

      ! get amplitude
      call runfile%get(sids(si), "VA", VA)

      ! get frequency
      call runfile%get(sids(si), "f0", f0)
      call runfile%get(sids(si), "f1", f1)
      call runfile%get(sids(si), "Nf", Nf)
      f = linspace(f0, f1, Nf)

      ! number of harmonics
      call runfile%get(sids(si), "NH", NH)

      ! number of time evaluation points
      call runfile%get(sids(si), "Nt", Nt)

      ! output file
      call runfile%get(sids(si), "output", output)

      ! init input source
      allocate (c(dev%par%nct,0:1))
      allocate (s(dev%par%nct,1))
      c(   :,0) = [(dev%volt(ict)%x, ict = 1, dev%par%nct)]
      c(   :,1) = 0
      s(   :,1) = 0
      s(isrc,1) = VA
      call input%init(0.0, c, s)

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      call ss%run(input = input, gummel = gummel)

      call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
      call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

      open (newunit = ofunit, file = output%s, status = "replace", action = "write")
      do i = 1, Nf
        call hb%select_harmonic(1, i, sine = .true.)
        power = 0.5 * VA * dev%curr(isrc)%x
        call hb%select_harmonic(0, i)
        curr  = dev%curr(idrn)%x
        resp = curr / power
        write (ofunit, "(4ES25.16E3)") denorm(f(i), "Hz"), denorm(power, "W/um"), denorm(dev%curr(idrn)%x, "A/um"), denorm(resp, "A/W")
      end do
      close (ofunit)

      deallocate (f, c, s)
    end do
  end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: it, ci, min_it, max_it
    logical            :: log
    real               :: error, err_pot, err_iref(2), atol
    real, allocatable  :: pot0(:), iref0(:,:)
    type(steady_state) :: ss_dd(2)

    if (.not. gummel_enabled) return
    if (gummel_once) gummel_enabled = .false.

    allocate (pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n,2))

    if (gummel_restart) then
      ! approximate imrefs
      do ci = dev%par%ci0, dev%par%ci1
        call approx_imref(dev%par, dev%iref(ci), dev%volt)
      end do

      ! approximate potential
      call approx_potential(dev%par, dev%pot, dev%iref)
    end if

    do ci = dev%par%ci0, dev%par%ci1
      call runfile%get("dd params", "log", log)
      call ss_dd(ci)%init(dev%sys_dd(ci), log = log, msg = CR_NAME(ci) // "DD: ")
      call ss_dd(ci)%input_newton_params(runfile, "dd params")
      call ss_dd(ci)%input_var_params(runfile, "dd params")
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

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()

        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

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

    integer                   :: it, nx, min_it, max_it
    real                      :: err, res0, res1, damping, dx0, atol, dx_lim
    real, allocatable         :: x0(:), f(:), dx(:)
    type(block_real), pointer :: dfdx

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
    call runfile%get("nlpe params", "atol", atol)
    call runfile%get("nlpe params", "dx_lim", dx_lim)
    call runfile%get("nlpe params", "min_it", min_it)
    call runfile%get("nlpe params", "max_it", max_it)

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
      call dfdx%factorize()
      call dfdx%solve_vec(-f, dx)
      call dfdx%reset(only_factorization = .true.)

      ! limit update
      dx0 = maxval(abs(dx) / dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! absolute error
      err = maxval(abs(dx))

      ! update variables
      call dev%sys_nlpe%set_x(x0 + dx)

      write (*, "(ES25.16E3)") denorm(err, "V")
    end do
  end subroutine

end program
