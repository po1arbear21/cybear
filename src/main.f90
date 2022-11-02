program main

  use approx_m,           only: approx_imref, approx_potential
  use device_m,           only: dev
  use error_m,            only: program_error
  use harmonic_balance_m, only: harmonic_balance
  use input_m,            only: input_file
  use input_src_m,        only: const_src, polygon_src, harmonic_src
  use math_m,             only: linspace, logspace, PI
  use newton_m,           only: newton_opt
  use normalization_m,    only: init_normconst, norm, denorm
  use output_file_m,      only: output_file
  use small_signal_m,     only: small_signal
  use steady_state_m,     only: steady_state
  use string_m,           only: string

  implicit none

  character(:), allocatable :: filename
  logical                   :: gummel_restart, gummel_once, gummel_enabled
  real                      :: T
  type(string)              :: dev_filename, ofilename
  type(input_file)          :: file
  type(newton_opt)          :: opt_nlpe, opt_dd, opt_gum, opt_full
  type(output_file)         :: ofile

  ! get input filename from command line arguments
  call parse_command_line_args()

  ! load input file
  call file%init(filename)

  ! get temperature and initialize normalization
  call file%get("", "temperature", T, normalize = .false.)
  call init_normconst(T)

  ! get device filename and initialize device
  call file%get("", "device", dev_filename)
  call dev%init(dev_filename%s)

  ! output file
  call file%get("", "output", ofilename)
  call ofile%init(ofilename%s)

  ! output grids
  call dev%par%gx%output(ofile, unit = "nm")
  call dev%par%gy%output(ofile, unit = "nm")
  call dev%par%g%output( ofile, unit = "nm")

  ! output variable data
  call dev%sys_full%output_info(ofile)

  ! get iteration options
  call load_iteration_params("nlpe params",        dev%sys_nlpe%n,  opt_nlpe)
  call load_iteration_params("dd params",          dev%sys_dd(1)%n, opt_dd  )
  call load_iteration_params("gummel params",      1,               opt_gum )
  call load_iteration_params("full newton params", dev%sys_full%n,  opt_full)

  ! solve
  call solve_steady_state()
  call solve_small_signal()
  call solve_harmonic_balance()
  call solve_responsivity()

  ! save and close output file
  call ofile%close()

contains

  subroutine parse_command_line_args()
    integer :: nargs, length, status

    ! get number of command line arguments
    nargs = command_argument_count()
    if (nargs < 1) call program_error("Command filename expected!")
    if (nargs > 1) print *, "Warning: Command-line arguments after first one are ignored"

    ! read device filename
    allocate (character(32) :: filename)
    do while (.true.)
      call get_command_argument(1, value = filename, length = length, status = status)
      if (status == -1) then
        deallocate (filename)
        allocate (character(2*len(filename)) :: filename)
      elseif (status == 0) then
        filename = filename(1:length)
        exit
      else
        call program_error("Error at reading command-line argument")
      end if
    end do
  end subroutine

  subroutine load_iteration_params(section_name, n, opt)
    character(*),     intent(in)  :: section_name
      !! section name in input file
    integer,          intent(in)  :: n
      !! system size
    type(newton_opt), intent(out) :: opt
      !! output newton options

    integer :: max_it
    logical :: log
    real    :: rtol, atol, lim

    call file%get(section_name, "max_it", max_it)
    call file%get(section_name, "rtol",   rtol  )
    call file%get(section_name, "atol",   atol  )
    call file%get(section_name, "lim",    lim   )
    call file%get(section_name, "log",    log   )

    call opt%init(n, atol = atol, rtol = rtol, dx_lim = lim, max_it = max_it, log = log)
  end subroutine

  subroutine solve_steady_state()
    integer              :: si, ict, ict_sweep, nsweep
    integer, allocatable :: sids(:)
    logical              :: status
    real,    allocatable :: Vi(:,:), Vbounds(:), t(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    allocate (Vi(size(dev%par%contacts),2))

    call file%get_sections("steady state", sids)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call file%get(sids(si), "name", name)

      ! voltage input source
      ict_sweep = 0
      do ict = 1, size(dev%par%contacts)
        call file%get(sids(si), "N_"//dev%par%contacts(ict)%name, nsweep, status = status)
        if (status) then
          if (ict_sweep /= 0) call program_error("only one voltage sweep per steady-state section allowed")
          ict_sweep = ict
          call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vbounds)
          Vi(ict,:) = Vbounds
        else
          call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vi(ict,1))
          Vi(ict,2) = Vi(ict,1)
        end if
      end do
      call input%init([0.0, 1.0], Vi)

      if (ict_sweep == 0) then
        t = [ 0.0 ]
      else
        t = linspace(0.0, 1.0, nsweep)
      end if

      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = t, gum = gummel, ofile = ofile, oname = name%s)
    end do
  end subroutine

  subroutine solve_small_signal()
    integer                   :: si, ict, Nf
    integer,      allocatable :: sids(:)
    logical                   :: flog
    real                      :: f0, f1
    real,         allocatable :: f(:)
    complex,      allocatable :: s(:)
    type(string)              :: name
    type(const_src)           :: input
    type(small_signal)        :: ac
    type(steady_state)        :: ss

    call file%get_sections("small signal", sids)
    do si = 1, size(sids)
      print "(A)", "small signal"

      call file%get(sids(si), "name", name)

      ! get DC voltages
      do ict = 1, size(dev%par%contacts)
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do
      call input%init([(dev%volt(ict)%x, ict = 1, size(dev%par%contacts))])

      ! get frequency
      call file%get(sids(si), "f0", f0)
      call file%get(sids(si), "f1", f1)
      call file%get(sids(si), "Nf", Nf)
      call file%get(sids(si), "flog", flog)
      if (flog) then
        f = logspace(f0, f1, Nf)
      else
        f = linspace(f0, f1, Nf)
      end if
      s = 2 * PI * (0.0, 1.0) * f

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, gum = gummel)

      ! run small-signal analysis
      call ac%run_analysis(dev%sys_full, s, ofile = ofile, oname = name%s)

      deallocate (f, s)
    end do
    print *
  end subroutine

  subroutine solve_harmonic_balance()
    integer                :: si, ict, Nf, NH, Nt
    integer, allocatable   :: sids(:)
    logical                :: flog
    real                   :: f0, f1
    real,    allocatable   :: volt(:), f(:), c(:,:), s(:,:)
    type(string)           :: name

    type(harmonic_src)     :: input
    type(harmonic_balance) :: hb
    type(steady_state)     :: ss
    type(newton_opt)       :: opt_hb

    allocate (c(size(dev%par%contacts),0:1), s(size(dev%par%contacts),1))

    call file%get_sections("harmonic balance", sids)
    do si = 1, size(sids)
      print "(A)", "harmonic balance"

      call file%get(sids(si), "name", name)

      ! get input source
      do ict = 1, size(dev%par%contacts)
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, volt)
        if (size(volt) /= 3) call program_error("3 values for voltage expected (c0, c1, s1 coefficients)")

        c(ict,0) = volt(1)
        c(ict,1) = volt(2)
        s(ict,1) = volt(3)
      end do
      call input%init(0.0, c, s)

      ! get frequency
      call file%get(sids(si), "f0", f0)
      call file%get(sids(si), "f1", f1)
      call file%get(sids(si), "Nf", Nf)
      call file%get(sids(si), "flog", flog)
      if (flog) then
        f = logspace(f0, f1, Nf)
      else
        f = linspace(f0, f1, Nf)
      end if

      ! number of harmonics
      call file%get(sids(si), "NH", NH)

      ! number of time evaluation points
      call file%get(sids(si), "Nt", Nt)

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, gum = gummel)

      ! run harmonic balance
      call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
      call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

      deallocate (f)
    end do
  end subroutine

  subroutine solve_responsivity()
    integer                   :: i, si, ict, isrc, idrn, Nf, NH, Nt, ofunit
    integer,      allocatable :: sids(:)
    real                      :: VA, f0, f1, power, curr, resp
    real,         allocatable :: f(:), c(:,:), s(:,:)
    type(string)              :: source, drain, ofile

    type(harmonic_src)        :: input
    type(harmonic_balance)    :: hb
    type(steady_state)        :: ss
    type(newton_opt)          :: opt_hb

    call file%get_sections("responsivity", sids)
    do si = 1, size(sids)
      ! get DC voltages
      do ict = 1, size(dev%par%contacts)
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do
      ! call input%init([(dev%volt(ict)%x, ict = 1, size(dev%par%contacts))])

      ! get source and drain contact ids
      call file%get(sids(si), "source", source)
      call file%get(sids(si), "drain",  drain )
      isrc = dev%par%contact_map%get(source)
      idrn = dev%par%contact_map%get(drain)

      ! get amplitude
      call file%get(sids(si), "VA", VA)

      ! get frequency
      call file%get(sids(si), "f0", f0)
      call file%get(sids(si), "f1", f1)
      call file%get(sids(si), "Nf", Nf)
      f = linspace(f0, f1, Nf)

      ! number of harmonics
      call file%get(sids(si), "NH", NH)

      ! number of time evaluation points
      call file%get(sids(si), "Nt", Nt)

      ! output file
      call file%get(sids(si), "output", ofile)

      ! init input source
      allocate (c(size(dev%par%contacts),0:1))
      allocate (s(size(dev%par%contacts),1))
      c(   :,0) = [(dev%volt(ict)%x, ict = 1, size(dev%par%contacts))]
      c(   :,1) = 0
      s(   :,1) = 0
      s(isrc,1) = VA
      call input%init(0.0, c, s)

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, gum = gummel)

      call load_iteration_params("harmonic balance params", dev%sys_full%n*(1+2*NH), opt_hb)
      call hb%run(dev%sys_full, NH, f, input, nopt = opt_hb, nt = Nt)

      open (newunit = ofunit, file = ofile%s, status = "replace", action = "write")
      do i = 1, Nf
        call hb%select_harmonic(1, i, sine = .true.)
        power = 0.5 * VA * dev%curr(isrc)%x
        call hb%select_harmonic(0, i)
        curr  = dev%curr(idrn)%x
        resp = curr / power
        write (ofunit, "(4ES24.16)") denorm(f(i), "Hz"), denorm(power, "W/um"), denorm(dev%curr(idrn)%x, "A/um"), denorm(resp, "A/W")
      end do

      close (ofunit)

      deallocate (f, c, s)
    end do
  end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: i, ci
    real               :: error, err_pot, err_iref(2)
    real, allocatable  :: pot0(:), iref0(:,:)
    type(steady_state) :: ss_nlpe, ss_dd(2)

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

    ! gummel iteration
    i = 0
    error = huge(1.0)
    do while ((error > opt_gum%atol(1)) .and. (i < opt_gum%max_it))
      i = i + 1

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call ss_nlpe%run(dev%sys_nlpe, nopt = opt_nlpe)
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()
        call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

      ! log
      if (opt_gum%log) then
        print "(A,I6,ES24.16)", "Gummel: ", i, denorm(error, "V")
      end if
    end do

    if ((i > opt_gum%max_it) .and. (error > opt_gum%atol(1))) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if
  end subroutine

end program
