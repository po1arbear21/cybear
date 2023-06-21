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
  ! integer                   :: idx_dir
  logical                   :: gummel_restart, gummel_once, gummel_enabled
  real                      :: T
  type(string)              :: dev_filename, ofilename
  type(input_file)          :: file
  type(newton_opt)          :: opt_nlpe, opt_dd(2), opt_gum, opt_full
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

  ! ! output grids
  ! do idx_dir = 1, dev%par%g%idx_dim
  !   call dev%par%g1D(idx_dir)%output(ofile, unit = "nm")
  ! end do

  ! output variable data
  call dev%sys_full%output_info(ofile)

  ! get iteration options
  call load_iteration_params("nlpe params",        dev%sys_nlpe%n,  opt_nlpe)
  call load_iteration_params("dd params",          dev%sys_dd(1)%n, opt_dd(1))
  call load_iteration_params("dd params",          dev%sys_dd(2)%n, opt_dd(2))
  call load_iteration_params("gummel params",      1,               opt_gum )
  call load_iteration_params("full newton params", dev%sys_full%n,  opt_full)

  opt_nlpe%msg  = "NLPE: "
  opt_dd(1)%msg = "NDD: "
  opt_dd(2)%msg = "PDD: "
  opt_gum%msg   = "Gummel: "
  opt_full%msg  = "Newton: "

  block
    integer :: ci, idens, iion, itab, ibl, i0, i1
    do ci = dev%par%ci0, dev%par%ci1
      idens = dev%sys_full%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_full%res2block(idens)%d)
        ibl = dev%sys_full%res2block(idens)%d(itab)
        i0 = dev%sys_full%i0(ibl)
        i1 = dev%sys_full%i1(ibl)
        opt_full%xmin(i0:i1) = norm(1e-100, "1/cm^3")
      end do

      idens = dev%sys_dd(ci)%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_dd(ci)%res2block(idens)%d)
        ibl = dev%sys_dd(ci)%res2block(idens)%d(itab)
        i0 = dev%sys_dd(ci)%i0(ibl)
        i1 = dev%sys_dd(ci)%i1(ibl)
        opt_dd(ci)%xmin(i0:i1) = norm(1e-100, "1/cm^3")
      end do

      iion = dev%sys_full%search_main_var(dev%ion(ci)%name)
      do itab = 1, size(dev%sys_full%res2block(iion)%d)
        ibl = dev%sys_full%res2block(iion)%d(itab)
        i0 = dev%sys_full%i0(ibl)
        i1 = dev%sys_full%i1(ibl)
        opt_full%xmin(i0:i1) = 0.0
        opt_full%xmax(i0:i1) = 1.0
      end do

      iion = dev%sys_dd(ci)%search_main_var(dev%ion(ci)%name)
      do itab = 1, size(dev%sys_dd(ci)%res2block(iion)%d)
        ibl = dev%sys_dd(ci)%res2block(iion)%d(itab)
        i0 = dev%sys_dd(ci)%i0(ibl)
        i1 = dev%sys_dd(ci)%i1(ibl)
        opt_dd(ci)%xmin(i0:i1) = 0.0
        opt_dd(ci)%xmax(i0:i1) = 1.0
      end do

opt_dd(ci)%error_if_not_converged = .false.
    end do
  end block


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

    allocate (Vi(dev%par%nct,2))

    call file%get_sections("steady state", sids)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call file%get(sids(si), "name", name)

      ! voltage input source
      ict_sweep = 0
      do ict = 1, dev%par%nct
        call file%get(sids(si), "N_"//dev%par%contacts(ict)%name, nsweep, status = status)
        if (status) status = (nsweep > 1)
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

block
  use normalization_m
  use grid_m

  ! integer :: i, funit

  ! ! call dev%calc_cdens(1,1)%test()
  ! ! call dev%calc_cdens(1,2)%test()
  ! ! call dev%calc_cdens(2,1)%test()
  ! ! call dev%calc_cdens(2,2)%test()

  ! open (newunit = funit, file = "IV.csv", status = "replace", action = "write")
  ! do i = 1, nsweep
  !   call ss%select(i)
  !   write (funit, "(2ES25.16E3)") denorm(dev%volt(2)%x, "V"), denorm(dev%curr(2)%x, "A")
  ! end do
  ! close (funit)
  ! stop

  integer :: i, j, i1, i2, funit, idx(2), idx1(2), idx2(2)
  ! ! real    :: p1(3), p2(3)

  call ss%select(1)

  ! open (newunit = funit, file = "pot_x.csv", status = "replace", action = "write")
  ! do i = 1, dev%par%g1D(1)%n
  !   write (funit, "(2ES25.16E3)") denorm(dev%par%g1D(1)%x(i), "nm"), denorm(dev%pot%x1(i), "V")
  ! end do
  ! close (funit)

  ! open (newunit = funit, file = "pot_xy.csv", status = "replace", action = "write")
  ! do i = 1, dev%par%g1D(2)%n
  !   write (funit, "(2ES25.16E3)") denorm(dev%par%g1D(2)%x(i), "nm"), denorm(dev%pot%x2(11,i), "V")
  ! end do
  ! close (funit)

  ! open (newunit = funit, file = "pot_z_14_11.csv", status = "replace", action = "write")
  ! do i = 1, dev%par%g1D(3)%n
  !   write (funit, "(2ES25.16E3)") denorm(dev%par%g1D(3)%x(i), "nm"), denorm(dev%pot%x3(14,11,i), "V")
  ! end do
  ! close (funit)
  ! open (newunit = funit, file = "pot_z_5_11.csv", status = "replace", action = "write")
  ! do i = 1, dev%par%g1D(3)%n
  !   write (funit, "(2ES25.16E3)") denorm(dev%par%g1D(3)%x(i), "nm"), denorm(dev%pot%x3(5,11,i), "V")
  ! end do
  ! close (funit)
  ! open (newunit = funit, file = "pot_z_14_5.csv", status = "replace", action = "write")
  ! do i = 1, dev%par%g1D(3)%n
  !   write (funit, "(2ES25.16E3)") denorm(dev%par%g1D(3)%x(i), "nm"), denorm(dev%pot%x3(14,5,i), "V")
  ! end do
  ! close (funit)

  open (newunit = funit, file = "pot_xy.csv", status = "replace", action = "write")
  write (funit, "(A)") "$ DATA=CURVE3D"
  do j = 1, dev%par%g1D(2)%n; do i = 1, dev%par%g1D(1)%n
    idx1 = [i,   j]
    idx2 = [i+1, j]
    if (i < dev%par%g1D(1)%n) then
      write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
        &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
        &                         denorm(dev%pot%get(idx1), "V")
      write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
        &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
        &                         denorm(dev%pot%get(idx2), "V")
      write (funit, *)
    end if

    idx1 = [i, j  ]
    idx2 = [i, j+1]
    if (j < dev%par%g1D(2)%n) then
      write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
        &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
        &                         denorm(dev%pot%get(idx1), "V")
      write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
        &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
        &                         denorm(dev%pot%get(idx2), "V")
      write (funit, *)
    end if
  end do; end do
  write (funit, "(A)") "$ END"
  close (funit)
  ! open (newunit = funit, file = "pot_xy_20.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "$ DATA=CURVE3D"
  ! do j = 1, dev%par%g1D(2)%n; do i = 1, dev%par%g1D(1)%n
  !   idx1 = [i,   j, 20]
  !   idx2 = [i+1, j, 20]
  !   if (i < dev%par%g1D(1)%n) then
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx1), "V")
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx2), "V")
  !     write (funit, *)
  !   end if

  !   idx1 = [i, j,   20]
  !   idx2 = [i, j+1, 20]
  !   if (j < dev%par%g1D(2)%n) then
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx1), "V")
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx2), "V")
  !     write (funit, *)
  !   end if
  ! end do; end do
  ! write (funit, "(A)") "$ END"
  ! close (funit)
  ! open (newunit = funit, file = "pot_xy_30.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "$ DATA=CURVE3D"
  ! do j = 1, dev%par%g1D(2)%n; do i = 1, dev%par%g1D(1)%n
  !   idx1 = [i,   j, 30]
  !   idx2 = [i+1, j, 30]
  !   if (i < dev%par%g1D(1)%n) then
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx1), "V")
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx2), "V")
  !     write (funit, *)
  !   end if

  !   idx1 = [i, j,   30]
  !   idx2 = [i, j+1, 30]
  !   if (j < dev%par%g1D(2)%n) then
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx1(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx1(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx1), "V")
  !     write (funit, "(3ES25.16E3)") denorm(dev%par%g1D(1)%x(idx2(1)), "nm"), &
  !       &                         denorm(dev%par%g1D(2)%x(idx2(2)), "nm"), &
  !       &                         denorm(dev%pot%get(idx2), "V")
  !     write (funit, *)
  !   end if
  ! end do; end do
  ! write (funit, "(A)") "$ END"
  ! close (funit)

  ! open (newunit = funit, file = "pot_tr.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "$ DATA=CURVE3D"
  ! do i = 1, dev%par%poisson(IDX_EDGE,1)%n
  !   idx(1:1) = dev%par%poisson(IDX_EDGE,1)%get_idx(i)
  !   call dev%par%g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(1:1), 1, idx1(1:1), status)
  !   call dev%par%g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(1:1), 2, idx2(1:1), status)
  !   call dev%par%g%get_vertex(idx1(1:1), p1)
  !   call dev%par%g%get_vertex(idx2(1:1), p2)
  !   write (funit, "(3ES25.16E3)") denorm(p1(1), "nm"), denorm(p1(2), "nm"), denorm(dev%pot%get(idx1(1:1)), "V")
  !   write (funit, "(3ES25.16E3)") denorm(p2(1), "nm"), denorm(p2(2), "nm"), denorm(dev%pot%get(idx2(1:1)), "V")
  !   write (funit, *)
  ! end do
  ! close (funit)

  ! open (newunit = funit, file = "pot_tr_z_10.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "$ DATA=CURVE3D"
  ! do i = 1, dev%par%gtr%nedge
  !   idx(1) = i
  !   idx(2) = 10
  !   call dev%par%g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(1:2), 1, idx1(1:2), status)
  !   call dev%par%g%get_neighb(IDX_EDGE, 1, IDX_VERTEX, 0, idx(1:2), 2, idx2(1:2), status)
  !   call dev%par%g%get_vertex(idx1(1:2), p1)
  !   call dev%par%g%get_vertex(idx2(1:2), p2)
  !   write (funit, "(3ES25.16E3)") denorm(p1(1), "nm"), denorm(p1(2), "nm"), denorm(dev%pot%get(idx1(1:2)), "V")
  !   write (funit, "(3ES25.16E3)") denorm(p2(1), "nm"), denorm(p2(2), "nm"), denorm(dev%pot%get(idx2(1:2)), "V")
  !   write (funit, *)
  ! end do
  ! close (funit)

  ! stop
end block
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
      do ict = 1, dev%par%nct
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do
      call input%init([(dev%volt(ict)%x, ict = 1, dev%par%nct)])

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

    allocate (c(dev%par%nct,0:1), s(dev%par%nct,1))

    call file%get_sections("harmonic balance", sids)
    do si = 1, size(sids)
      print "(A)", "harmonic balance"

      call file%get(sids(si), "name", name)

      ! get input source
      do ict = 1, dev%par%nct
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
      do ict = 1, dev%par%nct
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do

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
        write (ofunit, "(4ES25.16E3)") denorm(f(i), "Hz"), denorm(power, "W/um"), denorm(dev%curr(idrn)%x, "A/um"), denorm(resp, "A/W")
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

block
  use normalization_m
  integer :: ii

  do ii = 1, dev%par%g1D(2)%n
    print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(ii), "nm"), denorm(dev%dens(1)%x2(8,ii), "1/cm^3")
  end do
  print *
  print *
end block

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()

block
  use matrix_m

  integer :: ii, funit
  real, allocatable :: f(:), x0(:), x(:), dx(:)
  type(sparse_real) :: df

  allocate (f(dev%sys_dd(ci)%n), x0(dev%sys_dd(ci)%n), x(dev%sys_dd(ci)%n), dx(dev%sys_dd(ci)%n))

  x0 = dev%sys_dd(ci)%get_x()
  x  = x0
  where(x < opt_dd(ci)%xmin) x = opt_dd(ci)%xmin
  where(x > opt_dd(ci)%xmax) x = opt_dd(ci)%xmax

  call dev%sys_dd(ci)%set_x(x)

  call dev%sys_dd(ci)%eval(f = f, df = df)

  call df%output(file = "df")
  ! call dev%calc_genrec(ci)%test(10000)
  call df%factorize()
  call df%solve_vec(f, dx)

  open (newunit = funit, file = "t", status = "replace", action = "write")
  do ii = 1, dev%sys_dd(ci)%n
    write (funit, "(I6,4ES25.16E3)") ii, f(ii), x0(ii), x(ii), dx(ii)
  end do
  close (funit)

  call dev%sys_dd(ci)%print()
  stop
end block

! if (ci == 2) then
!   block
!     use matrix_m

!     type(sparse_real) :: df

!     call dev%sys_dd(ci)%eval(df = df)

!     call dev%contin(2)%test(100)
!     call dev%ion_contin(2)%test(100)
!     call dev%calc_rec(2)%test(100)

!     call df%output("df.txt")
!     call dev%sys_dd(ci)%print()
!     stop
!   end block
! end if
        call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd(ci))

! call dev%contin(ci)%test(10000)
! call dev%ion_contin(ci)%test(10000)
! call dev%calc_genrec(ci)%test(10000)
! call dev%calc_cdens(1,ci)%test(10000)
! call dev%calc_cdens(2,ci)%test(10000)
! call dev%calc_mob(1,ci)%test(10000)
! call dev%calc_mob(2,ci)%test(10000)

        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

! block
!   use matrix_m

!   integer :: ii, funit
!   real, allocatable :: f(:), x0(:), x(:), dx(:)
!   type(sparse_real) :: df

!   allocate (f(dev%sys_dd(2)%n), x0(dev%sys_dd(2)%n), x(dev%sys_dd(2)%n), dx(dev%sys_dd(2)%n))

!   x0 = dev%sys_dd(2)%get_x()
!   x  = x0
!   where(x < opt_dd(2)%xmin) x = opt_dd(2)%xmin
!   where(x > opt_dd(2)%xmax) x = opt_dd(2)%xmax

!   call dev%sys_dd(2)%set_x(x)

!   call dev%sys_dd(2)%eval(f = f, df = df)

!   call df%output(file = "df")
!   ! call dev%calc_genrec(2)%test(10000)
!   call df%factorize()
!   call df%solve_vec(f, dx)

!   open (newunit = funit, file = "t", status = "replace", action = "write")
!   do ii = 1, dev%sys_dd(2)%n
!     write (funit, "(I6,4ES25.16E3)") ii, f(ii), x0(ii), x(ii), dx(ii)
!   end do
!   close (funit)

!   call dev%sys_dd(2)%print()
!   stop
! end block

! call dev%contin(2)%test(10000)
! call dev%ion_contin(2)%test(10000)
call dev%calc_genrec(2)%test(10000)
! call dev%calc_cdens(1,2)%test(10000)
! call dev%calc_cdens(2,2)%test(10000)
! call dev%calc_mob(1,2)%test(10000)
! call dev%calc_mob(2,2)%test(10000)
stop

block
  use normalization_m
  integer :: ii

  do ii = 1, dev%par%g1D(2)%n
    print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(ii), "nm"), denorm(dev%dens(1)%x2(8,ii), "1/cm^3")
  end do
  print *
  print *
end block
stop

      ! log
      if (opt_gum%log) then
        print "(A,I6,ES25.16E3)", "Gummel: ", i, denorm(error, "V")
      end if
    end do

    if ((i > opt_gum%max_it) .and. (error > opt_gum%atol(1))) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if
  end subroutine

end program
