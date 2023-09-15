program main

  use approx_m,           only: approx_imref, approx_potential
  use contact_m,          only: CT_OHMIC
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
  use transient_m,        only: backward_euler, tr_bdf2

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
  call dev%init(dev_filename%s, T)

! block
!   use normalization_m
!   use grid_m
!   integer :: ii, idx1(1), idx2(1), funit1, funit2
!   logical, allocatable :: oxide(:), silicon(:)
!   real    :: p1(2), p2(2), x0, x1, y0, y1, Lx, Ly
!   real, allocatable :: x(:,:), y(:,:)

!   allocate (oxide(dev%par%gtr%nedge), silicon(dev%par%gtr%nedge), source = .false.)
!   do ii = 1, dev%par%oxide(IDX_CELL,0)%n
!     idx1 = dev%par%oxide(IDX_CELL,0)%get_idx(ii)
!     oxide(dev%par%gtr%cell2edge(1,idx1(1))) = .true.
!     oxide(dev%par%gtr%cell2edge(2,idx1(1))) = .true.
!     oxide(dev%par%gtr%cell2edge(3,idx1(1))) = .true.
!   end do
!   do ii = 1, dev%par%transport(IDX_CELL,0)%n
!     idx1 = dev%par%transport(IDX_CELL,0)%get_idx(ii)
!     silicon(dev%par%gtr%cell2edge(1,idx1(1))) = .true.
!     silicon(dev%par%gtr%cell2edge(2,idx1(1))) = .true.
!     silicon(dev%par%gtr%cell2edge(3,idx1(1))) = .true.
!   end do

!   allocate (x(2,dev%par%gtr%nedge), y(2,dev%par%gtr%nedge))

!   x0 = -8
!   x1 = 300
!   y0 = -215
!   y1 = 215

!   Lx = 12
!   Ly = 8

!   do ii = 1, dev%par%gtr%nedge
!     idx1 = [dev%par%gtr%edge2vert(1,ii)]
!     idx2 = [dev%par%gtr%edge2vert(2,ii)]
!     call dev%par%gtr%get_vertex(idx1, p1)
!     call dev%par%gtr%get_vertex(idx2, p2)
!     x(1,ii) = (denorm(p1(2), "nm") - y0) / (y1 - y0) * Lx
!     y(1,ii) = Ly - (denorm(p1(1), "nm") - x0) / (x1 - x0) * Ly
!     x(2,ii) = (denorm(p2(2), "nm") - y0) / (y1 - y0) * Lx
!     y(2,ii) = Ly - (denorm(p2(1), "nm") - x0) / (x1 - x0) * Ly
!   end do

!   open (newunit = funit1, file = "gtr_oxide.csv", status = "replace", action = "write")
!   open (newunit = funit2, file = "gtr_silicon.csv", status = "replace", action = "write")
!   do ii = 1, dev%par%gtr%nedge
!     if (oxide(ii)) then
!       write (funit1, "(2ES25.16E3)") x(1,ii), y(1,ii)
!       write (funit1, "(2ES25.16E3)") x(2,ii), y(2,ii)
!       write (funit1, *)
!     end if
!     if (silicon(ii)) then
!       write (funit2, "(2ES25.16E3)") x(1,ii), y(1,ii)
!       write (funit2, "(2ES25.16E3)") x(2,ii), y(2,ii)
!       write (funit2, *)
!     end if
!   end do
!   close (funit1)
!   close (funit2)
!   stop
! end block

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
    integer :: ci, icurr, idens, iion, itab, ibl, i0, i1
    icurr = dev%sys_full%search_main_var("currents")
    ibl = dev%sys_full%res2block(icurr)%d(1)
    i0  = dev%sys_full%i0(ibl)
    i1  = dev%sys_full%i1(ibl)
    opt_full%atol(i0:i1) = norm(1e-12, "A")

    do ci = dev%par%ci0, dev%par%ci1
      idens = dev%sys_full%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_full%res2block(idens)%d)
        ibl = dev%sys_full%res2block(idens)%d(itab)
        i0 = dev%sys_full%i0(ibl)
        i1 = dev%sys_full%i1(ibl)
        opt_full%xmin(i0:i1) = norm(1e-10, "1/cm^3")
      end do

      idens = dev%sys_dd(ci)%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_dd(ci)%res2block(idens)%d)
        ibl = dev%sys_dd(ci)%res2block(idens)%d(itab)
        i0 = dev%sys_dd(ci)%i0(ibl)
        i1 = dev%sys_dd(ci)%i1(ibl)
        opt_dd(ci)%xmin(i0:i1) = norm(1e-10, "1/cm^3")
      end do

      if (dev%par%smc%incomp_ion) then
        ! iion = dev%sys_full%search_main_var(dev%ion(ci)%name)
        ! do itab = 1, size(dev%sys_full%res2block(iion)%d)
        !   ibl = dev%sys_full%res2block(iion)%d(itab)
        !   i0 = dev%sys_full%i0(ibl)
        !   i1 = dev%sys_full%i1(ibl)
        !   opt_full%xmin(i0:i1) = norm(1e-10, "1/cm^3")
        !   opt_full%xmax(i0:i1) = huge(1.0)
        ! end do

        ! iion = dev%sys_dd(ci)%search_main_var(dev%ion(ci)%name)
        ! do itab = 1, size(dev%sys_dd(ci)%res2block(iion)%d)
        !   ibl = dev%sys_dd(ci)%res2block(iion)%d(itab)
        !   i0 = dev%sys_dd(ci)%i0(ibl)
        !   i1 = dev%sys_dd(ci)%i1(ibl)
        !   opt_dd(ci)%xmin(i0:i1) = norm(1e-10, "1/cm^3")
        !   opt_dd(ci)%xmax(i0:i1) = huge(1.0)
        ! end do
      end if
    end do
  end block


  ! solve
  call solve_steady_state()
  call solve_small_signal()
  call solve_harmonic_balance()
  call solve_transient()
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
  use grid_m, only: IDX_VERTEX
  integer :: ii, jj, funit
  integer :: nx, ny, ict, idx(2)
  real, allocatable :: x(:), y(:), pot(:,:), dens(:,:)

  nx = dev%par%g1D(1)%n
  ny = dev%par%g1D(2)%n

  allocate (x(nx), y(ny), pot(nx,ny), dens(nx,ny))
  x = denorm(dev%par%g1D(1)%x, "nm")
  y = denorm(dev%par%g1D(2)%x, "nm")

  open (newunit = funit, file = "x.csv", status = "replace", action = "write")
  do ii = 1, nx
    write (funit, "(ES25.16E3)") x(ii)
  end do
  close (funit)
  open (newunit = funit, file = "y.csv", status = "replace", action = "write")
  do ii = 1, ny
    write (funit, "(ES25.16E3)") y(ii)
  end do
  close (funit)

  call ss%select(1)
  do jj = 1, ny
    do ii = 1, nx
      idx = [ii, jj]
      ict = dev%par%ict%get(idx)
      if (ict == 0) then
        pot(ii,jj) = denorm(dev%pot%get(idx), "V")
      else
        pot(ii,jj) = denorm(dev%volt(ict)%x + dev%par%contacts(ict)%phims, "V")
      end if

      if (dev%par%transport(IDX_VERTEX,0)%flags%get(idx)) then
        dens(ii,jj) = denorm(dev%dens(1)%get(idx), "1/cm^3")
      end if
    end do
  end do
  open (newunit = funit, file = "pot0.csv", status = "replace", action = "write")
  do ii = 1, nx
    do jj = 1, ny
      write (funit, "(ES25.16E3)", advance = "no") pot(ii,jj)
    end do
    write (funit, *)
  end do
  close (funit)
  open (newunit = funit, file = "dens0.csv", status = "replace", action = "write")
  do ii = 1, nx
    do jj = 1, ny
      write (funit, "(ES25.16E3)", advance = "no") dens(ii,jj)
    end do
    write (funit, *)
  end do
  close (funit)
end block

! block
!   integer :: ii, funit
!   real, allocatable :: VD(:), ID(:)

!   allocate (VD(nsweep), ID(nsweep))

!   do ii = 1, nsweep
!     call ss%select(ii)
!     VD(ii) = denorm(dev%volt(2)%x, "V")
!     ID(ii) = denorm(dev%curr(2)%x, "A")
!   end do

!   open (newunit = funit, file = "iv.csv", status = "replace", action = "write")
!   write (funit, "(A)") "V I"
!   do ii = 1, nsweep
!     write (funit, "(2ES25.16E3)") VD(ii), ID(ii)
!   end do
!   close (funit)
! end block

! block
!   integer :: ii, funit, nx
!   real, allocatable :: x(:), E0(:,:), iref(:,:), dens(:,:)

!   nx = dev%par%g1D(1)%n
!   allocate (x(nx), E0(nx,4), iref(nx,4), dens(nx,4))

!   x = denorm(dev%par%g1D(1)%x, "um")
!   call ss%select(1)
!   call dev%sys_full%eval()
!   do ii = 1, nx
!     E0(  ii,1) = denorm(dev%par%smc%band_edge(1) - dev%pot%x1(ii), "eV")
!     E0(  ii,2) = denorm(dev%par%smc%band_edge(2) - dev%pot%x1(ii), "eV")
!     iref(ii,1) = denorm(- dev%iref(1)%x1(ii), "eV")
!     iref(ii,2) = denorm(- dev%iref(2)%x1(ii), "eV")
!     dens(ii,1) = denorm(dev%dens(1)%x1(ii), "1/cm^3")
!     dens(ii,2) = denorm(dev%dens(2)%x1(ii), "1/cm^3")
!   end do
!   call ss%select(2)
!   call dev%sys_full%eval()
!   do ii = 1, nx
!     E0(  ii,3) = denorm(dev%par%smc%band_edge(1) - dev%pot%x1(ii), "eV")
!     E0(  ii,4) = denorm(dev%par%smc%band_edge(2) - dev%pot%x1(ii), "eV")
!     iref(ii,3) = denorm(- dev%iref(1)%x1(ii), "eV")
!     iref(ii,4) = denorm(- dev%iref(2)%x1(ii), "eV")
!     dens(ii,3) = denorm(dev%dens(1)%x1(ii), "1/cm^3")
!     dens(ii,4) = denorm(dev%dens(2)%x1(ii), "1/cm^3")
!   end do

!   open (newunit = funit, file = "pn.csv", status = "replace", action = "write")
!   write (funit, "(A)") "x Ec0 Ev0 Ec1 Ev1 imrefn0 imrefp0 imrefn1 imrefp1 n0 p0 n1 p1"
!   do ii = 1, dev%par%g1D(1)%n
!     write (funit, "(13ES25.16E3)") x(ii), E0(ii,1:4), iref(ii,1:4), dens(ii,1:4)
!   end do
!   close (funit)
! end block

! block
!   integer :: ii, funit
!   character(32) :: fname

!   write (fname, "(A,F4.2,A)") "output_", denorm(Vi(3,1), "V"), ".csv"

!   open (newunit = funit, file = fname, status = "replace", action = "write")
!   do ii = 1, nsweep
!     call ss%select(ii)
!     write (funit, "(2ES25.16E3)") denorm(dev%volt(2)%x, "V"), denorm(dev%curr(2)%x, "A")
!   end do
!   close (funit)
! end block
    end do
  end subroutine

  subroutine solve_small_signal()
    integer                   :: si, ict, Nf, Nss
    integer,      allocatable :: sids(:)
    logical                   :: flog
    real                      :: f0, f1
    real,         allocatable :: f(:), Vi(:,:), tt(:)
    complex,      allocatable :: s(:)
    type(string)              :: name
    type(polygon_src)         :: input
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

      ! ramp up ohmic contact voltages
      call file%get(sids(si), "Nss", Nss)
      allocate (Vi(dev%par%nct, 2))
      do ict = 1, dev%par%nct
        Vi(ict,:) = dev%volt(ict)%x
        if (dev%par%contacts(ict)%type == CT_OHMIC) Vi(ict,1) = 0.0
      end do
      call input%init([0.0, 1.0], Vi)
      allocate (tt(Nss))
      tt = linspace(0.0, 1.0, Nss)

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = tt, gum = gummel)

      ! run small-signal analysis
      call ac%run_analysis(dev%sys_full, s, ofile = ofile, oname = name%s)

block
  integer :: funit, ii

  open (newunit = funit, file = "ac.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "f rY11 rY12 rY21 rY22"
  do ii = 1, Nf
    write (funit, "(ES25.16E3)", advance = "no") denorm(f(ii), "Hz")

    call ac%select_real(3, ii)
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(3)%x, "A/V")
    call ac%select_real(2, ii)
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(3)%x, "A/V")
    call ac%select_real(3, ii)
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(2)%x, "A/V")
    call ac%select_real(2, ii)
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(2)%x, "A/V")

    write (funit, *)
  end do
  close (funit)
end block

      deallocate (f, s, tt)
    end do
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

  subroutine solve_transient()
    integer              :: ci, i0, i1, ibl, ict, icurr, idens, itab, Nt, si
    integer, allocatable :: sids(:)
    real                 :: dt0
    real,    allocatable :: ti(:), Vtmp(:), Vi(:,:)
    type(string)         :: name
    type(newton_opt)     :: opt_tr
    type(steady_state)   :: ss
    type(tr_bdf2)        :: trbdf2
    type(polygon_src)    :: input

    ! iteration parameters
    call load_iteration_params("transient params", dev%sys_full%n, opt_tr)
    do ci = dev%par%ci0, dev%par%ci1
      idens = dev%sys_full%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_full%res2block(idens)%d)
        ibl = dev%sys_full%res2block(idens)%d(itab)
        i0 = dev%sys_full%i0(ibl)
        i1 = dev%sys_full%i1(ibl)
        opt_tr%xmin(i0:i1) = norm(1e-10, "1/cm^3")
      end do
    end do
    icurr = dev%sys_full%search_main_var("currents")
    ibl = dev%sys_full%res2block(icurr)%d(1)
    i0  = dev%sys_full%i0(ibl)
    i1  = dev%sys_full%i1(ibl)
    opt_tr%atol(i0:i1) = norm(1e-12, "A")
    opt_tr%msg = "Transient: "

    call file%get_sections("transient", sids)
    do si = 1, size(sids)
      print "(A)", "transient"

      call file%get(sids(si), "name", name)
      call file%get(sids(si), "dt0", dt0)
      call file%get(sids(si), "t", ti)
      Nt = size(ti)
      allocate (Vi(dev%par%nct, Nt))
      do ict = 1, dev%par%nct
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vtmp)
        Vi(ict,:) = Vtmp
      end do

      call input%init(ti, Vi)

      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.

      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = [0.0], gum = gummel)

      call trbdf2%init(dev%sys_full, ti(1), ti(Nt), delta_t = dt0, nopt = opt_tr, adaptive = .true., log = .true.)
      call trbdf2%run(input = input)

block
  integer, parameter :: Ntt = 201
  integer            :: funit, ii
  real, allocatable  :: tt(:)

  allocate (tt(Ntt))
  tt = logspace(ti(3), ti(Nt), Ntt)

  open (newunit = funit, file = "tr.csv", status = "replace", action = "write")
  do ii = 1, Ntt
    write (funit, "(ES25.16E3)", advance = "no") denorm(tt(ii), "s")

    call trbdf2%select(tt(ii))
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(1)%x, "A")
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(2)%x, "A")
    write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(3)%x, "A")

    write (funit, *)
  end do
  close (funit)
end block
    end do
    print *
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
      call ss_nlpe%select(1)
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()
        call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd(ci))
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

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
