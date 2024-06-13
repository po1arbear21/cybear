program dd

  use approx_m,           only: approx_imref, approx_potential
  use bin_search_m,       only: bin_search
  use contact_m,          only: CT_OHMIC
  use cl_options_m,       only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,           only: dev
  use error_m,            only: program_error
  use harmonic_balance_m, only: harmonic_balance
  use input_m,            only: input_file
  use input_src_m,        only: const_src, polygon_src, harmonic_src
  use math_m,             only: linspace, logspace, PI
  use matrix_m,           only: default_spsolver, SPSOLVER_MUMPS
  use newton_m,           only: newton_opt
  use normalization_m,    only: init_normconst, norm, denorm
  use small_signal_m,     only: small_signal
  use steady_state_m,     only: steady_state
  use string_m,           only: string
  use transient_m,        only: backward_euler, tr_bdf2
  use util_m,             only: int2str
  use variable_m,         only: variable_real

  implicit none

  logical          :: gummel_restart, gummel_once, gummel_enabled
  real             :: temperature
  type(string)     :: projectdir
  type(input_file) :: runfile
  type(newton_opt) :: opt_nlpe, opt_dd(2), opt_gum, opt_full

! block
!   use distributions_m, only: fermi_dirac_integral_1h
!   use fermi_m,         only: fermi_dirac_integral_1h_reg
!   use math_m,          only: linspace

!   integer, parameter :: N = 1001
!   integer            :: funit, i
!   real               :: tmp, y, yreg
!   real, allocatable  :: x(:)

!   allocate (x(N))
!   x = linspace(-10000.0, 100.0, N)
!   open (newunit = funit, file = "F.csv", status = "replace", action = "write")
!   do i = 1, N
!     call fermi_dirac_integral_1h(x(i), y, tmp)
!     call fermi_dirac_integral_1h_reg(x(i), yreg, tmp)

!     write (funit, "(3ES25.16E3)") x(i), y, yreg
!   end do
!   close (funit)
!   stop
! end block

  ! parse command line arguments
  call command_line()

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

  ! use mumps solver
  ! default_spsolver = SPSOLVER_MUMPS

  opt_nlpe%msg  = "NLPE: "
  opt_dd(1)%msg = "NDD: "
  opt_dd(2)%msg = "PDD: "
  opt_gum%msg   = "Gummel: "
  opt_full%msg  = "Newton: "

  block
    integer :: ci, icurr, idens, itab, ibl, i0, i1
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
        opt_full%xmin(i0:i1) = norm(1e-80, "1/cm^3")
      end do

      idens = dev%sys_dd(ci)%search_main_var(dev%dens(ci)%name)
      do itab = 1, size(dev%sys_dd(ci)%res2block(idens)%d)
        ibl = dev%sys_dd(ci)%res2block(idens)%d(itab)
        i0 = dev%sys_dd(ci)%i0(ibl)
        i1 = dev%sys_dd(ci)%i1(ibl)
        opt_dd(ci)%xmin(i0:i1) = norm(1e-80, "1/cm^3")
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

contains

  subroutine command_line()
    !! parse command line arguments, init normalization, device and run file

    integer                      :: idesc, i, j, funit
    integer,         allocatable :: iclopt(:), jclopt(:)
    type(cl_option), allocatable :: clopt(:)
    type(cl_option_descriptor)   :: desc(4) = [ &
      cl_option_descriptor('P', "project_dir", .true., .false., .true., .true.), &
      cl_option_descriptor('T', "temperature", .true., .false., .true., .true.), &
      cl_option_descriptor('d', "device",      .true., .false., .true., .true.), &
      cl_option_descriptor('r', "run",         .true., .false., .true., .true.)  &
    ]

    call get_cl_options(desc, clopt, iclopt, jclopt)
    do idesc = 1, size(desc)
      do i = iclopt(idesc), iclopt(idesc + 1) - 1
        j = jclopt(i)
        select case (clopt(j)%short)
        case ('P')
          projectdir%s = clopt(j)%arg

        case ('T')
          read (clopt(j)%arg, *) temperature
          call init_normconst(temperature)
          print "(A,ES25.16E3,A)", "T = ", temperature, " K"

        case ('d')
          call dev%init(projectdir%s, clopt(j)%arg, temperature)

        case ('r')
          call runfile%init(clopt(j)%arg)

        end select
      end do
    end do

    ! output grids
    open (newunit = funit, file = "x.csv", status = "replace", action = "write")
    do i = 1, dev%par%g1D(1)%n
      write (funit, "(ES25.16E3)") denorm(dev%par%g1D(1)%x(i), "um")
    end do
    close (funit)
    open (newunit = funit, file = "y.csv", status = "replace", action = "write")
    do i = 1, dev%par%g1D(2)%n
      write (funit, "(ES25.16E3)") denorm(dev%par%g1D(2)%x(i), "um")
    end do
    close (funit)
    open (newunit = funit, file = "z.csv", status = "replace", action = "write")
    do i = 1, dev%par%g1D(3)%n
      write (funit, "(ES25.16E3)") denorm(dev%par%g1D(3)%x(i), "um")
    end do
    close (funit)

    ! block
    !   use grid_m, only: IDX_VERTEX, IDX_EDGE, IDX_CELL
    !   integer :: idx(3), idx_bnd(2,3)


    !   call dev%par%g%get_idx_bnd(IDX_CELL, 0, idx_bnd)
    !   idx(3) = idx_bnd(1,3) - 1 + bin_search(dev%par%g1D(3)%x, 0.0)-1
    !   idx(1) = idx_bnd(2,1)
    !   do j = idx_bnd(1,2), idx_bnd(2,2)
    !     idx(2) = j
    !     print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(j), "um"), denorm(dev%par%eps(IDX_CELL,0)%get(idx), "eps0")
    !   end do
    !   print *
    !   idx(3) = idx_bnd(1,3) - 1 + bin_search(dev%par%g1D(3)%x, 0.0)
    !   idx(1) = idx_bnd(2,1)
    !   do j = idx_bnd(1,2), idx_bnd(2,2)
    !     idx(2) = j
    !     print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(j), "um"), denorm(dev%par%eps(IDX_CELL,0)%get(idx), "eps0")
    !   end do
    !   print *
    !   call dev%par%g%get_idx_bnd(IDX_EDGE, 1, idx_bnd)
    !   idx(3) = idx_bnd(1,3) - 1 + bin_search(dev%par%g1D(3)%x, 0.0)
    !   idx(1) = idx_bnd(2,1)
    !   do j = idx_bnd(1,2), idx_bnd(2,2)
    !     idx(2) = j
    !     print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(j), "um"), denorm(dev%par%eps(IDX_EDGE,1)%get(idx), "eps0")
    !   end do
    !   print *
    !   call dev%par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    !   idx(3) = idx_bnd(1,3) - 1 + bin_search(dev%par%g1D(3)%x, 0.0)
    !   idx(1) = idx_bnd(2,1)
    !   do j = idx_bnd(1,2), idx_bnd(2,2)
    !     idx(2) = j
    !     print "(2ES25.16E3)", denorm(dev%par%g1D(2)%x(j), "um"), denorm(dev%par%vol%get(idx), "um^3")
    !   end do
    !   stop
    ! end block
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

  subroutine solve_steady_state()
    integer              :: si, ict, ict_sweep, nsweep, isweep
    integer, allocatable :: sids(:)
    logical              :: status
    real,    allocatable :: Vi(:,:), Vbounds(:), t(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    allocate (Vi(dev%par%nct,2))

    call runfile%get_sections("steady state", sids)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call runfile%get(sids(si), "name", name)

      ! voltage input source
      ict_sweep = 0
      do ict = 1, dev%par%nct
        call runfile%get(sids(si), "N_"//dev%par%contacts(ict)%name, nsweep, status = status)
        if (status) status = (nsweep > 1)
        if (status) then
          if (ict_sweep /= 0) call program_error("only one voltage sweep per steady-state section allowed")
          ict_sweep = ict
          call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vbounds)
          Vi(ict,:) = Vbounds
        else
          call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, Vi(ict,1))
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
opt_full%error_if_not_converged = .false.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = t, gum = gummel)

! block
!   print *, "poisson"
!   call dev%poiss%test(500)
!   print *, "rho"
!   call dev%calc_rho%test(500, no_sign_change = [.true., .true., .true., .true.])
!   print *, "contin 1"
!   call dev%contin(1)%test(500, no_sign_change = [.true., .false., .false., .false., .true.])
!   print *, "contin 2"
!   call dev%contin(2)%test(500, no_sign_change = [.true., .false., .false., .false., .true.])
!   print *, "ion_contin 1"
!   call dev%ion_contin(1)%test(500, no_sign_change = [.true., .false.])
!   print *, "ion_contin 2"
!   call dev%ion_contin(2)%test(500, no_sign_change = [.true., .false.])
!   print *, "genrec 1"
!   call dev%calc_genrec(1)%test(500, no_sign_change = [.false., .false., .true.])
!   print *, "genrec 2"
!   call dev%calc_genrec(2)%test(500, no_sign_change = [.false., .false., .true.])
!   print *, "iref 1"
!   call dev%calc_iref(1)%test(500, no_sign_change = [.false., .true.])
!   print *, "iref 2"
!   call dev%calc_iref(2)%test(500, no_sign_change = [.false., .true.])
!   print *, "calc_cdens 1"
!   call dev%calc_cdens(1,1)%test(500, no_sign_change = [.false., .true., .true.])
!   call dev%calc_cdens(2,1)%test(500, no_sign_change = [.false., .true., .true.])
!   call dev%calc_cdens(3,1)%test(500, no_sign_change = [.false., .true., .true.])
!   print *, "calc_cdens 2"
!   call dev%calc_cdens(1,2)%test(500, no_sign_change = [.false., .true., .true.])
!   call dev%calc_cdens(2,2)%test(500, no_sign_change = [.false., .true., .true.])
!   call dev%calc_cdens(3,2)%test(500, no_sign_change = [.false., .true., .true.])
! end block

      if (nsweep < 1) then
        nsweep = 1
      end if
      do isweep = 1, nsweep
        call ss%select(isweep)

        call output_xy(dev%pot,     0.0, "potxy_" // int2str(isweep) // ".csv")
        call output_xz(dev%pot,     0.0, "potxz_" // int2str(isweep) // ".csv")
        call output_xy(dev%dens(1), 0.0,   "nxy_" // int2str(isweep) // ".csv")
        call output_xz(dev%dens(1), 0.0,   "nxz_" // int2str(isweep) // ".csv")
        call output_xy(dev%dens(2), 0.0,   "pxy_" // int2str(isweep) // ".csv")
        call output_xz(dev%dens(2), 0.0,   "pxz_" // int2str(isweep) // ".csv")
        call output_xy(dev%ion(1),  0.0, "ionxy_" // int2str(isweep) // ".csv")
        call output_xz(dev%ion(1),  0.0, "ionxz_" // int2str(isweep) // ".csv")
      end do


      ! block
      !   integer :: funit, i, j, k, k0
      !   real    :: dz, p

      !   k0 = bin_search(dev%par%g1D(3)%x, 0.0)

      !   open (newunit = funit, file = "pinv"//int2str(si)//".csv", status = "replace", action = "write")

      !   do ict_sweep = 1, nsweep
      !     call ss%select(ict_sweep)

      !     i = dev%par%g1D(1)%n
      !     j = dev%par%g1D(2)%n
      !     p = 0
      !     do k = 1, k0 - 1
      !       dz = dev%par%g1D(3)%x(k+1) - dev%par%g1D(3)%x(k)
      !       p = p + 0.5 * dz * (dev%dens(2)%get([i,j,k]) + dev%dens(2)%get([i,j,k+1]))
      !     end do

      !     write (funit, "(2ES25.16E3)") denorm(dev%volt(2)%x, "V"), denorm(p, "1/cm^2")
      !   end do

      !   close (funit)
      ! end block
    end do
  end subroutine

  subroutine solve_small_signal()
    integer                   :: si, ict, Nf, Nss, i
    integer,      allocatable :: sids(:)
    logical                   :: flog, status
    real                      :: f0, f1
    real,         allocatable :: f(:), Vi(:,:), tt(:)
    complex,      allocatable :: s(:)
    type(string)              :: name
    type(polygon_src)         :: input
    type(small_signal)        :: ac
    type(steady_state)        :: ss

    call runfile%get_sections("small signal", sids)
    do si = 1, size(sids)
      print "(A)", "small signal"

      call runfile%get(sids(si), "name", name)

      ! get DC voltages
      do ict = 1, dev%par%nct
        call runfile%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do

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
      s = 2 * PI * (0.0, 1.0) * f

      ! ramp up ohmic contact voltages
      call runfile%get(sids(si), "Nss", Nss, status = status)
      if (.not. status) Nss = 1
      if (status) status = (Nss > 1)
      allocate (Vi(dev%par%nct, 2))
      do ict = 1, dev%par%nct
        Vi(ict,:) = dev%volt(ict)%x
        if ((Nss > 1) .and. (dev%par%contacts(ict)%type == CT_OHMIC)) Vi(ict,1) = 0.0
      end do
      call input%init([0.0, 1.0], Vi)
      if (Nss == 0) then
        tt = [ 0.0 ]
      else
        print *, Nss
        tt = linspace(0.0, 1.0, Nss)
      end if

      block
        use grid_generator_m, only: DIR_NAME
        integer :: funit, i, dir

        do dir = 1, 3
          open (newunit = funit, file = DIR_NAME(dir)//".csv", status = "replace", action = "write")
          do i = 1, dev%par%g1D(dir)%n
            write (funit, "(ES25.16E3)") denorm(dev%par%g1D(dir)%x(i), "um")
          end do
          close (funit)
        end do
      end block

      ! solve steady-state
      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.
      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = tt, gum = gummel)

      ! run small-signal analysis
      call ac%run_analysis(dev%sys_full, s)

      call ac%select_abs(2, 1)
      call output_xy(dev%pot,     0.0, "ac1_potxy.csv")
      call output_xz(dev%pot,     0.0, "ac1_potxz.csv")
      call output_xy(dev%dens(1), 0.0, "ac1_nxy.csv")
      call output_xz(dev%dens(1), 0.0, "ac1_nxz.csv")
      call output_xy(dev%dens(2), 0.0, "ac1_pxy.csv")
      call output_xz(dev%dens(2), 0.0, "ac1_pxz.csv")

      call ac%select_abs(2, Nf)
      call output_xy(dev%pot,     0.0, "ac"//int2str(Nf)//"_potxy.csv")
      call output_xz(dev%pot,     0.0, "ac"//int2str(Nf)//"_potxz.csv")
      call output_xy(dev%dens(1), 0.0, "ac"//int2str(Nf)//"_nxy.csv")
      call output_xz(dev%dens(1), 0.0, "ac"//int2str(Nf)//"_nxz.csv")
      call output_xy(dev%dens(2), 0.0, "ac"//int2str(Nf)//"_pxy.csv")
      call output_xz(dev%dens(2), 0.0, "ac"//int2str(Nf)//"_pxz.csv")

      do i = 1, Nf
        print "(I6,ES25.16E3)", i, denorm(f(i), "Hz")
        call ac%select_abs(2, i)
        call output_xy(dev%dens(2), 0.0, "ac"//int2str(i)//"_pxy.csv")
      end do

      block
        integer :: funit, i

        open (newunit = funit, file = "ac_Y.csv", status = "replace", action = "write")
        do i = 1, Nf
          write (funit, "(ES25.16E3)", advance = "no") denorm(f(i), "Hz")

          call ac%select_real(2, i)
          write (funit, "(ES25.16E3)", advance = "no") denorm(dev%curr(2)%x, "A/V")
          call ac%select_imag(2, i)
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
        Vi(ict,:) = Vtmp
      end do

      call input%init(ti, Vi)

      gummel_restart = .true.
      gummel_once    = .true.
      gummel_enabled = .true.

      call ss%run(dev%sys_full, nopt = opt_full, input = input, t_input = [0.0], gum = gummel)

      call trbdf2%init(dev%sys_full, ti(1), ti(Nt), delta_t = dt0, nopt = opt_tr, adaptive = .true., log = .true.)
      call trbdf2%run(input = input)
    end do
  end subroutine

  subroutine solve_responsivity()
    integer                   :: i, si, ict, isrc, idrn, Nf, NH, Nt, ofunit
    integer,      allocatable :: sids(:)
    real                      :: VA, f0, f1, power, curr, resp
    real,         allocatable :: f(:), c(:,:), s(:,:)
    type(string)              :: source, drain, output

    type(harmonic_src)        :: input
    type(harmonic_balance)    :: hb
    type(steady_state)        :: ss
    type(newton_opt)          :: opt_hb

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
      call ss%run(dev%sys_full, nopt = opt_full, input = input, gum = gummel)

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

    integer            :: it, ci
    real               :: error, err_pot, err_iref(2)
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

!     ! no iteration necessary for equilibrium
!     if (dev%equilibrium()) then
!       call solve_nlpe()
!       do ci = dev%par%ci0, dev%par%ci1
!         ! call dev%sys_dd(ci)%eval()
!         call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd(ci))
!       end do
!       return
!     end if

    ! gummel iteration
    it = 0
    error = huge(1.0)
    do while (((error > opt_gum%atol(1)) .and. (it < opt_gum%max_it)) .or. (it < opt_gum%min_it))
      it = it + 1

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()

        call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd(ci))
        call ss_dd(ci)%select(1)
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

      ! log
      if (opt_gum%log) then
        print "(A,I6,ES25.16E3)", "Gummel: ", it, denorm(error, "V")
      end if
    end do

    if ((it > opt_gum%max_it) .and. (error > opt_gum%atol(1))) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if

! block
!   use steady_state_m, only: MYDEBUG, DEBUG_PROC

!   MYDEBUG = .true.
!   DEBUG_PROC => mydebug_proc
! end block
  end subroutine

  subroutine mydebug_proc()
    integer, save :: i = 0

    integer                       :: ix, iy, iz
    real                          :: dens
    class(variable_real), pointer :: vp

    i = i + 1
    print *, "mydebug_proc: ", i

    if (i == 2) then
      do iz = 1, dev%par%g1D(3)%n
        do iy = 1, dev%par%g1D(2)%n
          do ix = 1, dev%par%g1D(1)%n
            dens = denorm(dev%dens(1)%x3(ix,iy,iz), "1/cm^3")
            if ((dens > 0) .and. (dens < 1e-75)) then
              print "(3I6,ES25.16E3)", ix, iy, iz, dens
            end if
          end do
        end do
      end do
    end if
    ! 11     1    43


    call output_xz(dev%pot,       0.0,     "potxz_" // int2str(i) // ".csv")
    call output_xz(dev%dens(1),   0.0,       "nxz_" // int2str(i) // ".csv")
    call output_xz(dev%dens(2),   0.0,       "pxz_" // int2str(i) // ".csv")
    call output_xz(dev%iref(1),   0.0,  "nimrefxz_" // int2str(i) // ".csv")
    call output_xz(dev%iref(2),   0.0,  "pimrefxz_" // int2str(i) // ".csv")
    call output_xz(dev%ion(1),    0.0,    "nionxz_" // int2str(i) // ".csv")
    call output_xz(dev%ion(2),    0.0,    "pionxz_" // int2str(i) // ".csv")
    call output_xz(dev%genrec(1), 0.0, "ngenrecxz_" // int2str(i) // ".csv")
    call output_xz(dev%genrec(2), 0.0, "pgenrecxz_" // int2str(i) // ".csv")

    select type (p => dev%poiss%f%v(1)%p)
    class is (variable_real)
      vp => p
    end select
    call output_xz(vp, 0.0, "fpoissxz_" // int2str(i) // ".csv")

    select type (p => dev%contin(1)%f%v(1)%p)
    class is (variable_real)
      vp => p
    end select
    call output_xz(vp, 0.0, "fnxz_" // int2str(i) // ".csv")
    select type (p => dev%contin(2)%f%v(1)%p)
    class is (variable_real)
      vp => p
    end select
    call output_xz(vp, 0.0, "fpxz_" // int2str(i) // ".csv")

    select type (p => dev%ion_contin(1)%f%v(1)%p)
    class is (variable_real)
      vp => p
    end select
    call output_xz(vp, 0.0, "fnionxz_" // int2str(i) // ".csv")
    select type (p => dev%ion_contin(2)%f%v(1)%p)
    class is (variable_real)
      vp => p
    end select
    call output_xz(vp, 0.0, "fpionxz_" // int2str(i) // ".csv")
  end subroutine

  subroutine solve_nlpe()
    !! solve non-linear poisson equation

    use matrix_m, only: sparse_real

    integer           :: it, nx
    real              :: err, res0, res1, damping
    real, allocatable :: x0(:), f(:), dx(:)
    type(sparse_real) :: dfdx

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

    ! newton iteration
    do while (((err > opt_nlpe%atol(1)) .and. (it <= opt_nlpe%max_it)) .or. (it < opt_nlpe%min_it))
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
      call dfdx%reset()

      ! absolute error
      err = maxval(abs(dx))

      ! update variables
      call dev%sys_nlpe%set_x(x0 + dx)

      write (*, "(ES25.16E3)") denorm(err, "V")
    end do
  end subroutine

  subroutine output_xy(v, z, file)
    class(variable_real), intent(in) :: v
    real,                 intent(in) :: z
    character(*),         intent(in) :: file

    integer           :: i, idx(3), idx_bnd(2,3), j, k
    real, allocatable :: data(:,:)

    call dev%par%g%get_idx_bnd(v%idx_type, v%idx_dir, idx_bnd)
    allocate (data(idx_bnd(1,1):idx_bnd(2,1),idx_bnd(1,2):idx_bnd(2,2)), source = 0.0)
    k = idx_bnd(1,3) - 1 + bin_search(dev%par%g1D(3)%x, z)
    do j = idx_bnd(1,2), idx_bnd(2,2); do i = idx_bnd(1,1), idx_bnd(2,1)
      idx = [i, j, k]
      data(i,j) = denorm(v%get(idx), v%unit)
    end do; end do

    call write_to_file(data, file)
  end subroutine

  subroutine output_xz(v, y, file)
    class(variable_real), intent(in) :: v
    real,                 intent(in) :: y
    character(*),         intent(in) :: file

    integer           :: i, idx(3), idx_bnd(2,3), j, k
    real, allocatable :: data(:,:)

    call dev%par%g%get_idx_bnd(v%idx_type, v%idx_dir, idx_bnd)
    allocate (data(idx_bnd(1,1):idx_bnd(2,1),idx_bnd(1,3):idx_bnd(2,3)), source = 0.0)
    j = idx_bnd(1,2) - 1 + bin_search(dev%par%g1D(2)%x, y)
    do k = idx_bnd(1,3), idx_bnd(2,3); do i = idx_bnd(1,1), idx_bnd(2,1)
      idx = [i, j, k]
      data(i,k) = denorm(v%get(idx), v%unit)
    end do; end do

    call write_to_file(data, file)
  end subroutine

  subroutine write_to_file(data, file)
    real,         intent(in) :: data(:,:)
    character(*), intent(in) :: file

    integer :: i, j, funit

    open (newunit = funit, file = file, status = "replace", action = "write")
    do i = 1, size(data, 1)
      do j = 1, size(data, 2)
        write (funit, "(ES25.16E3)", advance = "no") data(i,j)
      end do
      write (funit, *)
    end do
    close (funit)
  end subroutine

end program
