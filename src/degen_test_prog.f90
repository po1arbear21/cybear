program degen_test_prog

  use cl_options_m, only: get_cl_options_simple
  use degen_m,      only: degen_init, degen_get, get_current, fd12, ifd12, weierstrass, DEGEN_DEBUG, DEGEN_TANH_SINH
  use gauss_table_m
  use fukushima_m, only: fd1h
  use math_m,      only: linspace, PI
  use omp_lib,     only: omp_get_num_threads, omp_get_thread_num
  use qsort_m,     only: qsort
  use string_m,    only: new_string, string
  use util_m,      only: int2str
  use vector_m,    only: vector_int, vector_real

  implicit none

  interface
    subroutine curve(t, deta, ddetadt, dpot, ddpotdt)
      real, intent(in)  :: t
      real, intent(out) :: deta
      real, intent(out) :: ddetadt
      real, intent(out) :: dpot
      real, intent(out) :: ddpotdt
    end subroutine
  end interface

  integer      :: ng_, ngexp_, i
  logical      :: output, use_weierstrass
  real         :: err, eta1, rad
  type(string) :: values(6)

  call get_cl_options_simple([new_string("NG"), new_string("_NGEXP"), new_string("output"), new_string("eta1"), new_string("rad"), new_string("weierstrass")], values)

  ! unescape minus
  do i = 1, 6
    if (values(i)%s(1:1) == "!") values(i)%s = values(i)%s(2:len(values(i)%s))
  end do

  read (values(1)%s,*) ng_
  read (values(2)%s,*) ngexp_
  read (values(3)%s,*) output
  read (values(4)%s,*) eta1
  read (values(5)%s,*) rad
  read (values(6)%s,*) use_weierstrass

  call degen_init(ng_, ngexp_)

  ! DEGEN_DEBUG = .false.
  ! call check_error(circle, 0.0, 2 * PI, output, err)
  ! print "(A,ES25.16E3)", "err = ", err
  ! stop

  ! DEGEN_DEBUG     = .false.
  ! DEGEN_TANH_SINH = .false.
  ! call check_continuity(circle, 0.0, 2 * PI, output, err)
  ! print "(A,ES25.16E3)", "err = ", err
  ! stop

  block
    real :: t, deta, dpot, eta2, n1, n2, j, j0, djdn1, djdn2, djddpot, nquad, dum1, dum2
    logical :: status

    DEGEN_DEBUG = .true.

    eta1 = -20
    deta = 98.9529
    eta2 = eta1 + deta
    dpot = 80.6667

    call fd12(eta1, F = n1)
    call fd12(eta2, F = n2)

    DEGEN_TANH_SINH = .true.
    call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
    print *, j
    j0 = j

    DEGEN_TANH_SINH = .false.
    call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
    print *, j

    print *, j0 - j
    print *, abs(j0 - j) / abs(j0)

    ! call weierstrass(fd12, ifd12, n1, n2, dpot, j, djdn1, djdn2, djddpot, status)
    ! print *, j

    stop
  end block

  ! block
  !   integer, parameter :: NN = 5001
  !   integer :: ii, funit
  !   real :: eta1, eta2, dpot0, n10, n20, j0, j, djdn1, djdn2, djddpot, nquad, tmp1, tmp2, tmp3
  !   real, allocatable :: n1(:), n2(:), dpot(:)

  !   dpot0 =   9.7245093720238190E-003
  !   eta1  =  -2.4999998108678277E+001
  !   eta2  =  -4.4408920985006271E-016

  !   ! DEGEN_DEBUG = .false.
  !   call fd12(eta1, F = n10)
  !   call fd12(eta2, F = n20)

  !   call get_current(fd12, ifd12, 0.0, sqrt(0.125), n10, n20, dpot0, j0, djdn1, djdn2, djddpot, nquad)

  !   n1 = linspace(n10 * (1-5e-1), n10 * (1+5e-1), NN)
  !   n2 = linspace(n20 * (1-5e-1), n20 * (1+5e-1), NN)
  !   dpot = linspace(dpot0 - 1e-3, dpot0 + 1e-3, NN)

  !   open (newunit = funit, file = "dj_n1.csv", status = "replace", action = "write")
  !   do ii = 1, NN
  !     call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1(ii), n20, dpot0, j, tmp1, tmp2, tmp3, nquad)
  !     write (funit, "(2ES25.16E3)") n1(ii), j - (j0 + djdn1 * (n1(ii) - n10))
  !   end do
  !   close (funit)

  !   open (newunit = funit, file = "dj_n2.csv", status = "replace", action = "write")
  !   do ii = 1, NN
  !     call get_current(fd12, ifd12, 0.0, sqrt(0.125), n10, n2(ii), dpot0, j, tmp1, tmp2, tmp3, nquad)
  !     write (funit, "(2ES25.16E3)") n2(ii), j - (j0 + djdn2 * (n2(ii) - n20))
  !   end do
  !   close (funit)

  !   open (newunit = funit, file = "dj_dpot.csv", status = "replace", action = "write")
  !   do ii = 1, NN
  !     call get_current(fd12, ifd12, 0.0, sqrt(0.125), n10, n20, dpot(ii), j, tmp1, tmp2, tmp3, nquad)
  !     write (funit, "(2ES25.16E3)") dpot(ii), j - (j0 + djddpot * (dpot(ii) - dpot0))
  !   end do
  !   close (funit)

  !   open (newunit = funit, file = "djdpot.csv", status = "replace", action = "write")
  !   do ii = 1, NN
  !     call get_current(fd12, ifd12, 0.0, sqrt(0.125), n10, n20, dpot(ii), j, tmp1, tmp2, djddpot, nquad)
  !     write (funit, "(2ES25.16E3)") dpot(ii), djddpot
  !   end do
  !   close (funit)

  !   stop
  ! end block

  ! block
  !   use math_m, only: logspace
  !   use reg_m,  only: dist => fermi_dirac_integral_1h_reg, idist => inv_fermi_dirac_integral_1h_reg

  !   integer, parameter :: NN = 1001
  !   integer :: ii, funit
  !   real :: deta, dpot, n1, j, djdn(2), djddpot, nquad
  !   real, allocatable :: n2(:)

  !   n1 = 2.0313141997174735E-044
  !   n2 = [2.0313564794884010E-044]
  !   dpot =   1.0406890996819129E-003

  !   ! n2 = logspace(2.0010501511193074E-014, 2.0010501511193071E-011, NN)
  !   ! n2 = linspace(0.9e-13, 1.0e-13, NN)

  !   DEGEN_TANH_SINH = .true.
  !   DEGEN_DEBUG = .true.
  !   open (newunit = funit, file = "j.csv", status = "replace", action = "write")
  !   do ii = 1, size(n2)
  !     call get_current(dist, idist, 0.0, sqrt(0.125), n1, n2(ii), dpot, j, djdn(1), djdn(2), djddpot, nquad)

  !     write (funit, "(2ES25.16E3)") n2(ii), j
  !   end do
  !   close (funit)

  !   stop
  ! end block

  ! block
  !   integer, parameter :: NN = 5001
  !   integer            :: ii, funit
  !   real               :: deta, dpot, dum, n1, n2, j, eta2
  !   real, allocatable  :: t(:)

  !   t = linspace(2.905, 2.906, NN)

  !   open (newunit = funit, file = "j.csv", status = "replace", action = "write")
  !   do ii = 1, NN
  !     call circle(t(ii), deta, dum, dpot, dum)
  !     eta2 = eta1 + deta
  !     call fd12(eta1, F = n1)
  !     call fd12(eta2, F = n2)
  !     call weierstrass(fd12, ifd12, n1, n2, dpot, j, dum, dum, dum)
  !     write (funit, "(2ES25.16E3)") t(ii), j
  !   end do
  !   close (funit)
  !   stop
  ! end block

  ! integer           :: i, funit
  ! real              :: n1, n2, dpot, j, djdn1, djdn2, djddpot, eta1, eta2, nquad, tt
  ! real, allocatable :: t(:), deta_(:), dpot_(:), j0(:), j_(:)

  ! call degen_test_init()

  ! t = linspace(0.0, 2*PI, N)
  ! dpot_ = RAD * cos(t)
  ! deta_ = RAD * sin(t)
  ! allocate (j0(N), j_(N))

  ! eta1 = 5.0
  ! call fd12(eta1, F = n1)

  ! USE_TANH_SINH = .true.
  ! open (newunit = funit, file = "j0.csv", status = "replace", action = "write")
  ! write (funit, "(A)") "%lc=4"
  ! do i = 1, N
  !   eta2 = eta1 + deta_(i)
  !   call fd12(eta2, F = n2)
  !   dpot = dpot_(i)

  !   print "(A,ES25.16E3)", "t = ", t(i) / PI
  !   call get_current(fd12, ifd12, -1.0, sqrt(0.125), n1, n2, dpot, j0(i), djdn1, djdn2, djddpot, nquad)

  !   if (.not. ieee_is_finite(j0(i))) then
  !     print "(A,ES25.16E3)", "eta1 = ", eta1
  !     print "(A,ES25.16E3)", "eta2 = ", eta2
  !     print "(A,ES25.16E3)", "dpot = ", dpot
  !     error stop "NaN"
  !   end if
  !   write (funit, "(2ES25.16E3)") t(i) / PI, j0(i)
  ! end do
  ! close (funit)

  ! USE_TANH_SINH = .false.
  ! open (newunit = funit, file = "j_"//int2str(NG)//".csv", status = "replace", action = "write")
  ! do i = 1, N
  !   eta2 = eta1 + deta_(i)
  !   call fd12(eta2, F = n2)
  !   dpot = dpot_(i)

  !   print "(A,ES25.16E3)", "t = ", t(i) / PI
  !   call get_current(fd12, ifd12, -1.0, sqrt(0.125), n1, n2, dpot, j_(i), djdn1, djdn2, djddpot, nquad)

  !   if (.not. ieee_is_finite(j_(i))) then
  !     print "(A,ES25.16E3)", "eta1 = ", eta1
  !     print "(A,ES25.16E3)", "eta2 = ", eta2
  !     print "(A,ES25.16E3)", "dpot = ", dpot
  !     error stop "NaN"
  !   end if
  !   write (funit, "(2ES25.16E3)") t(i) / PI, j_(i)
  ! end do
  ! close (funit)

  ! open (newunit = funit, file = "ej"//int2str(NG)//".csv", status = "replace", action = "write")
  ! do i = 1, N
  !   write (funit, "(2ES25.16E3)") t(i) / PI, abs(j_(i) - j0(i)) / abs(j_(i))
  ! end do
  ! close (funit)



  ! tt = 1.503
  ! dpot = RAD * cos(tt*pi)
  ! eta1 = 5.0
  ! eta2 = eta1 + RAD * sin(tt*pi)
  ! call fd12(eta1, F = n1)
  ! call fd12(eta2, F = n2)

  ! USE_TANH_SINH = .false.
  ! call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
  ! print "(A,ES25.16E3)", "j = ", j

  ! USE_TANH_SINH = .true.
  ! call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
  ! print "(A,ES25.16E3)", "j = ", j

  ! tt = 1.5004
  ! dpot = RAD * cos(tt*pi)
  ! eta1 = 5.0
  ! eta2 = eta1 + RAD * sin(tt*pi)
  ! call fd12(eta1, F = n1)
  ! call fd12(eta2, F = n2)

  ! USE_TANH_SINH = .false.
  ! call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
  ! print "(A,ES25.16E3)", "j = ", j

  ! USE_TANH_SINH = .true.
  ! call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j, djdn1, djdn2, djddpot, nquad)
  ! print "(A,ES25.16E3)", "j = ", j

contains

  subroutine circle(t, deta, ddetadt, dpot, ddpotdt)
    real, intent(in)  :: t
    real, intent(out) :: deta
    real, intent(out) :: ddetadt
    real, intent(out) :: dpot
    real, intent(out) :: ddpotdt

    deta    =   rad * sin(t)
    ddetadt =   rad * cos(t)
    dpot    =   rad * cos(t)
    ddpotdt = - rad * sin(t)
  end subroutine

  subroutine dpot_const(t, deta, ddetadt, dpot, ddpotdt)
    real, intent(in)  :: t
    real, intent(out) :: deta
    real, intent(out) :: ddetadt
    real, intent(out) :: dpot
    real, intent(out) :: ddpotdt

    real, parameter :: m     = 1.0
    real, parameter :: deta0 = -24.0
    real, parameter :: dpot0 = 3.4517979415322336E-002

    deta    = deta0 + m * t
    ddetadt = m
    dpot    = dpot0
    ddpotdt = 0.0
  end subroutine

  subroutine deta_const(t, deta, ddetadt, dpot, ddpotdt)
    real, intent(in)  :: t
    real, intent(out) :: deta
    real, intent(out) :: ddetadt
    real, intent(out) :: dpot
    real, intent(out) :: ddpotdt

    real, parameter :: m     = 1.0
    real, parameter :: deta0 = -24.0
    real, parameter :: dpot0 = 0.0

    deta    = deta0
    ddetadt = 0.0
    dpot    = dpot0 + m * t
    ddpotdt = m
  end subroutine

  subroutine check_error(c, t0, t1, output, max_err)
    procedure(curve)     :: c
    real,    intent(in)  :: t0
    real,    intent(in)  :: t1
    logical, intent(in)  :: output
    real,    intent(out) :: max_err

    integer, parameter :: NT = 64 * 256 + 1

    integer           :: i, funit
    logical           :: status
    real              :: eta2, n1, n2, deta, ddetadt, dpot, ddpotdt, djdn1, djdn2, djddpot, nquad
    real, allocatable :: t(:), j0(:), j(:), err(:)

    ! init grid
    allocate (t(NT), j(NT), j0(NT), err(NT))
    t = linspace(t0, t1, NT)

    call fd12(eta1, F = n1)

    DEGEN_TANH_SINH = .true.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i,deta,ddetadt,dpot,ddpotdt,eta2,n2,djdn1,djdn2,djddpot,nquad) &
    !$omp shared(eta1,n1,t,j0)
    do i = 1, NT
      call c(t(i), deta, ddetadt, dpot, ddpotdt)
      eta2 = eta1 + deta
      call fd12(eta2, F = n2)
      call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j0(i), djdn1, djdn2, djddpot, nquad)
    end do
    !$omp end parallel do

    DEGEN_TANH_SINH = .false.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i,deta,ddetadt,dpot,ddpotdt,eta2,n2,djdn1,djdn2,djddpot,status,nquad) &
    !$omp shared(eta1,n1,t,j0,j,err,use_weierstrass)
    do i = 1, NT
      call c(t(i), deta, ddetadt, dpot, ddpotdt)
      eta2 = eta1 + deta
      call fd12(eta2, F = n2)
      if (use_weierstrass) then
        call weierstrass(fd12, ifd12, n1, n2, dpot, j(i), djdn1, djdn2, djddpot, status)
        if (.not. status) then
          err(i) = 1e-16
          cycle
        end if
      else
        call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j(i), djdn1, djdn2, djddpot, nquad)
      end if

      err(i) = abs(j(i) - j0(i)) / (abs(j0(i)) + 1e-16)
    end do
    !$omp end parallel do

    max_err = maxval(err)

    if (output) then
      open (newunit = funit, file = "j0.csv", status = "replace", action = "write")
      do i = 1, NT
        write (funit, "(2ES25.16E3)") t(i), j0(i)
      end do
      close (funit)
      open (newunit = funit, file = "j.csv", status = "replace", action = "write")
      do i = 1, NT
        write (funit, "(2ES25.16E3)") t(i), j(i)
      end do
      close (funit)
      open (newunit = funit, file = "err.csv", status = "replace", action = "write")
      do i = 1, NT
        write (funit, "(2ES25.16E3)") t(i), max(err(i), 1e-16)
      end do
      close (funit)
    end if
  end subroutine

  subroutine check_continuity(c, t0, t1, output, max_err)
    procedure(curve)     :: c
    real,    intent(in)  :: t0
    real,    intent(in)  :: t1
    logical, intent(in)  :: output
    real,    intent(out) :: max_err

    integer, parameter :: NT0 = 64 * 64 + 1
    real,    parameter :: RTOL = 1e-12, ATOL = 5e-12, TTOL = 1e-12

    integer              :: i, i1, i2, i3, ithread, nthreads, ntot, funit
    integer, allocatable :: nth(:), perm(:)
    logical              :: status
    real                 :: eta2, n1, n2, dn2deta2, deta, ddetadt, dpot, ddpotdt, djdn1, djdn2, djddpot, nquad, dt, dj1, j0
    real,    allocatable :: t(:), j(:), dj(:), djdt(:), tmp(:), eps(:), err(:)
    type(vector_real)    :: vt, vj, vdjdt
    type(vector_int)     :: vi

    call fd12(eta1, F = n1)

    ! init coarse grid
    allocate (t(NT0), j(NT0), djdt(NT0))
    t = linspace(t0, t1, NT0)

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,ithread,eta2,n2,dn2deta2,deta,ddetadt,dpot,ddpotdt,djdn1,djdn2,djddpot,nquad,dt,j0,dj1,vt,vj,vdjdt,vi,status) &
    !$omp shared(nthreads,ntot,nth,eta1,n1,t,j,djdt,tmp,use_weierstrass)

    nthreads = omp_get_num_threads()
    ithread  = omp_get_thread_num() + 1

    !$omp single
    allocate (nth(0:nthreads + 1), source = 0)
    !$omp end single

    ! coarse grid
    !$omp do schedule(dynamic)
    do i = 1, NT0
      call c(t(i), deta, ddetadt, dpot, ddpotdt)
      eta2 = eta1 + deta
      call fd12(eta2, F = n2, dF1 = dn2deta2)
      if (use_weierstrass) then
        call weierstrass(fd12, ifd12, n1, n2, dpot, j(i), djdn1, djdn2, djddpot, status)
      else
        call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j(i), djdn1, djdn2, djddpot, nquad)
      end if
      djdt(i) = djdn2 * dn2deta2 * ddetadt + djddpot * ddpotdt
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    call vt%init(   NT0, c = 4 * NT0, x = t)
    call vj%init(   NT0, c = 4 * NT0, x = j)
    call vdjdt%init(NT0, c = 4 * NT0, x = djdt)

    ! interval stack
    call vi%init(0, c = 16)

    ! refinement
    !$omp do schedule(dynamic)
    do i = 1, NT0 - 1
      ! print *, i
      ! add interval
      call vi%push(i)
      call vi%push(i+1)

      do while (vi%n > 0)
        ! get interval from stack
        i1   = vi%d(vi%n - 1)
        i2   = vi%d(vi%n    )
        vi%n = vi%n - 2

        ! midpoint
        call vt%push(0.5 * (vt%d(i1) + vt%d(i2)))
        i3 = vt%n

        ! check for numerical resolution of t
        if ((abs(vt%d(i3) - vt%d(i1)) < TTOL * 0.5 * abs(vt%d(i3) + vt%d(i1))) .or. &
            (abs(vt%d(i3) - vt%d(i2)) < TTOL * 0.5 * abs(vt%d(i2) + vt%d(i3)))) then
          vt%n = vt%n - 1
          cycle
        end if

        ! make room for new point
        call vj%resize(i3)
        call vdjdt%resize(i3)

        ! evaluate current at new point
        call c(vt%d(i3), deta, ddetadt, dpot, ddpotdt)
        eta2 = eta1 + deta
        call fd12(eta2, F = n2, dF1 = dn2deta2)
        if (use_weierstrass) then
          call weierstrass(fd12, ifd12, n1, n2, dpot, vj%d(i3), djdn1, djdn2, djddpot, status)
        else
          call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, vj%d(i3), djdn1, djdn2, djddpot, nquad)
        end if
        vdjdt%d(i3) = djdn2 * dn2deta2 * ddetadt + djddpot * ddpotdt

        ! check for refinement
        dt = vt%d(i2) - vt%d(i1)
        j0 = abs(0.5 * (vj%d(i1) + vj%d(i3)))
        dj1 = abs((vj%d(i1) + 0.25 * dt * vdjdt%d(i1)) - (vj%d(i3) - 0.25 * dt * vdjdt%d(i3)))
        if (dj1 > max(RTOL * j0, ATOL)) then
          call vi%push(i1)
          call vi%push(i3)
        end if

        j0 = abs(0.5 * (vj%d(i3) + vj%d(i2)))
        dj1 = abs((vj%d(i3) + 0.25 * dt * vdjdt%d(i3)) - (vj%d(i2) - 0.25 * dt * vdjdt%d(i2)))
        if (dj1 > max(RTOL * j0, ATOL)) then
          call vi%push(i3)
          call vi%push(i2)
        end if
      end do
    end do
    !$omp end do

    ! save number of points generated by this thread
    nth(ithread+1) = vt%n - NT0

    !$omp barrier

    !$omp single

    ! count elements
    nth(0) = 1
    nth(1) = NT0
    do i = 1, nthreads + 1
      nth(i) = nth(i) + nth(i-1)
    end do
    ntot = nth(nthreads + 1) - 1

    ! reallocate global memory
    allocate (tmp(ntot))
    tmp(1:NT0) = t
    call move_alloc(tmp, t)
    allocate (tmp(ntot))
    tmp(1:NT0) = j
    call move_alloc(tmp, j)
    allocate (tmp(ntot))
    tmp(1:NT0) = djdt
    call move_alloc(tmp, djdt)

    !$omp end single

    ! copy local values to global memory
    t(   nth(ithread):nth(ithread+1)-1) = vt%d(   NT0+1:vt%n)
    j(   nth(ithread):nth(ithread+1)-1) = vj%d(   NT0+1:vt%n)
    djdt(nth(ithread):nth(ithread+1)-1) = vdjdt%d(NT0+1:vt%n)

    ! cleanup thread-local memory
    call vt%destruct()
    call vj%destruct()
    call vdjdt%destruct()

    !$omp end parallel

    ! sort by t
    allocate (perm(size(t)))
    call qsort(t, perm = perm)
    j    = j(perm)
    djdt = djdt(perm)

    allocate (dj(size(t)-1), eps(size(t)-1), err(size(t)-1), source = 0.0)
    do i = 1, size(t) - 1
      dt    = t(i+1) - t(i)
      dj(i) = abs((j(i) + 0.5 * dt * djdt(i)) - (j(i+1) - 0.5 * dt * djdt(i+1)))
      eps(i) = 0.5 * (max(RTOL * abs(j(i)), ATOL) + max(RTOL * abs(j(i+1)), ATOL))
      err(i) = max(dj(i) - eps(i), 0.0)
    end do
    max_err = maxval(err) + ATOL

    if (output) then
      open (newunit = funit, file = "j.csv", status = "replace", action = "write")
      do i = 1, size(j)
        write (funit, "(2ES25.16E3)") t(i), j(i)
      end do
      close (funit)

      open (newunit = funit, file = "eps.csv", status = "replace", action = "write")
      write (funit, "(A)") "%lc=4"
      do i = 1, size(eps)
        write (funit, "(2ES25.16E3)") 0.5 * (t(i) + t(i+1)), eps(i)
      end do
      close (funit)

      open (newunit = funit, file = "djdt.csv", status = "replace", action = "write")
      do i = 1, size(djdt)
        write (funit, "(2ES25.16E3)") t(i), djdt(i)
      end do
      close (funit)

      open (newunit = funit, file = "dj.csv", status = "replace", action = "write")
      do i = 1, size(dj)
        write (funit, "(2ES25.16E3)") 0.5 * (t(i) + t(i+1)), max(dj(i), ATOL)
      end do
      close (funit)

      open (newunit = funit, file = "err.csv", status = "replace", action = "write")
      do i = 1, size(err)
        write (funit, "(2ES25.16E3)") 0.5 * (t(i) + t(i+1)), err(i)
      end do
      close (funit)
    end if

  end subroutine

end program
