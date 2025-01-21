program paper2024

  use cl_options_m, only: cl_option_descriptor, cl_option, get_cl_options
  use degen_m,      only: degen_init, degen_get, get_current, fd12, ifd12, weierstrass, DEGEN_DEBUG, DEGEN_TANH_SINH
  use math_m,       only: linspace
  use string_m,     only: string
  use util_m,       only: int2str

  implicit none

  type(cl_option_descriptor) :: desc(7) = [ &
    cl_option_descriptor("n", "NG",           .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("a", "dpot0",        .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("b", "dpot1",        .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("c", "deta0",        .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("d", "deta1",        .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("e", "eta1",         .true.,  .false., .true.,  .true. ), &
    cl_option_descriptor("w", "weierstrass",  .false., .false., .false., .false.) &
  ]
  type(cl_option), allocatable :: opt(:)
  integer                      :: idesc, i, NG
  integer,         allocatable :: iopt(:), jopt(:)
  logical                      :: use_weierstrass
  real                         :: dpot0, dpot1, deta0, deta1, eta1, max_err

  ! command line
  call get_cl_options(desc, opt, iopt, jopt)
  use_weierstrass = .false.
  do idesc = 1, size(desc)
    do i = iopt(idesc), iopt(idesc+1)-1
      if (opt(jopt(i))%arg(1:1) == "M") opt(jopt(i))%arg(1:1) = "-"

      select case (opt(jopt(i))%short)
      case ("n")
        read (opt(jopt(i))%arg, *) NG
      case ("a")
        read (opt(jopt(i))%arg, *) dpot0
      case ("b")
        read (opt(jopt(i))%arg, *) dpot1
      case ("c")
        read (opt(jopt(i))%arg, *) deta0
      case ("d")
        read (opt(jopt(i))%arg, *) deta1
      case ("e")
        read (opt(jopt(i))%arg, *) eta1
      case ("w")
        use_weierstrass = .true.
      end select
    end do
  end do

  DEGEN_DEBUG = .false.

  call degen_init(NG, NG)

  call get_error2D(max_err)
  print "(A,ES25.16E3)", "err = ", max_err

contains

  subroutine get_error(max_err)
    real, intent(out) :: max_err

    integer, parameter :: NT = 64 * 256 + 1

    integer           :: it, funit
    logical           :: status
    real              :: n1, n2, dpot, deta, eta2, djdn1, djdn2, djddpot, nquad
    real, allocatable :: t(:), j(:), j0(:), err(:)

    ! init grid
    allocate (t(NT), j(NT), j0(NT), err(NT))
    t = linspace(0.0, 1.0, NT)

    call fd12(eta1, F = n1)

    DEGEN_TANH_SINH = .true.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(it,dpot,deta,eta2,n2,djdn1,djdn2,djddpot,nquad) &
    !$omp shared(eta1,n1,dpot0,dpot1,deta0,deta1,t,j0)
    do it = 1, NT
      dpot = dpot0 + (dpot1 - dpot0) * t(it)
      deta = deta0 + (deta1 - deta0) * t(it)
      eta2 = eta1 + deta
      call fd12(eta2, F = n2)
      call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j0(it), djdn1, djdn2, djddpot, nquad)
    end do
    !$omp end parallel do

    DEGEN_TANH_SINH = .false.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(it,dpot,deta,eta2,n2,djdn1,djdn2,djddpot,nquad,status) &
    !$omp shared(eta1,n1,dpot0,dpot1,deta0,deta1,t,j0,j,err,use_weierstrass)
    do it = 1, NT
      dpot = dpot0 + (dpot1 - dpot0) * t(it)
      deta = deta0 + (deta1 - deta0) * t(it)
      eta2 = eta1 + deta
      call fd12(eta2, F = n2)
      if (use_weierstrass) then
        call weierstrass(fd12, ifd12, n1, n2, dpot, j(it), djdn1, djdn2, djddpot, status)
        if (.not. status) then
          err(it) = 1e-16
          cycle
        end if
      else
        call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot, j(it), djdn1, djdn2, djddpot, nquad)
      end if

      err(it) = abs(j(it) - j0(it)) / (abs(j0(it)) + 1e-16)
    end do
    !$omp end parallel do

    max_err = maxval(err)

    ! output
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
  end subroutine

  subroutine get_error2D(max_err)
    real, intent(out) :: max_err

    integer, parameter :: NDETA = 192, NDPOT = 256

    integer           :: it, jt, funit
    logical           :: status
    real              :: n1, n2, eta2, djdn1, djdn2, djddpot, nquad
    real, allocatable :: deta(:), dpot(:), j(:,:), j0(:,:), err(:,:)

    ! init grid
    allocate (deta(NDETA), dpot(NDPOT), j(NDETA,NDPOT), j0(NDETA,NDPOT), err(NDETA,NDPOT))
    deta = linspace(deta0, deta1, NDETA)
    dpot = linspace(dpot0, dpot1, NDPOT)

    call fd12(eta1, F = n1)

    DEGEN_TANH_SINH = .true.
    !$omp parallel do schedule(dynamic) default(none) collapse(2) &
    !$omp private(it,jt,eta2,n2,djdn1,djdn2,djddpot,nquad) &
    !$omp shared(eta1,n1,deta,dpot,j0)
    do jt = 1, NDPOT
      do it = 1, NDETA
        eta2 = eta1 + deta(it)
        call fd12(eta2, F = n2)
        call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot(jt), j0(it,jt), djdn1, djdn2, djddpot, nquad)
      end do
    end do
    !$omp end parallel do

    DEGEN_TANH_SINH = .false.
    !$omp parallel do schedule(dynamic) default(none) collapse(2) &
    !$omp private(it,jt,eta2,n2,djdn1,djdn2,djddpot,nquad,status) &
    !$omp shared(eta1,n1,deta,dpot,j0,j,err,use_weierstrass)
    do jt = 1, NDPOT
      do it = 1, NDETA
        eta2 = eta1 + deta(it)
        call fd12(eta2, F = n2)
        if (use_weierstrass) then
          call weierstrass(fd12, ifd12, n1, n2, dpot(jt), j(it,jt), djdn1, djdn2, djddpot, status)
          if (.not. status) then
            err(it,jt) = 1e-16
            cycle
          end if
        else
          call get_current(fd12, ifd12, 0.0, sqrt(0.125), n1, n2, dpot(jt), j(it,jt), djdn1, djdn2, djddpot, nquad)
        end if

        err(it,jt) = abs(j(it,jt) - j0(it,jt)) / (abs(j0(it,jt)) + 1e-16)
      end do
    end do
    !$omp end parallel do

    max_err = maxval(err)

    ! output
    open (newunit = funit, file = "deta.csv", status = "replace", action = "write")
    do it = 1, NDETA
      write (funit, "(ES25.16E3)") deta(it)
    end do
    close (funit)
    open (newunit = funit, file = "dpot.csv", status = "replace", action = "write")
    do jt = 1, NDPOT
      write (funit, "(ES25.16E3)") dpot(jt)
    end do
    close (funit)

    open (newunit = funit, file = "j0.csv", status = "replace", action = "write")
    do it = 1, NDETA
      write (funit, "(" // int2str(NDPOT) // "ES25.16E3)") j0(it,:)
    end do
    close (funit)
    open (newunit = funit, file = "j.csv", status = "replace", action = "write")
    do it = 1, NDETA
      write (funit, "(" // int2str(NDPOT) // "ES25.16E3)") j(it,:)
    end do
    close (funit)
    open (newunit = funit, file = "err.csv", status = "replace", action = "write")
    do it = 1, NDETA
      write (funit, "(" // int2str(NDPOT) // "ES25.16E3)") max(err(it,:), 1e-16)
    end do
    close (funit)
  end subroutine

end program
