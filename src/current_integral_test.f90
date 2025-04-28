program current_integral_test

  use current_integral_m,   only: CURRENT_INTEGRAL_DEBUG, current_integral_get
  use distribution_table_m, only: distribution_table
  use error_m,              only: program_error
  use fukushima_m,          only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h, dfdm9h, ifd1h
  use ieee_arithmetic,      only: ieee_is_finite
  use math_m,               only: expm1, linspace, PI

  implicit none

  type(distribution_table) :: tab
  logical                  :: use_table, use_fermi_dirac

  use_table       = .true.
  use_fermi_dirac = .true.

  ! call tab%init("F12", parabolic_dos, fermi_dirac, -100.0, 500.0, 3)
  call tab%init("FGF", gauss_dos, fermi_dirac, -100.0, 500.0, 3)

  call test_rect(0.0, -20.0, 20.0, 201, -50.0, 50.0, 501)

contains

  subroutine test_single(eta, dpot)
    real, intent(in) :: eta(2)
    real, intent(in) :: dpot

    integer, parameter :: NN = 1001

    integer           :: i, funit
    real              :: n(2), dum, j, djdn(2), djddpot, j2, dum2(2)
    real, allocatable :: n1_ar(:), n2_ar(:), dpot_ar(:)

    call distribution(eta(1), 0, n(1), dum)
    call distribution(eta(2), 0, n(2), dum)

    CURRENT_INTEGRAL_DEBUG = .true.
    call current_integral_get(distribution, inv_distribution, n, dpot, j, djdn, djddpot)
    CURRENT_INTEGRAL_DEBUG = .false.
    print "(A,ES25.16E3)", "j       = ", j
    print "(A,ES25.16E3)", "djdn1   = ", djdn(1)
    print "(A,ES25.16E3)", "djdn2   = ", djdn(2)
    print "(A,ES25.16E3)", "djddpot = ", djddpot

    n1_ar   = linspace(n(1)*0.99, n(1)*1.01, NN)
    n2_ar   = linspace(n(2)*0.99, n(2)*1.01, NN)
    dpot_ar = linspace(dpot*0.99, dpot*1.01, NN)

    open (newunit = funit, file = "j_n1.csv", status = "replace", action = "write")
    do i = 1, NN
      call current_integral_get(distribution, inv_distribution, [n1_ar(i), n(2)], dpot, j2, dum2, dum)
      write (funit, "(2ES25.16E3)") n1_ar(i), j2 - (j + djdn(1) * (n1_ar(i) - n(1)))
    end do
    close (funit)

    open (newunit = funit, file = "j_n2.csv", status = "replace", action = "write")
    do i = 1, NN
      call current_integral_get(distribution, inv_distribution, [n(1), n2_ar(i)], dpot, j2, dum2, dum)
      write (funit, "(2ES25.16E3)") n2_ar(i), j2 - (j + djdn(2) * (n2_ar(i) - n(2)))
    end do
    close (funit)

    open (newunit = funit, file = "j_dpot.csv", status = "replace", action = "write")
    do i = 1, NN
      call current_integral_get(distribution, inv_distribution, n, dpot_ar(i), j2, dum2, dum)
      write (funit, "(2ES25.16E3)") dpot_ar(i), j2 - (j + djddpot * (dpot_ar(i) - dpot))
    end do
    close (funit)

  end subroutine

  subroutine test_line(eta1, deta1, deta2, dpot1, dpot2, NN)
    real,    intent(in) :: eta1
    real,    intent(in) :: deta1
    real,    intent(in) :: deta2
    real,    intent(in) :: dpot1
    real,    intent(in) :: dpot2
    integer, intent(in) :: NN

    integer           :: i, funit
    real              :: n(2), dum, dum2(2)
    real, allocatable :: deta(:), dpot(:), j(:)

    deta = linspace(deta1, deta2, NN)
    dpot = linspace(dpot1, dpot2, NN)

    call distribution(eta1, 0, n(1), dum)

    allocate (j(NN))
    CURRENT_INTEGRAL_DEBUG = .false.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i,dum,dum2) firstprivate(n) &
    !$omp shared(eta1,deta,dpot,NN,j)
    do i = 1, NN
      call distribution(eta1 + deta(i), 0, n(2), dum)
      call current_integral_get(distribution, inv_distribution, n, dpot(i), j(i), dum2, dum)
      print "(I6,A,I6)", i, " / ", NN
    end do
    !$omp end parallel do

    open (newunit = funit, file = "j.csv", status = "replace", action = "write")
    do i = 1, NN
      write (funit, "(3ES25.16E3)") deta(i), dpot(i), j(i)
    end do
    close (funit)
  end subroutine

  subroutine test_rect(eta1, deta1, deta2, Ndeta, dpot1, dpot2, Ndpot)
    real,    intent(in) :: eta1
    real,    intent(in) :: deta1
    real,    intent(in) :: deta2
    integer, intent(in) :: Ndeta
    real,    intent(in) :: dpot1
    real,    intent(in) :: dpot2
    integer, intent(in) :: Ndpot

    integer           :: ii, jj, funit
    real              :: n(2), dum, dum2(2)
    real, allocatable :: deta(:), dpot(:), j(:,:)

    deta = linspace(deta1, deta2, Ndeta)
    dpot = linspace(dpot1, dpot2, Ndpot)

    call distribution(eta1, 0, n(1), dum)

    open (newunit = funit, file = "deta.csv", status = "replace", action = "write")
    do ii = 1, Ndeta
      write (funit, "(1ES25.16E3)") deta(ii)
    end do
    close (funit)
    open (newunit = funit, file = "dpot.csv", status = "replace", action = "write")
    do ii = 1, Ndpot
      write (funit, "(1ES25.16E3)") dpot(ii)
    end do
    close (funit)

    allocate (j(Ndeta,Ndpot))
    CURRENT_INTEGRAL_DEBUG = .false.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(ii,jj,dum,dum2) firstprivate(n) &
    !$omp shared(eta1,deta,dpot,Ndeta,Ndpot,j)
    do ii = 1, Ndeta
      call distribution(eta1 + deta(ii), 0, n(2), dum)
      do jj = 1, Ndpot
        ! print "(2I6)", ii, jj
        call current_integral_get(distribution, inv_distribution, n, dpot(jj), j(ii,jj), dum2, dum)
      end do
      print "(I6,A,I6)", ii, " / ", Ndeta
    end do
    !$omp end parallel do

    open (newunit = funit, file = "j.csv", status = "replace", action = "write")
    do jj = 1, Ndpot
      do ii = 1, Ndeta
        write (funit, "(ES25.16E3)", advance = "no") j(ii,jj)
      end do
      write (funit, *)
    end do
    close (funit)
  end subroutine

  subroutine distribution(eta, k, F, dFdeta)
    real,    intent(in)  :: eta
    integer, intent(in)  :: k
    real,    intent(out) :: F
    real,    intent(out) :: dFdeta

    real, parameter :: G0 = 1.0 / gamma(1.5)
    real, parameter :: G1 = 1.0 / gamma(0.5)
    real, parameter :: G2 = 1.0 / gamma(-0.5)
    real, parameter :: G3 = 1.0 / gamma(-1.5)
    real, parameter :: G4 = 1.0 / gamma(-2.5)
    real, parameter :: G5 = 1.0 / gamma(-3.5)

    if (use_table) then
      ! Use distribution from table
      call tab%get(eta, k, F, dFdeta)
    elseif (use_fermi_dirac) then
      ! Fermi-Dirac
      select case (k)
      case (0)
        F      = fd1h(eta) * G0
        dFdeta = fdm1h(eta) * G1
      case (1)
        F      = fdm1h(eta) * G1
        dFdeta = fdm3h(eta) * G2
      case (2)
        F      = fdm3h(eta) * G2
        dFdeta = fdm5h(eta) * G3
      case (3)
        F      = fdm5h(eta) * G3
        dFdeta = fdm7h(eta) * G4
      case default
        F      = 0
        dFdeta = 0
        call program_error("not implemented for k outside of range 0..3")
      end select
    else
      ! Maxwell-Boltzmann
      F      = exp(-eta)
      dFdeta = F
    end if

  end subroutine

  subroutine inv_distribution(F, eta, detadF)
    real, intent(in)  :: F
    real, intent(out) :: eta
    real, intent(out) :: detadF

    real, parameter :: G = gamma(1.5)

    if (use_table) then
      call tab%inv(F, eta, detadF)
    else
      call ifd1h(F * G, eta, detadF)
      detadF = detadF * G
    end if
  end subroutine

  function parabolic_dos(t) result(Z)
    real(kind=16), intent(in) :: t
    real(kind=16)             :: Z

    if (t > 0) then
      Z = 2 * sqrt(t / PI)
    else
      Z = 0
    end if
  end function

  function gauss_dos(t) result(Z)
    real(kind=16), intent(in) :: t
    real(kind=16)             :: Z

    real, parameter :: SIGMA = 4.0

    Z = 1 / (sqrt(2 * PI) * SIGMA) * exp(-t**2 / (2 * SIGMA**2))
  end function

  function fermi_dirac(u, k) result(f)
    real(kind=16), intent(in) :: u
    integer,       intent(in) :: k
    real(kind=16)             :: f

    real(kind=16) :: e, f0

    e  = exp(u)
    f0 = 1 / (1 + e)
    f  = 0

    select case (k)
    case (0)
      f = f0
    case (1)
      if (ieee_is_finite(e)) f = (- e * f0) * f0
    case (2)
      if (ieee_is_finite(e)) f = ((e * f0) * (expm1(u) * f0)) * f0
    case (3)
      if (ieee_is_finite(e)) f = (e * f0) * (- 1 + 6 * (e * f0) * (1 - e * f0)) * f0
    case (4)
      if (ieee_is_finite(e)) f = (- (e * f0) + 14 * (e * f0)**2 - 36 * (e * f0)**3 + 24 * (e * f0)**4) * f0
    end select
  end function

end program
