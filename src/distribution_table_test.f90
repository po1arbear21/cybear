program distribution_table_test

  use distribution_table_m
  use ieee_arithmetic, only: ieee_is_finite
  use fukushima_m,     only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h
  use math_m,          only: linspace, expm1, PI

  implicit none

  type(distribution_table) :: tab
  real :: eta, F, dFdeta, eta2, deta2dF

  call tab%init("F12", parabolic_dos, 0.0, ieee_value(1.0, IEEE_POSITIVE_INF), fermi_dirac, -100.0, 500.0, 3)

  eta = 243.235225456
  call tab%get(eta, 0, F, dFdeta)
  call tab%inv(F, eta2, deta2dF)

  print *, eta, eta2

  ! call tab%output("F12.csv", -100.0, 500.0, 0, 10001)
  ! call tab%output("Fm12.csv", -100.0, 500.0, 1, 10001)
  ! call tab%output("Fm32.csv", -100.0, 500.0, 2, 10001)
  ! call tab%output("Fm52.csv", -100.0, 500.0, 3, 10001)

  ! do i = 1, size(tab%eta)
  !   print *, i, tab%eta(i), tab%lg(0,i)
  ! end do

  ! call tab%init("FMB", parabolic_dos, 0.0, ieee_value(1.0, IEEE_POSITIVE_INF), maxwell_boltzmann, -100.0, 100.0, 3)
  ! call tab%output("FMB.csv", -100.0, 100.0, 0, 10001)
  ! call tab%output("FMB1.csv", -100.0, 100.0, 1, 10001)
  ! call tab%output("FMB2.csv", -100.0, 100.0, 2, 10001)
  ! call tab%output("FMB3.csv", -100.0, 100.0, 3, 10001)

  ! call tab%init("FGF", gauss_dos, 0.0, ieee_value(1.0, IEEE_POSITIVE_INF), fermi_dirac, -100.0, 500.0, 3)
  ! call tab%output("FGF.csv", -100.0, 500.0, 0, 10001)
  ! call tab%output("FGF1.csv", -100.0, 500.0, 1, 10001)
  ! call tab%output("FGF2.csv", -100.0, 500.0, 2, 10001)
  ! call tab%output("FGF3.csv", -100.0, 500.0, 3, 10001)

contains

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

    real, parameter :: SIGMA = 5.0

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

  function maxwell_boltzmann(u, k) result(f)
    real(kind=16), intent(in) :: u
    integer,       intent(in) :: k
    real(kind=16)             :: f

    f = exp(-u)
    if (mod(k, 2) /= 0) f = - f
  end function

end program
