m4_include(util/macro.f90.inc)

module current_integral_m
  !! get edge-current by solving integral equation

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
  use, intrinsic :: iso_fortran_env, only: real128

  use error_m, only: program_error
  use math_m,  only: ber, log1p, expm1
  use quad_m,  only: quad
  use util_m,  only: int2str

  implicit none

  private
  public :: current_integral_get

  logical, public :: CURRENT_INTEGRAL_DEBUG   = .false.
  logical, public :: CURRENT_INTEGRAL_PLOTRES = .false.

  interface
    subroutine distribution(eta, k, F, dFdeta)
      !! k-th derivative of cumulative distribution function
      real,    intent(in)  :: eta
        !! normalized chemical potential
      integer, intent(in)  :: k
        !! derivative (can be 0)
      real,    intent(out) :: F
        !! output k-th derivative of distribution function
      real,    intent(out) :: dFdeta
        !! output derivative of F wrt eta (used for Newton iteration, might be slightly different from F with k + 1)
    end subroutine

    subroutine inv_distribution(F, eta, detadF)
      !! inverse of cumulative distribution function
      real, intent(in)  :: F
        !! value of cumulative distribution function
      real, intent(out) :: eta
        !! output normalized chemical potential corresponding to F
      real, intent(out) :: detadF
        !! output derivative of eta wrt F
    end subroutine

    subroutine int_distribution(eta, k, I, dIdeta)
      !! integrate cumulative distribution function
      real,       intent(in)  :: eta(2)
        !! integration bounds
      integer,    intent(in)  :: k
        !! exponent
      real,       intent(out) :: I
        !! output integral over dist^k
      real,       intent(out) :: dIdeta(2)
        !! output derivatives of I wrt eta
    end subroutine

  end interface

  integer,      parameter :: CASE0A  = 1
  integer,      parameter :: CASE0B  = 2
  integer,      parameter :: CASE0C  = 3
  integer,      parameter :: CASE1A  = 4
  integer,      parameter :: CASE1B  = 5
  integer,      parameter :: CASE2A  = 6
  integer,      parameter :: CASE2B  = 7
  character(2), parameter :: CASENAME(7) = ["0A", "0B", "0C", "1a", "1b", "2a", "2b"]
  real,         parameter :: CASETOL = 1e-3

  integer, parameter :: MIN_IT = 2
  integer, parameter :: MAX_IT = 20
  real,    parameter :: RTOL   = 2e-14
  real,    parameter :: ATOL   = 1e-24

  real,    parameter :: EPS = epsilon(1.0)

contains

  recursive subroutine current_integral_get(dist, inv_dist, int_dist, n, dpot, j, djdn, djddpot)
    !! get edge current
    procedure(distribution)     :: dist
      !! cumulative distribution function
    procedure(inv_distribution) :: inv_dist
      !! inverse of cumulative distribution function
    procedure(int_distribution) :: int_dist
      !! integral over cumulative distribution function
    real,           intent(in)  :: n(2)
      !! normalized densities at left/right edge end-point
    real,           intent(in)  :: dpot
      !! normalized potential drop
    real,           intent(out) :: j
      !! output normalized current density on edge
    real,           intent(out) :: djdn(2)
      !! output derivative of j wrt n
    real,           intent(out) :: djddpot
      !! output derivative of j wrt dpot

    integer :: cs, it
    real    :: eta(2), detadn(2), deta, g, tmp(2)
    real    :: jmin, jmax, t, tmin, tmax, told, dt, err
    real    :: r, drdt, drdeta(2), drddpot, djdt, djdeta(2)
    real    :: FL2, dFL2, epseta, depseta, sgn, A, dAdeta1, dAddpot, B, dBdeta1, C, dCdeta1, dCddpot

    ! flip edge direction if potential drop is negative
    if (dpot < 0) then
      call current_integral_get(dist, inv_dist, int_dist, [n(2), n(1)], - dpot, j, djdn, djddpot)

      ! flip current (dj/ddpot = d(-j)/d(-dpot) unchanged)
      j    = - j
      djdn = - djdn(2:1:-1)
      return
    end if

    ! get eta
    call inv_dist(n(1), eta(1), detadn(1))
    call inv_dist(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)
    sgn  = sign(1.0, deta)

    ! get case
    if (abs(dpot) > CASETOL) then
      if (deta - dpot > dpot * CASETOL) then
        cs = CASE0A
      elseif (deta - dpot >= - dpot * CASETOL) then
        cs = CASE1A
      elseif (deta > CASETOL) then
        cs = CASE0B
      elseif (deta >= - CASETOL) then
        cs = CASE1B
      else
        cs = CASE0C
      end if
    else
      if (abs(deta) > CASETOL) then
        cs = CASE2A
      else
        cs = CASE2B
      end if
    end if

    ! debug info
    if (CURRENT_INTEGRAL_DEBUG) then
      print *
      print "(A,ES24.16E3)", "n1   = ", n(1)
      print "(A,ES24.16E3)", "n2   = ", n(2)
      print "(A,ES24.16E3)", "eta1 = ", eta(1)
      print "(A,ES24.16E3)", "eta2 = ", eta(2)
      print "(A,ES24.16E3)", "dpot = ", dpot
      print "(A,ES24.16E3)", "deta = ", deta
      print "(2A)",          "case = ", CASENAME(cs)
    end if

    ! calculate current based on case
    select case (cs)
    case (CASE0A, CASE0B, CASE0C) ! general case
      ! initial range for j
      select case (cs)
      case (CASE0A)
        call int_dist([eta(2) - dpot, eta(2)], -1, jmin, tmp)
        jmin = dpot * (dpot - deta) / jmin
        jmax = n(1) * (dpot - deta)

      case (CASE0B)
        jmin = n(1) * (dpot - deta)
        jmax = min(n(1) * dpot, n(2) * (dpot - deta))

      case (CASE0C)
        jmin = max(n(1) * dpot, n(2) * (dpot - deta))
        jmax = n(1) * (dpot - deta)
      end select

      ! initial guess for j: enhanced diffusion approximation
      if (abs(deta) < 1e-6) then
        g = n(1) * detadn(1)
      else
        g = deta / log(n(2) / n(1))
      end if
      j = (ber(-dpot / g) * n(1) - ber(dpot / g) * n(2)) * g

      ! transform initial range and guess
      if (deta > 0) then ! CASE0A, CASE0B
        tmin = - log1p(-jmin / (dpot * n(1)))
        if (jmax < dpot * n(1)) then
          tmax = - log1p(-jmax / (dpot * n(1)))
        else
          tmax = 700.0
        end if
        if (j < dpot * n(1)) then
          t = - log1p(- j / (dpot * n(1)))
        else
          t = 700.0
        end if
      else ! CASE0C
        tmin = - log(jmax / (dpot * n(1)) - 1)
        if (jmin > dpot * n(1)) then
          tmax = - log(jmin / (dpot * n(1)) - 1)
        else
          tmax = 700.0
        end if
        if (j > dpot * n(1)) then
          t = - log(j / (dpot * n(1)) - 1)
        else
          t = 700.0
        end if
      end if

      ! adjust initial guess if not in range
      if ((t < tmin) .or. (t > tmax)) t = 0.5 * (tmin + tmax)

      ! get epseta such that F can be approximated by a linear function around eta(1) +/- epseta
      call dist(eta(1), 2, FL2, dFL2)                 ! use second derivative to estimate error of first order Taylor
      epseta  = 0.5 * sqrt(2 * EPS * n(1) / abs(FL2)) ! 0.5 is a safety factor
      depseta = (1.0 / detadn(1) - n(1) * dFL2 / FL2) * EPS / (4 * epseta * abs(FL2))

      ! get factors for analytical part of integration close to eta(1)
      A        = detadn(1) / dpot**2 ! = 1 / (dpot**2 * dndeta(1))
      dAdeta1  = - FL2 * (detadn(1) / dpot)**2
      dAddpot  = - 2 * detadn(1) / dpot**3
      B        = epseta / (detadn(1) * n(1))
      dBdeta1  = (depseta / detadn(1) + epseta * (FL2 - 1 / (detadn(1)**2 * n(1)))) / n(1)
      C        = sgn * epseta / dpot
      dCdeta1  = sgn * depseta / dpot
      dCddpot  = - sgn * epseta / dpot**2

      if (CURRENT_INTEGRAL_DEBUG) then
        print "(A,ES24.16E3)", "tmin = ", tmin
        print "(A,ES24.16E3)", "tmax = ", tmax
        print "(A,ES24.16E3)", "t    = ", t
      end if

      if (CURRENT_INTEGRAL_PLOTRES) then
        block
          use math_m, only: linspace
          integer, parameter :: NN = 101

          integer           :: ii, funit
          real, allocatable :: tt(:)

          ! output residual over total range
          tt = linspace(tmin, tmax, NN)
          open (newunit = funit, file = "res.csv", status = "replace", action = "write")
          do ii = 1, NN
            call residual(tt(ii), r, drdt, drdeta, drddpot, j, djdt, djdeta, djddpot)
            write (funit, "(2ES24.16E3)") tt(ii), r
          end do
          write (funit, *)

          ! output initial guess
          write (funit, "(A)") "% lt=0 mt=2 mc=4 ms=4"
          call residual(t, r, drdt, drdeta, drddpot, j, djdt, djdeta, djddpot)
          write (funit, "(2ES24.16E3)") t, r
          close (funit)
        end block
      end if

      ! Newton iteration
      err = huge(1.0)
      it  = 0
      do while (((it < MIN_IT) .or. (err > RTOL * abs(j) + ATOL)) .and. (tmin < tmax))
        it = it + 1
        if (it > MAX_IT) then
          print *
          print "(A,ES24.16E3)", "n1   = ", n(1)
          print "(A,ES24.16E3)", "n2   = ", n(2)
          print "(A,ES24.16E3)", "eta1 = ", eta(1)
          print "(A,ES24.16E3)", "eta2 = ", eta(2)
          print "(A,ES24.16E3)", "dpot = ", dpot
          print "(A,ES24.16E3)", "deta = ", deta
          print "(2A)",          "case = ", CASENAME(cs)
          call program_error("No convergence after " // int2str(MAX_IT) // " iterations")
        end if

        ! evaluate residual and get Newton update
        call residual(t, r, drdt, drdeta, drddpot, j, djdt, djdeta, djddpot)

        ! Newton update
        dt  = - r / drdt
        err = abs(dt)

        ! update bounds (assume monotonic behaviour)
        if (dt > 0) then
          tmin = min(t, tmax)
        else
          tmax = max(t, tmin)
        end if

        ! update solution
        told = t
        t    = t + dt

        ! bisection
        if ((t < tmin) .or. (t > tmax) .or. ((told == tmin) .and. (t == tmax))) then
          if (CURRENT_INTEGRAL_DEBUG) then
            print "(A)", "bisection"
            print "(A,ES24.16E3)", "  tmin = ", tmin
            print "(A,ES24.16E3)", "  tmax = ", tmax
            print "(A,ES24.16E3)", "  told = ", told
            print "(A,ES24.16E3)", "  t    = ", t
            print "(A,ES24.16E3)", "  tnew = ", 0.5 * (tmin + tmax)
          end if
          t   = 0.5 * (tmin + tmax)
          err = min(err, tmax - tmin)
        end if

        ! update current error
        err = abs(djdt) * err

        if (CURRENT_INTEGRAL_DEBUG) then
          print "(I6,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3)", it, ": t = ", t, ", j = ", j, ", err/j = ", err / abs(j)
        end if
      end do

      ! get current + derivatives with implicit differentiation
      call residual(t, r, drdt, drdeta, drddpot, j, djdt, djdeta, djddpot)
      djdeta  = djdeta  - drdeta  * djdt / drdt
      djddpot = djddpot - drddpot * djdt / drdt

    case (CASE1A)
      call case_1a(j, djdeta, djddpot)

    case (CASE1B)
      call case_1b(j, djdeta, djddpot)

    case (CASE2A)
      call case_2a(j, djdeta, djddpot)

    case (CASE2B)
      call case_2b(j, djdeta, djddpot)

    end select

    ! derivatives wrt densities
    djdn = djdeta * detadn

  contains

    subroutine case_1a(j, djdeta, djddpot)
      !! dpot >> 0 and deta ~ dpot
      real,       intent(out) :: j
        !! output normalized current density on edge
      real,       intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real,       intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer       :: k
      real          :: I(-5:0), dI(2,-5:0), jc(1:5), djc(2,1:5)
      real(real128) :: I_16(-5:0), jc_16(1:5), djc_16(2,1:5), tmp_16, dtmp_16(2), j_16

      ! get I_k = integral_eta1^eta2 F(eta)^k deta for k = -5 to 0
      do k = -5, -1
        call int_dist(eta, k, I(k), dI(:,k))
      end do
      I(0)    = deta
      dI(:,0) = [-1.0, 1.0]
      I_16 = I

      ! first Taylor coefficient
      jc_16(1)    = I_16(0) / I_16(-1)
      djc_16(:,1) = (dI(:,0) - dI(:,-1) * I_16(0) / I_16(-1)) / I_16(-1)

      ! second Taylor coefficient
      tmp_16   = I_16(-1)**2 - I_16(-2)*I_16(0)
      dtmp_16  = 2*I_16(-1)*dI(:,-1) - dI(:,-2)*I_16(0) - I_16(-2)*dI(:,0)
      jc_16(2) = tmp_16 / I_16(-1)**3
      djc_16(:,2) = (dtmp_16 - 3 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**3

      ! third Taylor coefficient
      tmp_16   = 2*I_16(-2)**2*I_16(0) - I_16(-1)**2*I_16(-2) - I_16(-1)*I_16(-3)*I_16(0)
      dtmp_16  = 4*I_16(-2)*dI(:,-2)*I_16(0) + 2*I_16(-2)**2*dI(:,0) - 2*I_16(-1)*dI(:,-1)*I_16(-2) - I_16(-1)**2*dI(:,-2) &
        &      - dI(:,-1)*I_16(-3)*I_16(0) - I_16(-1)*dI(:,-3)*I_16(0) - I_16(-1)*I_16(-3)*dI(:,0)
      jc_16(3) = tmp_16 / I_16(-1)**5
      djc_16(:,3) = (dtmp_16 - 5 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**5

      ! fourth Taylor coefficient
      tmp_16   = 2*I_16(-1)**2*I_16(-2)**2 - I_16(-1)**3*I_16(-3) - 5*I_16(-2)**3*I_16(0) &
        &      - I_16(-1)**2*I_16(-4)*I_16(0) + 5*I_16(-1)*I_16(-2)*I_16(-3)*I_16(0)
      dtmp_16  = 4*I_16(-1)*dI(:,-1)*I_16(-2)**2 + 4*I_16(-1)**2*I_16(-2)*dI(:,-2) - 3*I_16(-1)**2*dI(:,-1)*I_16(-3) &
        &      - I_16(-1)**3*dI(:,-3) - 15*I_16(-2)**2*dI(:,-2)*I_16(0) - 5*I_16(-2)**3*dI(:,0) &
        &      - 2*I_16(-1)*dI(:,-1)*I_16(-4)*I_16(0) - I_16(-1)**2*dI(:,-4)*I_16(0) - I_16(-1)**2*I_16(-4)*dI(:,0) &
        &      + 5*dI(:,-1)*I_16(-2)*I_16(-3)*I_16(0) + 5*I_16(-1)*dI(:,-2)*I_16(-3)*I_16(0) &
        &      + 5*I_16(-1)*I_16(-2)*dI(:,-3)*I_16(0) + 5*I_16(-1)*I_16(-2)*I_16(-3)*dI(:,0)
      jc_16(4) = tmp_16 / I_16(-1)**7
      djc_16(:,4) = (dtmp_16 - 7 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**7

      ! fifth Taylor coefficient
      tmp_16   = 14*I_16(0)*I_16(-2)**4 - I_16(-1)**4*I_16(-4) - 5*I_16(-1)**2*I_16(-2)**3 &
        &      - I_16(0)*I_16(-1)**3*I_16(-5) + 5*I_16(-1)**3*I_16(-2)*I_16(-3) + 3*I_16(0)*I_16(-1)**2*I_16(-3)**2 &
        &      - 21*I_16(0)*I_16(-1)*I_16(-2)**2*I_16(-3) + 6*I_16(0)*I_16(-1)**2*I_16(-2)*I_16(-4)
      dtmp_16  = 14*dI(:,0)*I_16(-2)**4 + 56*I_16(0)*I_16(-2)**3*dI(:,-2) - 4*I_16(-1)**3*dI(:,-1)*I_16(-4) &
        &      - I_16(-1)**4*dI(:,-4) - 10*I_16(-1)*dI(:,-1)*I_16(-2)**3 - 15*I_16(-1)**2*I_16(-2)**2*dI(:,-2) &
        &      - dI(:,0)*I_16(-1)**3*I_16(-5) - 3*I_16(0)*I_16(-1)**2*dI(:,-1)*I_16(-5) - I_16(0)*I_16(-1)**3*dI(:,-5) &
        &      + 15*I_16(-1)**2*dI(:,-1)*I_16(-2)*I_16(-3) + 5*I_16(-1)**3*dI(:,-2)*I_16(-3) &
        &      + 5*I_16(-1)**3*I_16(-2)*dI(:,-3) + 3*dI(:,0)*I_16(-1)**2*I_16(-3)**2 &
        &      + 6*I_16(0)*I_16(-1)*dI(:,-1)*I_16(-3)**2 + 6*I_16(0)*I_16(-1)**2*I_16(-3)*dI(:,-3) &
        &      - 21*dI(:,0)*I_16(-1)*I_16(-2)**2*I_16(-3) - 21*I_16(0)*dI(:,-1)*I_16(-2)**2*I_16(-3) &
        &      - 42*I_16(0)*I_16(-1)*I_16(-2)*dI(:,-2)*I_16(-3) - 21*I_16(0)*I_16(-1)*I_16(-2)**2*dI(:,-3) &
        &      + 6*dI(:,0)*I_16(-1)**2*I_16(-2)*I_16(-4) + 12*I_16(0)*I_16(-1)*dI(:,-1)*I_16(-2)*I_16(-4) &
        &      + 6*I_16(0)*I_16(-1)**2*dI(:,-2)*I_16(-4) + 6*I_16(0)*I_16(-1)**2*I_16(-2)*dI(:,-4)
      jc_16(5) = tmp_16 / I_16(-1)**9
      djc_16(:,5) = (dtmp_16 - 9 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**9

      ! evaluate Taylor series
      j_16     = 0
      djdeta  = 0
      djddpot = 0
      jc      = real(jc_16)
      djc     = real(djc_16)
      do k = 1, 5
        if (.not. ieee_is_finite(jc_16(k))) cycle
        j_16     = j_16 + jc_16(k) * (dpot - I_16(0))**k
        djdeta  = djdeta + (djc(:,k) * (dpot - I(0)) - k * jc(k) * dI(:,0)) * (dpot - I(0))**(k-1)
        djddpot = djddpot + jc(k) * k * (dpot - I(0))**(k-1)
      end do
      j = real(j_16)
    end subroutine

    subroutine case_1b(j, djdeta, djddpot)
      !! |dpot| >> 0 and deta ~ 0
      real,       intent(out) :: j
        !! output normalized current density on edge
      real,       intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real,       intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer :: k
      real    :: etam, F(0:3), dF(0:3), G(2), dG(2), e, dedetam, deddpot, u, dudetam, duddpot, v, dvdetam, dvddpot
      real    :: jc(0:3), djcdetam(0:3), djcddpot(0:3)

      ! evaluate distribution and derivatives at mean eta
      etam = 0.5 * (eta(1) + eta(2))
      do k = 0, 3
        call dist(etam, k, F(k), dF(k))
      end do

      ! F(1) and F(2) scaled by 1/F(0)
      do k = 1, 2
        G(k) = F(k) / F(0)
        dG(k) = (dF(k) - G(k) * dF(0)) / F(0)
      end do

      ! abbreviation
      e = 1 / tanh(0.5 * G(1) * dpot)
      u = 1 / sinh(0.5 * G(1) * dpot)**2
      dedetam = - 0.5 * dpot * dG(1) * u
      deddpot = - 0.5 * G(1) * u

      ! zero-th Taylor coefficient
      jc(0)       = dpot * F(0)
      djcdetam(0) = dpot * dF(0)
      djcddpot(0) = F(0)

      ! first Taylor coefficient
      jc(1)       = - 0.5 * F(1) * dpot * e
      djcdetam(1) = - 0.5 * dpot * (dF(1) * e + F(1) * dedetam)
      djcddpot(1) = - 0.5 * F(1) * (e + dpot * deddpot)

      ! second Taylor coefficient = 0.25 * dpot * (0.5 * F(2) + (F(2) - F(1)**2/F(0)) * (1 - e**2) * (0.5*dpot*F(1)/F(0)*e - 1))
      u           = (F(2) - F(1) * G(1)) * (1 - e**2)
      dudetam     = (dF(2) - dF(1) * G(1) - F(1) * dG(1)) * (1 - e**2) - 2 * (F(2) - F(1) * G(1)) * e * dedetam
      duddpot     = - 2 * (F(2) - F(1) * G(1)) * e * deddpot
      v           = 0.5 * dpot * G(1) * e - 1
      dvdetam     = 0.5 * dpot * (dG(1) * e + G(1) * dedetam)
      dvddpot     = 0.5 * G(1) * (e + dpot * deddpot)
      jc(2)       = 0.25 * dpot * (0.5 * F(2) + u * v)
      djcdetam(2) = 0.25 * dpot * (0.5 * dF(2) + dudetam * v + u * dvdetam)
      djcddpot(2) = 0.25 * (0.5 * F(2) + u * v + dpot * (duddpot * v + u * dvddpot))

      ! third Taylor coefficient
      u           = e * (1 - e**2)
      dudetam     = (1 - 3*e**2) * dedetam
      duddpot     = (1 - 3*e**2) * deddpot
      v           = G(1)**2 * F(1) + 0.25 * (F(2)/F(1))*F(2) - 1.5 * G(1) * F(2)
      dvdetam     = G(1) * (2*F(1)*dG(1) + G(1)*dF(1)) + 0.25 * (F(2)/F(1)) * (2*dF(2) - (F(2)/F(1)) * dF(1)) - 1.5 * (dG(1) * F(2) + G(1) * dF(2))
      jc(3)       = 0.25 * dpot * u * v
      djcdetam(3) = 0.25 * dpot * (dudetam * v + u * dvdetam)
      djcddpot(3) = 0.25 * (u + dpot * duddpot) * v
      dudetam     = dudetam * (1 - 2 * e**2) - 4 * u * e * dedetam
      duddpot     = duddpot * (1 - 2 * e**2) - 4 * u * e * deddpot
      u           = u * (1 - 2 * e**2)
      v           = G(1) * (G(1)**2*F(2) - 0.5*G(2)*F(2) - 0.5*G(1)**3*F(1))
      dvdetam     = G(1)**2*(3*F(2)*dG(1) + G(1)*dF(2)) - 0.5*(dG(1)*G(2)*F(2) + G(1)*dG(2)*F(2) + G(1)*G(2)*dF(2)) - 0.5*G(1)**3*(4*dG(1)*F(1) + G(1)*dF(1))
      jc(3)       = jc(3) + dpot**3/16 * u * v
      djcdetam(3) = djcdetam(3) + dpot**3/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot**2/16 * (3 * u + dpot * duddpot) * v
      u           = (1 - e**2) * (1 - 5*e**2)
      dudetam     = - (12 - 20*e**2) * e * dedetam
      duddpot     = - (12 - 20*e**2) * e * deddpot
      v           = G(2) * F(2)
      dvdetam     = dG(2) * F(2) + G(2) * dF(2)
      jc(3)       = jc(3) + dpot**2/32 * u * v
      djcdetam(3) = djcdetam(3) + dpot**2/32 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot/32 * (2 * u + dpot * duddpot) * v
      u           = (1 - e**2) * (1 - 4*e**2)
      dudetam     = - (10 - 16*e**2) * e * dedetam
      duddpot     = - (10 - 16*e**2) * e * deddpot
      v           = G(1)**3 * F(1)
      dvdetam     = G(1)**2 * (3 * F(1) * dG(1) + G(1) * dF(1))
      jc(3)       = jc(3) + dpot**2/16 * u * v
      djcdetam(3) = djcdetam(3) + dpot**2/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot/16 * (2 * u + dpot * duddpot) * v
      u           = e * (2 - 3*e**2)
      dudetam     = (2 - 9*e**2) * dedetam
      duddpot     = (2 - 9*e**2) * deddpot
      jc(3)       = jc(3) + dpot/48 * u * F(3)
      djcdetam(3) = djcdetam(3) + dpot/48 * (dudetam * F(3) + u * dF(3))
      djcddpot(3) = djcddpot(3) + (u + dpot * duddpot) * F(3) / 48
      u           = (1 - e**2) * (1.5 - 7*e**2)
      dudetam     = - (17 - 28*e**2) * e * dedetam
      duddpot     = - (17 - 28*e**2) * e * deddpot
      v           = G(1)**2 * F(2)
      dvdetam     = G(1) * (2 * F(2) * dG(1) + G(1) * dF(2))
      jc(3)       = jc(3) - dpot**2/16 * u * v
      djcdetam(3) = djcdetam(3) - dpot**2/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) - dpot/16 * (2 * u + dpot * duddpot) * v
      u           = e**2 * (1 - e**2)
      dudetam     = 2 * e * (1 - 2*e**2) * dedetam
      duddpot     = 2 * e * (1 - 2*e**2) * deddpot
      v           = G(1) * F(3)
      dvdetam     = dG(1) * F(3) + G(1) * dF(3)
      jc(3)       = jc(3) - dpot**2/32 * u * v
      djcdetam(3) = djcdetam(3) - dpot**2/32 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) - dpot/32 * (2 * u + dpot * duddpot) * v

      ! evaluate Taylor series
      j = jc(0) + jc(1) * deta + jc(2) * deta**2 + jc(3) * deta**3
      djdeta(1) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) - jc(1) - 2 * jc(2) * deta - 3 * jc(3) * deta**2
      djdeta(2) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) + jc(1) + 2 * jc(2) * deta + 3 * jc(3) * deta**2
      djddpot   = djcddpot(0) + djcddpot(1) * deta + djcddpot(2) * deta**2 + djcddpot(3) * deta**3
    end subroutine

    subroutine case_2a(j, djdeta, djddpot)
      !! dpot ~ 0 and |deta| >> 0
      real,       intent(out) :: j
        !! output normalized current density on edge
      real,       intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real,       intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer :: k
      real    :: I(5), dI(2,5), jc(0:4), djc(2,0:4), tmp, dtmp(2)

      ! get I_k = integral_eta1^eta2 F(eta)^k deta for k = 0 to 5
      do k = 1, 5
        call int_dist(eta, k, I(k), dI(:,k))
      end do

      ! normalize I_{2..5} wrt I_1 to improve numerical stability (avoid underflows etc.)
      do k = 2, 5
        I(   k) = I(k) / I(1)
        dI(:,k) = (dI(:,k) - I(k) * dI(:,1)) / I(1)
      end do

      ! zero-th Taylor coefficient
      jc(0)    = - I(1)
      djc(:,0) = - dI(:,1)

      ! first Taylor coefficient
      jc(1)    = I(2)
      djc(:,1) = dI(:,2)

      ! second Taylor coefficient
      tmp      = I(2)**2 - I(3)
      dtmp     = 2*I(2)*dI(:,2) - dI(:,3)
      jc(2)    = tmp / I(1)
      djc(:,2) = (dtmp - jc(2) * dI(:,1)) / I(1)

      ! third Taylor coefficient
      tmp      = I(4) + 2*I(2)**3 - 3*I(2)*I(3)
      dtmp     = dI(:,4) + 6*I(2)**2*dI(:,2) - 3*dI(:,2)*I(3) - 3*I(2)*dI(:,3)
      jc(3)    = tmp / I(1)**2
      djc(:,3) = (dtmp - 2 * tmp * dI(:,1) / I(1)) / I(1)**2

      ! fourth Taylor coefficient
      tmp      = 5*I(2)**4 - I(5) + 2*I(3)**2 - 10*I(2)**2*I(3) + 4*I(2)*I(4)
      dtmp     = 20*I(2)**3*dI(:,2) - dI(:,5) + 4*I(3)*dI(:,3) - 20*I(2)*dI(:,2)*I(3) - 10*I(2)**2*dI(:,3) + 4*dI(:,2)*I(4) + 4*I(2)*dI(:,4)
      jc(4)    = tmp / I(1)**3
      djc(:,4) = (dtmp - 3 * tmp * dI(:,1) / I(1)) / I(1)**3

      ! evaluate Taylor series
      j       = jc(0)
      djdeta  = djc(:,0)
      djddpot = 0
      do k = 1, 4
        j       = j + jc(k) * dpot**k
        djdeta  = djdeta + djc(:,k) * dpot**k
        djddpot = djddpot + jc(k) * k * dpot**(k-1)
      end do
    end subroutine

    subroutine case_2b(j, djdeta, djddpot)
      !! dpot ~ 0 and deta ~ 0
      real,       intent(out) :: j
        !! output normalized current density on edge
      real,       intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real,       intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer :: k
      real    :: etam, F(0:2), dF(0:2), tmp, dtmpdetam, dtmpddeta, dtmpddpot, djdetam, djddeta

      ! evaluate distribution and derivatives at mean eta
      etam = 0.5 * (eta(1) + eta(2))
      do k = 0, 2
        call dist(etam, k, F(k), dF(k))
      end do

      ! factor
      tmp       = F(0) - (F(1)/F(0)) * F(1) / 12 * deta * dpot + F(2) / 24 * deta**2
      dtmpdetam = dF(0) - (2 * dF(1) - F(1) * dF(0) / F(0)) * F(1) / (12 * F(0)) * deta * dpot
      dtmpddeta = - F(1)**2 / (12 * F(0)) * dpot + F(2) / 12 * deta
      dtmpddpot = - F(1)**2 / (12 * F(0)) * deta

      ! current
      j = (dpot - deta) * tmp
      djdetam = (dpot - deta) * dtmpdetam
      djddeta = - tmp + (dpot - deta) * dtmpddeta
      djdeta  = djdetam * [0.5, 0.5] + djddeta * [-1.0, 1.0]
      djddpot = tmp + (dpot - deta) * dtmpddpot
    end subroutine

    subroutine residual(t, r, drdt, drdeta, drddpot, j, djdt, djdeta, djddpot)
      real, intent(in)  :: t
      real, intent(out) :: r
      real, intent(out) :: drdt
      real, intent(out) :: drdeta(2)
      real, intent(out) :: drddpot
      real, intent(out) :: j
      real, intent(out) :: djdt
      real, intent(out) :: djdeta(2)
      real, intent(out) :: djddpot

      real :: e, e1, p(2), drdp(2), L, dLdB, dLdt

      ! transform: j = dpot * FL * (1 - sgn * exp(-t))
      e = exp(-t)
      if (sgn > 0) then ! CASE0A, CASE0B
        e1 = - expm1(-t)
      else ! CASE0C
        e1 = e + 1
      end if
      j         = dpot * n(1) * e1
      djdt      = dpot * n(1) * sgn * e
      djdeta(1) = dpot * e1 / detadn(1)
      djdeta(2) = 0
      djddpot   = n(1) * e1

      ! tanh-sinh quadrature, leave out part close to eta(1) (might have singularity)
      p = [j, dpot]
      call quad(integrand, eta(1) + sgn * epseta, eta(2), p, r, drdeta(1), drdeta(2), drdp, max_levels = 8)
      drdt      = drdp(1) * djdt
      drdeta(1) = drdeta(1) * (1 + sgn * depseta) + drdp(1) * djdeta(1)
      drddpot   = drdp(2) + drdp(1) * djddpot

      ! analytical integration for missing part (distribution function approximated by first order Taylor series)
      if (t <= - log(B)) then
        L = log1p(B / e)
      else
        L = t + log(e + B)
      end if
      dLdB      = 1 / (B + e)
      dLdt      = B * dLdB
      r         = r         + A * j * L + C
      drdt      = drdt      + A * (djdt * L + j * dLdt)
      drdeta(1) = drdeta(1) + (dAdeta1 * j + A * djdeta(1)) * L + A * j * dLdB * dBdeta1 + dCdeta1
      drddpot   = drddpot   + (dAddpot * j + A * djddpot) * L + dCddpot

      ! r = integral - 1 == 0
      r = r - 1
    end subroutine

    subroutine integrand(eta, p, u, dudeta, dudp)
      real, intent(in)  :: eta
      real, intent(in)  :: p(:)
      real, intent(out) :: u
      real, intent(out) :: dudeta
      real, intent(out) :: dudp(:)

      real :: q

      call dist(eta, 0, u, dudeta)
      q      = p(2) * u - p(1)
      u      = u / q
      dudeta = (1 - u * p(2)) / q * dudeta
      dudp(1) = u / q
      dudp(2) = - u**2
    end subroutine

  end subroutine

end module
