module quad_m
  !! general purpose integration module (tanh-sinh quadrature)

  use ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public quad

  interface
    subroutine integrand(x, p, f, dfdx, dfdp)
      real, intent(in)  :: x
        !! argument
      real, intent(in)  :: p(:)
        !! parameters
      real, intent(out) :: f
        !! output integrand
      real, intent(out) :: dfdx
        !! output derivative of f wrt x
      real, intent(out) :: dfdp(:)
        !! output derivatives of f wrt p
    end subroutine
  end interface

  integer, parameter :: MODE_TANH_SINH = 0
  integer, parameter :: MODE_EXP_SINH  = 1
  integer, parameter :: MODE_SINH_SINH = 2

contains

  subroutine quad(func, a, b, p, I, dIda, dIdb, dIdp, rtol, err, min_levels, max_levels, ncalls)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature (uses quad-precision internally)
    procedure(integrand)           :: func
      !! function to integrate
    real,              intent(in)  :: a
      !! lower integration bound
    real,              intent(in)  :: b
      !! upper integration bound
    real,              intent(in)  :: p(:)
      !! additional parameters passed to the integrand
    real,              intent(out) :: I
      !! output value of integral
    real,              intent(out) :: dIda
      !! output derivative of I wrt lower bound
    real,              intent(out) :: dIdb
      !! output derivative of I wrt upper bound
    real,              intent(out) :: dIdp(:)
      !! output derivatives of I wrt parameters
    real,    optional, intent(in)  :: rtol
      !! error tolerance, overestimate (default: 1e-15)
    real,    optional, intent(out) :: err
      !! output error estimate
    integer, optional, intent(in)  :: min_levels
      !! minimum number of levels (default: 2)
    integer, optional, intent(in)  :: max_levels
      !! maximal number of levels (default: 16)
    integer, optional, intent(out) :: ncalls
      !! output number of function calls

    integer       :: min_levels_, max_levels_, ncalls_, sgn, level, const_counter, mode
    real          :: bnd(2), rtol_, dIdbnd(2), scale, h, dpartdbnd(2), dpartdp(size(p))
    real          :: dsumdbnd(2), dsumdp(size(p)), sum, sum_err, x1, x1_old, x2, x2_old, dx1dbnd(2), dx2dbnd(2), xref
    real          :: f, dfdx, dfdp(size(p)), df1dbnd(2), df2dbnd(2), df1dp(size(p)), df2dp(size(p))
    real(kind=16) :: sum_16, part_16, part_16_old, exp_h_16, exp_nh_16, w1_16, w2_16, f1_16, f2_16, r_16, s_16

    ! optional arguments
    min_levels_ = 2
    max_levels_ = 16
    ncalls_     = 0
    rtol_       = 1e-15
    if (present(min_levels)) min_levels_ = min_levels
    if (present(max_levels)) max_levels_ = max_levels
    if (present(rtol)) rtol_ = rtol

    ! sort integration bounds
    if (a < b) then
      bnd = [a, b]
      sgn = 1
    else
      bnd = [b, a]
      sgn = -1
    end if

    ! determine integration mode and first node
    if (ieee_is_finite(bnd(1)) .and. ieee_is_finite(bnd(2))) then
      mode    = MODE_TANH_SINH
      scale   = 0.5 * (bnd(2) - bnd(1))
      x1      = 0.5 * (bnd(1) + bnd(2)) ! x1 = (a + b) / 2 + (b - a) / 2 * tanh(sinh(0))
      dx1dbnd = [0.5, 0.5]
    elseif (ieee_is_finite(bnd(1))) then
      mode    = MODE_EXP_SINH
      scale   = 1
      xref    = bnd(1)
      x1      = xref + scale ! x1 = a + exp(sinh(0))
      dx1dbnd = [1.0, 0.0]
    elseif (ieee_is_finite(bnd(2))) then
      mode    = MODE_EXP_SINH
      scale   = -1
      xref    = bnd(2)
      x1      = xref + scale ! x1 = b - exp(sinh(0))
      dx1dbnd = [0.0, 1.0]
    else
      mode    = MODE_SINH_SINH
      scale   = 1
      x1      = 0 ! x1 = sinh(sinh(0))
      dx1dbnd = 0
    end if
    dx2dbnd = dx1dbnd

    ! start sum with first point, w(t=0) = 1 regardless of integration mode
    call func(x1, p, f, dfdx, dfdp)
    ncalls_  = ncalls_ + 1
    sum_16   = f
    dsumdbnd = dfdx * dx1dbnd
    dsumdp   = dfdp
    sum_err  = huge(1.0)

    ! loop over levels
    level   = 0
    h       = 2 ! twice initial step size
    do while ((level < min_levels_) .or. ((sum_err > rtol_ * abs(sum_16)) .and. (level < max_levels_)))
      h         = 0.5 * h
      exp_h_16  = exp(real(h, kind = 16)) ! exp(h)
      exp_nh_16 = exp_h_16                ! exp(n*h)

      ! jump over every second point for higher levels (these are already included in the previous level)
      if (level > 0) exp_h_16 = exp_h_16 * exp_h_16 ! exp(2*h)

      ! reset temporary values
      part_16       = 0
      dpartdbnd     = 0
      dpartdp       = 0
      const_counter = 0
      x1            = 0
      x2            = 0
      f1_16         = 0
      f2_16         = 0
      df1dbnd       = 0
      df2dbnd       = 0
      df1dp         = 0
      df2dp         = 0

      ! loop over n, until part_16 does not change anymore (check twice for safety)
      do while (const_counter < 2)
        x1_old = x1
        x2_old = x2

        ! get nodes and weights depending on integration mode
        select case (mode)
        case (MODE_TANH_SINH)
          r_16    = exp(1 / exp_nh_16 - exp_nh_16)                  ! r  = exp(-2 * sinh(n * h))
          s_16    = 2 * r_16 / (1 + r_16)                           ! s  = 1 - tanh(sinh(n * h))
          w1_16   = (exp_nh_16 + 1 / exp_nh_16) * s_16 / (1 + r_16) ! w1 = cosh(n * h) / cosh(sinh(n * h))**2
          w2_16   = w1_16                                           ! w2 = cosh(n * h) / cosh(sinh(n * h))**2
          x1      = real(bnd(1) + scale * s_16)                     ! x1 = (a+b)/2 + (b-a)/2 * tanh(sinh(-n*h))
          dx1dbnd = [1.0, 0.0] + [-0.5, 0.5] * real(s_16)
          x2      = real(bnd(2) - scale * s_16)                     ! x2 = (a+b)/2 + (b-a)/2 * tanh(sinh( n*h))
          dx2dbnd = [0.0, 1.0] + [0.5, -0.5] * real(s_16)
        case (MODE_EXP_SINH)
          r_16    = exp(0.5*(exp_nh_16 - 1 / exp_nh_16))     ! r  = exp(sinh(n * h))
          w1_16   = 0.5 * (exp_nh_16 + 1 / exp_nh_16) / r_16 ! w1 = cosh(n * h) * exp(-sinh(n * h))
          w2_16   = 0.5 * (exp_nh_16 + 1 / exp_nh_16) * r_16 ! w2 = cosh(n * h) * exp( sinh(n * h))
          x1      = real(xref + scale / r_16)                ! x1 = a[b] +[-] exp(-sinh(n * h))
          x2      = real(xref + scale * r_16)                ! x2 = a[b] +[-] exp( sinh(n * h))
        case (MODE_SINH_SINH)
          r_16    = exp(0.5*(exp_nh_16 - 1 / exp_nh_16))     ! r  = exp(sinh(n * h))
          s_16    = 0.5 * (r_16 - 1 / r_16)                  ! s  = sinh(sinh(n * h))
          r_16    = 0.5 * (r_16 + 1 / r_16)                  ! r  = cosh(sinh(n * h))
          w1_16   = 0.5 * (exp_nh_16 + 1 / exp_nh_16) * r_16 ! w1 = cosh(n * h) * cosh(sinh(n * h))
          w2_16   = w1_16                                    ! w2 = cosh(n * h) * cosh(sinh(n * h))
          x1      = - real(s_16)                             ! x1 = - sinh(sinh(n * h))
          x2      = - x1                                     ! x2 =   sinh(sinh(n * h))
        end select

        ! evaluate function at first point (reuse value if too close to prev. x)
        if (x1 /= x1_old) then
          call func(x1, p, f, dfdx, dfdp)
          ncalls_ = ncalls_ + 1
          if (ieee_is_finite(f)) then
            f1_16   = f
            df1dbnd = dfdx * dx1dbnd
            df1dp   = dfdp
          end if
        end if

        ! evaluate function at second point (reuse value if too close to prev. x)
        if (x2 /= x2_old) then
          call func(x2, p, f, dfdx, dfdp)
          ncalls_ = ncalls_ + 1
          if (ieee_is_finite(f)) then
            f2_16   = f
            df2dbnd = dfdx * dx2dbnd
            df2dp   = dfdp
          end if
        end if

        ! update partial sum
        part_16_old = part_16
        part_16     = part_16 + w1_16 * f1_16 + w2_16 * f2_16
        if (part_16 == part_16_old) const_counter = const_counter + 1
        dpartdbnd = dpartdbnd + real(w1_16) * df1dbnd + real(w2_16) * df2dbnd
        dpartdp   = dpartdp   + real(w1_16) * df1dp   + real(w2_16) * df2dp

        ! exp(n * h), n = 1, 2, 3, ... for level = 0 and n = 1, 3, 5, ... for level > 0
        exp_nh_16 = exp_nh_16 * exp_h_16
      end do

      ! update sum and estimate error (error is based on difference to prev. level => overestimate)
      sum_err  = abs(real(sum_16 - part_16)) ! corresponds to |I_{level} - I_{level-1}| >= rtol_ |I_{level}|
      sum_16   = sum_16   + part_16
      dsumdbnd = dsumdbnd + dpartdbnd
      dsumdp   = dsumdp   + dpartdp

      ! go to next level
      level = level + 1
    end do

    ! calculate integral from raw sum
    sum    = real(sum_16)
    I      = sgn * scale * sum * h
    dIdbnd = sgn * scale * dsumdbnd * h
    if (mode == MODE_TANH_SINH) then
      dIdbnd = dIdbnd + sgn * [-0.5, 0.5] * sum * h
    end if
    dIdp   = sgn * scale * dsumdp * h
    if (sgn == 1) then
      dIda = dIdbnd(1)
      dIdb = dIdbnd(2)
    else
      dIda = dIdbnd(2)
      dIdb = dIdbnd(1)
    end if

    ! optional output
    if (present(err)) err = sum_err / (abs(sum) + 1e-16)
    if (present(ncalls)) ncalls = ncalls_
  end subroutine

end module
