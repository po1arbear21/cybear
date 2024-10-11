m4_include(macro.f90.inc)

module quad_m
  !! general purpose integration module (tanh-sinh quadrature)
  !! double precision and multiprecision (requires MPFR library)

  m4_ifdef({m4_mpfr},{
  use mpfr_m, only: abs, add, div, exp, fma, fms, mpfr, mul, neg, sqr, sub, operator(<), operator(<=), operator(>), operator(==)
  })
  use ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public quad
  m4_ifdef({m4_mpfr},{
  public quad_mpfr
  })

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
    m4_ifdef({m4_mpfr},{
    subroutine integrand_mpfr(x, f)
      import mpfr
      type(mpfr), intent(in)    :: x
        !! argument
      type(mpfr), intent(inout) :: f
        !! output integrand (already initialized)
    end subroutine
    })
  end interface

  integer, parameter :: MODE_TANH_SINH = 0
  integer, parameter :: MODE_EXP_SINH  = 1
  integer, parameter :: MODE_SINH_SINH = 2

contains

  subroutine quad(func, a, b, p, I, dIda, dIdb, dIdp, eps, err, nlevels, ncalls)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature with double precision (see qthsh)
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
    real,    optional, intent(in)  :: eps
      !! error tolerance (default: 1e-14)
    real,    optional, intent(out) :: err
      !! output error estimate
    integer, optional, intent(in)  :: nlevels
      !! maximal number of levels (default: 8)
    integer, optional, intent(out) :: ncalls
      !! output number of function calls

    integer :: lvl, mode, sgn, nlevels_, ncalls_
    real    :: bnd(2), c, d, eh, eps_, f, fm, fp, h, pp, qq, r, s, t, u, tol, v, w, x
    real    :: dfdp(size(p)), dfpdp(size(p)), dfmdp(size(p)), dppdp(size(p)), dqqdp(size(p)), dsdp(size(p))
    real    :: dxdbnd(2), dfdx, dfpdbnd(2), dfmdbnd(2), dppdbnd(2), dqqdbnd(2), dsdbnd(2), dIdbnd(2)

    bnd      = [a, b]
    nlevels_ = 8
    ncalls_  = 0
    if (present(nlevels)) nlevels_ = nlevels
    eps_ = 1e-14
    if (present(eps)) eps_ = eps

    tol = 10 * eps_
    sgn = 1
    h   = 2
    lvl = 0

    ! swap bounds if necessary
    if (bnd(2) < bnd(1)) then
      v      = bnd(2)
      bnd(2) = bnd(1)
      bnd(1) = v
      sgn    = -1
    end if

    ! determine mode
    if (ieee_is_finite(bnd(1)) .and. ieee_is_finite(bnd(2))) then
      mode   = MODE_TANH_SINH
      c      = 0.5 * (bnd(1) + bnd(2))
      d      = 0.5 * (bnd(2) - bnd(1))
      x      = c
      dxdbnd = [0.5, 0.5]
    elseif (ieee_is_finite(bnd(1))) then
      mode   = MODE_EXP_SINH
      c      = bnd(1)
      d      = 1
      x      = c + d
      dxdbnd = [1.0, 0.0]
    elseif (ieee_is_finite(bnd(2))) then
      mode   = MODE_EXP_SINH
      c      = bnd(2)
      d      = -1
      sgn    = -sgn
      x      = c + d
      dxdbnd = [0.0, 1.0]
    else
      mode   = MODE_SINH_SINH
      c      = 0
      d      = 1
      x      = 0
      dxdbnd = 0
    end if

    ! start sum with first point
    call func(x, p, f, dfdx, dfdp)
    ncalls_ = ncalls_ + 1
    s       = f
    dsdp    = dfdp
    dsdbnd  = dfdx * dxdbnd

    do
      pp      = 0
      dppdp   = 0
      dppdbnd = 0
      fp      = 0
      dfpdp   = 0
      dfpdbnd = 0
      fm      = 0
      dfmdp   = 0
      dfmdbnd = 0
      h       = h / 2
      eh      = exp(h)
      t       = eh
      if (lvl > 0) eh = eh * eh

      if (mode == MODE_TANH_SINH) then ! tanh-sinh
        do
          u      = exp(1 / t - t)            ! = exp(-2*sinh(j*h)) = 1 / exp(sinh(j*h))**2
          r      = 2 * u / (1 + u)           ! = 1 - tanh(sinh(j*h))
          w      = (t + 1 / t) * r / (1 + u) ! = cosh(j*h) / cosh(sinh(j*h))**2
          x      = d * r
          dxdbnd = r * [-0.5, 0.5]

          ! only recalculate integrand if not too close to previous value
          if (bnd(1) + x > bnd(1)) then
            call func(bnd(1) + x, p, f, dfdx, dfdp)
            ncalls_ = ncalls_ + 1
            if (ieee_is_finite(f)) then
              fp      = f
              dfpdp   = dfdp
              dfpdbnd = dfdx * ([1.0, 0.0] + dxdbnd)
            end if
          end if
          if (bnd(2) - x < bnd(2)) then
            call func(bnd(2) - x, p, f, dfdx, dfdp)
            ncalls_ = ncalls_ + 1
            if (ieee_is_finite(f)) then
              fm      = f
              dfmdp   = dfdp
              dfmdbnd = dfdx * ([0.0, 1.0] - dxdbnd)
            end if
          end if

          qq      = w * (fp + fm)
          dqqdp   = w * (dfpdp + dfmdp)
          dqqdbnd = w * (dfpdbnd + dfmdbnd)
          pp      = pp + qq
          dppdp   = dppdp + dqqdp
          dppdbnd = dppdbnd + dqqdbnd
          t       = t * eh
          if (abs(qq) <= eps_ * abs(pp)) exit
        end do
      else
        t = t / 2
        do
          r       = exp(t - 0.25 / t) ! = exp(sinh(j*h))
          w       = r
          qq      = 0
          dqqdp   = 0
          dqqdbnd = 0
          if (mode == MODE_EXP_SINH) then
            x = c + d / r
            if (x == c) exit
            call func(x, p, f, dfdx, dfdp)
            ncalls_ = ncalls_ + 1
            if (ieee_is_finite(f)) then
              qq      = qq + f / w
              dqqdp   = dqqdp + dfdp / w
              dqqdbnd = dqqdbnd + dfdx * dxdbnd / w
            end if
          else ! MODE_SINH_SINH
            r = 0.5 * (r - 1 / r) ! = sinh(sinh(j*h))
            w = 0.5 * (w + 1 / w) ! = cosh(sinh(j*h))
            x = c - d * r
            call func(x, p, f, dfdx, dfdp)
            if (ieee_is_finite(f)) then
              qq    = qq + f * w
              dqqdp = dqqdp + dfdp * w
            end if
          end if
          x = c + d * r
          call func(x, p, f, dfdx, dfdp)
          ncalls_ = ncalls_ + 1
          if (ieee_is_finite(f)) then
            qq      = qq + f * w
            dqqdp   = dqqdp + dfdp * w
            dqqdbnd = dqqdbnd + dfdx * dxdbnd * w
          end if
          qq      = qq * (t + 0.25 / t) ! qq = qq * cosh(j*h)
          dqqdp   = dqqdp * (t + 0.25 / t)
          dqqdbnd = dqqdbnd * (t + 0.25 / t)
          pp      = pp + qq
          dppdp   = dppdp + dqqdp
          dppdbnd = dppdbnd + dqqdbnd
          t       = t * eh
          if (abs(qq) <= eps_ * abs(pp)) exit
        end do
      end if
      v      = s - pp
      s      = s + pp
      dsdp   = dsdp + dppdp
      dsdbnd = dsdbnd + dppdbnd
      lvl = lvl + 1
      if ((abs(v) <= tol * abs(s)) .or. (lvl > nlevels_)) exit
    end do
    I      = sgn * d * s * h
    dIdp   = sgn * d * dsdp * h
    dIdbnd = sgn * d * dsdbnd * h
    if (mode == MODE_TANH_SINH) then
      dIdbnd = dIdbnd + sgn * [-0.5, 0.5] * s * h
    end if
    if (sgn == 1) then
      dIda = dIdbnd(1)
      dIdb = dIdbnd(2)
    else
      dIda = dIdbnd(2)
      dIdb = dIdbnd(1)
    end if
    if (present(err)) err = abs(v) / (abs(s) + eps_)
    if (present(ncalls)) ncalls = ncalls_
  end subroutine

  m4_divert(m4_ifdef({m4_mpfr},0,-1))
  subroutine quad_mpfr(func, a, b, I, eps, err, nlevels, ncalls)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature with multiprecision
    procedure(integrand_mpfr)           :: func
      !! integrand
    real,                 intent(in)    :: a
      !! lower integration bound
    real,                 intent(in)    :: b
      !! upper integration bound
    type(mpfr),           intent(inout) :: I
      !! output integral value (must be already initialized)
    real,       optional, intent(in)    :: eps
      !! error tolerance (default: 1e-32)
    type(mpfr), optional, intent(inout) :: err
      !! output error (must be already initialized)
    integer,    optional, intent(in)    :: nlevels
      !! maximal number of levels (default: 16)
    integer,    optional, intent(out)   :: ncalls
      !! output number of function calls

    integer    :: lvl, mode, sgn, nlevels_, ncalls_
    real       :: a_, b_, tmp, tol, eps_
    type(mpfr) :: c, d, eh, fm, fp, h, p, q, r, s, t, t1, t2, u, v, w, x, y

    a_   = a
    b_   = b
    nlevels_   = 16
    ncalls_ = 0
    if (present(nlevels)) nlevels_ = nlevels
    eps_ = 1e-32
    if (present(eps)) eps_ = eps

    ! initialize memory
    call c%init()
    call d%init()
    call eh%init()
    call fm%init()
    call fp%init()
    call h%init()
    call p%init()
    call q%init()
    call r%init()
    call s%init()
    call t%init()
    call t1%init()
    call t2%init()
    call u%init()
    call v%init()
    call w%init()
    call x%init()
    call y%init()

    tol = 10 * eps_
    call c%set(0)
    call d%set(1)
    sgn = 1
    call h%set(2)
    lvl = 0

    ! swap bounds if necessary
    if (b_ < a_) then
      tmp = b_
      b_  = a_
      a_  = tmp
      sgn = -1
    end if

    ! determine mode
    if (ieee_is_finite(a_) .and. ieee_is_finite(b_)) then
      mode = MODE_TANH_SINH
      call c%set(a_)       ! c = a_
      call add(c, c, b_)   ! c = a_ + b_
      call mul(c, c, 0.5)  ! c = 0.5 * (a_ + b_)
      call d%set(b_)       ! d = b_
      call sub(d, d, a_)   ! d = b_ - a_
      call mul(d, d, 0.5)  ! d = 0.5 * (b_ - a_)
      call v%set(c)        ! v = c
    elseif (ieee_is_finite(a_)) then
      mode = MODE_EXP_SINH
      call c%set(a_)       ! c = a_
      call add(v, c, d)    ! v = a_ + d
    elseif (ieee_is_finite(b_)) then
      mode = MODE_EXP_SINH
      call neg(d, d)       ! d = - d
      sgn  = -sgn
      call c%set(b_)       ! c = b_
      call add(v, c, d)    ! v = b_ + d
    else
      mode = MODE_SINH_SINH
      call v%set(0)        ! v = 0
    end if

    ! start sum with first point
    call func(v, s)
    ncalls_ = ncalls_ + 1

    do
      call p%set(0)     ! p  = 0
      call fp%set(0)    ! fp = 0
      call fm%set(0)    ! fm = 0
      call div(h, h, 2) ! h  = h / 2
      call exp(eh, h)   ! eh = exp(h)
      call t%set(eh)    ! t  = eh
      if (lvl > 0) call sqr(eh, eh) ! eh = eh * eh

      if (mode == MODE_TANH_SINH) then
        do
          call div(t1, 1, t)  ! t1 = 1 / t
          call sub(u, t1, t)  ! u = 1 / t - t
          call exp(u, u)      ! u = exp(1 / t - t) = exp(-2*sinh(j*h)) = 1 / exp(sinh(j*h))**2

          call add(t2, u, 1)  ! t2 = 1 + u
          call div(r, u, t2)  ! r = u / (1 + u)
          call add(r, r, r)   ! r = 2 * u / (1 + u) = 1 - tanh(sinh(j*h))

          call div(w, r, t2)  ! w = r / (1 + u)
          call add(t1, t, t1) ! t1 = t + 1 / t
          call mul(w, t1, w)  ! w = (t + 1 / t) * r / (1 + u) = cosh(j*h) / cosh(sinh(j*h))**2

          call mul(x, d, r)   ! x = d * r

          ! only recalculate integrand if not too close to previous value
          call add(t1, x, a_) ! t1 = a_ + x
          if (t1 > a_) then
            call func(t1, fp)
            ncalls_ = ncalls_ + 1
          end if
          call sub(t1, b_, x) ! t1 = b_ - x
          if (t1 < b_) then
            call func(t1, fm)
            ncalls_ = ncalls_ + 1
          end if

          call add(q, fp, fm)    ! q = fp + fm
          call mul(q, w, q)      ! q = w * (fp + fm)
          call add(p, p, q)      ! p = p + q
          call mul(t, t, eh)     ! t = t * eh

          call abs(t1, q)        ! t1 = abs(q)
          call abs(t2, p)        ! t2 = abs(p)
          call mul(t2, t2, eps_) ! t2 = eps_ * abs(p)
          if (t1 <= t2) exit
        end do
      else
        call div(t, t, 2) ! t = t / 2
        do
          call div(r, 0.25, t) ! r = 0.25 / t
          call sub(r, t, r)    ! r = t - 0.25 / t
          call exp(r, r)       ! r = exp(t - 0.25 / t) = exp(sinh(j*h))
          call w%set(r)        ! w = r
          call q%set(0)        ! q = 0
          if (mode == MODE_EXP_SINH) then
            call div(x, d, r)  ! x = d / r
            call add(x, c, x)  ! x = c + d / r
            if (x == c) exit
            call func(x, y)
            ncalls_ = ncalls_ + 1
            call div(y, y, w)  ! y = y / w
            call add(q, q, y)  ! q = q + y / w
          else ! MODE_SINH_SINH
            call div(t1, 1, r)   ! t1 = 1 / r
            call sub(t1, r, t1)  ! t1 = r - 1 / r
            call mul(r, t1, 0.5) ! r = 0.5 * (r - 1 / r) = sinh(sinh(j*h))
            call div(t1, 1, w)   ! t1 = 1 / w
            call add(t1, w, t1)  ! t1 = w + 1 / w
            call mul(w, t1, 0.5) ! w = 0.5 * (w + 1 / w) = cosh(sinh(j*h))
            call fms(x, d, r, c) ! x = d * r - c
            call neg(x, x)       ! x = c - d * r
            call func(x, y)
            ncalls_ = ncalls_ + 1
            call fma(q, y, w, q) ! q = q + y * w
          end if
          call fma(x, d, r, c) ! x = c + d * r
          call func(x, y)
          ncalls_ = ncalls_ + 1
          call fma(q, y, w, q) ! q = q + y * w

          call div(t1, 0.25, t)  ! t1 = 0.25 / t
          call add(t1, t, t1)    ! t1 = t + 0.25 / t
          call mul(q, q, t1)     ! q  = q * cosh(j*h)

          call add(p, p, q)      ! p = p + q
          call mul(t, t, eh)     ! t = t * eh

          call abs(t1, q)        ! t1 = abs(q)
          call abs(t2, p)        ! t2 = abs(p)
          call mul(t2, t2, eps_) ! t2 = eps_ * abs(p)
          if (t1 <= t2) exit
        end do
      end if
      call sub(v, s, p) ! v = s - p
      call add(s, s, p) ! s = s + p

      lvl = lvl + 1
      if (lvl > nlevels_) exit

      call abs(t1, v)       ! t1 = abs(v)
      call abs(t2, s)       ! t2 = abs(s)
      call mul(t2, t2, tol) ! t2 = tol * abs(s)
      if (t1 <= t2) exit
    end do

    ! integral
    call mul(I, d, sgn)      ! I = sgn * d
    call mul(I, I, s)        ! I = sgn * d * s
    call mul(I, I, h)        ! I = sgn * d * s * h

    ! error estimate
    if (present(err)) then
      call abs(t1, v)        ! t1 = abs(v)
      call abs(t2, s)        ! t2 = abs(s)
      call add(t2, t2, eps_) ! t2 = abs(s) + eps_
      call div(err, t1, t2)  ! err = abs(v) / (abs(s) + eps_)
    end if

    if (present(ncalls)) ncalls = ncalls_

    ! cleanup
    call c%destruct()
    call d%destruct()
    call eh%destruct()
    call fm%destruct()
    call fp%destruct()
    call h%destruct()
    call p%destruct()
    call q%destruct()
    call r%destruct()
    call s%destruct()
    call t%destruct()
    call t1%destruct()
    call t2%destruct()
    call u%destruct()
    call v%destruct()
    call w%destruct()
    call x%destruct()
    call y%destruct()
  end subroutine

  m4_divert(0)

end module
