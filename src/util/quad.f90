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
    subroutine integrand(x, f)
      real, intent(in)  :: x
        !! argument
      real, intent(out) :: f
        !! output integrand
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

contains

  subroutine quad(func, a, b, I, n, eps, e, num)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature with double precision (see qthsh)
    procedure(integrand)           :: func
      !! function to integrate
    real,              intent(in)  :: a
      !! lower integration bound
    real,              intent(in)  :: b
      !! upper integration bound
    real,              intent(out) :: I
      !! output value of integral
    integer, optional, intent(in)  :: n
      !! maximal number of levels (default: 8)
    real,    optional, intent(in)  :: eps
      !! error tolerance (default: 1e-14)
    real,    optional, intent(out) :: e
      !! output error estimate
    integer, optional, intent(out) :: num
      !! output number of function calls

    integer :: k, mode, sign, n_, num_
    real    :: a_, b_, c, d, eh, eps_, fm, fp, h, p, q, r, s, t, u, tol, v, w, x, y

    a_   = a
    b_   = b
    n_   = 8
    num_ = 0
    if (present(n)) n_ = n
    eps_ = 1e-14
    if (present(eps)) eps_ = eps

    tol  = 10 * eps_
    c    = 0
    d    = 1
    sign = 1
    h    = 2
    k    = 0

    ! swap bounds if necessary
    if (b_ < a_) then
      v    = b_
      b_   = a_
      a_   = v
      sign = -1
    end if

    ! determine mode
    if (ieee_is_finite(a_) .and. ieee_is_finite(b_)) then
      mode = 0 ! tanh-sinh
      c    = 0.5 * (a_ + b_)
      d    = 0.5 * (b_ - a_)
      v    = c
    elseif (ieee_is_finite(a_)) then
      mode = 1 ! exp-sinh
      c    = a_
      v    = a_ + d
    elseif (ieee_is_finite(b_)) then
      mode = 1 ! exp-sinh
      d    = -d
      sign = -sign
      c    = b_
      v    = b_ + d
    else
      mode = 2 ! sinh-sinh
      v    = 0
    end if

    ! start sum with first point
    call func(v, s)
    num_ = num_ + 1

    do
      p  = 0
      fp = 0
      fm = 0
      h  = h / 2
      eh = exp(h)
      t  = eh
      if (k > 0) eh = eh * eh

      if (mode == 0) then ! tanh-sinh
        do
          u = exp(1 / t - t)            ! = exp(-2*sinh(j*h)) = 1 / exp(sinh(j*h))**2
          r = 2 * u / (1 + u)           ! = 1 - tanh(sinh(j*h))
          w = (t + 1 / t) * r / (1 + u) ! = cosh(j*h) / cosh(sinh(j*h))**2
          x = d * r

          ! only recalculate integrand if not too close to previous value
          if (a_ + x > a_) then
            call func(a_ + x, y)
            num_ = num_ + 1
            if (ieee_is_finite(y)) fp = y
          end if
          if (b_ - x < b_) then
            call func(b_ - x, y)
            num_ = num_ + 1
            if (ieee_is_finite(y)) fm = y
          end if

          q = w * (fp + fm)
          p = p + q
          t = t * eh
          if (abs(q) <= eps_ * abs(p)) exit
        end do
      else
        t = t / 2
        do
          r = exp(t - 0.25 / t) ! = exp(sinh(j*h))
          w = r
          q = 0
          if (mode == 1) then
            ! exp-sinh
            x = c + d / r
            if (x == c) exit
            call func(x, y)
            num_ = num_ + 1
            if (ieee_is_finite(y)) q = q + y / w
          else
            ! sinh-sinh
            r = 0.5 * (r - 1 / r) ! = sinh(sinh(j*h))
            w = 0.5 * (w + 1 / w) ! = cosh(sinh(j*h))
            x = c - d * r
            call func(x, y)
            if (ieee_is_finite(y)) q = q + y * w
          end if
          x = c + d * r
          call func(x, y)
          num_ = num_ + 1
          if (ieee_is_finite(y)) q = q + y * w
          q = q * (t + 0.25 / t) ! q = q * cosh(j*h)
          p = p + q
          t = t * eh
          if (abs(q) <= eps_ * abs(p)) exit
        end do
      end if
      v = s - p
      s = s + p
      k = k + 1
      if ((abs(v) <= tol * abs(s)) .or. (k > n_)) exit
    end do
    I = sign * d * s * h
    if (present(e)) e = abs(v) / (abs(s) + eps_)
    if (present(num)) num = num_
  end subroutine

  m4_divert(m4_ifdef({m4_mpfr},0,-1))

  subroutine quad_mpfr(func, a, b, I, n, eps, e, num)
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
    type(mpfr), optional, intent(inout) :: e
      !! output error (must be already initialized)
    integer,    optional, intent(in)    :: n
      !! maximal number of levels (default: 16)
    real,       optional, intent(in)    :: eps
      !! error tolerance (default: 1e-32)
    integer,    optional, intent(out)   :: num
      !! output number of function calls

    integer    :: k, mode, sign, n_, num_
    real       :: a_, b_, tmp, tol, eps_
    type(mpfr) :: c, d, eh, fm, fp, h, p, q, r, s, t, t1, t2, u, v, w, x, y

    a_   = a
    b_   = b
    n_   = 16
    num_ = 0
    if (present(n)) n_ = n
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

    tol  = 10 * eps_
    call c%set(0)
    call d%set(1)
    sign = 1
    call h%set(2)
    k    = 0

    ! swap bounds if necessary
    if (b_ < a_) then
      tmp  = b_
      b_   = a_
      a_   = tmp
      sign = -1
    end if

    ! determine mode
    if (ieee_is_finite(a_) .and. ieee_is_finite(b_)) then
      mode = 0 ! tanh-sinh
      call c%set(a_)       ! c = a_
      call add(c, c, b_)   ! c = a_ + b_
      call mul(c, c, 0.5)  ! c = 0.5 * (a_ + b_)
      call d%set(b_)       ! d = b_
      call sub(d, d, a_)   ! d = b_ - a_
      call mul(d, d, 0.5)  ! d = 0.5 * (b_ - a_)
      call v%set(c)        ! v = c
    elseif (ieee_is_finite(a_)) then
      mode = 1 ! exp-sinh
      call c%set(a_)       ! c = a_
      call add(v, c, d)    ! v = a_ + d
    elseif (ieee_is_finite(b_)) then
      mode = 1 ! exp-sinh
      call neg(d, d)       ! d = - d
      sign = -sign
      call c%set(b_)       ! c = b_
      call add(v, c, d)    ! v = b_ + d
    else
      mode = 2 ! sinh-sinh
      call v%set(0)        ! v = 0
    end if

    ! start sum with first point
    call func(v, s)
    num_ = num_ + 1

    do
      call p%set(0)               ! p  = 0
      call fp%set(0)              ! fp = 0
      call fm%set(0)              ! fm = 0
      call div(h, h, 2)           ! h  = h / 2
      call exp(eh, h)             ! eh = exp(h)
      call t%set(eh)              ! t  = eh
      if (k > 0) call sqr(eh, eh) ! eh = eh * eh

      if (mode == 0) then ! tanh-sinh
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
            num_ = num_ + 1
          end if
          call sub(t1, b_, x) ! t1 = b_ - x
          if (t1 < b_) then
            call func(t1, fm)
            num_ = num_ + 1
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
          if (mode == 1) then
            ! exp-sinh
            call div(x, d, r)  ! x = d / r
            call add(x, c, x)  ! x = c + d / r
            if (x == c) exit
            call func(x, y)
            num_ = num_ + 1
            call div(y, y, w)  ! y = y / w
            call add(q, q, y)  ! q = q + y / w
          else
            ! sinh-sinh
            call div(t1, 1, r)   ! t1 = 1 / r
            call sub(t1, r, t1)  ! t1 = r - 1 / r
            call mul(r, t1, 0.5) ! r = 0.5 * (r - 1 / r) = sinh(sinh(j*h))
            call div(t1, 1, w)   ! t1 = 1 / w
            call add(t1, w, t1)  ! t1 = w + 1 / w
            call mul(w, t1, 0.5) ! w = 0.5 * (w + 1 / w) = cosh(sinh(j*h))
            call fms(x, d, r, c) ! x = d * r - c
            call neg(x, x)       ! x = c - d * r
            call func(x, y)
            num_ = num_ + 1
            call fma(q, y, w, q) ! q = q + y * w
          end if
          call fma(x, d, r, c) ! x = c + d * r
          call func(x, y)
          num_ = num_ + 1
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

      k = k + 1
      if (k > n_) exit

      call abs(t1, v)       ! t1 = abs(v)
      call abs(t2, s)       ! t2 = abs(s)
      call mul(t2, t2, tol) ! t2 = tol * abs(s)
      if (t1 <= t2) exit
    end do

    ! integral
    call mul(I, d, sign)     ! I = sign * d
    call mul(I, I, s)        ! I = sign * d * s
    call mul(I, I, h)        ! I = sign * d * s * h

    ! error estimate
    if (present(e)) then
      call abs(t1, v)        ! t1 = abs(v)
      call abs(t2, s)        ! t2 = abs(s)
      call add(t2, t2, eps_) ! t2 = abs(s) + eps_
      call div(e, t1, t2)    ! e = abs(v) / (abs(s) + eps_)
    end if

    if (present(num)) num = num_

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
