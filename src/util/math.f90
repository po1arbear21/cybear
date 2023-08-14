m4_include(macro.f90.inc)

module math_m

  use bin_search_m,    only: bin_search, BS_LESS
  use error_m,         only: assert_failed
  use ieee_arithmetic
  use lapack95
  use qsort_m,         only: qsort

  implicit none

  private
  public PI
  public check_lin_dep
  public cross_product, cross_product_2d
  public heaviside
  public interp1
  public isinf
  public ber, dberdx
  public eye_int, eye_real
  public expm1, log1p
  public linspace, logspace
  public norm_inf
  public phi1, dphi1dx, phi2, dphi2dx
  public roots
  public polyg_area_2d

  real, parameter :: PI = 3.141592653589793238462643

  interface expm1
    module procedure :: m_expm1
  end interface
  interface log1p
    module procedure :: m_log1p
  end interface

contains

  pure function cross_product(a, b) result(axb)
    !! calculates the cross product between a and b
    real, intent(in) :: a(3), b(3)
      !! 3-dimensional vectors
    real             :: axb(3)

    axb = [a(2)*b(3) - a(3)*b(2), &
      &    a(3)*b(1) - a(1)*b(3), &
      &    a(1)*b(2) - a(2)*b(1)  ]
  end function

  pure function cross_product_2d(a, b) result(axb)
    !! calculates the 2D cross product (x,y) -> z
    real, intent(in) :: a(2), b(2)
      !! 2-dimensional vectors
    real             :: axb

    axb = a(1)*b(2) - a(2)*b(1)
  end function

  elemental function heaviside(x) result(h)
    !! heaviside step function
    real, intent(in) :: x
    real             :: h

    if      (x < 0) then
      h = 0
    else if (x > 0) then
      h = 1
    else
      h = 0.5
    end if
  end function

  elemental function isinf(x) result (r)
    !! return true if argument is positive of negative infinity
    real, intent(in) :: x
      !! real argument
    logical          :: r

    r = (ieee_class(x) == IEEE_POSITIVE_INF) .or. (ieee_class(x) == IEEE_NEGATIVE_INF)
  end function

  elemental function ber(x) result(b)
    !! bernoulli function ber(x) = x / (exp(x) - 1)
    real, intent(in) :: x
    real             :: b

    if (abs(x) > 1e-6) then
      b = 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x)
    else
      b = 1.0 + x * (-0.5 + x / 12.0)
    end if
  end function

  elemental function dberdx(x) result(dbdx)
    !! derivative of bernoulli function (max rel. error ~ 3e-15 at x = 0.075; FIXME: improve)
    real, intent(in) :: x
    real             :: dbdx

    real :: x_

    x_ = x
    if (x < 0) x_ = -x

    if (x_ > 1e2) then
      dbdx = (1 - x_) * exp(-x_)
    elseif (x_ > 0.075) then
      dbdx = (2.0 * exp(-0.5 * x_) * sinh(0.5 * x_) - x_) / (4.0 * sinh(0.5 * x_)**2)
    elseif (x_ > 0.025) then
      dbdx = 0.5 * (tanh(x_*(1.0/3 + x_**2*(1.0/810 - x_**2*(1.0/68040 + x_**2/6123600)))) - 1)
    elseif (x_ > 1e-3) then
      dbdx = 0.5 * (tanh(x_*(1.0/3 + x_**2*(1.0/810 - x_**2/68040))) - 1)
    elseif (x_ > 1e-5) then
      dbdx = 0.5 * (tanh(x_*(1.0/3 + x_**2/810)) - 1)
    else
      dbdx = 0.5 * (tanh(x_/3) - 1)
    end if

    if (x < 0) dbdx = - dbdx - 1
  end function

  elemental function phi1(x) result(phi)
    !! phi1(x) = (exp(x) - 1) / x = 1 / ber(x)
    real, intent(in) :: x
    real             :: phi

    phi = 1.0 / ber(x)
  end function

  elemental function dphi1dx(x) result(dphidx)
    !! derivative of phi1(x)
    real, intent(in) :: x
    real             :: dphidx

    dphidx = - dberdx(x) / (ber(x)**2)
  end function

  elemental function phi2(x) result(phi)
    !! phi2(x) = (exp(x) - 1 - x) / x**2
    real, intent(in) :: x
    real             :: phi

    ! local variables
    real,    parameter :: xr = 0.2
    integer, parameter :: m = 10
    integer            :: i
    real               :: fact

    if (abs(x) > xr) then
      phi = (exp(x) - 1 - x) / x**2
    else ! argument x is too small => use taylor series
      phi  = 0
      fact = 1
      do i = 0, m
        fact = fact * (i + 2)
        phi  = phi + x**i / fact
      end do
    end if
  end function

  elemental function dphi2dx(x) result(dphidx)
    !! derivative of phi2(x)
    real, intent(in) :: x
    real             :: dphidx

    ! local variables
    real,    parameter :: xr = 0.2
    integer, parameter :: m = 10
    integer            :: i
    real               :: fact

    if (abs(x) > xr) then
      dphidx = (exp(x) * (x - 2) + x + 2) / x**3
    else ! argument x is too small => use taylor series
      dphidx = 0
      fact   = 2
      do i = 1, m
        fact = fact * (i + 2)
        dphidx = dphidx + i * x**(i-1) / fact
      end do
    end if
  end function

  elemental function m_expm1(x) result(e)
    !! exp(x) - 1; accurate even for x close to 0
    real, intent(in) :: x
    real             :: e

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      e = x
    else
      e = exp(x)

      if (e == 1.0) then
        e = x
      else if (e - 1.0 == - 1.0) then
        e = -1
      else
        e = (e - 1.0) * x / log(e)
      end if
    end if
  end function

  elemental function m_log1p(x) result(l)
    !! log(1 + x); accurate even for x close to 0
    real, intent(in) :: x
    real             :: l

    real :: u, d

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      l = x
    else
      u = 1.0 + x
      d = u - 1.0

      if (d == 0) then
        l = x
      else
        l = log(u) * x / d
      end if
    end if
  end function

  pure function linspace(x0, x1, nx) result(x)
    !! create array of linear spaced values

    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    integer :: i
    real    :: dx

    ! spacing between values
    dx = (x1 - x0) / (nx - 1)

    x(1) = x0
    do i = 2, nx - 1
      x(i) = x0 + (i - 1) * dx
    end do
    x(nx) = x1
  end function

  function logspace(x0, x1, nx) result(x)
    !! create array of logarithmic spaced values

    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    integer :: i
    real    :: e0, e1, de

    m4_assert(x0 > 0.0)
    m4_assert(x1 > 0.0)

    e0 = log(x0)
    e1 = log(x1)
    de = (e1 - e0) / (nx - 1)

    x(1) = x0
    do i = 2, nx - 1
      x(i) = exp(e0 + (i - 1) * de)
    end do
    x(nx) = x1
  end function

  pure function eye_int(n) result(e)
    !! integer identity matrix

    integer, intent(in) :: n
      !! dimension
    integer             :: e(n,n)
      !! return integer identity matrix

    integer :: i

    e = 0
    do i = 1, n
      e(i,i) = 1
    end do
  end function

  pure function eye_real(n) result(e)
    !! real identity matrix

    integer, intent(in) :: n
      !! dimension
    real                :: e(n,n)
      !! return real identity matrix

    integer :: i

    e = 0.0
    do i = 1, n
      e(i,i) = 1.0
    end do
  end function

  pure function norm_inf(arr)
    !! calculates the infinity norm $$ || \ \cdot \  ||_\infty $$ of the given array.

    real, intent(in) :: arr(:)
    real             :: norm_inf

    norm_inf = maxval(abs(arr))
  end function

  pure function check_lin_dep(M, rtol) result(l)
    !! check if row vectors of M are linearly dependent
    real,           intent(in) :: M(:,:)
      !! input matrix
    real, optional, intent(in) :: rtol
      !! optional: relative tolerance (default: 1e-14)
    logical                    :: l
      !! return true if rows of M are lin. dep.; otherwise false

    integer :: nrows, ncols, ipiv(size(M,1)), row, row2, col, col2, pv, tmp
    real    :: rtol_, T(size(M,1),size(M,2)), nrm

    ! system size
    nrows = size(M,1)
    ncols = size(M,2)

    ! tolerance
    rtol_ = 1e-14
    if (present(rtol)) rtol_ = rtol

    ! init pivot indices
    do row = 1, nrows
      ipiv(row) = row
    end do

    ! work with copy of M
    T = M

    ! gauss elimination
    col = 0
    do row = 1, nrows
      col = col + 1

      ! find pivot
      do while (.true.)
        pv = row
        do row2 = row+1, nrows
          if (abs(T(ipiv(row2),col)) > abs(T(ipiv(pv),col))) pv = row2
        end do
        if (T(ipiv(pv),col) /= 0.0) exit
        col = col + 1
      end do

      ! virtually exchange current row with pivot
      if (pv /= row) then
        tmp       = ipiv(row)
        ipiv(row) = ipiv(pv)
        ipiv(pv ) = tmp
      end if

      ! gauss elimination for one column
      do row2 = row+1, nrows
        ! get row norm for zero check
        nrm = norm_inf(T(ipiv(row2),col:ncols))

        ! update row
        do col2 = col+1, ncols
          T(ipiv(row2),col2) = T(ipiv(row2),col2) - T(ipiv(row),col2) * T(ipiv(row2),col) / T(ipiv(row),col)
        end do

        ! check if row is zero => linearly dependent
        if (maxval(abs(T(ipiv(row2),col+1:ncols))) < rtol_ * nrm) then
          l = .true.
          return
        end if
      end do
    end do

    l = .false.
  end function

  function interp1(x, v, xq) result(vq)
    !! interpolates data piecewise linearly
    !!
    !! similar to MATLAB interp1: https://de.mathworks.com/help/matlab/ref/interp1.html
    real, intent(in) :: x(:)
      !! x-data vector
    real, intent(in) :: v(:)
      !! y-data vector
    real, intent(in) :: xq
      !! x query point
    real             :: vq
      !! y query value

    integer :: i

    m4_assert(size(x) == size(v))

    if      (xq <= x(1      )) then
      vq = v(1)

    else if (xq >= x(size(x))) then
      vq = v(size(x))

    else
      i  = bin_search(x, xq, mode=BS_LESS)
      vq = (v(i) * (x(i+1) - xq) + v(i+1) * (xq - x(i))) / (x(i+1) - x(i))
    end if
  end function

  function roots(p) result(rts)
    !! Find complex roots of polynomial
    real, intent(in) :: p(:)
      !! coefficients of polynomial in reduced form: f(x) = p(1) + p(2) * x + p(3) * x**2 + ... + x**n
    complex          :: rts(size(p))
      !! return roots of polynomial

    integer :: n, nz, i, perm(size(p))
    real    :: real_rts(size(p))

    ! get degree
    n = size(p)
    m4_assert(n > 0)

    ! every leading zero in p is a zero in rts
    do nz = 0, n-1
      if (p(nz+1) /= 0) exit
    end do
    rts(1:nz) = 0
    n = n - nz

    ! calculate non-zero roots
    select case (n)
    case (0) ! constant: f(x) = 1 == 0

    case (1) ! linear: f(x) = p(nz+1) + x == 0
      rts(nz + 1) = - p(nz + 1)

    case (2) ! quadratic: f(x) = p(nz+1) + p(nz+2) * x + x**2 == 0
      block
        complex :: d

        ! discriminant
        d = p(nz + 2)**2 - 4*p(nz + 1)

        ! add numbers of similar sign to prevent loss of accuracy
        if (p(nz + 2) >= 0) then
          rts(nz + 1) = - 0.5 * (sqrt(d) + p(nz + 2))
          rts(nz + 2) = p(nz + 1) / rts(nz + 1)
        else
          rts(nz + 2) = 0.5 * (sqrt(d) - p(nz + 2))
          rts(nz + 1) = p(nz + 1) / rts(nz + 2)
        end if
      end block

    case (3) ! cubic: f(x) = p(nz+1) + p(nz+2) * x + p(nz+3) * x**2 + x**3 == 0
      block
        real :: q, r, d, A, t1, x2, y2, th, ph1, ph2, ph3

        q = p(nz+2) / 3.0 - p(nz+3)**2 / 9.0
        r = (p(nz+2) * p(nz+3) - 3 * p(nz+1)) / 6.0 - p(nz+3)**3 / 27.0

        ! discriminant
        d = r**2 + q**3

        if (d > 0) then ! one real solution
          A = (abs(r) + sqrt(d)) ** (1.0 / 3.0)

          if (r >= 0) then
            t1 = A - q/A
          else
            t1 = q/A - A
          end if

          x2 = - 0.5 * t1 - p(nz+3) / 3.0
          y2 = 0.5 * sqrt(3.0) * (A + q/A)

          rts(nz + 1) = t1 - p(nz+3) / 3.0
          rts(nz + 2) = cmplx(x2,  y2)
          rts(nz + 3) = cmplx(x2, -y2)
        else ! three real solutions
          if (q == 0) then
            th = 0
          else
            th = acos(r / (-q) ** 1.5)
          end if

          ph1 = th / 3.0
          ph2 = ph1 - 2 * PI / 3.0
          ph3 = ph1 + 2 * PI / 3.0

          rts(nz + 1) = 2 * sqrt(-q) * cos(ph1) - p(nz+3) / 3.0
          rts(nz + 2) = 2 * sqrt(-q) * cos(ph2) - p(nz+3) / 3.0
          rts(nz + 3) = 2 * sqrt(-q) * cos(ph3) - p(nz+3) / 3.0
        end if
      end block

    case default ! default: f(x) = p(nz+1) + p(nz+2) * x + p(nz+3) * x**2 + ... + p(nz+n) * x**(n-1) + x**n == 0
      block
        real :: a(n,n), wr(n), wi(n)

        ! construct companion matrix
        a = 0
        do i = 1, n
          a(1,i) = -p(nz + n - i + 1)
        end do
        do i = 2, n
          a(i,i-1) = 1
        end do

        ! get eigenvalues of companion matrix
        call geev(a, wr, wi)

        ! construct complex roots
        do i = 1, n
          rts(nz + i) = cmplx(wr(i), wi(i))
        end do
      end block

    end select

    ! sort by real part
    real_rts = real(rts)
    call qsort(real_rts, perm = perm)
    rts = rts(perm)
  end function

  function polyg_area_2d(p) result(A)
    !! calculate the area of an arbitrary polygon in 2D
    real, intent(in) :: p(:,:)
      !! points of polygon
    real             :: A
      !! area of polygon

    integer :: i, n

    m4_assert(size(p, 1) == 2)
    n = size(p, 2)

    A = 0.0
    do i=1, n-1
      A = A + cross_product_2d(p(:,i), p(:,i+1))
    end do
    A = 0.5*abs(A + cross_product_2d(p(:,n), p(:,1)))
  end function

end module
