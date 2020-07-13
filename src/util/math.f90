#include "macro.f90.inc"

module math_m
  use error_m
  use ieee_arithmetic
  implicit none

  real, parameter :: PI = 3.141592653589793238462643

contains

  pure function cross_product(a, b) result(axb)
    real, intent(in) :: a(3), b(3)
    real             :: axb(3)

    axb = [a(2)*b(3) - a(3)*b(2), &
      &    a(3)*b(1) - a(1)*b(3), &
      &    a(1)*b(2) - a(2)*b(1)  ]
  end function

  pure real function cross_product_2d(a, b)
    real, intent(in) :: a(2), b(2)

    cross_product_2d = a(1)*b(2) - a(2)*b(1)
  end function

  elemental function heaviside(x) result(h)
    !! heaviside step function
    real, intent(in) :: x
    real             :: h

    h = 0.5 * (sign(1.0, x) + 1.0)
  end function

  elemental function isinf(x) result (r)
    !! return true if argument is positive of negative infinity
    real, intent(in) :: x
    logical          :: r

    r = (ieee_class(x) == ieee_positive_inf) .or. (ieee_class(x) == ieee_negative_inf)
  end function

  elemental function ber(x) result(b)
    !! bernoulli function
    real, intent(in) :: x
    real             :: b

    if (abs(x) > 1e-6) then
      b = 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x)
    else
      b = 1.0 + x * (-0.5 + x / 12.0)
    end if
  end function

  elemental function dberdx(x) result(dbdx)
    !! derivative of bernoulli function
    real, intent(in) :: x
    real             :: dbdx

    if (abs(x) > 1e-6) then
      dbdx = (2.0 * exp(-0.5 * x) * sinh(0.5 * x) - x) / (4.0 * sinh(0.5 * x)**2)
    else
      dbdx = -0.5 + x / 6.0
    end if
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

    ASSERT(x0 > 0.0)
    ASSERT(x1 > 0.0)

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
    !! Calculates the infinity norm $$ || \ \cdot \  ||_\infty $$ of the given array.

    real, intent(in) :: arr(:)
    real             :: norm_inf

    norm_inf = maxval(abs(arr))
  end function

  pure function check_lin_dep(M, rtol) result(l)
    !! Check if row vectors of M are linearly dependent
    real,           intent(in) :: M(:,:)
      !! input matrix
    real, optional, intent(in) :: rtol
      !! optional: relative tolerance (default: 1e-14)
    logical                    :: l
      !! return true if rows of M are lin. dep.; otherwise false

    integer :: nrows, ncols, i, j, k, ipiv(size(M,1))
    real    :: T(size(M,1),size(M,2)), tmp, rtol_

    ! tolerance
    rtol_ = 1e-14
    if (present(rtol)) rtol_ = rtol

    ! system size
    nrows = size(M,1)
    ncols = size(M,2)

    ! init pivot indices
    do i = 1, nrows
      ipiv(i) = i
    end do

    ! work with copy of M
    T = M

    ! gauss elimination
    do i = 1, nrows
      ! find pivot
      k = i
      do j = i+1, nrows
        if (abs(T(ipiv(j),i)) > abs(T(ipiv(k),i))) k = j
      end do

      ! virtually exchange row i with row k
      if (k /= i) then
        j       = ipiv(i)
        ipiv(i) = ipiv(k)
        ipiv(k) = j
      end if

      ! divide by pivot
      T(ipiv(i),i+1:ncols) = T(ipiv(i),i+1:ncols) / T(ipiv(i),i)

      ! gauss elimination for one column
      do j = i+1, nrows
        ! get norm of row for zero check
        tmp = norm_inf(T(ipiv(j),:))

        ! update row
        do k = i+1, ncols
          T(ipiv(j),k) = T(ipiv(j),k) - T(ipiv(i),k) * T(ipiv(j),i)
        end do

        ! check if row is zero => linearly dependent
        if (maxval(abs(T(ipiv(j),i+1:ncols))) < 1e-14 * tmp) then
          l = .true.
          return
        end if
      end do
    end do

    l = .false.
  end function

end module
