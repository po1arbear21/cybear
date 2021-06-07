module test_math_m

  use ieee_arithmetic, only: ieee_next_after, ieee_value, ieee_positive_inf, ieee_negative_inf
  use math_m
  use util_m,          only: int2str
  use test_case_m,     only: test_case

  implicit none

  private
  public test_math

contains

  subroutine test_math()
    type(test_case) :: tc

    call tc%init("math")

    ! cross product
    ! 1. check: (1,2,3)x(1,2,3)=(0,0,0)
    ! 2. check: (1,2,3)x(3,2,1)=(-4,8,-4)
    block
      real :: a(3), b(3), axa(3), axa_exp(3), axb(3), axb_exp(3)

      a = [1, 2, 3]
      b = [3, 2, 1]

      axa = cross_product(a, a)
      axb = cross_product(a, b)

      axa_exp = 0
      axb_exp = [-4, 8, -4]

      call tc%assert_eq(axa_exp, axa, 1e-15, "cross product 1")
      call tc%assert_eq(axb_exp, axb, 1e-15, "cross product 2")
    end block

    ! cross product 2d
    ! 1. check: (1,2)x(1,2)=0
    ! 2. check: (1,2)x(2,1)=-3
    block
      real :: a(2), b(2), axa, axb, axa_exp, axb_exp

      a = [1, 2]
      b = [2, 1]

      axa = cross_product_2d(a, a)
      axb = cross_product_2d(a, b)

      axa_exp =  0.0
      axb_exp = -3.0

      call tc%assert_eq(axa_exp, axa, 1e-15, "cross product 2d 1")
      call tc%assert_eq(axb_exp, axb, 1e-15, "cross product 2d 2")
    end block

    ! heaviside
    ! testing single values from x
    block
      real :: x(7), y(7), y_exp(7)

      x     = [-100.0, -10.0, -epsilon(1.0), 0.0, epsilon(1.0), 10.0, 100.0 ]
      y_exp = [   0.0,   0.0,  0.0,          0.5, 1.0,           1.0,   1.0]

      y     = heaviside(x)
      call tc%assert_eq(y_exp, y, 1e-15, "heavyside")
    end block

    ! isinf(x)
    ! testing single values from x
    block
      logical :: y(5), y_exp(5)
      real    :: x(5)

      x     = [ieee_value(1.0, ieee_negative_inf),  ieee_value(1.0, ieee_positive_inf), 0.0,    100.0, 1.0/epsilon(1.0)]
      y_exp = [  .true.,                               .true.,                        .false., .false.,     .false.]

      y = isinf(x)

      call tc%assert(y_exp .eqv. y, "is inf")
    end block

    ! ber(x) = x / (exp(x) - 1)
    ! testing single values from x and smoothness around |x| = 1e-6
    block
      real :: x(6), y(6), y_exp(6)

      x     = [0.0, -5.0,                 -100.0, 15.0,              -1e-6,              1e-6]
      y_exp = [1.0,  5.0339182745315121,   100.0,  0.000004588536211, 1.000000500000083, 0.999999500000083]

      y     = ber(x)

      call tc%assert_eq(y_exp, y, 1e-13, "ber(x) values")
      call tc%assert_eq(y(5), ber(ieee_next_after(-1e-6,-1.0)), 1e-15, "ber(x) smoothness negative")
      call tc%assert_eq(y(6), ber(ieee_next_after( 1e-6, 1.0)), 1e-15, "ber(x) smoothness positive")
    end block

    ! dberdx(x)
    ! testing single values from x and smoothness around |x| = 1e-6
    block
      real :: x(6), y(6), y_exp(6)

      x     = [ 0.0, -5.0,                -100.0, 15.0,             -1e-6,                1e-6]
      y_exp = [-0.5, -0.9726352905053439,  - 1.0, -0.00000428263520,-0.5000001666666667, -0.4999998333333333]

      y     = dberdx(x)
      call tc%assert_eq(y_exp, y,                                   1e-15, "dberdx(x) values"             )
      call tc%assert_eq(y(5),  dberdx(ieee_next_after(-1e-6,-1.0)), 1e-8,  "dberdx(x) smoothness negative")
      call tc%assert_eq(y(6),  dberdx(ieee_next_after( 1e-6, 1.0)), 1e-8,  "dberdx(x) smoothness positive")
    end block

    ! phi1(x)
    ! testing single values from x and smoothness around |x| = 1e-6
    block
      real :: x(5), y(5), y_exp(5)

      x     = [0.0, -5.0,              -100.0, -1e-6,              1e-6]
      y_exp = [1.0,  0.198652410600182,   0.01, 0.999999500000166, 1.000000500000166]

      y     = phi1(x)

      call tc%assert_eq(y_exp, y, 1e-15, "phi1(x) values")
      call tc%assert_eq(y(4), phi1(ieee_next_after(-1e-6,-1.0)), 1e-15, "phi1(x) smoothness negative")
      call tc%assert_eq(y(5), phi1(ieee_next_after( 1e-6, 1.0)), 1e-15, "phi1(x) smoothness positive")
    end block

    ! dphi1dx
    ! testing single values from x and smoothness around |x| = 1e-6
    block
      real :: x(5), y(5), y_exp(5)

      x     = [0.0, -5.0,                 -100.0,   -1e-6,                  1e-6]
      y_exp = [0.5,  0.03838289272021949,    0.0001, 0.4999996666667916666, 0.5000003333334583334]

      y     = dphi1dx(x)

      call tc%assert_eq(y_exp, y, 1e-15, "dphi1dx(x) values")
      call tc%assert_eq(y(4), dphi1dx(ieee_next_after(-1e-6,-1.0)), 1e-9, "dphi1dx(x) smoothness negative")
      call tc%assert_eq(y(5), dphi1dx(ieee_next_after( 1e-6, 1.0)), 1e-9, "dphi1dx(x) smoothness positive")
    end block

    ! phi2(x)
    ! testing single values from x and smoothness around |x| = 0.2
    block
      real :: x(5), y(5), y_exp(5)

      x     = [0.0, -5.0,               -100.0,   -0.2,                 0.2]
      y_exp = [0.5,  0.1602695178799634,   0.0099, 0.46826882694954647, 0.53506895400424585]

      y     = phi2(x)

      call tc%assert_eq(y_exp, y, 1e-15, "phi2(x) values")
      call tc%assert_eq(y(4), phi2(ieee_next_after(-0.2,-1.0)), 1e-15, "phi2(x) smoothness negative")
      call tc%assert_eq(y(5), phi2(ieee_next_after( 0.2, 1.0)), 1e-15, "phi2(x) smoothness positive")
    end block

    ! dphi2dx
    ! testing single values from x and smoothness around |x| = 0.2
    block
      real :: x(5), y(5), y_exp(5)

      x     = [    0.0,                -5.0,   -100.0,                -0.2,               0.2]
      y_exp = [1.0/6.0, 0.02437732503194879, 0.000098, 0.1509570964450111, 0.1843794139617874]

      y     = dphi2dx(x)

      call tc%assert_eq(y_exp, y, 1e-15, "dphi2dx(x) values")
      call tc%assert_eq(y(4), dphi2dx(ieee_next_after(-0.2,-1.0)), 1e-13, "dphi2dx(x) smoothness negative")
      call tc%assert_eq(y(5), dphi2dx(ieee_next_after( 0.2, 1.0)), 1e-13, "dphi2dx(x) smoothness positive")
    end block

    ! expm1 -> exp(x) - 1
    ! testing single values from x
    block
      real :: x(6), y(6), y_exp(6)

      x     = [epsilon(1.0), -epsilon(1.0), 0.1,                -1.0,                1e-5,                   -0.2143]
      y_exp = [epsilon(1.0), -epsilon(1.0), 0.1051709180756476, -0.6321205588285577, 0.00001000005000016667, -0.1928937831657808]

      y     = expm1(x)

      call tc%assert_eq(y_exp, y, 1e-15, "expm1 (exp(x) - 1)")
    end block

    ! log1p(x) -> log(1+x)
    ! testing single values from x
    block
      real :: x(5), y(5), y_exp(5)

      x     = [epsilon(1.0), -epsilon(1.0), 0.1,                 1e-5,             -0.2143]
      y_exp = [epsilon(1.0), -epsilon(1.0), 0.09531017980432486, 0.00000999995000, -0.2411802388003611]

      y     = log1p(x)

      call tc%assert_eq(y_exp, y, 1e-15,                                  "log1p (log(1+x)")
      call tc%assert((ieee_value(1.0, ieee_negative_inf) == log1p(-1.0)), "log1p x=-1 y=-inf")
    end block

    ! linspace
    ! testing the spacing of the different lists
    block
      integer           :: n(4), i, j
      real              :: a(4), b(4)
      real, allocatable :: list(:), dif(:), dif_exp(:)

      a = [  0.0,  10.0,  431.23140,   0.0]
      b = [100.0, -10.0, 1233.42576,  1e-5]
      n = [    5,    20,         41, 12352]

      do i = 1, 4
        list    = linspace(a(i), b(i), n(i))
        dif_exp = [((b(i)-a(i))/(n(i)-1), j = 1, n(i)-1)]

        dif     = [(list(j+1) - list(j), j = 1, n(i)-1)]

        call tc%assert_eq(dif_exp, dif, 1e-15*abs(b(i)-a(i)), "linspace "//int2str(i))
      end do
    end block

    ! logspace
    ! testing the spacing of the different lists
    block
      real                :: a(4), b(4)
      integer             :: n(4), i, j
      real , allocatable  :: list(:), quot(:), quot_exp(:)

      a = [1.0,   10.0,  431.2314,  1.0]
      b = [100.0,  3.0, 1233.42576, 1e-5]
      n = [5,     20,     41,       12352]

      do i = 1, 4
        list     = logspace(a(i), b(i), n(i))
        quot_exp = [((b(i)/a(i))**(1.0/(n(i)-1)), j = 1, n(i)-1)]

        quot     = [(list(j+1)/list(j), j = 1, n(i)-1)] ! taking the difference between the j-th and (j+1)th element

        call tc%assert_eq(quot_exp, quot, 1e-12, "logspace "//int2str(i))
      end do

    end block

    ! eye_int
    ! testing the spacing of the different lists
    block
      integer :: e(5,5), a(5), b(5)

      a = [0, 10, 431, -531, 3]
      e = eye_int(5)

      b = matmul(e,a)

      call tc%assert_eq(a, b, "eye_int")
    end block

    ! eye_real
    ! testing the spacing of the different lists
    block
      real :: e(5,5), a(5), b(5)

      a = [0.0, 10.24535, 431.2354, -4563.7895446, epsilon(1.0)]
      e = eye_real(5)

      b = matmul(e, a)

      call tc%assert_eq(a, b, 1e-15, "eye_int")
    end block

    ! norm_inf
    ! testing the biggest value of a list
    block
      real :: a(5), b, b_exp

      a = [0.0, -10.24535, 431.2354, -4563.7895446, epsilon(1.0)]
      b_exp = 4563.7895446

      b = norm_inf(a)

      call tc%assert_eq(b_exp, b, 1e-15, "norm_inf")
    end block

    ! check linear dependence
    block
      real    :: M(3,8)
      logical :: l

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      M(2,:) = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 1")

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      M(2,:) = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 2")

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      M(2,:) = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(l, "check linear dependence 3")

      M(1,:) = [0.3371, 0.6020, 0.0838, 0.9961, 0.1622, 0.2630, 0.2290, 0.0782]
      M(2,:) = [0.7943, 0.6541, 0.9133, 0.4427, 0.3112, 0.6892, 0.1524, 0.1067]
      M(3,:) = [2.1427, 3.0621, 1.2485, 4.4271, 0.9600, 1.7412, 1.0684, 0.4195]
      l = check_lin_dep(M)
      call tc%assert(l, "check linear dependence 4")

      M(2,8) = M(2,8) + 0.1
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 5")
    end block

    ! interp1
    ! testing an f(x)=x interpolation
    block
      integer :: i
      real    :: x(10), y(10), xq(5), yq(5), y_exp(5)

      x = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
      y = [-1, 2, 2, 3, 5, 4, 6, 7, 7, 10]

      xq    = [-10.0,  0.0, epsilon(1.0),        4.5, 10.0]
      y_exp = [ -1.0, -1.0, -1 + 3*epsilon(1.0), 4.5, 10.0]

      yq = [(interp1(x, y, xq(i)), i = 1, 5)]

      call tc%assert_eq(y_exp, yq, 1e-15, "interp1")
    end block

    call tc%finish()
  end subroutine


end module
