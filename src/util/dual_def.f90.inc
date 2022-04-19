type dual_{}N
  !! dual number type for automatic differentiation (forward mode)

  real :: x
    !! value of expression
  real :: dx(N)
    !! derivatives wrt. base variables
contains
  procedure :: init => dual_{}N{}_init
end type

interface operator(+)
  module procedure :: dual_add_dual_{}N
  module procedure :: dual_add_real_{}N
  module procedure :: dual_add_int_{}N
  module procedure :: real_add_dual_{}N
  module procedure :: int_add_dual_{}N
end interface

interface operator(-)
  module procedure :: dual_neg_{}N
  module procedure :: dual_sub_dual_{}N
  module procedure :: dual_sub_real_{}N
  module procedure :: dual_sub_int_{}N
  module procedure :: real_sub_dual_{}N
  module procedure :: int_sub_dual_{}N
end interface

interface operator(*)
  module procedure :: dual_mul_dual_{}N
  module procedure :: dual_mul_real_{}N
  module procedure :: dual_mul_int_{}N
  module procedure :: real_mul_dual_{}N
  module procedure :: int_mul_dual_{}N
end interface

interface operator(/)
  module procedure :: dual_div_dual_{}N
  module procedure :: dual_div_real_{}N
  module procedure :: dual_div_int_{}N
  module procedure :: real_div_dual_{}N
  module procedure :: int_div_dual_{}N
end interface

interface operator(**)
  module procedure :: dual_pow_dual_{}N
  module procedure :: dual_pow_real_{}N
  module procedure :: dual_pow_int_{}N
  module procedure :: real_pow_dual_{}N
  module procedure :: int_pow_dual_{}N
end interface

interface dot_product
  module procedure :: dual_dot_dual_{}N
  module procedure :: dual_dot_real_{}N
  module procedure :: real_dot_dual_{}N
end interface

interface abs
  module procedure :: dual_abs_{}N
end interface

interface sum
  module procedure :: dual_sum_{}N
end interface

interface sqrt
  module procedure :: dual_sqrt_{}N
end interface

interface exp
  module procedure :: dual_exp_{}N
end interface

interface log
  module procedure :: dual_log_{}N
end interface

interface sin
  module procedure :: dual_sin_{}N
end interface

interface cos
  module procedure :: dual_cos_{}N
end interface

interface tan
  module procedure :: dual_tan_{}N
end interface

m4_undefine({N})
