#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

#define DUAL_TYPE            CONCAT(dual,N)
#define DUAL_INIT            CONCAT(dual_init,N)
#define DUAL_ADD_DUAL        CONCAT(dual_add_dual,N)
#define DUAL_ADD_REAL        CONCAT(dual_add_real,N)
#define REAL_ADD_DUAL        CONCAT(real_add_dual,N)
#define DUAL_NEG             CONCAT(dual_neg,N)
#define DUAL_SUB_DUAL        CONCAT(dual_sub_dual,N)
#define DUAL_SUB_REAL        CONCAT(dual_sub_real,N)
#define REAL_SUB_DUAL        CONCAT(real_sub_dual,N)
#define DUAL_MUL_DUAL        CONCAT(dual_mul_dual,N)
#define DUAL_MUL_REAL        CONCAT(dual_mul_real,N)
#define REAL_MUL_DUAL        CONCAT(real_mul_dual,N)
#define DUAL_DIV_DUAL        CONCAT(dual_div_dual,N)
#define DUAL_DIV_REAL        CONCAT(dual_div_real,N)
#define REAL_DIV_DUAL        CONCAT(real_div_dual,N)
#define DUAL_POW_DUAL        CONCAT(dual_pow_dual,N)
#define DUAL_POW_REAL        CONCAT(dual_pow_real,N)
#define DUAL_POW_INT         CONCAT(dual_pow_int,N)
#define REAL_POW_DUAL        CONCAT(real_pow_dual,N)
#define DUAL_DOT_PRODUCT     CONCAT(dual_dot_product,N)
#define DUAL_ABS             CONCAT(dual_abs,N)
#define DUAL_SQRT            CONCAT(dual_sqrt,N)
#define DUAL_EXP             CONCAT(dual_exp,N)
#define DUAL_LOG             CONCAT(dual_log,N)
#define DUAL_SIN             CONCAT(dual_sin,N)
#define DUAL_COS             CONCAT(dual_cos,N)
#define DUAL_TAN             CONCAT(dual_tan,N)

type DUAL_TYPE
  !! dual number type for automatic differentiation (forward mode)

  real :: x
    !! value of expression
  real :: dx(N)
    !! derivatives wrt. base variables
contains
  procedure :: init => DUAL_INIT
end type

interface operator(+)
  module procedure :: DUAL_ADD_DUAL
  module procedure :: DUAL_ADD_REAL
  module procedure :: REAL_ADD_DUAL
end interface

interface operator(-)
  module procedure :: DUAL_NEG
  module procedure :: DUAL_SUB_DUAL
  module procedure :: DUAL_SUB_REAL
  module procedure :: REAL_SUB_DUAL
end interface

interface operator(*)
  module procedure :: DUAL_MUL_DUAL
  module procedure :: DUAL_MUL_REAL
  module procedure :: REAL_MUL_DUAL
end interface

interface operator(/)
  module procedure :: DUAL_DIV_DUAL
  module procedure :: DUAL_DIV_REAL
  module procedure :: REAL_DIV_DUAL
end interface

interface operator(**)
  module procedure :: DUAL_POW_DUAL
  module procedure :: DUAL_POW_REAL
  module procedure :: DUAL_POW_INT
  module procedure :: REAL_POW_DUAL
end interface

interface operator(.dot.)
  module procedure :: DUAL_DOT_PRODUCT
end interface

interface abs
  module procedure :: DUAL_ABS
end interface

interface sqrt
  module procedure :: DUAL_SQRT
end interface

interface exp
  module procedure :: DUAL_EXP
end interface

interface log
  module procedure :: DUAL_LOG
end interface

interface sin
  module procedure :: DUAL_SIN
end interface

interface cos
  module procedure :: DUAL_COS
end interface

interface tan
  module procedure :: DUAL_TAN
end interface

#undef N

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef DUAL_TYPE
#undef DUAL_INIT
#undef DUAL_ADD_DUAL
#undef DUAL_ADD_REAL
#undef REAL_ADD_DUAL
#undef DUAL_NEG
#undef DUAL_SUB_DUAL
#undef DUAL_SUB_REAL
#undef REAL_SUB_DUAL
#undef DUAL_MUL_DUAL
#undef DUAL_MUL_REAL
#undef REAL_MUL_DUAL
#undef DUAL_DIV_DUAL
#undef DUAL_DIV_REAL
#undef REAL_DIV_DUAL
#undef DUAL_POW_DUAL
#undef DUAL_POW_REAL
#undef DUAL_POW_INT
#undef REAL_POW_DUAL
#undef DUAL_DOT_PRODUCT
#undef DUAL_ABS
#undef DUAL_SQRT
#undef DUAL_EXP
#undef DUAL_LOG
#undef DUAL_SIN
#undef DUAL_COS
#undef DUAL_TAN