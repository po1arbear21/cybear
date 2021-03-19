#ifndef MACRO_F90_INC
#define MACRO_F90_INC

! macro for suppressing unused argument warnings
#define IGNORE(x) associate( x => x ); end associate

#ifdef DEBUG
#ifdef __GFORTRAN__
#define ASSERT(expr) if (.not. (expr)) call assert_failed("expr", __FILE__, __LINE__)
#else
#define ASSERT(expr) if (.not. (expr)) call assert_failed(#expr, __FILE__, __LINE__)
#endif
#else
#define ASSERT(expr)
#endif

#endif
