#ifndef ASSERT_F90_INC
#define ASSERT_F90_INC

#ifdef DEBUG
#define ASSERT(expr) if (.not. (expr)) call assert_failed(#expr, __FILE__, __LINE__)
#else
#define ASSERT(expr)
#endif

#endif