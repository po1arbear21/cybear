#ifndef MACRO_F90_INC
#define MACRO_F90_INC

! macro for suppressing unused argument warnings
#define IGNORE(x) associate( x => x ); end associate

#endif