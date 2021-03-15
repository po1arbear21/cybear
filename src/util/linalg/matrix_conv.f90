#include "../macro.f90.inc"

submodule (matrix_m) matrix_conv_m
  !! defines conversions of matrix types.
  !! submodule: will result in faster compilation and better code separation.

contains

#define T real
#define TT real
#include "matrix_conv_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_conv_imp.f90.inc"

end submodule
