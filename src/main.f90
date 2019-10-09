#include "src/util/assert.f90.inc"

program main
  use error_m
  implicit none
  integer :: x

  x = 8

  ASSERT(x == 7)
end program main
