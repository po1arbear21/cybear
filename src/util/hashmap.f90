module hashmap_m
  use util_m, only : hash
  use array_m
  use vector_m
  implicit none

  private
  public :: int2, int3, int4
  public :: hashmap_int_int, hashmap_int2_int, hashmap_int3_int, hashmap_int4_int

  type int2
    integer :: i(2)
  end type

  type int3
    integer :: i(3)
  end type

  type int4
    integer :: i(4)
  end type

  interface hash
    module procedure :: hash_int2, hash_int3, hash_int4
  end interface

  interface key_compare
    module procedure :: key_compare_int, key_compare_int2, key_compare_int3, key_compare_int4
  end interface

#define T int2
#define TT type(int2)
#include "vector_def.f90.inc"

#define T int3
#define TT type(int3)
#include "vector_def.f90.inc"

#define T int4
#define TT type(int4)
#include "vector_def.f90.inc"

#define TKEY int
#define TTKEY integer
#define TVALUE int
#define TTVALUE integer
#include "hashmap_def.f90.inc"

#define TKEY int2
#define TTKEY type(int2)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_def.f90.inc"

#define TKEY int3
#define TTKEY type(int3)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_def.f90.inc"

#define TKEY int4
#define TTKEY type(int4)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_def.f90.inc"

contains

#define T int2
#define TT type(int2)
#include "vector_imp.f90.inc"

#define T int3
#define TT type(int3)
#include "vector_imp.f90.inc"

#define T int4
#define TT type(int4)
#include "vector_imp.f90.inc"

#define TKEY int
#define TTKEY integer
#define TVALUE int
#define TTVALUE integer
#include "hashmap_imp.f90.inc"

#define TKEY int2
#define TTKEY type(int2)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_imp.f90.inc"

#define TKEY int3
#define TTKEY type(int3)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_imp.f90.inc"

#define TKEY int4
#define TTKEY type(int4)
#define TVALUE int
#define TTVALUE integer
#include "hashmap_imp.f90.inc"

  function hash_int2(i2) result(h)
    type(int2), intent(in) :: i2
    integer                :: h

    h = hash(i2%i)
  end function

  function hash_int3(i3) result(h)
    type(int3), intent(in) :: i3
    integer                :: h

    h = hash(i3%i)
  end function

  function hash_int4(i4) result(h)
    type(int4), intent(in) :: i4
    integer                :: h

    h = hash(i4%i)
  end function

  function key_compare_int(k1, k2) result(equal)
    integer, intent(in) :: k1
    integer, intent(in) :: k2
    logical             :: equal

    equal = (k1 == k2)
  end function

  function key_compare_int2(k1, k2) result(equal)
    type(int2), intent(in) :: k1
    type(int2), intent(in) :: k2
    logical                :: equal

    equal = all(k1%i == k2%i)
  end function

  function key_compare_int3(k1, k2) result(equal)
    type(int3), intent(in) :: k1
    type(int3), intent(in) :: k2
    logical                :: equal

    equal = all(k1%i == k2%i)
  end function

  function key_compare_int4(k1, k2) result(equal)
    type(int4), intent(in) :: k1
    type(int4), intent(in) :: k2
    logical                :: equal

    equal = all(k1%i == k2%i)
  end function

end module
