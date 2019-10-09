program test
  use test_matrix_m
  use test_qsort_m
  implicit none

#define all 1

#define matrix 2
#if ((USE_TEST==matrix) .or. (USE_TEST==all))
  call test_matrix()
#endif
#undef matrix

#define qsort 3
#if ((USE_TEST==qsort) .or. (USE_TEST==all))
  call test_qsort()
#endif
#undef qsort

#undef all

end program test
