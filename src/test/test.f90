program test

  use test_analysis_m
  use test_arnoldi_m
  use test_bin_search_m
  use test_deque_m
  use test_distributions_m
  use test_dual_m
  use test_esystem_m
  use test_esystem_prec_m
  use test_expm_m
  use test_gmres_m
  use test_grid_m
  use test_grid_table_m
  use test_hashmap_m
  use test_high_precision_m
  use test_input_m
  use test_map_m
  use test_math_m
  use test_matop_m
  use test_matrix_m
  use test_mkl_ilu_m
  use test_newton_m
  use test_normalization_m
  use test_plotmtv_m
  use test_poly_m
  use test_qsort_m
  use test_radau5_m
  use test_random_m
  use test_schur_m
  use test_sqrtm_m
  use test_string_m
  use test_util_m
  use test_variable_m
  use test_vector_m

#ifdef USE_FEAST
  use test_feast_m
#endif
#ifdef USE_ILUPACK
  use test_ilupack_m
#endif
#ifdef USE_MUMPS
  use test_mumps_m
#endif
#ifdef USE_QUADPACK
  use test_quadpack_m
#endif

  implicit none

  call test_analysis()
  call test_arnoldi()
  call test_bin_search()
  call test_deque()
  call test_distributions()
  call test_dual()
  call test_esystem()
  call test_esystem_prec()
  call test_expm()
  call test_gmres()
  call test_grid()
  call test_grid_table()
  call test_hashmap()
  call test_high_precision()
  call test_input()
  call test_map()
  call test_math()
  call test_matop()
  call test_matrix()
  call test_mkl_ilu()
  call test_newton()
  call test_normalization()
  call test_plotmtv()
  call test_poly()
  call test_qsort()
  call test_radau5()
  call test_random()
  call test_schur()
  call test_sqrtm()
  call test_string()
  call test_util()
  call test_variable()
  call test_vector()

#ifdef USE_FEAST
  call test_feast()
#endif
#ifdef USE_ILUPACK
  call test_ilupack()
#endif
#ifdef USE_MUMPS
  call test_mumps()
#endif
#ifdef USE_QUADPACK
  call test_quadpack()
#endif

end program
