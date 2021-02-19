program test

  use test_arnoldi_m
  use test_deque_m
  use test_dual_m
  use test_esystem_m
  use test_expm_m
  use test_gmres_m
  use test_grid_m
  use test_grid_table_m
  use test_high_precision_m
  use test_input_m
  use test_math_m
  use test_matop_m
  use test_matrix_m
  use test_newton_m
  use test_normalization_m
  use test_plotmtv_m
  use test_poly_m
  use test_qsort_m
  use test_radau5_m
  use test_random_m
  use test_schur_m
  use test_sqrtm_m
  use test_util_m
  use test_vector_m

  implicit none

  call test_gmres()
  call test_matop()
  call test_poly()
  call test_random()
  call test_high_precision()
  call test_math()
  call test_deque()
  call test_input()
  call test_expm()
  call test_normalization()
  call test_dual()
  call test_radau5()
  call test_newton()
  call test_matrix()
  call test_qsort()
  call test_vector()
  call test_sqrtm()
  call test_arnoldi()
  call test_schur()
  call test_plotmtv()
  call test_util()
  call test_grid()
  call test_grid_table()
  call test_esystem()

#ifdef USE_ILUPACK
  block
    use test_ilupack_m

    call test_ilupack()
  end block
#endif

#ifdef USE_MUMPS
  block
    use test_mumps_m

    call test_mumps()
  end block
#endif

#ifdef USE_QUADPACK
  block
    use quadpack_m

    call test_quadpack()
  end block
#endif

end program
