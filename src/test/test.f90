m4_include(../util/macro.f90.inc)

program test

  use test_analysis_m
  use test_bin_search_m
  use test_circuit_m
  use test_deque_m
  use test_distributions_m
  use test_dual_m
  use test_esystem_m
  use test_esystem_depgraph_m
  use test_esystem_prec_m
  use test_expm_m
  use test_gmres_m
  use test_grid_m
  use test_grid_table_m
  use test_hashmap_m
  use test_high_precision_m
  use test_input_m
  use test_input_src_m
  use test_json_m
  use test_logging_m
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
  use test_string_m
  use test_util_m
  use test_vector_m
  m4_ifdef({m4_feast},{
  use test_feast_m
  })
  m4_ifdef({m4_ilupack},{
  use test_ilupack_m
  })
  m4_ifdef({m4_mumps},{
  use test_mumps_m
  })
  m4_ifdef({m4_quadpack},{
  use test_quadpack_m
  })
  m4_ifdef({m4_mpfr},{
  use test_mpfr_m
  })
  m4_ifdef({m4_spike},{
  use test_spike_m
  })
  m4_ifdef({m4_triangle},{
  use test_triangle_m
  })

  implicit none

  call test_analysis()
  call test_bin_search()
  call test_circuit()
  call test_deque()
  call test_distributions()
  call test_dual()
  call test_esystem()
  call test_eval_list1()
  call test_eval_list2()
  call test_esystem_prec()
  call test_expm()
  call test_gmres()
  call test_grid()
  call test_grid_table()
  call test_hashmap()
  call test_high_precision()
  call test_input()
  call test_input_src()
  call test_json()
  call test_logging()
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
  call test_string()
  call test_util()
  call test_vector()

  m4_ifdef({m4_feast},{call test_feast()})
  m4_ifdef({m4_ilupack},{call test_ilupack()})
  m4_ifdef({m4_mumps},{call test_mumps()})
  m4_ifdef({m4_quadpack},{call test_quadpack()})
  m4_ifdef({m4_mpfr},{call test_mpfr()})
  m4_ifdef({m4_spike},{call test_spike()})
  m4_ifdef({m4_triangle},{call test_triangle()})

end program
