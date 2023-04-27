m4_include(macro.f90.inc)

module dual_m

  use error_m, only: assert_failed
  use math_m,  only: expm1, log1p

  implicit none

  private

  public operator(+), operator(-), operator(*), operator(/), operator(**)
  public abs, cos, dot_product, exp, expm1, log, log1p, sin, sqrt, sum, tan

  ! maximum supported dimension
  m4_define({m4_max_dim},{8})

  ! get dimension list (1 to max_dim)
  m4_define({m4_list_help},{
    m4_ifelse($1,0,,{m4_list_help(m4_decr($1))})
    m4_X($1)
  })
  m4_define({m4_list},{m4_list_help(m4_max_dim)})

  ! make defined types public
  m4_define({m4_X},{public dual_$1})
  m4_list

  ! include type definitions
  m4_define({m4_X},{
    m4_define({N},$1)
    m4_include(dual_def.f90.inc)
  })
  m4_list

contains

  ! include type procedure implementations
  m4_define({m4_X},{
    m4_define({N},$1)
    m4_include(dual_imp.f90.inc)
  })
  m4_list

end module
