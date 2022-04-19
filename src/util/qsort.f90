m4_include(macro.f90.inc)

module qsort_m

  use error_m,  only: assert_failed
  use string_m, only: string, operator(<), operator(<=), operator(>)

  implicit none

  private
  public qsort

  interface qsort
    module procedure :: qsort_int, qsort_string, qsort_real
  end interface

contains

  m4_define({m4_list},{
    m4_X(int)
    m4_X(real)
    m4_X(string)
  })
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(qsort_imp.f90.inc)
  })
  m4_list

end module
