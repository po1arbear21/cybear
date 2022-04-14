m4_include(macro.f90.inc)

module hashmap_m

  use error_m,  only: assert_failed
  use string_m, only: string
  use util_m,   only: hash
  use vector_m, only: vector_int, vector_log, vector_string, vector_real, vector_cmplx

  implicit none

  private

  ! list of built-in types
  m4_define({m4_list},{
    m4_X(int)
    m4_X(log)
    m4_X(real)
    m4_X(cmplx)
    m4_X(string)
  })

  ! make defined types public
  m4_define({m4_X},{public hashmap_$1})
  m4_list

  ! include type definitions
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(hashmap_def.f90.inc)
  })
  m4_list

contains

  ! include type procedure implementations
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(hashmap_imp.f90.inc)
  })
  m4_list

end module
