m4_include(macro.f90.inc)

module map_m

  use error_m, only: assert_failed
  use string_m

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
  m4_define({m4_X},{public map_string_$1, mapnode_string_$1})
  m4_list

  ! include type definitions
  m4_define({m4_X},{
    m4_define({T},string)
    m4_define({U},$1)
    m4_include(map_def.f90.inc)
  })
  m4_list

contains

  ! include type procedure implementations
  m4_define({m4_X},{
    m4_define({T},string)
    m4_define({U},$1)
    m4_include(map_imp.f90.inc)
  })
  m4_list

end module
