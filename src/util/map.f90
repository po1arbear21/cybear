m4_include(macro.f90.inc)

module map_m

  use error_m,         only: assert_failed
  use iso_fortran_env, only: int64
  use string_m

  implicit none

  private

  ! list of built-in key types
  m4_define({m4_keylist},{
    m4_X(int)
    m4_X(string)
  })

  ! list of built-in value types
  m4_define({m4_valuelist},{
    m4_Y($1,int)
    m4_Y($1,int64)
    m4_Y($1,log)
    m4_Y($1,real)
    m4_Y($1,cmplx)
    m4_Y($1,string)
  })

  ! combine key-list with value-list to create full key-value list
  m4_define({m4_X},{m4_valuelist($1)})
  m4_define({m4_list},m4_keylist)

  ! make defined types public
  m4_define({m4_Y},{public map_$1_$2, mapnode_$1_$2})
  m4_list

  ! include type definitions
  m4_define({m4_Y},{
    m4_define({T},$1)
    m4_define({U},$2)
    m4_include(map_def.f90.inc)
  })
  m4_list

contains

  ! include type procedure implementations
  m4_define({m4_Y},{
    m4_define({T},$1)
    m4_define({U},$2)
    m4_include(map_imp.f90.inc)
  })
  m4_list

end module
