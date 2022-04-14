m4_include(macro.f90.inc)

module array_m

  use string_m, only: string

  implicit none

  private

  ! maximum supported array dimension
  m4_define({m4_max_dim},{8})

  ! combine dimensions 0 to max_dim with built-in type list to get a full list
  m4_define({m4_list_help},{
    m4_ifelse($1,0,,{m4_list_help(m4_decr($1),$2)})
    m4_X($1,$2)
  })
  m4_define({m4_list},{
    m4_list_help(m4_max_dim,int)
    m4_list_help(m4_max_dim,log)
    m4_list_help(m4_max_dim,real)
    m4_list_help(m4_max_dim,cmplx)
    m4_list_help(m4_max_dim,string)
  })

  ! make defined types public
  m4_define({m4_X},{
    m4_ifelse($1,0,{
      public :: ptr_$2
    },{
      public :: array$1_$2, ptr$1_$2
    })
  })
  m4_list

  ! include type definitions (N: array dimension, T: short typename)
  m4_define({m4_X},{
    m4_define({N},$1)
    m4_define({T},$2)
    m4_include(array_def.f90.inc)
  })
  m4_list

end module
