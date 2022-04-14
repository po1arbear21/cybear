m4_include(../macro.f90.inc)

submodule (matrix_m) band_m

contains

  m4_define({m4_list},{
    m4_X(real)
    m4_X(cmplx)
  })
  m4_define({m4_X},{
    m4_define({T},{$1})
    m4_include(band_imp.f90.inc)
    m4_undefine({T})
  })
  m4_list

end submodule
