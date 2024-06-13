program test_galene

  use galene_m
  use normalization_m, only: init_normconst

  implicit none

  type(gal_file)           :: gat, sil
  type(gal_block), pointer :: b

  call init_normconst(300.0)

  call gat%init("GATE.OSV")
  call sil%init("SILIZ.GEO")

  b => sil%get_block("norm fac")
  print *, b%rdata(1)

  b => sil%get_block("x-coord")
  print *, b%rdata
end program
