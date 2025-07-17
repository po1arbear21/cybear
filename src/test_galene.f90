program test_galene

  use galene_m
  use normalization_m, only: init_normconst, denorm

  implicit none

  type(gal_file)           :: gat, sil
  type(gal_block), pointer :: b

  call init_normconst(300.0)

  call gat%init("DD.OSV")
  call sil%init("DGMOS.GEO")

  b => sil%get_block("norm fac")
  print *, b%rdata(1)

  b => sil%get_block("x-coord")
  print *, b%rdata

  b => sil%get_block("y-coord")
  print *, denorm(b%rdata, "um")
  print *, size(b%rdata)

  ! b => gat%get_block("n-density")
  ! print *, b%rdata
end program
