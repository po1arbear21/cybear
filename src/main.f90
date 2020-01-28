#include "util/assert.f90.inc"

program main
  use error_m
  use matrix_m
  use input_m
  use arpack_m
  use dual_m
  implicit none

  type(input_file)          :: f
  integer                   :: i, Nx
  integer,      allocatable :: sid(:)
  real,         allocatable :: y(:)
  real                      :: y0, T
  character(:), allocatable :: ctname
  type(dual)                    :: dl

  call f%init("run/example.inp")

  call f%get("", "temperature", T, normalize = .false.)
  call init_normconst(T)

  call f%get("grid", "Nx", Nx)
  print *, Nx
  call f%get("grid", "y", y)
  print *, y(1)
  print *, y(2)

  call f%get_sections("region", sid)
  do i = 1, size(sid)
    call f%get(sid(i), "y0", y0)
    print "(1A, 1I0)", "region ", i
    print *, "y0 = ", denorm(y0, "um")
  end do

  call f%get_sections("contact", sid)
  do i = 1, size(sid)
    call f%get(sid(i), "name", ctname)
    print "(1A, 1I0)", "contact ", i
    print *, "name = ", ctname
  end do

  call dl%init(3, -1.0, i = 1)
  dl = abs(dl)
  dl = cos(dl)

  print *, dl%x
  print *, dl%dx

  print *, abs(dl%x)

  dl = dl ** dl
  print *, dl%x
  print *, dl%dx

end program main
