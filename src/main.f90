#include "util/assert.f90.inc"

program main
  use error_m
  use matrix_m
  use input_m
  !use arpack_m
  implicit none

  type(input_file)              :: f
  integer                       :: i, Nx
  integer,          allocatable :: sid(:)
  real,             allocatable :: y(:)
  real                          :: y0, T
  character(len=:), allocatable :: ctname

  call f%init("run/example.inp")

  call f%get("", "temperature", T, normalize = .false.)
  call init_normconst(T)

  call f%get("grid", "Nx", Nx)
  write(*,*) Nx
  call f%get("grid", "y", y)
  write(*,*) y(1)
  write(*,*) y(2)

  call f%get_sections("region", sid)
  do i = 1, size(sid)
    call f%get(sid(i), "y0", y0)
    write(*,"(1A, 1I0)") "region ", i
    write(*,*) "y0 = ", denorm(y0, "um")
  end do

  call f%get_sections("contact", sid)
  do i = 1, size(sid)
    call f%get(sid(i), "name", ctname)
    write(*,"(1A, 1I0)") "contact ", i
    write(*,*) "name = ", ctname
  end do

end program main
