m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_blosc},0,-1))

module test_blosc_m

  use iso_fortran_env, only: real64
  use, intrinsic :: iso_c_binding

  use string_m,    only: string, new_string
  use test_case_m, only: test_case
  use blosc_m
  use util_m

  implicit none

  private
  public test_blosc

  integer, parameter :: SIZE = 1000
  integer, parameter :: MAXSIZE = 1002

contains

  subroutine test_blosc()
    type(test_case)           :: tc
    integer(kind=c_size_t)    :: isize, csize, dsize, osize, i
    real(kind=real64), target :: data(SIZE), data_out(MAXSIZE), data_dest(SIZE)
    
    call tc%init("blosc")
    
    isize = sizeof(data)
    osize = sizeof(data_out)
    do i = 1, SIZE
      data(i) = 30.0 + i * 1e-2
    end do

    call blosc_init()

    csize = blosc_compress(int(9, kind=c_int), BLOSC_BITSHUFFLE, int(8, kind=c_size_t), isize, c_loc(data), c_loc(data_out), osize)

    call tc%assert(csize /= 0, "Buffer is incompressible")
    call tc%assert(csize > 0, "Compression error: " // int2str(int(csize, kind=4)))

    dsize = blosc_decompress(c_loc(data_out), c_loc(data_dest), isize)
    call tc%assert(dsize > 0 , "Decompression error: " // int2str(int(dsize, kind=4)))

    call blosc_destroy()

    call tc%assert(all(data == data_dest), "Decompressed data not equal to input data.")

    call tc%finish()
  end subroutine

end module

m4_divert(0)
