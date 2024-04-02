m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_zlib},0,-1))

module test_zlib_m

  use, intrinsic :: iso_c_binding

  use string_m,    only: string, new_string
  use test_case_m, only: test_case
  use zlib_m

  implicit none

  private
  public test_zlib

contains

  subroutine test_zlib()
    type(test_case)        :: tc
    type(z_stream)         :: defstream, infstream
    character(100), target :: a, b, c
    integer                :: rc

    call tc%init("zlib")

    a = "Lorem ipsum dolor lorem ipsum dolor lorem ipsum dolor lorem ipsum dolor lorem ipsum dolor lorem ipsu"
    defstream%avail_in  = len_trim(a)
    defstream%next_in   = c_loc(a)
    defstream%avail_out = len(b)
    defstream%next_out  = c_loc(b)

    rc = deflate_init(defstream, Z_BEST_COMPRESSION)
    rc = deflate(defstream, Z_FINISH)
    rc = deflate_end(defstream)

    infstream%avail_in  = defstream%total_out
    infstream%next_in   = c_loc(b)
    infstream%avail_out = len(c)
    infstream%next_out  = c_loc(c)

    rc = inflate_init(infstream)
    rc = inflate(infstream, Z_NO_FLUSH)
    rc = inflate_end(infstream)

    call tc%assert_eq(new_string(trim(a)), new_string(trim(c)), "deflate and inflate")

    call tc%finish()
  end subroutine

end module

m4_divert(0)
