m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_blosc},0,-1))

module blosc_m

  use, intrinsic :: iso_c_binding

  implicit none

  private

  integer(kind=c_int), parameter, public :: BLOSC_MIN_HEADER_LENGTH = 16
    !! Minimum header length
  integer(kind=c_int), parameter, public :: BLOSC_MAX_OVERHEAD = BLOSC_MIN_HEADER_LENGTH
    !! The maximum overhead during compression in bytes. This equals to BLOSC_MIN_HEADER_LENGTH now, but can be higher in future implementations
  
  integer(kind=c_int), parameter, public :: BLOSC_NOSHUFFLE  =  0 
    !! no shuffle
  integer(kind=c_int), parameter, public :: BLOSC_SHUFFLE    =  1 
    !! byte-wise shuffle
  integer(kind=c_int), parameter, public :: BLOSC_BITSHUFFLE =  2
    !! bit-wise shuffle
  
  character(*), parameter, public :: BLOSC_BLOSCLZ_COMPNAME = "blosclz"
  character(*), parameter, public :: BLOSC_LZ4_COMPNAME     = "lz4"
  character(*), parameter, public :: BLOSC_LZ4HC_COMPNAME   = "lz4hc"
  character(*), parameter, public :: BLOSC_SNAPPY_COMPNAME  = "snappy"
  character(*), parameter, public :: BLOSC_ZLIB_COMPNAME    = "zlib"
  character(*), parameter, public :: BLOSC_ZSTD_COMPNAME    = "zstd"

  public :: blosc_init
  public :: blosc_set_compressor
  public :: blosc_compress
  public :: blosc_decompress
  public :: blosc_destroy

  interface
  
  subroutine blosc_init() bind(c, name="blosc_init")
    implicit none
  end subroutine

  function blosc_set_compressor(compname) bind(c, name="blosc_set_compressor")
    import :: c_char, c_int
    implicit none

    character(*, kind=c_char) :: compname

    integer(kind=c_int)       :: blosc_set_compressor
  end function

  function blosc_compress(clevel, doshuffle, typesize, nbytes, src, dest, destsize) bind(c, name="blosc_compress")
    import :: c_int, c_size_t, c_ptr
    implicit none

    integer(kind=c_int),    value :: clevel
    integer(kind=c_int),    value :: doshuffle
    integer(kind=c_size_t), value :: typesize
    integer(kind=c_size_t), value :: nbytes
    type(c_ptr),            value :: src
    type(c_ptr),            value :: dest
    integer(kind=c_size_t), value :: destsize

    integer(kind=c_int)    :: blosc_compress
  end function

  function blosc_decompress(src, dest, destsize) bind(c, name="blosc_decompress")
    import :: c_int, c_size_t, c_ptr
    implicit none

    type(c_ptr),            value :: src
    type(c_ptr),            value :: dest
    integer(kind=c_size_t), value :: destsize

    integer(kind=c_int)    :: blosc_decompress
  end function

  subroutine blosc_destroy() bind(c, name="blosc_destroy")
    implicit none
  end subroutine

  end interface

end module

m4_divert(0)