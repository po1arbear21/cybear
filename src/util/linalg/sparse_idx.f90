module sparse_idx_m
!  use iso_fortran_env ! gives warnings when real size is 8
  implicit none

  private
  public :: SPARSE_IDX

  ! 32 or 64 bit integers for sparse matrix indices
#ifdef IDXSIZE32
  integer, parameter :: SPARSE_IDX = 4 ! int32
#endif
#ifdef IDXSIZE64
  integer, parameter :: SPARSE_IDX = 8 ! int64
#endif
end module