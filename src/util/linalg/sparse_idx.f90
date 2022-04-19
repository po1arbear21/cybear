m4_include(../macro.f90.inc)

module sparse_idx_m

  use iso_fortran_env, only: int32, int64

  implicit none

  private
  public SPARSE_IDX

  ! 32 or 64 bit integers for sparse matrix indices
  integer, parameter :: SPARSE_IDX = int{}m4_idxsize

end module
