#ifdef USE_MUMPS

module mumps_m
  use error_m
  use sparse_idx_m
  use vector_m
  implicit none

  private
  public :: create_mumps_handle_r, create_mumps_handle_c
  public :: destruct_mumps_handle_r, destruct_mumps_handle_c
  public :: mumps_factorize
  public :: mumps_solve

  include 'dmumps_struc.h'
  include 'zmumps_struc.h'

#define T dmumps_struc
#define TT type(dmumps_struc)
#include "../util/vector_def.f90.inc"
#define T zmumps_struc
#define TT type(zmumps_struc)
#include "../util/vector_def.f90.inc"

  type(vector_dmumps_struc) :: dmumps_handles
  type(vector_zmumps_struc) :: zmumps_handles
  type(vector_int)          :: free_dmumps_handles
  type(vector_int)          :: free_zmumps_handles

  interface
    subroutine dmumps(d)
      import dmumps_struc
      type(dmumps_struc), intent(inout) :: d
    end subroutine
    subroutine zmumps(z)
      import zmumps_struc
      type(zmumps_struc), intent(inout) :: z
    end subroutine
  end interface

  interface mumps_factorize
    module procedure :: mumps_factorize_r
    module procedure :: mumps_factorize_c
  end interface

  interface mumps_solve
    module procedure :: mumps_solve_r
    module procedure :: mumps_solve_c
  end interface

contains

#define T dmumps_struc
#define TT type(dmumps_struc)
#include "../util/vector_imp.f90.inc"
#define T zmumps_struc
#define TT type(zmumps_struc)
#include "../util/vector_imp.f90.inc"

  subroutine exec_dmumps(m, job)
    !! set job and execute dmumps(m); handle potential error
    type(dmumps_struc), intent(inout) :: m
      !! dmumps structure
    integer,            intent(in)    :: job
      !! job to execute

    m%JOB = job
    call dmumps(m)
    if (m%INFOG(1) < 0) then
      print "(A,I0)"     , "JOB   = ", m%JOB
      print "(A,I0,A,I0)", "INFOG = ", m%INFOG(1), ", ", m%INFOG(2)
      call program_error("DMUMPS error")
    end if
  end subroutine

  subroutine exec_zmumps(m, job)
    !! set job and execute zmumps(m); handle potential error
    type(zmumps_struc), intent(inout) :: m
      !! zmumps structure
    integer,            intent(in)    :: job
      !! job to execute

    m%JOB = job
    call zmumps(m)
    if (m%INFOG(1) < 0) then
      print "(A,I0)"     , "JOB   = ", m%JOB
      print "(A,I0,A,I0)", "INFOG = ", m%INFOG(1), ", ", m%INFOG(2)
      call program_error("DMUMPS error")
    end if
  end subroutine

  function create_mumps_handle_r() result(h)
    !! create real mumps handle
    integer :: h
      !! return real mumps handle (index)

    if (.not. allocated(dmumps_handles%d)) then
      call dmumps_handles%init(0, c = 4)
    end if

    if (free_dmumps_handles%n > 0) then
      h = free_dmumps_handles%d(free_dmumps_handles%n)
      call free_dmumps_handles%resize(free_dmumps_handles%n-1)
    else
      block
        type(dmumps_struc) :: m
        call dmumps_handles%push(m)
      end block
      h = dmumps_handles%n
    end if

    associate (m => dmumps_handles%d(h))
      m%SYM  = 0              ! unsymmetric
      m%PAR  = 1              ! run calculation on host cpu
      call exec_dmumps(m, -1) ! initialize

      m%ICNTL(4) = 1
    end associate
  end function

  function create_mumps_handle_c() result(h)
    !! create complex mumps handle
    integer :: h
      !! return complex mumps handle (index)

    if (.not. allocated(zmumps_handles%d)) then
      call zmumps_handles%init(0, c = 4)
    end if

    if (free_zmumps_handles%n > 0) then
      h = free_zmumps_handles%d(free_zmumps_handles%n)
      call free_zmumps_handles%resize(free_zmumps_handles%n-1)
    else
      block
        type(zmumps_struc) :: m
        call zmumps_handles%push(m)
      end block
      h = zmumps_handles%n
    end if

    associate (m => zmumps_handles%d(h))
      m%SYM  = 0              ! unsymmetric
      m%PAR  = 1              ! run calculation on host cpu
      call exec_zmumps(m, -1) ! initialize

      m%ICNTL(4) = 1
    end associate
  end function

  subroutine destruct_mumps_handle_r(h)
    !! destruct real mumps hanlde
    integer, intent(inout) :: h
      !! real mumps handle (index)

    associate (m => dmumps_handles%d(h))
      ! deallocate user data
      if (associated(m%IRN)) deallocate (m%IRN)
      if (associated(m%JCN)) deallocate (m%JCN)
      if (associated(m%A  )) deallocate (m%A  )
      if (associated(m%RHS)) deallocate (m%RHS)

      ! deallocate internal data
      call exec_dmumps(m, -2)
    end associate

    if (.not. allocated(free_dmumps_handles%d)) then
      call free_dmumps_handles%init(0, c = 4)
    end if
    call free_dmumps_handles%push(h)
    h = 0
  end subroutine

  subroutine destruct_mumps_handle_c(h)
    !! destruct complex mumps hanlde
    integer, intent(inout) :: h
      !! complex mumps handle (index)

    associate (m => zmumps_handles%d(h))
      ! deallocate user data
      if (associated(m%IRN)) deallocate (m%IRN)
      if (associated(m%JCN)) deallocate (m%JCN)
      if (associated(m%A  )) deallocate (m%A  )
      if (associated(m%RHS)) deallocate (m%RHS)

      ! deallocate internal data
      call exec_zmumps(m, -2)
    end associate

    if (.not. allocated(free_zmumps_handles%d)) then
      call free_zmumps_handles%init(0, c = 4)
    end if
    call free_zmumps_handles%push(h)
    h = 0
  end subroutine

  subroutine mumps_factorize_r(h, ia, ja, a)
    !! factorize real sparse matrix using MUMPS
    integer,             intent(in) :: h
      !! real mumps handle (index)
    integer(SPARSE_IDX), intent(in) :: ia(:)
    integer,             intent(in) :: ja(:)
    real,                intent(in) :: a(:)

    integer             :: i, n
    integer(SPARSE_IDX) :: j, nnz

    associate (m => dmumps_handles%d(h))
      ! matrix size and number of zeros
      n   = size(ia) - 1
      nnz = size(a, kind = SPARSE_IDX)

      ! specify mumps matrix
      m%n   = n
      m%nnz = nnz
      allocate (m%IRN(nnz), m%JCN(nnz), m%A(nnz))
      do i = 1, n
        do j = ia(i), ia(i+1)-1
          m%IRN(j) = i     ! row number
          m%JCN(j) = ja(j) ! col number
          m%A(  j) = a(j)  ! value
        end do
      end do

      ! perform analysis and factorization (JOB=4 is the combination of JOB=1 and JOB=2)
      call exec_dmumps(m, 4)
    end associate
  end subroutine

  subroutine mumps_factorize_c(h, ia, ja, a)
    !! factorize complex sparse matrix using MUMPS
    integer,             intent(in) :: h
      !! complex mumps handle (index)
    integer(SPARSE_IDX), intent(in) :: ia(:)
    integer,             intent(in) :: ja(:)
    complex,             intent(in) :: a(:)

    integer             :: i, n
    integer(SPARSE_IDX) :: j, nnz

    associate (m => zmumps_handles%d(h))
      ! matrix size and number of zeros
      n   = size(ia) - 1
      nnz = size(a, kind = SPARSE_IDX)

      ! specify mumps matrix
      m%N   = n
      m%NNZ = nnz
      allocate (m%IRN(nnz), m%JCN(nnz), m%A(nnz))
      do i = 1, n
        do j = ia(i), ia(i+1)-1
          m%IRN(j) = i     ! row number
          m%JCN(j) = ja(j) ! col number
          m%A(  j) = a(j)  ! value
        end do
      end do

      ! perform analysis and factorization (JOB=4 is the combination of JOB=1 and JOB=2)
      call exec_zmumps(m, 4)
    end associate
  end subroutine

  subroutine mumps_solve_r(h, b, x)
    !! solve real system with MUMPS
    integer, intent(in)  :: h
    real,    intent(in)  :: b(:)
    real,    intent(out) :: x(:)

    associate (m => dmumps_handles%d(h))
      ! set rhs
      if (.not. associated(m%RHS)) allocate (m%RHS(size(b)))
      m%RHS = b

      ! number of right-hand sides
      m%NRHS = size(b)/m%N

      ! solve
      call exec_dmumps(m, 3)

      ! extract solution
      x = m%RHS
    end associate
  end subroutine

  subroutine mumps_solve_c(h, b, x)
    !! solve complex system with MUMPS
    integer, intent(in)  :: h
    complex, intent(in)  :: b(:)
    complex, intent(out) :: x(:)

    associate (m => zmumps_handles%d(h))
      ! set rhs
      if (.not. associated(m%RHS)) allocate (m%RHS(size(b)))
      m%RHS = b

      ! number of right-hand sides
      m%NRHS = size(b)/m%N

      ! solve
      call exec_zmumps(m, 3)

      ! extract solution
      x = m%RHS
    end associate
  end subroutine

end module

#endif