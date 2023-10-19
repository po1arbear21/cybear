m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mumps},0,-1))

module mumps_m

  use deque_m,      only: deque_int
  use error_m,      only: program_error
  use sparse_idx_m, only: sparse_idx

  implicit none

  ! FIXME: not thread-safe

  private
  public create_mumps_handle_r, create_mumps_handle_c
  public destruct_mumps_handle_r, destruct_mumps_handle_c
  public mumps_factorize
  public mumps_solve

  include 'dmumps_struc.h'
  include 'zmumps_struc.h'

  integer, parameter :: MUMPS_NUM_HANDLES = 128
  type(dmumps_struc) :: dmumps_handles(MUMPS_NUM_HANDLES)
  type(zmumps_struc) :: zmumps_handles(MUMPS_NUM_HANDLES)
  type(deque_int)    :: dmumps_free_handles
  type(deque_int)    :: zmumps_free_handles

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

    integer :: i

    ! get free handle
    !$omp critical (omp_dmumps_handles)
    if (.not. allocated(dmumps_free_handles%d)) then
      call dmumps_free_handles%init(MUMPS_NUM_HANDLES, x = [(i, i=1, MUMPS_NUM_HANDLES)])
    end if
    if (dmumps_free_handles%n < 1) call program_error("No free dmumps handles!")
    h = dmumps_free_handles%front()
    call dmumps_free_handles%pop_front()
    !$omp end critical (omp_dmumps_handles)

    associate (m => dmumps_handles(h))
      m%SYM  = 0              ! unsymmetric
      m%PAR  = 1              ! run calculation on host cpu
      call exec_dmumps(m, -1) ! initialize

      m%ICNTL( 4) = 1
      m%ICNTL(14) = 100
    end associate
  end function

  function create_mumps_handle_c() result(h)
    !! create complex mumps handle
    integer :: h
      !! return complex mumps handle (index)

    integer :: i

    ! get free handle
    !$omp critical (omp_zmumps_handles)
    if (.not. allocated(zmumps_free_handles%d)) then
      call zmumps_free_handles%init(MUMPS_NUM_HANDLES, x = [(i, i=1, MUMPS_NUM_HANDLES)])
    end if
    if (zmumps_free_handles%n < 1) call program_error("No free zmumps handles!")
    h = zmumps_free_handles%front()
    call zmumps_free_handles%pop_front()
    !$omp end critical (omp_zmumps_handles)

    associate (m => zmumps_handles(h))
      m%SYM  = 0              ! unsymmetric
      m%PAR  = 1              ! run calculation on host cpu
      call exec_zmumps(m, -1) ! initialize

      m%ICNTL( 4) = 1
      m%ICNTL(14) = 100       ! increase memory in steps of 100%
    end associate
  end function

  subroutine destruct_mumps_handle_r(h)
    !! destruct real mumps hanlde
    integer, intent(inout) :: h
      !! real mumps handle (index)

    associate (m => dmumps_handles(h))
      ! deallocate user data
      if (associated(m%IRN)) deallocate (m%IRN)
      if (associated(m%JCN)) deallocate (m%JCN)
      if (associated(m%A  )) deallocate (m%A  )
      if (associated(m%RHS)) deallocate (m%RHS)

      ! deallocate internal data
      call exec_dmumps(m, -2)
    end associate

    !$omp critical (omp_dmumps_handles)
    call dmumps_free_handles%push_back(h)
    !$omp end critical (omp_dmumps_handles)

    h = 0
  end subroutine

  subroutine destruct_mumps_handle_c(h)
    !! destruct complex mumps hanlde
    integer, intent(inout) :: h
      !! complex mumps handle (index)

    associate (m => zmumps_handles(h))
      ! deallocate user data
      if (associated(m%IRN)) deallocate (m%IRN)
      if (associated(m%JCN)) deallocate (m%JCN)
      if (associated(m%A  )) deallocate (m%A  )
      if (associated(m%RHS)) deallocate (m%RHS)

      ! deallocate internal data
      call exec_zmumps(m, -2)
    end associate

    !$omp critical (omp_zmumps_handles)
    call zmumps_free_handles%push_back(h)
    !$omp end critical (omp_zmumps_handles)

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

    associate (m => dmumps_handles(h))
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

    associate (m => zmumps_handles(h))
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

    associate (m => dmumps_handles(h))
      ! set rhs
      if (.not. associated(m%RHS)) allocate (m%RHS(size(b)))
      m%RHS = b

      ! number of right-hand sides
      m%NRHS = size(b)/m%N
      m%LRHS = m%N

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

    associate (m => zmumps_handles(h))
      ! set rhs
      if (.not. associated(m%RHS)) allocate (m%RHS(size(b)))
      m%RHS = b

      ! number of right-hand sides
      m%NRHS = size(b)/m%N
      m%LRHS = m%N

      ! solve
      call exec_zmumps(m, 3)

      ! extract solution
      x = m%RHS
    end associate
  end subroutine

end module

m4_divert(0)
