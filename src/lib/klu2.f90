m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_klu2},0,-1))

module klu2_m

  use, intrinsic :: iso_c_binding

  use deque_m,      only: deque_int
  use error_m,      only: program_error
  use sparse_idx_m, only: SPARSE_IDX

  implicit none

  private
  public create_klu2_handle
  public destruct_klu2_handle
  public klu2_factorize
  public klu2_solve

  m4_define({m4_c_int},{c_int{}m4_idxsize{}_t})

  type, bind(c) :: klu2_handle
    type(c_ptr)       :: Common   = c_null_ptr
    type(c_ptr)       :: Symbolic = c_null_ptr
    type(c_ptr)       :: Numeric  = c_null_ptr
    integer(m4_c_int) :: n
    integer(m4_c_int) :: nnz
    type(c_ptr)       :: irow = c_null_ptr
    type(c_ptr)       :: row  = c_null_ptr
    type(c_ptr)       :: val  = c_null_ptr
  end type

  integer, parameter :: KLU2_NUM_HANDLES = 128
  type(klu2_handle)  :: klu2_handles(KLU2_NUM_HANDLES)
  type(deque_int)    :: klu2_free_handles

  interface
    ! extern "C" void klu2_init(klu2_handle *h)
  subroutine klu2_init_lib(h) bind(c, name="klu2_init")
    import klu2_handle, m4_c_int
    type(klu2_handle), intent(out) :: h
  end subroutine

  ! extern "C" void klu2_factorize(klu2_handle *h, const Int n, const Int nnz, const Int *ia, const Int *ja, const __float128 *a)
  subroutine klu2_factorize_lib(h, n, nnz, ia, ja, a) bind(c, name="klu2_factorize")
    import klu2_handle, m4_c_int
    type(klu2_handle),        intent(inout) :: h
    integer(m4_c_int), value, intent(in)    :: n
    integer(m4_c_int), value, intent(in)    :: nnz
    integer(m4_c_int),        intent(in)    :: ia(*)
    integer(m4_c_int),        intent(in)    :: ja(*)
    real(kind=16),            intent(in)    :: a(*)
  end subroutine

  ! extern "C" void klu2_solve(const klu2_handle *h, Int nrhs, const Real *b, Real *x)
  subroutine klu2_solve_lib(h, nrhs, b, x) bind(c, name="klu2_solve")
    import klu2_handle, m4_c_int
    type(klu2_handle),        intent(in)  :: h
    integer(m4_c_int), value, intent(in)  :: nrhs
    real(kind=16),            intent(in)  :: b(*)
    real(kind=16),            intent(out) :: x(*)
  end subroutine

  ! extern "C" void klu2_cleanup(klu2_handle *h)
  subroutine klu2_cleanup_lib(h) bind(c, name="klu2_cleanup")
    import klu2_handle, m4_c_int
    type(klu2_handle), intent(inout) :: h
  end subroutine
  end interface

contains

  function create_klu2_handle() result(h)
    !! create klu2 handle
    integer :: h
      !! return klu2 handle index

    integer :: i

    ! get free handle
    !$omp critical (omp_klu2_handles)
    if (.not. allocated(klu2_free_handles%d)) then
      call klu2_free_handles%init(KLU2_NUM_HANDLES, x = [(i, i=1, KLU2_NUM_HANDLES)])
    end if
    if (klu2_free_handles%n < 1) call program_error("No free KLU2 handles!")
    h = klu2_free_handles%front()
    call klu2_free_handles%pop_front()
    !$omp end critical (omp_klu2_handles)

    ! initialize
    call klu2_init_lib(klu2_handles(h))
  end function

  subroutine destruct_klu2_handle(h)
    !! destruct klu2 handle (free internal memory)
    integer, intent(inout) :: h
      !! klu2 handle index

    call klu2_cleanup_lib(klu2_handles(h))
    h = 0
  end subroutine

  subroutine klu2_factorize(h, ia, ja, a)
    !! factorize sparse matrix using KLU2
    integer,             intent(in) :: h
      !! klu2 handle index
    integer(SPARSE_IDX), intent(in) :: ia(:)
      !! column pointers
    integer,             intent(in) :: ja(:)
      !! columns
    real(kind=16),       intent(in) :: a(:)
      !! values

    integer(SPARSE_IDX) :: n, nnz

    m4_ifelse(m4_idxsize,m4_intsize,,{
    integer(SPARSE_IDX), allocatable :: ja_(:)

    allocate (ja_(size(ja)), source = int(ja, kind=SPARSE_IDX))
    })

    ! matrix size and number of non-zeros
    n   = size(ia) - 1
    nnz = size(a, kind = SPARSE_IDX)

    ! perform factorization
    call klu2_factorize_lib(klu2_handles(h), n, nnz, ia, ja{}m4_ifelse(m4_idxsize,m4_intsize,{},{_}), a)
  end subroutine

  subroutine klu2_solve(h, b, x)
    !! solve sparse system using KLU2
    integer,       intent(in)  :: h
      !! klu2 handle index
    real(kind=16), intent(in)  :: b(:)
      !! right-hand side(s)
    real(kind=16), intent(out) :: x(:)
      !! solution(s)

    integer(SPARSE_IDX) :: nrhs

    ! number of right-hand sides
    nrhs = size(b, kind=SPARSE_IDX) / klu2_handles(h)%n

    call klu2_solve_lib(klu2_handles(h), nrhs, b, x)
  end subroutine

end module

m4_divert(0)
