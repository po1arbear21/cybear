#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_PTR_TYPE        CONCAT(band_ptr,T)
#define BLOCK_TYPE           CONCAT(block,T)
#define BLOCK_PTR_TYPE       CONCAT(block_ptr,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define DENSE_PTR_TYPE       CONCAT(dense_ptr,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define HESSENBERG_PTR_TYPE  CONCAT(hessenberg_ptr,T)
#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_PTR_TYPE      CONCAT(matrix_ptr,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPARSE_PTR_TYPE      CONCAT(sparse_ptr,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define TRIANG_TYPE          CONCAT(triang,T)
#define TRIANG_PTR_TYPE      CONCAT(triang_ptr,T)

#define SPARSE_INIT          CONCAT3(sparse,T,init)
#define SPARSE_INIT_SOLVER   CONCAT3(sparse,T,init_solver)
#define SPARSE_DESTRUCT      CONCAT3(sparse,T,destruct)
#define SPARSE_RESET         CONCAT3(sparse,T,reset)
#define SPARSE_SCALE         CONCAT3(sparse,T,scale)
#define SPARSE_MUL_VEC       CONCAT3(sparse,T,mul_vec)
#define SPARSE_MUL_VEC_SLICE CONCAT3(sparse,T,mul_vec_slice)
#define SPARSE_MUL_MAT       CONCAT3(sparse,T,mul_mat)
#define SPARSE_FACTORIZE     CONCAT3(sparse,T,factorize)
#define SPARSE_SOLVE_VEC     CONCAT3(sparse,T,solve_vec)
#define SPARSE_SOLVE_MAT     CONCAT3(sparse,T,solve_mat)
#define SPARSE_TRANSPOSE     CONCAT3(sparse,T,transpose)
#define SPARSE_TRANSPOSE2    CONCAT3(sparse,T,transpose2)
#define SPARSE_MUL_SPARSE    CONCAT3(sparse,T,mul_sparse)
#define SPARSE_INPUT         CONCAT3(sparse,T,input)
#define SPARSE_OUTPUT        CONCAT3(sparse,T,output)
#define SPARSE_NNZ           CONCAT3(sparse,T,nnz)
#define SPARSE_IS_EMPTY      CONCAT3(sparse,T,is_empty)
#define SPARSE_EYE           CONCAT3(sparse,eye,T)
#define SPARSE_ZERO          CONCAT3(sparse,zero,T)

#define SPBUILD_INIT         CONCAT3(spbuild,T,init)
#define SPBUILD_DESTRUCT     CONCAT3(spbuild,T,destruct)
#define SPBUILD_RESET        CONCAT3(spbuild,T,reset)
#define SPBUILD_RESET_ROW    CONCAT3(spbuild,T,reset_row)
#define SPBUILD_SET          CONCAT3(spbuild,T,set)
#define SPBUILD_SET_ROW      CONCAT3(spbuild,T,set_row)
#define SPBUILD_ADD          CONCAT3(spbuild,T,add)
#define SPBUILD_SORT         CONCAT3(spbuild,T,sort)
#define SPBUILD_LOAD         CONCAT3(spbuild,T,load)
#define SPBUILD_SAVE         CONCAT3(spbuild,T,save)

module subroutine SPARSE_INIT(this, nrows, ncols)
  !! Initialize sparse matrix.

  class(SPARSE_TYPE), intent(out) :: this
    !! Sparse matrix
  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional: Number of columns (default: ncols = nrows)

  ! init base
  call this%init("", nrows, ncols = ncols)
  this%solver = default_solver
end subroutine

module subroutine SPARSE_DESTRUCT(this)
  !! Destruct sparse matrix by deallocating memory and deleting factorization.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  call this%reset()
end subroutine

module subroutine SPARSE_RESET(this, only_factorization)
  !! Reset sparse matrix to empty state (optionally only the factorization).

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  logical, optional,  intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  ! deallocate CSR data
  if (.not. only_fact) then
    if (allocated(this%ia)) deallocate (this%ia)
    if (allocated(this%ja)) deallocate (this%ja)
    if (allocated(this%a )) deallocate (this%a )
  end if

  ! destruct factorization
  if (this%solver_handle > 0) then
    select case (this%solver)
      case (SOLVER_PARDISO)
        call destruct_pardiso_handle(this%solver_handle)

#ifdef USE_MUMPS
      case (SOLVER_MUMPS)
#ifdef TCMPLX
        call destruct_mumps_handle_c(this%solver_handle)
#else
        call destruct_mumps_handle_r(this%solver_handle)
#endif
#endif

#ifdef USE_ILUPACK
      case (SOLVER_ILUPACK)
        call destruct_ilupack_handle(this%solver_handle)
#endif
    end select
  end if

  this%factorized = .false.
end subroutine

module subroutine SPARSE_SCALE(this, fact)
  !! Scale sparse matrix by a factor.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  TT,                 intent(in)    :: fact
    !! Scaling factor

  ! h == 0 is default value for unfactorized
  ASSERT(this%solver_handle == 0)

  ! scale values
  if (.not. this%is_empty()) this%a = fact * this%a
end subroutine

module subroutine SPARSE_MUL_VEC_SLICE(this, x, y, fact_y, i0, i1)
  !! Multiply sparse matrix by a raw vector where only some rows of y are computed.
  !! Performs operation: y = this * x + fact_y * y

  class(SPARSE_TYPE), intent(in)    :: this
    !! Sparse matrix
  TT,                 intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                 intent(inout) :: y(:)
    !! Store result in y
  TT,      optional,  intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  integer, optional,  intent(in)    :: i0
    !! first row index. default: 1
  integer, optional,  intent(in)    :: i1
    !! last row index. default: this%nrows

  integer             :: i0_, i1_, i
  integer(SPARSE_IDX) :: j

  i0_ = 1
  if (present(i0)) i0_ = i0
  i1_ = this%nrows
  if (present(i1)) i1_ = i1

  ASSERT(size(y) == i1_-i0_+1)
  ASSERT((i0_ >= 1) .and. (i0_ <= this%nrows))
  ASSERT((i1_ >= 1) .and. (i1_ <= this%nrows))
  ASSERT(i0_ <= i1_)

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! matrix vector multiplication
  !$omp parallel do default(none) private(i, j) shared(this, x, y, i0_, i1_) schedule(dynamic)
  do i = 1, i1_-i0_+1
    do j = this%ia(i+i0_-1), this%ia(i+i0_)-1
      y(i) = y(i) + this%a(j) * x(this%ja(j))
    end do
  end do
  !$omp end parallel do
end subroutine

module subroutine SPARSE_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(SPARSE_TYPE),  intent(in)    :: this
    !! Sparse matrix
  TT,                  intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                  intent(inout) :: y(:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transpose flag, must be 'N', 'T' or 'C' (default: 'N')
    !! 'N':   y =                     this  * x + fact_y * y
    !! 'T':   y =           transpose(this) * x + fact_y * y
    !! 'C':   y = conjugate_transpose(this) * x + fact_y * y

  character           :: trans_
  integer             :: i
  integer(SPARSE_IDX) :: j

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! transpose flag
  trans_ = 'N'
  if (present(trans)) trans_ = trans
#ifdef TCMPLX
  ASSERT((trans_ == 'N') .or. (trans_ == 'T') .or. (trans_ == 'C'))
#else
  ASSERT((trans_ == 'N') .or. (trans_ == 'T'))
#endif

  if (trans_ == 'N') then
    ! matrix vector multiplication
    !$omp parallel do default(none) private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        y(i) = y(i) + this%a(j) * x(this%ja(j))
      end do
    end do
    !$omp end parallel do

  else
    ! fixme y_tmp is unnecessary! use directly y whereever there is y_tmp now. (in reduction and assignments). yields catastrophic compiler error for ifort ...
    block
      TT, allocatable :: y_tmp(:)

      allocate (y_tmp(size(y)))
      y_tmp = 0

#ifdef TCMPLX
      if (trans_ == 'C') then
        ! matrix vector multiplication with complex conjugated transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j)) = y_tmp(this%ja(j)) + conjg(this%a(j)) * x(i)
          end do
        end do
        !$omp end parallel do

      else
#endif
        ! matrix vector multiplication with transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j)) = y_tmp(this%ja(j)) + this%a(j) * x(i)
          end do
        end do
        !$omp end parallel do
#ifdef TCMPLX
      end if
#endif
      y = y + y_tmp
    end block
  end if
end subroutine

module subroutine SPARSE_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(SPARSE_TYPE),  intent(in)    :: this
    !! Sparse matrix
  TT,                  intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                  intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transpose flag for this matrix, must be 'N', 'T' or 'C' (default: 'N')
    !! 'N':   y =                     this  * x + fact_y * y
    !! 'T':   y =           transpose(this) * x + fact_y * y
    !! 'C':   y = conjugate_transpose(this) * x + fact_y * y

  character           :: trans_
  integer             :: i
  integer(SPARSE_IDX) :: j

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! transpose flag
  trans_ = 'N'
  if (present(trans)) trans_ = trans
#ifdef TCMPLX
  ASSERT((trans_ == 'N') .or. (trans_ == 'T') .or. (trans_ == 'C'))
#else
  ASSERT((trans_ == 'N') .or. (trans_ == 'T'))
#endif

  if (trans_ == 'N') then
    ! matrix matrix multiplication
    !$omp parallel do default(none) private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        y(i,:) = y(i,:) + this%a(j) * x(this%ja(j),:)
      end do
    end do
    !$omp end parallel do

  else
    ! fixme y_tmp is unnecessary! use directly y whereever there is y_tmp now. (in reduction and assignments). yields catastrophic compiler error for ifort ...
    block
      TT, allocatable :: y_tmp(:,:)

      allocate (y_tmp(size(y, dim=1), size(y, dim=2)))
      y_tmp = 0

#ifdef TCMPLX
      if (trans_ == 'C') then
        ! matrix matrix multiplication with complex conjugated transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j),:) = y_tmp(this%ja(j),:) + conjg(this%a(j)) * x(i,:)
          end do
        end do
        !$omp end parallel do

      else
#endif
        ! matrix matrix multiplication with transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j),:) = y_tmp(this%ja(j),:) + this%a(j) * x(i,:)
          end do
        end do
        !$omp end parallel do
#ifdef TCMPLX
      end if
#endif
      y = y + y_tmp
    end block
  end if
end subroutine

module subroutine SPARSE_INIT_SOLVER(this)
  !! creates solver handle
  !!
  !! this routine exists s.t. one can init ilupack, then get its handler and change parameters.

  class(SPARSE_TYPE), intent(inout) :: this

  ! make sure matrix is unfactorized and not empty
  ASSERT(this%solver_handle == 0)
  ASSERT(.not. this%is_empty())

  select case (this%solver)
    case (SOLVER_PARDISO)
#ifdef TCMPLX
      this%solver_handle = create_pardiso_handle(this%nrows, .true. )
#else
      this%solver_handle = create_pardiso_handle(this%nrows, .false.)
#endif

#ifdef USE_MUMPS
    case (SOLVER_MUMPS)
#ifdef TCMPLX
      this%solver_handle = create_mumps_handle_c()
#else
      this%solver_handle = create_mumps_handle_r()
#endif
#endif

#ifdef USE_ILUPACK
    case (SOLVER_ILUPACK)
      this%solver_handle = create_ilupack_handle(this%ia, this%ja, this%a)
#endif
  end select
end subroutine

module subroutine SPARSE_FACTORIZE(this)
  !! Factorize Matrix using PARDISO or ILUPACK.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  if (this%solver_handle == 0) call this%init_solver

  select case (this%solver)
    case (SOLVER_PARDISO)
      call pardiso_factorize(this%solver_handle, this%ia, this%ja, this%a)

#ifdef USE_MUMPS
    case (SOLVER_MUMPS)
      call mumps_factorize(this%solver_handle, this%ia, this%ja, this%a)
#endif

#ifdef USE_ILUPACK
    case (SOLVER_ILUPACK)
      call ilupack_factorize(this%solver_handle)
#endif
  end select

  this%factorized = .true.
end subroutine

module subroutine SPARSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using PARDISO.
  !!
  !! Call this%factorize() first!

  class(SPARSE_TYPE),  intent(in)  :: this
    !! Sparse matrix
  TT,                  intent(in)  :: rhs(:)
    !! Right hand side
  TT,                  intent(out) :: x(:)
    !! Solution vector
  character, optional, intent(in)  :: trans
    !! Not implemented for 'T' or 'C' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  ! matrix must be factorized before solving
  ASSERT(this%factorized)

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed sparse matrix not implemented!")
  end if
  if (present(iter_refine)) then
    if (.not. iter_refine) call program_error("Solving with sparse matrix without iterative refinement not implemented!")
  end if

  ! solve system with pardiso/ilupack
  select case (this%solver)
    case (SOLVER_PARDISO)
      call pardiso_solve(this%solver_handle, this%ia, this%ja, this%a, rhs, x)

#ifdef USE_MUMPS
    case (SOLVER_MUMPS)
      call mumps_solve(this%solver_handle, rhs, x)
#endif

#ifdef USE_ILUPACK
    case (SOLVER_ILUPACK)
      call ilupack_solve(this%solver_handle, rhs, x)
#endif
  end select
end subroutine

module subroutine SPARSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using PARDISO.
  !!
  !! Call this%factorize() first.

  class(SPARSE_TYPE),  intent(in)  :: this
    !! Sparse matrix
  TT,                  intent(in)  :: rhs(:,:)
    !! Right hand side
  TT,                  intent(out) :: x(:,:)
    !! Solution vector
  character, optional, intent(in)  :: trans
    !! Not implemented for 'T' or 'C' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  TT, allocatable :: rhs_(:), x_(:)

  ! matrix must be factorized before solving
  ASSERT(this%factorized)

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed sparse matrix not implemented!")
  end if
  if (present(iter_refine)) then
    if (.not. iter_refine) call program_error("Solving with sparse matrix without iterative refinement not implemented!")
  end if

  ! allocate space for temporary copy
  allocate (rhs_(size(rhs,1) * size(rhs,2)))
  allocate (  x_(size(  x,1) * size(  x,2)))

  ! solve with pardiso/ilupack, pass right-hand sides one after another in 1D array
  rhs_ = reshape(rhs, [size(rhs,1) * size(rhs,2)])
  select case (this%solver)
    case (SOLVER_PARDISO)
      call pardiso_solve(this%solver_handle, this%ia, this%ja, this%a, rhs_, x_)

#ifdef USE_MUMPS
    case (SOLVER_MUMPS)
      call mumps_solve(this%solver_handle, rhs_, x_)
#endif

#ifdef USE_ILUPACK
    case (SOLVER_ILUPACK)
      call ilupack_solve(this%solver_handle, rhs_, x_)
#endif
  end select
  x = reshape(x_, [size(x,1), size(x,2)])
end subroutine

module subroutine SPARSE_TRANSPOSE(this)
  !! Transpose sparse matrix.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  type(SPARSE_TYPE) :: s2

  ! h == 0 is default value for unfactorized
  ASSERT(this%solver_handle == 0)

  ! transpose matrix
  call this%transpose(s2)

  ! copy transposed matrix
  this%nrows = s2%nrows
  this%ncols = s2%ncols
  if (.not. s2%is_empty()) then
    this%ia = s2%ia
    this%ja = s2%ja
    this%a  = s2%a
  end if
end subroutine

module subroutine SPARSE_TRANSPOSE2(this, s2)
  !! Transpose sparse matrix, store result in second matrix.

  class(SPARSE_TYPE), intent(in)  :: this
    !! First sparse matrix
  type(SPARSE_TYPE),  intent(out) :: s2
    !! Output transposed sparse matrix

  integer                        :: i
  integer(SPARSE_IDX)            :: j, k
  type(vector_int),  allocatable :: col(:)
  type(VECTOR_TYPE), allocatable :: val(:)

  ! init s2
  call s2%init(this%ncols, ncols = this%nrows)

  ! do not copy values of an empty matrix
  if (this%is_empty()) return

  ! allocate memory
  allocate (col(this%ncols))
  allocate (val(this%ncols))
  do i = 1, this%ncols
    call col(i)%init(n = 0, c = 8)
    call val(i)%init(n = 0, c = 8)
  end do

  ! transpose matrix
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      call col(this%ja(j))%push(i)
      call val(this%ja(j))%push(this%a(j))
    end do
  end do

  ! allocate memory for csr
  allocate (s2%ia(s2%nrows + 1))
  allocate (s2%ja(size(this%ja,kind=SPARSE_IDX)))
  allocate (s2%a( size(this%a ,kind=SPARSE_IDX)))

  ! construct csr
  k        = 1
  s2%ia(1) = k
  do i = 1, size(col)
    do j = 1, col(i)%n
      s2%ja(k) = col(i)%d(j)
      s2%a( k) = val(i)%d(j)
      k        = k + 1
    end do
    s2%ia(i+1) = k
  end do
end subroutine

module subroutine SPARSE_MUL_SPARSE(this, s2, s3, sort)
  !! Multiply two sparse matrices, store result in sparse third matrix.
  !!
  !! operation: s3 <- this * s2

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse first matrix (unchanged)
  type(SPARSE_TYPE),  intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE),  intent(out) :: s3
    !! Output sparse third matrix
  logical, optional,  intent(in)  :: sort
    !! Optional: sort column indices (default: true)

  integer                          :: i, nthreads, ithread
  integer(SPARSE_IDX)              :: j, k, icol0, icol1, ncol
  integer(SPARSE_IDX), allocatable :: ops(:)
  integer,             allocatable :: i1(:), col(:)
  logical                          :: sort_
  logical,             allocatable :: fil(:)
  TT                               :: next_val
  TT,                  allocatable :: val(:), val1(:)

  ! check dimensions
  ASSERT(this%ncols == s2%nrows)

  ! init s3
  call s3%init(this%nrows, ncols = s2%ncols)
  if (this%is_empty() .or. s2%is_empty()) then
    ! return empty matrix if either this or s2 is empty
    return
  end if
  allocate (s3%ia(s3%nrows + 1), source = int(0, kind = SPARSE_IDX))

  ! optional args
  sort_ = .true.
  if (present(sort)) sort_ = sort

  ! multiplication operation counters for each result row
  allocate (ops(0:this%nrows), source = int(0, kind = SPARSE_IDX))

  ! work in parallel
  !$omp parallel default(none) &
  !$omp private(i, j, k, ithread, col, icol0, icol1, ncol, val, val1, next_val, fil) &
  !$omp shared(this, s2, s3, sort_, nthreads, i1, ops)
    !$omp single
      ! get actual number of threads
      nthreads = omp_get_num_threads()

      ! allocate upper work chunk index bounds
      allocate (i1(0:nthreads), source = 0)
    !$omp end single nowait

    ! get current thread id
    ithread = omp_get_thread_num() + 1

    ! count number of multiplications per row
    !$omp do schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        ops(i) = ops(i) + s2%ia(this%ja(j)+1) - s2%ia(this%ja(j))
      end do
    end do
    !$omp end do

    ! prefix sum (parallelization here probably not worth it)
    !$omp single
      do i = 2, this%nrows
        ops(i) = ops(i) + ops(i-1)
      end do
    !$omp end single

    ! give every thread roughly the same number of ops
    if (ithread < nthreads) then
      ! equidistant index
      j = nint(ops(this%nrows) * (real(ithread) / real(nthreads)), kind = SPARSE_IDX)

      ! find row index
      i1(ithread) = bin_search(ops, j) - 1
    else
      i1(ithread) = this%nrows
    end if
    !$omp barrier

    ! get upper limit on memory required
    j = ops(i1(ithread)) - ops(i1(ithread-1))
    j = max(0, int(min(int(j, kind=int64), int(i1(ithread) - i1(ithread-1), kind=int64) * int(s2%ncols, kind=int64)), kind=int32))

    ! allocate thread local memory
    allocate (col(j))
    allocate (val(j))
    allocate (val1(s2%ncols))
    allocate (fil( s2%ncols), source = .false.)

    ! reset column counters
    icol1 = 0
    ncol  = 0

    ! calculate result rows
    do i = i1(ithread-1) + 1, i1(ithread)
      ! loop over columns of first matrix
      do j = this%ia(i), this%ia(i+1)-1
        ! loop over columns of second matrix where the row is equal to column of the first matrix
        do k = s2%ia(this%ja(j)), s2%ia(this%ja(j)+1)-1
          ! calculate new value
          next_val = this%a(j) * s2%a(k)

          ! check whether element already exists
          if (fil(s2%ja(k))) then
            ! update existing value
            val1(s2%ja(k)) = val1(s2%ja(k)) + next_val
          else
            ! set new value and fill flag
            val1(s2%ja(k)) = next_val
            fil( s2%ja(k)) = .true.

            ! save column index
            ncol = ncol + 1
            col(ncol) = s2%ja(k)
          end if
        end do
      end do

      ! update column pointers
      icol0 = icol1 + 1
      icol1 = ncol

      ! sort column indices
      if (sort_) call qsort(col(icol0:icol1))

      ! save column pointers (at this point: indices to thread-local memory)
      s3%ia(i+1) = icol1 + 1

      ! save values
      val(icol0:icol1) = val1(col(icol0:icol1))

      ! reset fill flag (where it was set)
      if (i < i1(ithread)) fil(col(icol0:icol1)) = .false.
    end do
    !$omp barrier

    !$omp single
      ! update column pointers (change into indices to global memory)
      s3%ia(1) = 1
      do i = 1, nthreads
        s3%ia(i1(i-1)+2:i1(i)+1) = s3%ia(i1(i-1)+2:i1(i)+1) + (s3%ia(i1(i-1)+1) - 1)
      end do

      ! allocate global arrays
      allocate (s3%ja(s3%ia(s3%nrows+1)-1))
      allocate (s3%a (s3%ia(s3%nrows+1)-1))
    !$omp end single

    ! copy thread-local data to global arrays
    s3%ja(s3%ia(i1(ithread-1)+1):s3%ia(i1(ithread)+1)-1) = col(1:ncol)
    s3%a( s3%ia(i1(ithread-1)+1):s3%ia(i1(ithread)+1)-1) = val(1:ncol)

    ! deallocate thread-local memory (avoid memory leak due to openmp bug)
    deallocate (val1, fil, col, val)
  !$omp end parallel
end subroutine

module subroutine SPARSE_INPUT(this, file, funit)
  !! Read sparse matrix from file or stream (does not work for large matrices with more than 2^31 non-zero entries)

  class(SPARSE_TYPE),     intent(out) :: this
    !! Sparse matrix
  character(*), optional, intent(in)  :: file
    !! File name (alternative: provide funit)
  integer,      optional, intent(in)  :: funit
    !! File unit (alternative: provide file name)

  integer                   :: i, row, col, nrows, ncols, status, funit_
  character(:), allocatable :: fmt_
  real                      :: rval
#ifdef TCMPLX
  real                      :: ival
#endif
  type(vector_int)          :: rows, cols
  type(VECTOR_TYPE)         :: vals
  type(SPBUILD_TYPE)        :: sb

  allocate (character(0) :: fmt_) ! remove gfortran warning

  ! open file if needed
  ASSERT(.not. (present(file) .and. present(funit)))
  ASSERT(present(file) .or. present(funit))
  if (present(file)) then
    open (newunit = funit_, file = file, status = "old", action = "read", iostat = status)
    if (status /= 0) call program_error("File could not be opened for reading")
  else
    funit_ = funit
  end if

  ! read and save into vectors
  call rows%init(0, c = 1024)
  call cols%init(0, c = 1024)
  call vals%init(0, c = 1024)
  do
#ifdef TCMPLX
    read (funit_, *, iostat = status) row, col, rval, ival ! FIXME
#else
    read (funit_, *, iostat = status) row, col, rval
#endif

    if (status < 0) then
      ! end of file reached
      close (funit_)
      exit
    elseif (status > 0) then
      ! IO error
      close(funit_)
      call program_error("IO Error")
    end if

    call rows%push(row)
    call cols%push(col)
#ifdef TCMPLX
    call vals%push(cmplx(rval, ival))
#else
    call vals%push(rval)
#endif
  end do

  ! close file if openend here
  if (present(file)) close(funit_)

  ! get dimensions
  nrows = maxval(rows%d(1:rows%n))
  ncols = maxval(cols%d(1:cols%n))

  ! create matrix
  call this%init(nrows, ncols = ncols)
  call sb%init(this)
  do i = 1, rows%n
    call sb%add(rows%d(i), cols%d(i), vals%d(i))
  end do
  call sb%save()
end subroutine

module subroutine SPARSE_OUTPUT(this, file, funit, fmt)
  !! Write sparse matrix to file or stream

  class(SPARSE_TYPE),     intent(in) :: this
    !! Sparse matrix
  character(*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,      optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(*), optional, intent(in) :: fmt
    !! Value format string

  integer                   :: i, status, funit_
  integer(SPARSE_IDX)       :: j
  character(:), allocatable :: fmt_

  allocate (character(0) :: fmt_) ! remove gfortran warning
  if (present(fmt)) then
    fmt_ = fmt
  else
#ifdef TCMPLX
    fmt_ = '(ES24.16,SP,ES24.16,"i")'
#else
    fmt_ = '(ES24.16)'
#endif
  end if

  ! open file if needed
  ASSERT(.not. (present(file) .and. present(funit)))
  ASSERT(present(file) .or. present(funit))
  if (present(file)) then
    open (newunit = funit_, file = file, status = "replace", action = "write", iostat = status)
    if (status /= 0) call program_error("File could not be opened")
  else
    funit_ = funit
  end if

  ! output csr
  if (.not. this%is_empty()) then
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        ! write row and column
        write (funit_, "(I0,A,I0,A)", advance = "no") i, " ", this%ja(j), " "

        ! write value
        write (funit_, fmt_) this%a(j)
      end do
    end do
  end if

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

module function SPARSE_NNZ(this, only_nonzeros) result(nnz)
  !! computes number of non-zero elements

  class(SPARSE_TYPE), intent(in) :: this
  logical, optional,  intent(in) :: only_nonzeros
    !! assume that all entries of sparse matrix are non-zero. thus, is much faster.
    !! default: false
  integer(SPARSE_IDX)            :: nnz

  logical :: only_nonzeros_

  if (this%is_empty()) then
    nnz = 0
  else
    only_nonzeros_ = .false.
    if (present(only_nonzeros)) only_nonzeros_ = only_nonzeros

    if (only_nonzeros_) then
      nnz = this%ia(this%nrows+1)-1
    else
      nnz = count(abs(this%a) > 0.0, kind=SPARSE_IDX)
    end if
  end if
end function

module function SPARSE_IS_EMPTY(this) result(empty)
  !! check if sparse matrix is empty (CSR not allocated)

  class(SPARSE_TYPE), intent(in) :: this
  logical                        :: empty
    !! return true if CSR is not allocated, false otherwise

  if ((this%nrows == 0) .or. &
    & (this%ncols == 0) .or. &
    & (.not. allocated(this%ia)) .or. &
    & (.not. allocated(this%ja)) .or. &
    & (.not. allocated(this%a))) then
    empty = .true.
  else
    empty = .false.
  end if
end function

module function SPARSE_EYE(nrows, ncols) result(eye)
  !! sparse identity matrix

  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(SPARSE_TYPE)               :: eye
    !! return sparse identity matrix

  integer :: i, n

  call eye%init(nrows, ncols = ncols)
  n = min(eye%nrows, eye%ncols)

  ! allocate csr memory
  allocate (eye%ia(eye%nrows + 1))
  allocate (eye%ja(n))
  allocate (eye%a( n))

  eye%ia(1) = 1
  do i = 1, n
    eye%ja(i  ) = i
    eye%a( i  ) = 1
    eye%ia(i+1) = int(i+1, kind=SPARSE_IDX)
  end do
  do i = n + 1, eye%nrows
    eye%ia(i+1) = eye%ia(n+1)
  end do
end function

module function SPARSE_ZERO(nrows, ncols) result(zero)
  !! create empty sparse matrix

  integer,           intent(in) :: nrows
    !! Number of rows
  integer, optional, intent(in) :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(SPARSE_TYPE)             :: zero
    !! return empty sparse matrix

  call zero%init(nrows, ncols = ncols)

  ! allocate csr memory
  allocate (zero%ia(nrows+1), source = int(1, kind=SPARSE_IDX))
  allocate (zero%ja(0))
  allocate (zero%a( 0))
end function

module subroutine SPBUILD_INIT(this, s, n_per_row, assume_unique)
  !! Initialize sparse matrix builder.

  class(SPBUILD_TYPE),       intent(out) :: this
    !! Sparse matrix builder
  type(SPARSE_TYPE), target, intent(in)  :: s
    !! Sparse matrix this builder is tied to
  integer, optional,         intent(in)  :: n_per_row
    !! Optional: provide hint for number of elements per row (default: 8)
  logical, optional,         intent(in)  :: assume_unique
    !! Optional: user ensures no duplicate entries (default: .false.)

  integer :: i, c

  ! number of elements per row
  c = 8
  if (present(n_per_row)) c = n_per_row

  ! assume unique
  this%assume_unique = .false.
  if (present(assume_unique)) this%assume_unique = assume_unique

  ! save pointer to matrix
  this%sp => s

  ! allocate data
  allocate (this%cols(  s%nrows))
  allocate (this%vals(  s%nrows))
  allocate (this%filled(s%nrows))
  do i = 1, s%nrows
    call this%cols(i)%init(  n = 0, c = c)
    call this%vals(i)%init(  n = 0, c = c)
    call this%filled(i)%init(n = 0, c = c)
  end do

  ! allocate row flags
  allocate (this%sorted(s%nrows),      source = .true.)
  allocate (this%keep_struct(s%nrows), source = .true.)

  ! load existing data from sparse matrix
  call this%load()
end subroutine

module subroutine SPBUILD_DESTRUCT(this)
  !! Deallocate memory.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  nullify (this%sp)
  if (allocated(this%cols       )) deallocate (this%cols       )
  if (allocated(this%vals       )) deallocate (this%vals       )
  if (allocated(this%filled     )) deallocate (this%filled     )
  if (allocated(this%sorted     )) deallocate (this%sorted     )
  if (allocated(this%keep_struct)) deallocate (this%keep_struct)
end subroutine

module subroutine SPBUILD_RESET(this)
  !! Delete all values of each row

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer :: i

  do i = 1, this%sp%nrows
    call this%reset_row(i)
  end do
end subroutine

module subroutine SPBUILD_RESET_ROW(this, row)
  !! Delete all values of a row

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder
  integer,             intent(in)    :: row
    !! row index

  if (this%keep_struct(row)) then
    ! reset data while keeping structure
    this%vals(  row)%d(1:this%vals(  row)%n) = 0
    this%filled(row)%d(1:this%filled(row)%n) = .false.

  else
    ! reset data without keeping structure
    call this%cols(  row)%resize(0)
    call this%vals(  row)%resize(0)
    call this%filled(row)%resize(0)
    this%sorted(row)      = .true.
    this%keep_struct(row) = .true.
  end if
end subroutine

module subroutine SPBUILD_SET(this, row, col, val, fact, search)
  !! Insert/Update value.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder
  integer,             intent(in)    :: row
    !! Row index
  integer,             intent(in)    :: col
    !! Column index
  TT,                  intent(in)    :: val
    !! New value
  TT,      optional,   intent(in)    :: fact
    !! Scale existing element before adding new value (default: 0)
  logical, optional,   intent(in)    :: search
    !! Search for existing value instead of appending to end of row (default: true)

  TT      :: fact_
  logical :: search_
  integer :: i, l, r, m

  ASSERT(row > 0)
  ASSERT(row <= this%sp%nrows)
  ASSERT(col > 0)
  ASSERT(col <= this%sp%ncols)

  ! optional args
  search_ = .true.
  if (present(search)) search_ = search
  fact_ = 0
  if (present(fact)) fact_ = fact

  ! get index i for new value
  if (search_) then
    ! decide if binary or linear search
    if (this%sorted(row)) then
      ! binary search in sorted row
      i = -1
      l = 1
      r = this%cols(row)%n
      do while (l <= r)
        m = (l + r) / 2
        if (this%cols(row)%d(m) < col) then
          l = m + 1
        elseif (this%cols(row)%d(m) > col) then
          r = m - 1
        else
          i = m
          exit
        end if
      end do

      ! not found, append to end
      if (i == -1) i = this%cols(row)%n + 1
    else
      ! linear search
      do i = 1, this%cols(row)%n
        if (this%cols(row)%d(i) == col) exit
      end do
    end if
  else
    ! append element without searching
    i = this%cols(row)%n + 1

    ! do not keep structure when resetting
    this%keep_struct(row) = .false.
  end if

  ! Insert/Update value
  if (i > this%cols(row)%n) then ! append to end
    ! update sorted flag
    if (this%sorted(row) .and. (this%cols(row)%n > 0)) then
      if (col <= this%cols(row)%back()) this%sorted(row) = .false.
    end if

    ! make sure that appending without searching happens in order
    ! user can disable this check if he ensures no duplicates
    if (.not. this%assume_unique) then
      ASSERT(search_ .or. this%sorted(row))
    end if

    ! append new value
    call this%cols(row)%push(col)
    call this%vals(row)%push(val)
    call this%filled(row)%push(.true.)
  else
    ! update existing value
    if (this%filled(row)%d(i)) then
      this%vals(row)%d(i) = fact_ * this%vals(row)%d(i) + val
    else
      this%vals(row)%d(i) = val
      this%filled(row)%d(i) = .true.
    end if
  end if
end subroutine

module subroutine SPBUILD_SET_ROW(this, row, vals, j0, j1)
  !! Insert continuous block of values.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder
  integer,             intent(in)    :: row
    !! Row index
  TT,                  intent(in)    :: vals(:)
    !! Values to insert (j1 - j0 + 1 elements)
  integer, optional,   intent(in)    :: j0
    !! Start column index (default: 1)
  integer, optional,   intent(in)    :: j1
    !! End column index (default: this%sp%ncols)

  integer :: j, j0_, j1_

  ! optional arguments
  j0_ = 1
  if (present(j0)) j0_ = j0
  j1_ = this%sp%ncols
  if (present(j1)) j1_ = j1
  ASSERT(j0_ > 0)
  ASSERT(j0_ <= this%sp%ncols)
  ASSERT(j1_ > 0)
  ASSERT(j1_ <= this%sp%ncols)

  ! setting whole row -> reset previous data of that row.
  ! thus, making sure row is sorted.
  if ((1 == j0_) .and. (this%sp%ncols == j1_)) then
    call this%cols(  row)%reset
    call this%vals(  row)%reset
    call this%filled(row)%reset
    this%sorted(row)      = .true.
    this%keep_struct(row) = .true.
  end if

  ! insert values
  do j = j0_, j1_
    call this%set(row, j, vals(j-j0_+1), search = .false.)
  end do
end subroutine

module subroutine SPBUILD_ADD(this, row, col, val)
  !! Insert new element or add to existing value.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder
  integer,             intent(in)    :: row
    !! Row index
  integer,             intent(in)    :: col
    !! Column index
  TT,                  intent(in)    :: val
    !! New value

  TT :: fact

  ! set scaling factor to 1 both for real and complex values
  fact = 1

  ! add value
  call this%set(row, col, val, fact = fact)
end subroutine

module subroutine SPBUILD_SORT(this)
  !! Sort column indices.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer              :: i, n
  integer, allocatable :: perm(:)

  ! sort rows in parallel
  !$omp parallel default(none) private(i, n, perm) shared(this)

  ! thread local memory: allocate enough memory for permutation array
  allocate (perm(maxval(this%cols%n)))

  !$omp do schedule(dynamic)
  do i = 1, this%sp%nrows
    ! do nothing if row is already sorted
    if (this%sorted(i)) cycle
    this%sorted(i) = .true.

    ! number of elements in row
    n = this%cols(i)%n

    ! sort column indices and save permutation
    call qsort(this%cols(i)%d(1:n), perm = perm(1:n))

    ! apply permutation to values
    this%vals(  i)%d(1:n) = this%vals(  i)%d(perm(1:n))
    this%filled(i)%d(1:n) = this%filled(i)%d(perm(1:n))
  end do
  !$omp end do

  ! deallocate thread local memory (otherwise: memory leak!!!)
  deallocate (perm)

  !$omp end parallel
end subroutine

module subroutine SPBUILD_LOAD(this)
  !! Load existing values from sparse matrix into builder.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer             :: i
  integer(SPARSE_IDX) :: j

  ! load all rows
  do i = 1, this%sp%nrows
    ! delete values
    call this%cols(  i)%resize(0)
    call this%vals(  i)%resize(0)
    call this%filled(i)%resize(0)

    ! reset row flags
    this%sorted(i)      = .true.
    this%keep_struct(i) = .true.

    ! add values from csr (if not empty)
    if (allocated(this%sp%ia)) then
      do j = this%sp%ia(i), this%sp%ia(i+1)-1
        call this%cols(  i)%push(this%sp%ja(j))
        call this%vals(  i)%push(this%sp%a( j))
        call this%filled(i)%push(.true.       )
      end do
    end if
  end do
end subroutine

module subroutine SPBUILD_SAVE(this)
  !! Save values in sparse matrix.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer             :: i, j
  integer(SPARSE_IDX) :: k

  ! reset sparse matrix. needed in case it has old values
  call this%sp%reset()

  ! sort column indices
  call this%sort()

  ! get total number of elements
  k = 0
  do i = 1, this%sp%nrows
    do j = 1, this%cols(i)%n
      if (.not. this%filled(i)%d(j)) cycle
      k = k + 1
    end do
  end do

  ! allocate memory
  allocate (this%sp%ia(this%sp%nrows + 1))
  allocate (this%sp%ja(k))
  allocate (this%sp%a( k))

  ! init csr
  this%sp%ia(1) = 1
  k = 0
  do i = 1, this%sp%nrows
    do j = 1, this%cols(i)%n
      if (.not. this%filled(i)%d(j)) cycle
      k             = k + 1
      this%sp%ja(k) = this%cols(i)%d(j)
      this%sp%a( k) = this%vals(i)%d(j)
    end do

    ! next row pointer
    this%sp%ia(i+1) = k + 1
  end do
end subroutine

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef ARRAY_TYPE
#undef VECTOR_TYPE

#undef BAND_TYPE
#undef BAND_PTR_TYPE
#undef BLOCK_TYPE
#undef BLOCK_PTR_TYPE
#undef DENSE_TYPE
#undef DENSE_PTR_TYPE
#undef HESSENBERG_TYPE
#undef HESSENBERG_PTR_TYPE
#undef MATRIX_TYPE
#undef MATRIX_PTR_TYPE
#undef SPARSE_TYPE
#undef SPARSE_PTR_TYPE
#undef SPBUILD_TYPE
#undef TRIANG_TYPE
#undef TRIANG_PTR_TYPE

#undef SPARSE_INIT
#undef SPARSE_INIT_SOLVER
#undef SPARSE_DESTRUCT
#undef SPARSE_RESET
#undef SPARSE_SCALE
#undef SPARSE_MUL_VEC
#undef SPARSE_MUL_VEC_SLICE
#undef SPARSE_MUL_MAT
#undef SPARSE_FACTORIZE
#undef SPARSE_SOLVE_VEC
#undef SPARSE_SOLVE_MAT
#undef SPARSE_TRANSPOSE
#undef SPARSE_TRANSPOSE2
#undef SPARSE_MUL_SPARSE
#undef SPARSE_INPUT
#undef SPARSE_OUTPUT
#undef SPARSE_NNZ
#undef SPARSE_IS_EMPTY
#undef SPARSE_EYE
#undef SPARSE_ZERO

#undef SPBUILD_INIT
#undef SPBUILD_DESTRUCT
#undef SPBUILD_RESET
#undef SPBUILD_RESET_ROW
#undef SPBUILD_SET
#undef SPBUILD_SET_ROW
#undef SPBUILD_ADD
#undef SPBUILD_SORT
#undef SPBUILD_LOAD
#undef SPBUILD_SAVE
