#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define BAND_TYPE            CONCAT(band,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define TRIANG_TYPE          CONCAT(triang,T)
#define BLOCK_TYPE           CONCAT(block,T)

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
#define SPARSE_TO_BAND       CONCAT3(sparse,T,to_band)
#define SPARSE_TO_DENSE      CONCAT3(sparse,T,to_dense)
#define SPARSE_TO_SPARSE     CONCAT3(sparse,T,to_sparse)
#ifdef TCMPLX
#define SPARSE_TO_REAL       CONCAT3(sparse,T,to_real)
#else
#define SPARSE_TO_CMPLX      CONCAT3(sparse,T,to_cmplx)
#endif
#define SPARSE_TRANSPOSE     CONCAT3(sparse,T,transpose)
#define SPARSE_TRANSPOSE2    CONCAT3(sparse,T,transpose2)
#define SPARSE_ADD_SPARSE    CONCAT3(sparse,T,add_sparse)
#define SPARSE_ADD_SPARSE3   CONCAT3(sparse,T,add_sparse3)
#define SPARSE_ADD_BAND      CONCAT3(sparse,T,add_band)
#define SPARSE_ADD_BAND3     CONCAT3(sparse,T,add_band3)
#define SPARSE_MUL_SPARSE    CONCAT3(sparse,T,mul_sparse)
#define SPARSE_OUTPUT        CONCAT3(sparse,T,output)
#define SPARSE_FROM_DIAG     CONCAT3(sparse,T,from_diag)
#define SPARSE_TO_DIAG       CONCAT3(sparse,T,to_diag)
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

subroutine SPARSE_INIT(this, nrows, ncols)
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

subroutine SPARSE_DESTRUCT(this)
  !! Destruct sparse matrix by deallocating memory and deleting factorization.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  call this%reset()
end subroutine

subroutine SPARSE_RESET(this, only_factorization)
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

subroutine SPARSE_SCALE(this, fact)
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

subroutine SPARSE_MUL_VEC_SLICE(this, x, y, fact_y, i0, i1)
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

subroutine SPARSE_MUL_VEC(this, x, y, fact_y, trans)
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

subroutine SPARSE_MUL_MAT(this, x, y, fact_y, trans)
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

subroutine SPARSE_INIT_SOLVER(this)
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

subroutine SPARSE_FACTORIZE(this)
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

subroutine SPARSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
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

subroutine SPARSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
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

subroutine SPARSE_TO_BAND(this, b, i0, j0)
  !! Insert sparse matrix into a band matrix.

  class(SPARSE_TYPE), intent(in)    :: this
    !! Sparse matrix (unchanged)
  type(BAND_TYPE),    intent(inout) :: b
    !! Band matrix (updated)
  integer, optional,  intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional,  intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j

  ASSERT(.not. this%is_empty())
  ASSERT(allocated(b%d))

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values into band matrix
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      call b%set_elem(i0_+i-1, j0_+this%ja(j)-1, this%a(j))
    end do
  end do
end subroutine

subroutine SPARSE_TO_DENSE(this, d, i0, j0)
  !! Insert sparse matrix into a dense matrix.

  class(SPARSE_TYPE), intent(in)    :: this
    !! Sparse matrix (unchanged)
  type(DENSE_TYPE),   intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional,  intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional,  intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j

  ASSERT(.not. this%is_empty())
  ASSERT(allocated(d%d))

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values into dense matrix
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      d%d(i0_+i-1,j0_+this%ja(j)-1) = this%a(j)
    end do
  end do
end subroutine

subroutine SPARSE_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert sparse matrix into another sparse matrix.

  class(SPARSE_TYPE), intent(in)    :: this
    !! First sparse matrix (unchanged)
  type(SPBUILD_TYPE), intent(inout) :: sb
    !! Second sparse matrix represented by sparse matrix builder (updated)
  integer, optional,  intent(in)    :: i0
    !! optional: row offset in sparse matrix
  integer, optional,  intent(in)    :: j0
    !! optional: col offset in sparse matrix
  logical, optional,  intent(in)    :: struct(:,:)
    !! optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional,  intent(in)    :: drop_zeros
    !! optional: determine sparsity structure by checking for zero entries

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j
  logical             :: drop_zeros_

  ! do nothing if empty
  if (this%is_empty()) return

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values into second sparse matrix
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,this%ja(j))) cycle
      else if (drop_zeros_) then
        if (this%a(j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+this%ja(j)-1, this%a(j), search = .false.)
    end do
  end do
end subroutine

#ifdef TCMPLX
subroutine SPARSE_TO_REAL(this, s)
  !! Convert sparse complex matrix to real by ignoring imaginary part

  class(sparse_cmplx), intent(in)  :: this
    !! Sparse complex matrix
  class(sparse_real),  intent(out) :: s
    !! Output real(this)

  ! init s
  call s%init(this%nrows, ncols = this%ncols)

  ! copy csr
  if (.not. this%is_empty()) then
    allocate (s%ia(size(this%ia,kind=SPARSE_IDX)), source = this%ia     )
    allocate (s%ja(size(this%ja,kind=SPARSE_IDX)), source = this%ja     )
    allocate (s%a( size(this%a ,kind=SPARSE_IDX)), source = real(this%a))
  end if
end subroutine
#else
subroutine SPARSE_TO_CMPLX(this, s)
  !! Convert sparse real matrix to complex

  class(sparse_real),  intent(in)  :: this
    !! Sparse real matrix
  class(sparse_cmplx), intent(out) :: s
    !! Output cmplx(this)

  ! init s
  call s%init(this%nrows, ncols = this%ncols)

  ! copy csr
  if (.not. this%is_empty()) then
    allocate (s%ia(size(this%ia,kind=SPARSE_IDX)), source = this%ia      )
    allocate (s%ja(size(this%ja,kind=SPARSE_IDX)), source = this%ja      )
    allocate (s%a( size(this%a ,kind=SPARSE_IDX)), source = cmplx(this%a))
  end if
end subroutine
#endif

subroutine SPARSE_TRANSPOSE(this)
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

subroutine SPARSE_TRANSPOSE2(this, s2)
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

subroutine SPARSE_ADD_SPARSE(this, s, fact)
  !! Add two sparse matrices, store result in first matrix.
  !!
  !! operation: this <- this + fact * s

  class(SPARSE_TYPE), intent(inout) :: this
    !! First sparse matrix (updated)
  type(SPARSE_TYPE),  intent(in)    :: s
    !! Second sparse matrix (unchanged)
  TT, optional,       intent(in)    :: fact
    !! Second sparse matrix scaling factor (default: 1)

  integer                          :: i
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, j2, k
  integer(SPARSE_IDX), allocatable :: ia(:)
  TT                               :: fact_
  TT,                  allocatable :: a(:)

  ASSERT(this%nrows == s%nrows)
  ASSERT(this%ncols == s%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(this%solver_handle == 0)

  ! do nothing if s is empty
  if (s%is_empty()) return

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! simply copy and scale s if this matrix is empty
  if (this%is_empty()) then
    this%ia = s%ia
    this%ja = s%ja
    this%a  = s%a * fact_
    return
  end if

  ! allocate memory
  allocate (ia(this%nrows + 1))
  allocate (ja(size(this%ja,kind=SPARSE_IDX) + size(s%ja,kind=SPARSE_IDX)))
  allocate ( a(size(this%ja,kind=SPARSE_IDX) + size(s%ja,kind=SPARSE_IDX)))

  ! add matrices
  k     = 1
  ia(1) = k
  do i = 1, this%nrows
    j1 = this%ia(i)
    j2 = s%ia(i)
    do while ((j1 < this%ia(i+1)) .and. (j2 < s%ia(i+1)))
      if (this%ja(j1) < s%ja(j2)) then
        ja(k) = this%ja(j1)
        a( k) = this%a( j1)
        k     = k  + 1
        j1    = j1 + 1
      elseif (this%ja(j1) > s%ja(j2)) then
        ja(k) = s%ja(j2)
        a( k) = fact_ * s%a(j2)
        k     = k  + 1
        j2    = j2 + 1
      else
        ja(k) = this%ja(j1)
        a( k) = this%a( j1) + fact_ * s%a(j2)
        k     = k  + 1
        j1    = j1 + 1
        j2    = j2 + 1
      end if
    end do
    do while (j1 < this%ia(i+1))
      ja(k) = this%ja(j1)
      a( k) = this%a( j1)
      k     = k  + 1
      j1    = j1 + 1
    end do
    do while (j2 < s%ia(i+1))
      ja(k) = s%ja(j2)
      a( k) = fact_ * s%a(j2)
      k     = k  + 1
      j2    = j2 + 1
    end do
    ia(i+1) = k
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, this%ia)
  deallocate (this%ja)
  deallocate (this%a )
  allocate (this%ja(k-1), source = ja(1:k-1))
  allocate (this%a (k-1), source =  a(1:k-1))
end subroutine

subroutine SPARSE_ADD_SPARSE3(this, s2, s3, fact1, fact2)
  !! Add two sparse matrices, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * this + fact2 * s2

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse first matrix (unchanged)
  type(SPARSE_TYPE),  intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE),  intent(out) :: s3
    !! Output sparse third matrix
  TT, optional,       intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,       intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ASSERT(this%nrows == s2%nrows)
  ASSERT(this%ncols == s2%ncols)

  ! init result matrix
  call s3%init(this%nrows, ncols = this%ncols)

  ! s3 = fact1 * this + fact2 * s2
  if (.not. this%is_empty()) then
    s3%ia = this%ia
    s3%ja = this%ja
    if (present(fact1)) then
      s3%a = fact1 * this%a
    else
      s3%a = this%a
    end if
  end if
  call s3%add_sparse(s2, fact = fact2)
end subroutine

subroutine SPARSE_ADD_BAND(this, b, fact)
  !! Add sparse and band matrix, store result in sparse matrix.
  !!
  !! operation: this <- this + fact * b

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix (updated)
  type(BAND_TYPE),    intent(in)    :: b
    !! Band matrix (unchanged)
  TT, optional,       intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  integer                          :: i, j2
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, k, b_nrows, b_nlower, b_nupper
  integer(SPARSE_IDX), allocatable :: ia(:)
  TT                               :: fact_
  TT,                  allocatable :: a(:)

  ASSERT(this%nrows == b%nrows)
  ASSERT(this%ncols == b%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(this%solver_handle == 0)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! if this is empty, set to zero matrix first
  if (this%is_empty()) then
    allocate (this%ia(this%nrows+1), source = int(1, kind = SPARSE_IDX))
    allocate (this%ja(0))
    allocate (this%a( 0))
  end if

  ! allocate memory
  allocate (ia(this%nrows + 1))
  b_nrows  = int(b%nrows ,kind=SPARSE_IDX)
  b_nlower = int(b%nlower,kind=SPARSE_IDX)
  b_nupper = int(b%nupper,kind=SPARSE_IDX)
  k = b_nrows + b_nlower * b_nrows - b_nlower * (b_nlower + 1) / 2 + b_nupper * b_nrows - b_nupper * (b_nupper + 1) / 2
  k = k + size(this%ja,kind=SPARSE_IDX)
  allocate (ja(k))
  allocate ( a(k))

  ! add matrices
  k     = 1
  ia(1) = k
  do i = 1, this%nrows
    j1 = this%ia(i)
    j2 = max(i - b%nlower, 1)
    do while ((j1 < this%ia(i+1)) .and. (j2 <= min(i + b%nupper, b%ncols)))
      if (this%ja(j1) < j2) then
        ja(k) = this%ja(j1)
        a( k) = this%a( j1)
        k     = k  + 1
        j1    = j1 + 1
      elseif (this%ja(j1) > j2) then
        ja(k) = j2
        a( k) = fact_ * b%d(i-j2,j2)
        k     = k  + 1
        j2    = j2 + 1
      else
        ja(k) = j2
        a( k) = this%a(j1) + fact_ * b%d(i-j2,j2)
        k     = k  + 1
        j1    = j1 + 1
        j2    = j2 + 1
      end if
    end do
    do while (j1 < this%ia(i+1))
      ja(k) = this%ja(j1)
      a( k) = this%a( j1)
      k     = k  + 1
      j1    = j1 + 1
    end do
    do while (j2 <= min(i + b%nupper, b%ncols))
      ja(k) = j2
      a( k) = fact_ * b%d(i-j2,j2)
      k     = k  + 1
      j2    = j2 + 1
    end do
    ia(i+1) = k
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, this%ia)
  deallocate (this%ja)
  deallocate (this%a )
  allocate (this%ja(k-1), source = ja(1:k-1))
  allocate (this%a( k-1), source =  a(1:k-1))
end subroutine

subroutine SPARSE_ADD_BAND3(this, b2, s3, fact1, fact2)
  !! Add sparse and band matrix, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * this + fact2 * b2

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse first matrix (unchanged)
  type(BAND_TYPE),    intent(in)  :: b2
    !! Band second matrix (unchanged)
  type(SPARSE_TYPE),  intent(out) :: s3
    !! Output sparse third matrix
  TT, optional,       intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,       intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ASSERT(this%nrows == b2%nrows)
  ASSERT(this%ncols == b2%ncols)

  ! s3 = fact1 * this + fact2 * s2
  call s3%init(this%nrows, ncols = this%ncols)
  if (.not. this%is_empty()) then
    s3%ia = this%ia
    s3%ja = this%ja
    if (present(fact1)) then
      s3%a = fact1 * this%a
    else
      s3%a = this%a
    end if
  end if
  call s3%add_band(b2, fact = fact2)
end subroutine

subroutine SPARSE_MUL_SPARSE(this, s2, s3, sort)
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

  integer              :: i, nthreads, ithread
  integer(SPARSE_IDX)  :: j, k, nnz1, nnz2, nnz3, icol0, icol1, ncol, nval
  integer, allocatable :: i0(:), i1(:), col(:), tmp(:)
  logical              :: sort_
  logical, allocatable :: fil(:)
  real                 :: weight
  TT                   :: next_val
  TT,      allocatable :: val(:), val1(:), vtmp(:)

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

  ! number of non-zero entries
  nnz1 = this%ia(this%nrows+1)
  nnz2 = s2%ia(    s2%nrows+1)

  ! work in parallel
  !$omp parallel default(none) &
  !$omp private(i, j, k, ithread, col, icol0, icol1, ncol, val, nval, val1, next_val, fil, tmp, vtmp) &
  !$omp shared(this, s2, s3, sort_, nthreads, nnz1, nnz2, nnz3, i0, i1, weight)

  ! manual work partitioning
  !$omp single

  ! get actual number of threads
  nthreads = omp_get_num_threads()

  ! allocate work chunk index bounds
  allocate (i0(nthreads), source = 1)
  allocate (i1(nthreads), source = 0)

  ! first thread
  ithread = 1
  weight  = 0.0

  do i = 1, this%nrows
    ! weight estimate work per result row
    weight = weight + real(this%ia(i+1)-this%ia(i)) / real(nnz1)

    if (weight >= 1.0 / real(nthreads)) then
      ! split work here
      i1(ithread) = i

      ! next thread
      ithread     = ithread + 1
      i0(ithread) = i1(ithread-1) + 1
      if (ithread == nthreads) then
        i1(ithread) = this%nrows
        exit
      end if

      ! reset weight for next thread
      weight = 0.0
    end if
  end do
  !$omp end single

  ! get current thread
  ithread = omp_get_thread_num() + 1

  ! allocate thread local memory
  allocate (val1(s2%ncols))
  allocate (fil( s2%ncols), source = .false.)
  allocate (col(4*(nnz1+nnz2))) ! try to overallocate, but not too much
  allocate (val(4*(nnz1+nnz2)))

  icol1 = 0
  ncol  = 0
  nval  = 0

  ! loop over chunk of rows
  do i = i0(ithread), i1(ithread)
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
          if (ncol == size(col, kind = SPARSE_IDX)) then
            ! increase capacity
            allocate (tmp(ncol * 2))
            tmp(1:ncol) = col
            call move_alloc(tmp, col)
          end if
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

    ! save column pointers (here: indices to thread-local memory)
    s3%ia(i+1) = icol1 + 1

    ! save values
    if (size(val, kind = SPARSE_IDX) < size(col, kind = SPARSE_IDX)) then
      ! increase capacity
      allocate (vtmp(size(col, kind = SPARSE_IDX)))
      vtmp(1:nval) = val
      call move_alloc(vtmp, val)
    end if
    val(icol0:icol1) = val1(col(icol0:icol1))
    nval = ncol

    ! reset fill flag (where it was set)
    fil(col(icol0:icol1)) = .false.
  end do
  !$omp barrier

  !$omp single

  ! shift column indices to global
  s3%ia(1) = 1
  do i = 1, nthreads
    s3%ia(i0(i)+1:i1(i)+1) = s3%ia(i0(i)+1:i1(i)+1) + (s3%ia(i0(i)) - 1)
  end do

  ! allocate global memory
  nnz3 = s3%ia(s3%nrows+1)-1
  allocate (s3%ja(nnz3))
  allocate (s3%a (nnz3))

  !$omp end single

  ! copy data to global array
  s3%ja(s3%ia(i0(ithread)):s3%ia(i1(ithread)+1)-1) = col(1:ncol)
  s3%a( s3%ia(i0(ithread)):s3%ia(i1(ithread)+1)-1) = val(1:nval)

  ! deallocate thread-local memory (avoid memory leak due to openmp bug)
  deallocate (val1, fil, col, val)

  !$omp end parallel
end subroutine

subroutine SPARSE_OUTPUT(this, file, funit, fmt)
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
    ! FIXME: check this
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

subroutine SPARSE_FROM_DIAG(this, diag)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  class(SPARSE_TYPE), intent(out) :: this
    !! sparse diagonal matrix
  TT,                 intent(in)  :: diag(:)
    !! diagonal array

  integer :: i, n

  n = size(diag)

  call this%init(n)
  allocate (this%ia(n+1), this%ja(n), this%a(n))

  this%a = diag
  do i = 1, n
    this%ia(i) = int(i, kind = SPARSE_IDX)
    this%ja(i) = i
  end do
  this%ia(n+1) = n+1
end subroutine

subroutine SPARSE_TO_DIAG(this, diag)
  !! extract main diagonal of a sparse matrix

  class(SPARSE_TYPE), intent(in)  :: this
  TT                , intent(out) :: diag(:)
    !! diagonal elements of matrix

  integer             :: i, n
  integer(SPARSE_IDX) :: j

  n = min(this%nrows, this%ncols)
  ASSERT(n == size(diag))

  ! empty matrix => zero diagonal
  if (this%is_empty()) then
    diag = 0
    return
  end if

  !$omp parallel do   &
  !$omp default(none) &
  !$omp private(i,j)  &
  !$omp shared(n,diag,this)
  do i = 1, n
    ! get value closest to main diagonal
    j = this%ia(i) + (bin_search(this%ja(this%ia(i):this%ia(i+1)-1), i) - 1)

    ! output value only, if found exactly
    if (this%ja(j) == i) then
      diag(i) = this%a(j)
    else
      diag(i) = 0
    end if
  end do
  !$omp end parallel do
end subroutine

function SPARSE_NNZ(this, only_nonzeros) result(nnz)
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

function SPARSE_IS_EMPTY(this) result(empty)
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

function SPARSE_EYE(nrows, ncols) result(eye)
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

function SPARSE_ZERO(nrows, ncols) result(zero)
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

subroutine SPBUILD_INIT(this, s, n_per_row)
  !! Initialize sparse matrix builder.

  class(SPBUILD_TYPE),       intent(out) :: this
    !! Sparse matrix builder
  type(SPARSE_TYPE), target, intent(in)  :: s
    !! Sparse matrix this builder is tied to
  integer, optional,         intent(in)  :: n_per_row
    !! Optional: provide hint for number of elements per row (default: 8)

  integer :: i, c

  ! number of elements per row
  c = 8
  if (present(n_per_row)) c = n_per_row

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

subroutine SPBUILD_DESTRUCT(this)
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

subroutine SPBUILD_RESET(this)
  !! Delete all values of each row

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer :: i

  do i = 1, this%sp%nrows
    call this%reset_row(i)
  end do
end subroutine

subroutine SPBUILD_RESET_ROW(this, row)
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

subroutine SPBUILD_SET(this, row, col, val, fact, search)
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
    ASSERT(search_ .or. this%sorted(row))

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

subroutine SPBUILD_SET_ROW(this, row, vals, j0, j1)
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

subroutine SPBUILD_ADD(this, row, col, val)
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

subroutine SPBUILD_SORT(this)
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

subroutine SPBUILD_LOAD(this)
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

subroutine SPBUILD_SAVE(this)
  !! Save values in sparse matrix.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer             :: i, j
  integer(SPARSE_IDX) :: k

  ! reset sparse matrix. needed in case it has old values
  call this%sp%reset

  ! sort column indices
  call this%sort

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
#undef DENSE_TYPE
#undef SPARSE_TYPE
#undef SPBUILD_TYPE
#undef BAND_TYPE
#undef HESSENBERG_TYPE
#undef TRIANG_TYPE
#undef BLOCK_TYPE

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
#undef SPARSE_TO_BAND
#undef SPARSE_TO_DENSE
#undef SPARSE_TO_SPARSE
#undef SPARSE_TO_REAL
#undef SPARSE_TO_CMPLX
#undef SPARSE_TRANSPOSE
#undef SPARSE_TRANSPOSE2
#undef SPARSE_ADD_SPARSE
#undef SPARSE_ADD_SPARSE3
#undef SPARSE_ADD_BAND
#undef SPARSE_ADD_BAND3
#undef SPARSE_MUL_SPARSE
#undef SPARSE_OUTPUT
#undef SPARSE_FROM_DIAG
#undef SPARSE_TO_DIAG
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
