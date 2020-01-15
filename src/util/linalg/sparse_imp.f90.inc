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
#define SPARSE_DESTRUCT      CONCAT3(sparse,T,destruct)
#define SPARSE_RESET         CONCAT3(sparse,T,reset)
#define SPARSE_SCALE         CONCAT3(sparse,T,scale)
#define SPARSE_MUL_VEC       CONCAT3(sparse,T,mul_vec)
#define SPARSE_MUL_MAT       CONCAT3(sparse,T,mul_mat)
#define SPARSE_FACTORIZE     CONCAT3(sparse,T,factorize)
#define SPARSE_SOLVE_VEC     CONCAT3(sparse,T,solve_vec)
#define SPARSE_SOLVE_MAT     CONCAT3(sparse,T,solve_mat)
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
#define SPARSE_EYE           CONCAT3(sparse,T,eye)

#define SPBUILD_INIT         CONCAT3(spbuild,T,init)
#define SPBUILD_DESTRUCT     CONCAT3(spbuild,T,destruct)
#define SPBUILD_RESET        CONCAT3(spbuild,T,reset)
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
end subroutine

subroutine SPARSE_DESTRUCT(this)
  !! Destruct sparse matrix by deallocating memory and deleting factorization.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! deallocate CSR data
  if (allocated(this%ia)) deallocate(this%ia)
  if (allocated(this%ja)) deallocate(this%ja)
  if (allocated(this%a )) deallocate(this%a )

  ! destruct pardiso handle
  if (this%pard > 0) call destruct_pardiso_handle(this%pard)
end subroutine

subroutine SPARSE_RESET(this, only_factorization)
  !! Reset dense matrix to zero (optionally only the factorization).

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  logical, optional,  intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  ! local variables
  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  ! deallocate CSR data
  if (.not. only_fact) then
    if (allocated(this%ia)) deallocate(this%ia)
    if (allocated(this%ja)) deallocate(this%ja)
    if (allocated(this%a )) deallocate(this%a )
  end if

  ! destruct factorization
  if (this%pard > 0) call destruct_pardiso_handle(this%pard)
end subroutine

subroutine SPARSE_SCALE(this, fact)
  !! Scale sparse matrix by a factor.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  TT,                 intent(in)    :: fact
    !! Scaling factor

  ! h == 0 is default value for unfactorized
  ASSERT(this%pard == 0)

  ! scale values
  this%a = fact * this%a
end subroutine

subroutine SPARSE_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(SPARSE_TYPE),         intent(in)    :: this
    !! Sparse matrix
  TT,                         intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                         intent(inout) :: y(:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transpose flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer          :: i, j
  TT               :: fact_y_
  character(len=1) :: trans_

  ! y scaling factor
  fact_y_ = 0
  if (present(fact_y)) fact_y_ = fact_y
  y = fact_y_ * y

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
    !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
      y(i) = y(i) + this%a(j) * x(this%ja(j))
      end do
    end do
    !$omp end parallel do
  else
#ifdef TCMPLX
    if (trans_ == 'C') then
      ! matrix vector multiplication with complex conjugated transposed matrix
      !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
      do i = 1, this%nrows
        do j = this%ia(i), this%ia(i+1)-1
          y(this%ja(j)) = y(this%ja(j)) + conjg(this%a(j)) * x(i)
        end do
      end do
      !$omp end parallel do
    else
#endif
      ! matrix vector multiplication with transposed matrix
      !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
      do i = 1, this%nrows
        do j = this%ia(i), this%ia(i+1)-1
          y(this%ja(j)) = y(this%ja(j)) + this%a(j) * x(i)
        end do
      end do
      !$omp end parallel do
#ifdef TCMPLX
    end if
#endif
  end if
end subroutine

subroutine SPARSE_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(SPARSE_TYPE),         intent(in)    :: this
    !! Sparse matrix
  TT,                         intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                         intent(inout) :: y(:,:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transpose flag for this matrix, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer          :: i, j
  TT               :: fact_y_
  character(len=1) :: trans_

  ! y scaling factor
  fact_y_ = 0
  if (present(fact_y)) fact_y_ = fact_y
  y = fact_y_ * y

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
    !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        y(i,:) = y(i,:) + this%a(j) * x(this%ja(j),:)
      end do
    end do
    !$omp end parallel do
  else
#ifdef TCMPLX
    if (trans_ == 'C') then
      ! matrix matrix multiplication with complex conjugated transposed matrix
      !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
      do i = 1, this%nrows
        do j = this%ia(i), this%ia(i+1)-1
          y(this%ja(j),:) = y(this%ja(j),:) + conjg(this%a(j)) * x(i,:)
        end do
      end do
      !$omp end parallel do
    else
#endif
      ! matrix matrix multiplication with transposed matrix
      !$omp parallel do private(i, j) shared(this, x, y) schedule(dynamic)
      do i = 1, this%nrows
        do j = this%ia(i), this%ia(i+1)-1
          y(this%ja(j),:) = y(this%ja(j),:) + this%a(j) * x(i,:)
        end do
      end do
      !$omp end parallel do
#ifdef TCMPLX
    end if
#endif
  end if
end subroutine

subroutine SPARSE_FACTORIZE(this)
  !! Factorize Matrix using PARDISO.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! h == 0 is default value for unfactorized
  ASSERT(this%pard == 0)

  ! init pardiso
#ifdef TCMPLX
  this%pard = create_pardiso_handle(this%nrows, .true.)
#else
  this%pard = create_pardiso_handle(this%nrows, .false.)
#endif

  ! factorize with pardiso
  call pardiso_factorize(this%pard, this%ia, this%ja, this%a)
end subroutine

subroutine SPARSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using PARDISO.
  !!
  !! Call this%factorize() first!

  class(SPARSE_TYPE),         intent(in)  :: this
    !! Sparse matrix
  TT,                         intent(in)  :: rhs(:)
    !! Right hand side
  TT,                         intent(out) :: x(:)
    !! Solution vector
  character(len=1), optional, intent(in)  :: trans
    !! Not implemented for 'T' or 'C' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed sparse matrix not implemented!")
  end if
  if (present(iter_refine)) then
    if (.not. iter_refine) call program_error("Solving with sparse matrix without iterative refinement not implemented!")
  end if

  ! solve system with pardiso
  call pardiso_solve(this%pard, this%ia, this%ja, this%a, rhs, x)
end subroutine

subroutine SPARSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using PARDISO.
  !!
  !! Call this%factorize() first.

  class(SPARSE_TYPE),         intent(in)  :: this
    !! Sparse matrix
  TT,                         intent(in)  :: rhs(:,:)
    !! Right hand side
  TT,                         intent(out) :: x(:,:)
    !! Solution vector
  character(len=1), optional, intent(in)  :: trans
    !! Not implemented for 'T' or 'C' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  ! local variables
  TT, allocatable :: rhs_(:), x_(:)

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

  ! solve with pardiso, pass right-hand sides one after another in 1D array
  rhs_ = reshape(rhs, [size(rhs,1) * size(rhs,2)])
  call pardiso_solve(this%pard, this%ia, this%ja, this%a, rhs_, x_)
  x = reshape(x_, [size(x,1), size(x,2)])
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

  ! local variables
  integer          :: i, j, i0_, j0_

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

  ! local variables
  integer          :: i, j, i0_, j0_
  logical          :: drop_zeros_

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
      elseif (drop_zeros_) then
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
  allocate (s%ia(size(this%ia)), source = this%ia     )
  allocate (s%ja(size(this%ja)), source = this%ja     )
  allocate (s%a( size(this%a )), source = real(this%a))
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
  allocate (s%ia(size(this%ia)), source = this%ia      )
  allocate (s%ja(size(this%ja)), source = this%ja      )
  allocate (s%a( size(this%a )), source = cmplx(this%a))
end subroutine
#endif

subroutine SPARSE_TRANSPOSE(this)
  !! Transpose sparse matrix.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! local variables
  type(SPARSE_TYPE) :: s2

  ! h == 0 is default value for unfactorized
  ASSERT(this%pard == 0)

  ! transpose matrix
  call this%transpose(s2)

  ! copy transposed matrix
  this%nrows = s2%nrows
  this%ncols = s2%ncols
  this%ia    = s2%ia
  this%ja    = s2%ja
  this%a     = s2%a
end subroutine

subroutine SPARSE_TRANSPOSE2(this, s2)
  !! Transpose sparse matrix, store result in second matrix.

  class(SPARSE_TYPE), intent(in)  :: this
    !! First sparse matrix
  type(SPARSE_TYPE),  intent(out) :: s2
    !! Output transposed sparse matrix

  ! local variables
  integer                        :: i, j, k
  type(vector_int),  allocatable :: col(:)
  type(VECTOR_TYPE), allocatable :: val(:)

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

  ! init s2
  call s2%init(this%ncols, ncols = this%nrows)

  ! allocate memory for csr
  allocate (s2%ia(s2%nrows + 1))
  allocate (s2%ja(size(this%ja)))
  allocate (s2%a( size(this%a)))

  ! construct csr
  s2%ia(1) = 1
  k = 0
  do i = 1, size(col)
    do j = 1, col(i)%n
      k = k + 1
      s2%ja(k) = col(i)%d(j)
      s2%a( k) = val(i)%d(j)
    end do
    s2%ia(i+1) = k + 1
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

  ! local variables
  integer              :: i, j1, j2
  TT                   :: fact_
  integer, allocatable :: ia(:)
  type(vector_int)     :: ja
  type(VECTOR_TYPE)    :: a

  ASSERT(this%nrows == s%nrows)
  ASSERT(this%ncols == s%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(this%pard == 0)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! allocate memory
  allocate (ia(this%nrows + 1))
  call ja%init(0, c = size(this%ja) + size(s%ja))
  call a%init( 0, c = size(this%ja) + size(s%ja))

  ! add matrices
  ia(1) = 1
  do i = 1, this%nrows
    j1 = this%ia(i)
    j2 = s%ia(i)
    do while ((j1 < this%ia(i+1)) .and. (j2 < s%ia(i+1)))
      if (this%ja(j1) < s%ja(j2)) then
        call ja%push(this%ja(j1))
        call a%push(this%a(j1))
        j1 = j1 + 1
      elseif (this%ja(j1) > s%ja(j2)) then
        call ja%push(s%ja(j2))
        call a%push(fact_ * s%a(j2))
        j2 = j2 + 1
      else
        call ja%push(this%ja(j1))
        call a%push(this%a(j1) + fact_ * s%a(j2))
        j1 = j1 + 1
        j2 = j2 + 1
      end if
    end do
    do while (j1 < this%ia(i+1))
      call ja%push(this%ja(j1))
      call a%push(this%a(j1))
      j1 = j1 + 1
    end do
    do while (j2 < s%ia(i+1))
      call ja%push(s%ja(j2))
      call a%push(fact_ * s%a(j2))
      j2 = j2 + 1
    end do
    ia(i+1) = ja%n + 1
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, this%ia)
  deallocate (this%ja)
  deallocate (this%a )
  allocate (this%ja(ja%n), source = ja%d(1:ja%n))
  allocate (this%a(  a%n), source =  a%d(1: a%n))
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

  ! s3 = fact1 * this + fact2 * s2
  call s3%init(this%nrows, ncols = this%ncols)
  s3%ia = this%ia
  s3%ja = this%ja
  if (present(fact1)) then
    s3%a = fact1 * this%a
  else
    s3%a = this%a
  end if
  call s3%add_sparse(s2, fact = fact2)
end subroutine

subroutine SPARSE_ADD_BAND(this, b, fact)
  !! Add sparse and band matrix, store result in sparse matrix.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix (updated)
  type(BAND_TYPE),    intent(in)    :: b
    !! Band matrix (unchanged)
  TT, optional,       intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  ! local variables
  integer              :: i, j1, j2
  integer, allocatable :: ia(:)
  TT                   :: fact_
  type(vector_int)     :: ja
  type(VECTOR_TYPE)    :: a

  ASSERT(this%nrows == b%nrows)
  ASSERT(this%ncols == b%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(this%pard == 0)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! add matrices
  do i = 1, this%nrows
    j1 = this%ia(i)
    j2 = max(i - b%nupper, 1)
    do while ((j1 < this%ia(i+1)) .and. (j2 <= min(i + b%nlower, b%ncols)))
      if (this%ja(j1) < j2) then
        call ja%push(this%ja(j1))
        call a%push(this%a(j1))
        j1 = j1 + 1
      elseif (this%ja(j1) > j2) then
        call ja%push(j2)
        call a%push(b%d(i-j2,j2))
        j2 = j2 + 1
      else
        call ja%push(j2)
        call a%push(this%a(j1) + fact_ * b%d(i-j2,j2))
        j1 = j1 + 1
        j2 = j2 + 1
      end if
    end do
    do while (j1 < this%ia(i+1))
      call ja%push(this%ja(j1))
      call a%push(this%a(j1))
      j1 = j1 + 1
    end do
    do while (j2 <= min(i + b%nlower, b%ncols))
      call ja%push(j2)
      call a%push(b%d(i-j2,j2))
      j2 = j2 + 1
    end do
    ia(i+1) = ja%n
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, this%ia)
  deallocate (this%ja)
  deallocate (this%a )
  allocate (this%ja(ja%n), source = ja%d(1:ja%n))
  allocate (this%a(  a%n), source =  a%d(1: a%n))
end subroutine

subroutine SPARSE_ADD_BAND3(this, b2, s3, fact1, fact2)
  !! Add sparse and band matrix, store result in sparse third matrix.

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
  s3%ia = this%ia
  s3%ja = this%ja
  if (present(fact1)) then
    s3%a = fact1 * this%a
  else
    s3%a = this%a
  end if
  call s3%add_band(b2, fact = fact2)
end subroutine

subroutine SPARSE_MUL_SPARSE(this, s2, s3)
  !! Multiply two sparse matrices, store result in sparse third matrix.
  !!
  !! operation: s3 <- this * s2

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse first matrix (unchanged)
  type(SPARSE_TYPE),  intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE),  intent(out) :: s3
    !! Output sparse third matrix

  ! local variables
  integer              :: i, j, k, num_threads, thread_i, icol0, icol1
  TT                   :: next_val
  integer, allocatable :: ia0(:), ia1(:), ja(:), thread_n(:), thread_i0(:), thread_i1(:)
  TT,      allocatable :: a(:), row(:)
  logical, allocatable :: fil(:)
  type(vector_int)     :: col, row_idx, colptr0, colptr1
  type(VECTOR_TYPE)    :: val

  ! check dimensions
  ASSERT(this%ncols == s2%nrows)

  ! init s3
  call s3%init(this%nrows, ncols = s2%ncols)

  ! allocate memory for global 4-array csr row pointers
  allocate (ia0(s3%nrows))
  allocate (ia1(s3%nrows))

  ! allocate memory for thread element counters and start stop indices
  num_threads = omp_get_max_threads()
  allocate (thread_n( num_threads), source = 0)
  allocate (thread_i0(num_threads), source = 0)
  allocate (thread_i1(num_threads), source = 0)

  ! work in parallel
  !$omp parallel private(i, j, k, thread_i, icol0, icol1, next_val, col, row_idx, colptr0, colptr1, val, row, fil) &
                 shared(this, s2, s3, ia0, ia1, ja, a, num_threads, thread_n, thread_i0, thread_i1)

  ! get current thread id
  thread_i = omp_get_thread_num() + 1

  ! allocate thread local memory
  call col%init(    0, c = 8*(s3%nrows+s3%ncols)) ! hopefully overallocated
  call val%init(    0, c = 8*(s3%nrows+s3%ncols))
  call row_idx%init(0, c =    s3%nrows          )
  call colptr0%init(0, c =    s3%nrows          )
  call colptr1%init(0, c =    s3%nrows          )
  allocate (row(s2%ncols)                  )
  allocate (fil(s2%ncols), source = .false.)

  ! init first column pointer
  icol1 = 0

  ! loop in parallel over rows of the first matrix
  !$omp do schedule(dynamic)
  do i = 1, this%nrows
    ! loop over columns of first matrix
    do j = this%ia(i), this%ia(i+1)-1
      ! loop over columns of second matrix where the row is equal to column of the first matrix
      do k = s2%ia(this%ja(j)), s2%ia(this%ja(j)+1)-1
        ! calculate new value
        next_val = this%a(j) * s2%a(k)

        ! check whether element already exists
        if (fil(s2%ja(k))) then
          ! update value
          row(s2%ja(k)) = row(s2%ja(k)) + next_val
        else
          ! set new value and fill flag, save column index
          row(s2%ja(k)) = next_val
          fil(s2%ja(k)) = .true.
          call col%push(s2%ja(k))
        end if
      end do
    end do

    ! update column pointers and thread element counter
    icol0 = icol1 + 1
    icol1 = col%n

    ! save row index and column pointers
    call row_idx%push(i)
    call colptr0%push(icol0)
    call colptr1%push(icol1)

    ! sort newly added column indices
    call qsort(col%d(icol0:icol1))

    ! save values
    call val%reserve(col%c)
    call val%resize( col%n)
    val%d(icol0:icol1) = row(col%d(icol0:icol1))

    ! reset fill flag where it was set
    fil(col%d(icol0:icol1)) = .false.
  end do
  !$omp end do nowait

  ! set thread element counter in global array
  thread_n(thread_i) = col%n
  !$omp barrier

  !$omp single
  ! set thread start/stop indices in global csr memory
  j = 0
  do i = 1, num_threads
    thread_i0(i) = j + 1
    thread_i1(i) = j + thread_n(i)
    j            = thread_i1(i)
  end do

  ! allocate global 4-array csr memory
  allocate (ja(thread_i1(num_threads)))
  allocate (a( thread_i1(num_threads)))
  !$omp end single

  ! copy data from each thread to global 4-array csr
  ja(thread_i0(thread_i):thread_i1(thread_i)) = col%d(1:col%n)
  a (thread_i0(thread_i):thread_i1(thread_i)) = val%d(1:val%n)

  ! set row pointers
  ia0(row_idx%d(1:row_idx%n)) = thread_i0(thread_i) + colptr0%d(1:row_idx%n) - 1
  ia1(row_idx%d(1:row_idx%n)) = thread_i0(thread_i) + colptr1%d(1:row_idx%n)

  !$omp end parallel

  ! convert 4-array csr to 3-array csr
  allocate (s3%ia(s3%nrows + 1))
  allocate (s3%ja(size(ja)))
  allocate (s3%a( size(a)))
  j = 0
  do i = 1, s3%nrows
    k = j + 1
    j = j + ia1(i) - ia0(i)
    s3%ja(k:j) = ja(ia0(i):ia1(i)-1)
    s3%a( k:j) = a( ia0(i):ia1(i)-1)
    s3%ia(  i) = k
  end do
  s3%ia(s3%nrows+1) = j + 1
end subroutine

subroutine SPARSE_OUTPUT(this, file, funit, fmt)
  !! Write sparse matrix to file or stream

  class(SPARSE_TYPE),         intent(in) :: this
    !! Sparse matrix
  character(len=*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,          optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(len=*), optional, intent(in) :: fmt
    !! Value format string

  ! local variables
  integer                       :: i, j, status, funit_
  character(len=:), allocatable :: fmt_

  if (present(fmt)) then
    fmt_ = fmt
  else
#ifdef TCMPLX
    ! FIXME: check this
    fmt_ = '(1A,F0.0,SP,F0.0,"i")'
#else
    fmt_ = "(1A,1F0.0)"
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
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      ! write row and column
      write(funit_, "(1I0,1A,1I0,1A)", advance = "no") i, " ", this%ja(j), " "

      ! write value
      write(funit_, fmt_) this%a(j)
    end do
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

function SPARSE_EYE(nrows, ncols) result(eye)
  !! sparse identity matrix

  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(SPARSE_TYPE)               :: eye

  ! local variables
  integer :: i, n

  call eye%init(nrows, ncols)
  n = min(eye%nrows, eye%ncols)

  ! allocate csr memory
  allocate (eye%ia(eye%nrows + 1))
  allocate (eye%ja(n))
  allocate (eye%a( n))

  eye%ia(1) = 1
  do i = 1, n
    eye%ja(i  ) = i
    eye%a( i  ) = 1
    eye%ia(i+1) = i+1
  end do
  do i = n + 1, eye%nrows
    eye%ia(i+1) = eye%ia(n+1)
  end do
end function

subroutine SPBUILD_INIT(this, s, n_per_row)
  !! Initialize sparse matrix builder.

  class(SPBUILD_TYPE),       intent(out) :: this
    !! Sparse matrix builder
  type(SPARSE_TYPE), target, intent(in)  :: s
    !! Sparse matrix this builder is tied to
  integer, optional,         intent(in)  :: n_per_row
    !! Optional: provide hint for number of elements per row (default: 8)

  ! local variables
  integer :: i, c

  ! number of elements per row
  c = 8
  if (present(n_per_row)) c = n_per_row

  ! save pointer to matrix
  this%sp => s

  ! allocate data
  allocate(this%cols(  s%nrows))
  allocate(this%vals(  s%nrows))
  allocate(this%filled(s%nrows))
  do i = 1, s%nrows
    call this%cols(i)%init(n = 0, c = c)
    call this%vals(i)%init(n = 0, c = c)
    call this%filled(i)%init(n = 0, c = c)
  end do

  ! allocate row flags
  allocate(this%sorted(s%nrows), source = .true.)
  allocate(this%keep_struct(s%nrows), source = .true.)

  ! load existing data from sparse matrix
  call this%load()
end subroutine

subroutine SPBUILD_DESTRUCT(this)
  !! Deallocate memory.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  nullify(this%sp)
  if (allocated(this%cols)) deallocate(this%cols)
  if (allocated(this%vals)) deallocate(this%vals)
  if (allocated(this%filled)) deallocate(this%filled)
  if (allocated(this%sorted)) deallocate(this%sorted)
  if (allocated(this%keep_struct)) deallocate(this%keep_struct)
end subroutine

subroutine SPBUILD_RESET(this)
  !! Delete all values

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  ! local variables
  integer :: i

  ! reset all rows
  do i = 1, this%sp%nrows
    if (this%keep_struct(i)) then
      ! reset data while keeping structure
      this%vals(  i)%d(1:this%vals(  i)%n) = 0
      this%filled(i)%d(1:this%filled(i)%n) = .false.
    else
      ! reset data without keeping structure
      call this%cols(i)%resize(0)
      call this%vals(i)%resize(0)
      call this%filled(i)%resize(0)
      this%sorted(i)      = .true.
      this%keep_struct(i) = .true.
    end if
  end do
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

  ! local variables
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
      if (col <= this%cols(row)%last()) this%sorted(row) = .false.
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

  ! local variables
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

  ! local variables
  TT :: fact

  ! set scaling factor to 1 both for real and complex values
  fact = 1

  ! add value
  call this%set(row, col ,val, fact = fact)
end subroutine

subroutine SPBUILD_SORT(this)
  !! Sort column indices.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  ! local variables
  integer              :: i, n
  integer, allocatable :: perm(:)

  ! sort rows in parallel
  !$omp parallel private(i, n, perm) shared(this)

  ! allocate enough memory for permutation array
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
  !$omp end parallel
end subroutine

subroutine SPBUILD_LOAD(this)
  !! Load existing values from sparse matrix into builder.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  ! local variables
  integer :: i, j

  ! load all rows
  do i = 1, this%sp%nrows
    ! delete values
    call this%cols(i)%resize(0)
    call this%vals(i)%resize(0)
    call this%filled(i)%resize(0)

    ! reset row flags
    this%sorted(i)      = .true.
    this%keep_struct(i) = .true.

    ! add values from csr (if not empty)
    if (allocated(this%sp%ia)) then
      do j = this%sp%ia(i), this%sp%ia(i+1)-1
        call this%cols(i)%push(this%sp%ja(j))
        call this%vals(i)%push(this%sp%a(j))
        call this%filled(i)%push(.true.)
      end do
    end if
  end do
end subroutine

subroutine SPBUILD_SAVE(this)
  !! Save values in sparse matrix.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  integer :: i, j, k

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
#undef SPARSE_DESTRUCT
#undef SPARSE_RESET
#undef SPARSE_SCALE
#undef SPARSE_MUL_VEC
#undef SPARSE_MUL_MAT
#undef SPARSE_FACTORIZE
#undef SPARSE_SOLVE_VEC
#undef SPARSE_SOLVE_MAT
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
#undef SPARSE_EYE

#undef SPBUILD_INIT
#undef SPBUILD_DESTRUCT
#undef SPBUILD_RESET
#undef SPBUILD_SET
#undef SPBUILD_SET_ROW
#undef SPBUILD_ADD
#undef SPBUILD_SORT
#undef SPBUILD_LOAD
#undef SPBUILD_SAVE
