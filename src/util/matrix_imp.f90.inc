#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_TYPE          CONCAT(vector,T)

#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_ALC_TYPE      CONCAT(matrix_alc,T)
#define MATRIX_INIT          CONCAT3(matrix,T,init)
#define MATRIX_DESTRUCT      CONCAT3(matrix,T,destruct)
#define MATRIX_RESET         CONCAT3(matrix,T,reset)
#define MATRIX_SCALE         CONCAT3(matrix,T,scale)
#define MATRIX_MUL_VEC       CONCAT3(matrix,T,mul_vec)
#define MATRIX_MUL_MAT       CONCAT3(matrix,T,mul_mat)
#define MATRIX_FACTORIZE     CONCAT3(matrix,T,factorize)
#define MATRIX_SOLVE_VEC     CONCAT3(matrix,T,solve_vec)
#define MATRIX_SOLVE_MAT     CONCAT3(matrix,T,solve_mat)
#define MATRIX_TO_DENSE      CONCAT3(matrix,T,to_dense)
#define MATRIX_TO_SPARSE     CONCAT3(matrix,T,to_sparse)

#define DENSE_TYPE           CONCAT(dense,T)
#define DENSE_INIT           CONCAT3(dense,T,init)
#define DENSE_DESTRUCT       CONCAT3(dense,T,destruct)
#define DENSE_RESET          CONCAT3(dense,T,reset)
#define DENSE_SCALE          CONCAT3(dense,T,scale)
#define DENSE_MUL_VEC        CONCAT3(dense,T,mul_vec)
#define DENSE_MUL_MAT        CONCAT3(dense,T,mul_mat)
#define DENSE_FACTORIZE      CONCAT3(dense,T,factorize)
#define DENSE_SOLVE_VEC      CONCAT3(dense,T,solve_vec)
#define DENSE_SOLVE_MAT      CONCAT3(dense,T,solve_mat)
#define DENSE_TO_DENSE       CONCAT3(dense,T,to_dense)
#define DENSE_TO_SPARSE      CONCAT3(dense,T,to_sparse)
#ifdef TCMPLX
#define DENSE_TO_REAL        CONCAT3(dense,T,to_real)
#else
#define DENSE_TO_CMPLX       CONCAT3(dense,T,to_cmplx)
#endif
#define DENSE_ADD_DENSE      CONCAT3(dense,T,add_dense)
#define DENSE_ADD_DENSE3     CONCAT3(dense,T,add_dense3)
#define DENSE_ADD_SPARSE     CONCAT3(dense,T,add_sparse)
#define DENSE_ADD_SPARSE3    CONCAT3(dense,T,add_sparse3)
#define DENSE_ADD_BAND       CONCAT3(dense,T,add_band)
#define DENSE_ADD_BAND3      CONCAT3(dense,T,add_band3)
#define DENSE_MUL_DENSE      CONCAT3(dense,T,mul_dense)
#define DENSE_EIG            CONCAT3(dense,T,eig)
#define DENSE_OUTPUT         CONCAT3(dense,T,output)

#define SPARSE_TYPE          CONCAT(sparse,T)
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
#define SPARSE_ADD_SPARSE    CONCAT3(sparse,T,add_sparse)
#define SPARSE_ADD_SPARSE3   CONCAT3(sparse,T,add_sparse3)
#define SPARSE_ADD_BAND      CONCAT3(sparse,T,add_band)
#define SPARSE_ADD_BAND3     CONCAT3(sparse,T,add_band3)
#define SPARSE_MUL_SPARSE    CONCAT3(sparse,T,mul_sparse)
#define SPARSE_GET_CSR       CONCAT3(sparse,T,get_csr)
#define SPARSE_SET_CSR       CONCAT3(sparse,T,set_csr)
#define SPARSE_DESTRUCT_CSR  CONCAT3(sparse,T,destruct_csr)
#define SPARSE_SET_MV_HINT   CONCAT3(sparse,T,set_mv_hint)
#define SPARSE_SET_MM_HINT   CONCAT3(sparse,T,set_mm_hint)
#define SPARSE_OUTPUT        CONCAT3(sparse,T,output)

#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define SPBUILD_INIT         CONCAT3(spbuild,T,init)
#define SPBUILD_DESTRUCT     CONCAT3(spbuild,T,destruct)
#define SPBUILD_LOAD         CONCAT3(spbuild,T,load)
#define SPBUILD_SET          CONCAT3(spbuild,T,set)
#define SPBUILD_SET_ROW      CONCAT3(spbuild,T,set_row)
#define SPBUILD_ADD          CONCAT3(spbuild,T,add)
#define SPBUILD_SORT         CONCAT3(spbuild,T,sort)
#define SPBUILD_FINISH       CONCAT3(spbuild,T,finish)
#define SPBUILD_RESET        CONCAT3(spbuild,T,reset)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_INIT            CONCAT3(band,T,init)
#define BAND_DESTRUCT        CONCAT3(band,T,destruct)
#define BAND_RESET           CONCAT3(band,T,reset)
#define BAND_SCALE           CONCAT3(band,T,scale)
#define BAND_MUL_VEC         CONCAT3(band,T,mul_vec)
#define BAND_MUL_MAT         CONCAT3(band,T,mul_mat)
#define BAND_FACTORIZE       CONCAT3(band,T,factorize)
#define BAND_SOLVE_VEC       CONCAT3(band,T,solve_vec)
#define BAND_SOLVE_MAT       CONCAT3(band,T,solve_mat)
#define BAND_TO_DENSE        CONCAT3(band,T,to_dense)
#define BAND_TO_SPARSE       CONCAT3(band,T,to_sparse)
#ifdef TCMPLX
#define BAND_TO_REAL         CONCAT3(band,T,to_real)
#else
#define BAND_TO_CMPLX        CONCAT3(band,T,to_cmplx)
#endif
#define BAND_ADD_BAND        CONCAT3(band,T,add_band)
#define BAND_ADD_BAND3       CONCAT3(band,T,add_band3)
#define BAND_OUTPUT          CONCAT3(band,T,output)

#define BLOCK_TYPE           CONCAT(block,T)
#define BLOCK_INIT           CONCAT3(block,T,init)
#define BLOCK_DESTRUCT       CONCAT3(block,T,destruct)
#define BLOCK_RESET          CONCAT3(block,T,reset)
#define BLOCK_SCALE          CONCAT3(block,T,scale)
#define BLOCK_MUL_VEC        CONCAT3(block,T,mul_vec)
#define BLOCK_MUL_MAT        CONCAT3(block,T,mul_mat)
#define BLOCK_FACTORIZE      CONCAT3(block,T,factorize)
#define BLOCK_SOLVE_VEC      CONCAT3(block,T,solve_vec)
#define BLOCK_SOLVE_MAT      CONCAT3(block,T,solve_mat)
#define BLOCK_TO_DENSE       CONCAT3(block,T,to_dense)
#define BLOCK_TO_SPARSE      CONCAT3(block,T,to_sparse)
#define BLOCK_GET_GLOBAL_IDX CONCAT3(block,T,get_global_idx)
#define BLOCK_GET_LOCAL_IDX  CONCAT3(block,T,get_local_idx)
#define BLOCK_SET_BLOCK_TYPE CONCAT3(block,T,set_block_type)
#define BLOCK_GET_DENSE      CONCAT3(block,T,get_dense)
#define BLOCK_GET_SPARSE     CONCAT3(block,T,get_sparse)
#define BLOCK_GET_BAND       CONCAT3(block,T,get_band)
#define BLOCK_GET_BLOCK      CONCAT3(block,T,get_block)

subroutine MATRIX_INIT(this, dum, nrows, ncols, symm)
  !! Initialize matrix.

  class(MATRIX_TYPE), intent(out) :: this
    !! Matrix
  character(len=*),   intent(in)  :: dum
    !! Dummy, unused
  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)
  logical, optional,  intent(in)  :: symm
    !! Optional symmetry flag (default: false)

  IGNORE(dum)

  this%nrows = nrows
  this%ncols = nrows
  if (present(ncols)) this%ncols = ncols
  this%symm = .false.
  if (present(symm)) this%symm = symm

  ASSERT(this%nrows > 0)
  ASSERT(this%ncols > 0)
end subroutine

subroutine DENSE_INIT(this, nrows, ncols, d0, symm)
  !! Initialize dense matrix.

  class(DENSE_TYPE), intent(out) :: this
    !! Dense matrix
  integer,           intent(in)  :: nrows
    !! Number of rows
  integer, optional, intent(in)  :: ncols
    !! Optional: Number of columns. (Default: ncols = nrows)
  TT,      optional, intent(in)  :: d0(:,:)
    !! Optional: Provide initial data. (Default: d0 = 0)
  logical, optional,  intent(in)  :: symm
    !! Optional symmetry flag (default: false)

  ! init base
  call this%MATRIX_INIT("", nrows, ncols = ncols, symm = symm)

  ! allocate data
  if (present(d0)) then
    allocate(this%d(this%nrows,this%ncols), source = d0)
  else
    allocate(this%d(this%nrows,this%ncols))
    this%d = 0
  end if
end subroutine

subroutine DENSE_DESTRUCT(this)
  !! Destruct dense matrix; i.e. deallocate its members.

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix

  if (allocated(this%d   )) deallocate(this%d   )
  if (allocated(this%lu  )) deallocate(this%lu  )
  if (allocated(this%ipiv)) deallocate(this%ipiv)
end subroutine

subroutine DENSE_RESET(this, only_factorization)
  !! Reset dense matrix to zero (optionally only the factorization).

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix
  logical, optional, intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  ! local variables
  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  ! reset data
  if (.not. only_fact) then
    this%d = 0
  end if

  ! reset factorization
  if (allocated(this%lu  )) deallocate(this%lu  )
  if (allocated(this%ipiv)) deallocate(this%ipiv)
end subroutine

subroutine DENSE_SCALE(this, fact)
  !! Scale dense matrix by a factor.

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix
  TT,                intent(in)    :: fact
    !! Scaling factor

  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  ! scale data
  this%d = fact * this%d
end subroutine

subroutine DENSE_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(DENSE_TYPE),          intent(in)    :: this
    !! Dense matrix
  TT,                         intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                         intent(inout) :: y(:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! Ff present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  call gemv(this%d, x, y, beta = fact_y, trans = trans)
end subroutine

subroutine DENSE_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(DENSE_TYPE),          intent(in)    :: this
    !! Dense matrix
  TT,                         intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                         intent(inout) :: y(:,:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! Scale existing data in y and add result to it (default: 0)
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag for this matrix, must be 'N', 'C' or 'T' (default: 'N')

  call gemm(this%d, x, y, transa = trans, beta = fact_y)
end subroutine

subroutine DENSE_FACTORIZE(this)
  !! Calculate LU decomposition.
  !! Call this subroutine before solving!

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix

  ! local variables
  integer :: info

  ! allocate memory for factorization
  if (.not. allocated(this%lu  )) allocate(this%lu(  this%nrows,this%ncols))
  if (.not. allocated(this%ipiv)) allocate(this%ipiv(this%nrows))

  ! factorization with lapack
  this%lu = this%d
  call getrf(this%lu, ipiv = this%ipiv, info = info)
  if (info /= 0) call program_error("Dense Factorization error in LAPACK", code = info)
end subroutine

subroutine DENSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for a single right-hand side.
  !! Solves: this * x = rhs, where x and rhs are raw vectors. Call this%factorize() first!

  class(DENSE_TYPE),          intent(in)  :: this
    !! Dense matrix
  TT,                         intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                         intent(out) :: x(:)
    !! Store result in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  TT :: rhs_(size(rhs),1), x_(size(x),1)

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

subroutine DENSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for multiple right-hand sides.
  !! Solves: this * x = rhs, where x and rhs are raw matrices. Call this%factorize() first!

  class(DENSE_TYPE),          intent(in)  :: this
    !! Dense matrix
  TT,                         intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                         intent(out) :: x(:,:)
    !! Store results in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  integer :: info
  logical :: iter_refine_

  ! optional arguments
  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine

  ! solve with lapack
  x = rhs
  call getrs(this%lu, this%ipiv, x, trans = trans, info = info)
  if (info /= 0) call program_error("Dense solving error in LAPACK", code = info)

  ! iterative refinement
  if (iter_refine_) then
    call gerfs(this%d, this%lu, this%ipiv, rhs, x, trans = trans, info = info)
    if (info /= 0) call program_error("Dense iterative refinement error in LAPACK", code = info)
  end if
end subroutine

subroutine DENSE_TO_DENSE(this, d, i0, j0)
  !! Insert dense matrix into another dense matrix.

  class(DENSE_TYPE), intent(in)    :: this
    !! First dense matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Second dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, i1, j1

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! end indices
  i1 = i0_ + this%nrows - 1
  j1 = j0_ + this%ncols - 1

  d%d(i0_:i1, j0_:j1) = this%d
end subroutine

subroutine DENSE_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert dense matrix into a sparse matrix.

  class(DENSE_TYPE),  intent(in)    :: this
    !! Dense matrix (unchanged)
  type(SPBUILD_TYPE), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional,  intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional,  intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional,  intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional,  intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, i, j
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  ASSERT(i0_ > 0)
  ASSERT(i0_ + this%nrows -1 <= sb%sp%nrows)
  j0_ = 1
  if (present(j0)) j0_ = j0
  ASSERT(j0_ > 0)
  ASSERT(j0_ + this%ncols -1 <= sb%sp%ncols)

  ! insert values
  do i = 1, this%nrows
    do j = 1, this%ncols
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,j)) cycle
      elseif (drop_zeros_) then
        if (this%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+j-1, this%d(i,j), search = .false.)
    end do
  end do
end subroutine

#ifdef TCMPLX
subroutine DENSE_TO_REAL(this, d)
  !! Convert dense complex matrix to real by ignoring imaginary part.

  class(dense_cmplx), intent(in)  :: this
    !! Dense complex matrix
  class(dense_real),  intent(out) :: d
    !! Output real(this)

  call d%init(this%nrows, ncols = this%ncols, d0 = real(this%d))
end subroutine
#else
subroutine DENSE_TO_CMPLX(this, d)
  !! Convert dense real matrix to complex with zero imaginary part.

  class(dense_real),  intent(in)  :: this
    !! Dense real matrix
  class(dense_cmplx), intent(out) :: d
    !! Output cmplx(this)

  call d%init(this%nrows, ncols = this%ncols, d0 = cmplx(this%d))
end subroutine
#endif

subroutine DENSE_ADD_DENSE(this, d, fact)
  !! Add two dense matrices, store the result in first matrix.

  class(DENSE_TYPE), intent(inout) :: this
    !! First dense matrix (updated)
  type(DENSE_TYPE),  intent(in)    :: d
    !! Second dense matrix (unchanged)
  TT, optional,      intent(in)    :: fact
    !! Second dense matrix scaling factor (default: 1)

  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  if (present(fact)) then
    this%d = this%d + fact * d%d
  else
    this%d = this%d + d%d
  end if
end subroutine

subroutine DENSE_ADD_DENSE3(this, d2, d3, fact1, fact2)
  !! Add two dense matrices, store the result in third dense matrix.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix (unchanged)
  type(DENSE_TYPE),  intent(in)  :: d2
    !! Dense second matrix (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ! local variables
  TT :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! init output matrix
  call d3%init(this%nrows, ncols = this%ncols, d0 = fact1_ * this%d + fact2_ * d2%d)
end subroutine

subroutine DENSE_ADD_SPARSE(this, s, fact)
  !! Add dense and sparse matrix, store result in dense matrix.

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix (updated)
  type(SPARSE_TYPE), intent(in)    :: s
    !! Sparse matrix (unchanged)
  TT, optional,      intent(in)    :: fact
    !! Sparse matrix scaling factor (default: 1)

  ! local variables
  integer          :: i, j
  TT               :: fact_
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)

  ASSERT(this%nrows == s%nrows)
  ASSERT(this%ncols == s%ncols)
  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! get csr
  call s%get_csr(ia0, ia1, ja, a)

  ! add contents of sparse matrix
  do i = 1, s%nrows
    do j = ia0(i), ia1(i)-1
      this%d(i,ja(j)) = this%d(i,ja(j)) + fact_ * a(j)
    end do
  end do
end subroutine

subroutine DENSE_ADD_SPARSE3(this, s2, d3, fact1, fact2)
  !! Add dense and sparse matrix, store result in another dense matrix.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix (unchanged)
  type(SPARSE_TYPE), intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ! local variables
  TT :: fact1_

  ASSERT(this%nrows == s2%nrows)
  ASSERT(this%ncols == s2%ncols)

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1

  ! init output matrix
  call d3%init(this%nrows, ncols = this%ncols, d0 = fact1_ * this%d)

  ! add sparse matrix to d3
  call d3%add_sparse(s2, fact = fact2)
end subroutine

subroutine DENSE_ADD_BAND(this, b, fact)
  !! Add dense and band matrix, store result in dense matrix.

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix (updated)
  type(BAND_TYPE),   intent(in)    :: b
    !! Band matrix (unchanged)
  TT, optional,      intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  ! local variables
  integer :: i, j, k
  TT      :: fact_

  ASSERT(this%nrows == b%nrows)
  ASSERT(this%ncols == b%ncols)
  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! add contents of band matrix
  do j = 1, b%ncols
    do i = b%i0, b%i1
      ! get actual column
      k = j + i
      if ((k < 1) .or. (k > b%nrows)) cycle

      this%d(k,j) = this%d(k,j) + fact_ * b%d(i,j)
    end do
  end do
end subroutine

subroutine DENSE_ADD_BAND3(this, b2, d3, fact1, fact2)
  !! Add dense and band matrix, store result in a third dense matrix.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix (unchanged)
  type(BAND_TYPE),   intent(in)  :: b2
    !! Band second matrix (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ! local variables
  TT :: fact1_

  ASSERT(this%nrows == b2%nrows)
  ASSERT(this%ncols == b2%ncols)

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1

  ! init output matrix
  call d3%init(this%nrows, ncols = this%ncols, d0 = fact1_ * this%d)

  ! add band matrix to d3
  call d3%add_band(b2, fact = fact2)
end subroutine

subroutine DENSE_MUL_DENSE(this, d2, d3)
  !! Multiply two dense matrices, store result in a third dense matrix.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix (unchanged)
  type(DENSE_TYPE),  intent(in)  :: d2
    !! Dense second matrix (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix

  ! init output matrix
  call d3%init(this%nrows, ncols = d2%ncols)

  ! multiply matrices
  call gemm(this%d, d2%d, d3%d)
end subroutine

subroutine DENSE_EIG(this, e, R, L, sort)
  !! Perform a dense eigenvalue decomposition.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense matrix
  complex,           intent(out) :: e(:)
    !! Output eigenvalues (complex in general)
  complex, optional, intent(out) :: R(:,:)
    !! Output right eigenvectors in columns (complex in general)
  complex, optional, intent(out) :: L(:,:)
    !! Output left  eigenvectors in columns (complex in general)
  logical, optional, intent(in)  :: sort
    !! Sort by magnitude of eigenvalues (default: false)

  ! local variables
  integer :: info
  logical :: sort_
#ifdef TCMPLX
  complex, allocatable :: a(:,:)

  ! allocate memory
  allocate (a(this%nrows,this%nrows), source = this%d)

  ! get eigenvalue decomposition
  call geev(a, e, vl = L, vr = R, info = info)
  if (info /= 0) call program_error("error in complex LAPACK geev", code = info)

#else
  integer           :: i
  real, allocatable :: ar(:,:), wr(:), wi(:), vl(:,:), vr(:,:)

  ! allocate memory
  allocate (ar(this%nrows,this%nrows), source = this%d)
  allocate (wr(this%nrows))
  allocate (wi(this%nrows))

  if (present(L) .and. present(R)) then
    allocate (vl(this%nrows,this%nrows))
    allocate (vr(this%nrows,this%nrows))
    call geev(ar, wr, wi, vl = vl, vr = vr, info = info)
  elseif (present(L)) then
    allocate (vl(this%nrows,this%nrows))
    call geev(ar, wr, wi, vl = vl, info = info)
  elseif (present(R)) then
    allocate (vr(this%nrows,this%nrows))
    call geev(ar, wr, wi, vr = vr, info = info)
  else
    call geev(ar, wr, wi, info = info)
  end if
  if (info /= 0) call program_error("error in real LAPACK geev", code = info)

  ! construct eigenvalues and eigenvectors
  i = 1
  do while (i <= this%nrows)
    if (wi(i) == 0) then
      e(i) = wr(i)
      if (present(L)) L(:,i) = vl(:,i)
      if (present(R)) R(:,i) = vr(:,i)
      i = i + 1
    else
      e(i  ) = cmplx(wr(i  ), wi(i  ))
      e(i+1) = cmplx(wr(i+1), wi(i+1))
      if (present(L)) then
        L(:,i  ) = cmplx(vl(:,i), + vl(:,i+1))
        L(:,i+1) = cmplx(vl(:,i), - vl(:,i+1))
      end if
      if (present(R)) then
        R(:,i  ) = cmplx(vr(:,i), + vr(:,i+1))
        R(:,i+1) = cmplx(vr(:,i), - vr(:,i+1))
      end if
      i = i + 2
    end if
  end do

#endif

  ! sort
  sort_ = .false.
  if (present(sort)) sort_ = sort
  if (sort) then
    block
      real,    allocatable :: abs_e(:)
      integer, allocatable :: perm(:)

      allocate (abs_e(this%nrows), source = abs(e))
      allocate (perm(this%nrows))

      ! sort eigenvalues
      call qsort(abs_e, perm = perm)

      e = e(perm)
      if (present(L)) L = L(:,perm)
      if (present(R)) R = R(:,perm)
    end block
  end if
end subroutine

subroutine DENSE_OUTPUT(this, file, funit, fmt)
  !! Write dense matrix to file or stream.

  class(DENSE_TYPE),          intent(in) :: this
    !! Dense matrix
  character(len=*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,          optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(len=*), optional, intent(in) :: fmt
    !! Value format string

  ! local variables
  integer                       :: i, j, funit_, status
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

  do i = 1, this%nrows
    do j = 1, this%ncols
      if (j > 1) then
        write(funit_, fmt_, advance = "no") ", ", this%d(i,j)
      else
        write(funit_, fmt_, advance = "no") "", this%d(i,j)
      end if
    end do

    ! next line
    write(funit_,*)
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

subroutine SPARSE_INIT(this, nrows, ncols, symm)
  !! Initialize sparse matrix.

  class(SPARSE_TYPE), intent(out) :: this
    !! Sparse matrix
  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional: Number of columns (default: ncols = nrows)
  logical, optional,  intent(in)  :: symm
    !! Optional symmetry flag (default: false)

  ! init base
  call this%MATRIX_INIT("", nrows, ncols = ncols, symm = symm)

  ! empty mkl handle
  this%mkl_spmat%ptr = 0

  ! set matrix descriptor type
  this%mkl_descr%type = SPARSE_MATRIX_TYPE_GENERAL
end subroutine

subroutine SPARSE_DESTRUCT(this)
  !! Destruct sparse matrix by deallocating memory and deleting factorization.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! free mkl handle
  call this%destruct_csr()

  ! destruct pardiso handle
  call this%pard%destruct()
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

  ! free mkl handle
  if (.not. only_fact) then
    call this%destruct_csr()
  end if

  ! destruct factorization
  call this%pard%destruct()
end subroutine

subroutine SPARSE_SCALE(this, fact)
  !! Scale sparse matrix by a factor.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  TT,                 intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer          :: i
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)

  ASSERT(.not. this%pard%factorized)

  ! get csr
  call this%get_csr(ia0, ia1, ja, a)

  ! scale values
  do i = 1, this%nrows
    a(ia0(i):ia1(i)-1) = fact * a(ia0(i):ia1(i)-1)
  end do
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
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  TT                  :: fact_y_
  integer             :: info
  integer(kind=c_int) :: mkl_trans

  ! y scaling factor
  fact_y_ = 0
  if (present(fact_y)) fact_y_ = fact_y

  ! matrix operation
  mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
  if (present(trans)) then
    if (trans == 'N') then
      mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
    elseif (trans == 'T') then
      mkl_trans = SPARSE_OPERATION_TRANSPOSE
    elseif (trans == 'C') then
      mkl_trans = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
    else
      print *, trans
      call program_error("trans must be 'N', 'T' or 'C'")
    end if
  end if

#ifdef TCMPLX
  info = mkl_sparse_z_mv(mkl_trans, (1.0,0.0), this%mkl_spmat, this%mkl_descr, x, fact_y_, y)
#else
  info = mkl_sparse_d_mv(mkl_trans, 1.0, this%mkl_spmat, this%mkl_descr, x, fact_y_, y)
#endif
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix vector multiplication error", code = info)
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
    !! Transposed flag for this matrix, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  TT                  :: fact_y_
  integer             :: info
  integer(kind=c_int) :: mkl_trans

  ! y scaling factor
  fact_y_ = 0
  if (present(fact_y)) fact_y_ = fact_y

  ! matrix operation
  mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
  if (present(trans)) then
    if (trans == 'N') then
      mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
    elseif (trans == 'T') then
      mkl_trans = SPARSE_OPERATION_TRANSPOSE
    elseif (trans == 'C') then
      mkl_trans = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
    else
      print *, trans
      call program_error("trans must be 'N', 'T' or 'C'")
    end if
  end if

#ifdef TCMPLX
  info = mkl_sparse_z_mm(mkl_trans, (1.0,0.0), this%mkl_spmat, this%mkl_descr, SPARSE_LAYOUT_COLUMN_MAJOR, &
                         x, size(y,2), size(x,1), fact_y_, y, size(y,1))
#else
  info = mkl_sparse_d_mm(mkl_trans, 1.0, this%mkl_spmat, this%mkl_descr, SPARSE_LAYOUT_COLUMN_MAJOR, &
                         x, size(y,2), size(x,1), fact_y_, y, size(y,1))
#endif
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix matrix multiplication error", code = info)
end subroutine

subroutine SPARSE_FACTORIZE(this)
  !! Factorize Matrix using PARDISO.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! local variables
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)

  ASSERT(.not. this%pard%factorized)

#ifdef TCMPLX
  call this%pard%init(this%nrows, .true.)
#else
  call this%pard%init(this%nrows, .false.)
#endif

  ! get csr
  call this%get_csr(ia0, ia1, ja, a)

  ! factorize with pardiso
  call this%pard%factorize(ia0, ia1, ja, a)
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
  call this%pard%solve(rhs, x)
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
  call this%pard%solve(rhs_, x_)
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
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! get csr
  call this%get_csr(ia0, ia1, ja, a)

  ! insert values into dense matrix
  do i = 1, this%nrows
    do j = ia0(i), ia1(i)-1
      d%d(i0_+i-1,j0_+ja(j)-1) = a(j)
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
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)
  logical          :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! get csr
  call this%get_csr(ia0, ia1, ja, a)

  ! insert values into second sparse matrix
  do i = 1, this%nrows
    do j = ia0(i), ia1(i)-1
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,ja(j))) cycle
      elseif (drop_zeros_) then
        if (a(j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+ja(j)-1, a(j), search = .false.)
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

  ! local variables
  integer, pointer     :: ia0(:), ia1(:), ja(:)
  complex, pointer     :: ac(:)
  real,    allocatable :: a(:)

  ! get csr
  call this%get_csr(ia0, ia1, ja, ac)

  ! convert values to real
  allocate (a(size(ac)), source = real(ac))

  ! init s
  call s%init(this%nrows, ncols = this%ncols)
  call s%set_csr(ia0, ia1, ja, a)
end subroutine
#else
subroutine SPARSE_TO_CMPLX(this, s)
  !! Convert sparse real matrix to complex

  class(sparse_real),  intent(in)  :: this
    !! Sparse real matrix
  class(sparse_cmplx), intent(out) :: s
    !! Output cmplx(this)

  ! local variables
  integer, pointer     :: ia0(:), ia1(:), ja(:)
  real,    pointer     :: a(:)
  complex, allocatable :: ac(:)

  ! get csr
  call this%get_csr(ia0, ia1, ja, a)

  ! convert values to complex
  allocate (ac(size(a)), source = cmplx(a))

  ! init s
  call s%init(this%nrows, ncols = this%ncols)
  call s%set_csr(ia0, ia1, ja, ac)
end subroutine
#endif

subroutine SPARSE_ADD_SPARSE(this, s, fact)
    !! Add two sparse matrices, store result in first matrix.

  class(SPARSE_TYPE), intent(inout) :: this
    !! First sparse matrix (updated)
  type(SPARSE_TYPE),  intent(in)    :: s
    !! Second sparse matrix (unchanged)
  TT, optional,       intent(in)    :: fact
    !! Second sparse matrix scaling factor (default: 1)

  ! local variables
  integer               :: info
  TT                    :: fact_
  type(sparse_matrix_t) :: mkl_c

  ASSERT(.not. this%pard%factorized)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

#ifdef TCMPLX
  info = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, s%mkl_spmat, fact_, this%mkl_spmat, mkl_c)
#else
  info = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, s%mkl_spmat, fact_, this%mkl_spmat, mkl_c)
#endif
  if (info /= 0) call program_error("MKL sparse matrix addition error", code = info)

  ! delete old matrix
  call this%destruct_csr()

  ! save result
  this%mkl_spmat = mkl_c
end subroutine

subroutine SPARSE_ADD_SPARSE3(this, s2, s3, fact1, fact2)
  !! Add two sparse matrices, store result in sparse third matrix.

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

  ! local variables
  integer :: info

  ! init s3
  call s3%init(this%nrows, ncols = this%ncols)

  ! copy this matrix to s3
  info = mkl_sparse_copy(this%mkl_spmat, this%mkl_descr, s3%mkl_spmat)
  if (info /= 0) call program_error("MKL sparse matrix copy error", code = info)

  ! scale
  if (present(fact1)) then
    call s3%scale(fact1)
  end if

  ! add s2 to s3
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
  integer               :: i, j, k, n, info
  TT                    :: fact_
  integer, allocatable  :: ia0(:), ia1(:), ja(:)
  TT,      allocatable  :: a(:)
  type(sparse_matrix_t) :: mkl_b, mkl_c

  ASSERT(.not. this%pard%factorized)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! total number of elements in band matrix
  n = b%ndiags * b%nrows - (b%ndiags * b%ndiags - 1) / 4

  ! convert band matrix to csr
  allocate (ia0(b%nrows))
  allocate (ia1(b%nrows))
  allocate (ja(n))
  allocate (a( n))
  k = 0
  do i = 1, b%nrows
    ia0(i) = k + 1
    do j = i + b%i0, i + b%i1
      if ((j < 1) .or. (j > b%ncols)) cycle
      k     = k + 1
      ja(k) = j
      a( k) = b%d(i-j,j)
    end do
    ia1(i) = k + 1
  end do

  ! create MKL handle for band matrix
#ifdef TCMPLX
  info = mkl_sparse_z_create_csr(mkl_b, SPARSE_INDEX_BASE_ONE, b%nrows, b%ncols, ia0, ia1, ja, a)
#else
  info = mkl_sparse_d_create_csr(mkl_b, SPARSE_INDEX_BASE_ONE, b%nrows, b%ncols, ia0, ia1, ja, a)
#endif
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix create error", code = info)

  ! add this and b
#ifdef TCMPLX
  info = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, mkl_b, fact_, this%mkl_spmat, mkl_c)
#else
  info = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, mkl_b, fact_, this%mkl_spmat, mkl_c)
#endif
  if (info /= 0) call program_error("MKL sparse matrix addition error", code = info)

  ! delete old matrix and b MKL handle
  call this%destruct_csr()
  info = mkl_sparse_destroy(mkl_b)
  if (info /= 0) call program_error("MKL sparse matrix destroy error", code = info)

  ! save result
  this%mkl_spmat = mkl_c
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

  ! local variables
  integer :: info

  ! init s3
  call s3%init(this%nrows, ncols = this%ncols)

  ! copy this matrix to s3
  info = mkl_sparse_copy(this%mkl_spmat, this%mkl_descr, s3%mkl_spmat)
  if (info /= 0) call program_error("MKL sparse matrix copy error", code = info)

  ! scale
  if (present(fact1)) then
    call s3%scale(fact1)
  end if

  ! add b2 to s3
  call s3%add_band(b2, fact = fact2)
end subroutine

subroutine SPARSE_MUL_SPARSE(this, s2, s3)
  !! Multiply two sparse matrices, store result in sparse third matrix.

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse first matrix (unchanged)
  type(SPARSE_TYPE),  intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE),  intent(out) :: s3
    !! Output sparse third matrix

  ! local variables
  integer :: info

  ! init s3
  call s3%init(this%nrows, ncols = s2%ncols)

  ! s3 = this * s2
  info = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, this%mkl_spmat, s2%mkl_spmat, s3%mkl_spmat)
  if (info /= 0) call program_error("MKL sparse matrix multiplication error", code = info)
end subroutine

subroutine SPARSE_GET_CSR(this, ia0, ia1, ja, a)
  !! Get pointers to CSR.

  class(SPARSE_TYPE), intent(in)  :: this
    !! Sparse matrix
  integer, pointer,   intent(out) :: ia0(:)
    !! Output pointer to row start pointers
  integer, pointer,   intent(out) :: ia1(:)
    !! Output pointer to row end pointers (exclusive)
  integer, pointer,   intent(out) :: ja(:)
    !! Output pointer to column indices
  TT,      pointer,   intent(out) :: a(:)
    !! Output pointer to values

  ! local variables
  integer             :: info
  integer(kind=c_int) :: nrows, ncols, indexing

  ! init pointers
  nullify(ia0, ia1, ja, a)

  ! return null pointers if empty
  if (this%mkl_spmat%ptr == 0) return

#ifdef TCMPLX
  info = mkl_sparse_z_export_csr(this%mkl_spmat, indexing, nrows, ncols, ia0, ia1, ja, a)
#else
  info = mkl_sparse_d_export_csr(this%mkl_spmat, indexing, nrows, ncols, ia0, ia1, ja, a)
#endif
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix get CSR error", code = info)

  ! FIXME: Remove if not necessary
  ASSERT(indexing == SPARSE_INDEX_BASE_ONE)
end subroutine

subroutine SPARSE_SET_CSR(this, ia0, ia1, ja, a)
  !! Set CSR data.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix
  integer,            intent(in)    :: ia0(:)
    !! Row start pointers
  integer,            intent(in)    :: ia1(:)
    !! Row end pointers (exclusive)
  integer,            intent(in)    :: ja(:)
    !! Column indices
  TT,                 intent(in)    :: a(:)
    !! Values

  ! local variables
  integer :: info

  ! delete old values
  call this%destruct_csr()

  ! create MKL handle
#ifdef TCMPLX
  info = mkl_sparse_z_create_csr(this%mkl_spmat, SPARSE_INDEX_BASE_ONE, this%nrows, this%ncols, ia0, ia1, ja, a)
#else
  info = mkl_sparse_d_create_csr(this%mkl_spmat, SPARSE_INDEX_BASE_ONE, this%nrows, this%ncols, ia0, ia1, ja, a)
#endif
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix create error", code = info)
end subroutine

subroutine SPARSE_DESTRUCT_CSR(this)
  !! Free memory of MKL handle.

  class(SPARSE_TYPE), intent(inout) :: this
    !! Sparse matrix

  ! local variables
  integer :: info

  if (this%mkl_spmat%ptr /= 0) then
    info = mkl_sparse_destroy(this%mkl_spmat)
    if (info /= 0) call program_error("MKL sparse matrix destroy error", code = info)
    this%mkl_spmat%ptr = 0 ! FIXME: check if unnecessary
  end if
end subroutine

subroutine SPARSE_SET_MV_HINT(this, expected_calls, trans)
  !! provide hint to MKL about matrix vector multiplications

  class(SPARSE_TYPE),         intent(inout) :: this
    !! Sparse matrix
  integer,                    intent(in)    :: expected_calls
    !! Expected number of matrix vector multiplications
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer             :: info
  integer(kind=c_int) :: mkl_trans

  ! matrix operation
  mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
  if (present(trans)) then
    if (trans == 'N') then
      mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
    elseif (trans == 'T') then
      mkl_trans = SPARSE_OPERATION_TRANSPOSE
    elseif (trans == 'C') then
      mkl_trans = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
    else
      print *, trans
      call program_error("trans must be 'N', 'T' or 'C'")
    end if
  end if

  ! provide hint
  info = mkl_sparse_set_mv_hint(this%mkl_spmat, mkl_trans, this%mkl_descr, expected_calls)
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix set mv hint error", code = info)

  ! optimize sparse matrix
  info = mkl_sparse_optimize(this%mkl_spmat)
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix optimize error", code = info)
end subroutine

subroutine SPARSE_SET_MM_HINT(this, expected_calls, ncols, trans)
  !! Provide hint to MKL about matrix matrix multiplications

  class(SPARSE_TYPE),         intent(inout) :: this
    !! Sparse matrix
  integer,                    intent(in)    :: expected_calls
    !! Expected number of matrix vector multiplications
  integer,                    intent(in)    :: ncols
    !! Number of columns of second matrix
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer             :: info
  integer(kind=c_int) :: mkl_trans

  ! matrix operation
  mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
  if (present(trans)) then
    if (trans == 'N') then
      mkl_trans = SPARSE_OPERATION_NON_TRANSPOSE
    elseif (trans == 'T') then
      mkl_trans = SPARSE_OPERATION_TRANSPOSE
    elseif (trans == 'C') then
      mkl_trans = SPARSE_OPERATION_CONJUGATE_TRANSPOSE
    else
      print *, trans
      call program_error("trans must be 'N', 'T' or 'C'")
    end if
  end if

  ! provide hint
  info = mkl_sparse_set_mm_hint(this%mkl_spmat, mkl_trans, this%mkl_descr, SPARSE_LAYOUT_COLUMN_MAJOR, &
                                ncols, expected_calls)
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix set mm hint error", code = info)

  ! optimize sparse matrix
  info = mkl_sparse_optimize(this%mkl_spmat)
  if (info /= SPARSE_STATUS_SUCCESS) call program_error("MKL sparse matrix optimize error", code = info)
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
  integer, pointer              :: ia0(:), ia1(:), ja(:)
  TT,      pointer              :: a(:)
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
  call this%get_csr(ia0, ia1, ja, a)
  do i = 1, this%nrows
    do j = ia0(i), ia1(i)-1
      ! write row and column
      write(funit_, "(1I0,1A,1I0,1A)", advance = "no") i, " ", ja(j), " "

      ! write value
      write(funit_, fmt_) a(j)
    end do
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

subroutine SPBUILD_INIT(this, s)
  !! Initialize sparse matrix builder.

  class(SPBUILD_TYPE),       intent(out) :: this
    !! Sparse matrix builder
  type(SPARSE_TYPE), target, intent(in)  :: s
    !! Sparse matrix this builder is tied to

  ! local variables
  integer :: i

  ! save pointer to matrix
  this%sp => s

  ! allocate data
  allocate(this%cols(  s%nrows))
  allocate(this%vals(  s%nrows))
  allocate(this%filled(s%nrows))
  do i = 1, s%nrows
    call this%cols(i)%init(n = 0, c = 8)
    call this%vals(i)%init(n = 0, c = 8)
    call this%filled(i)%init(n = 0, c = 8)
  end do

  ! allocate row flags
  allocate(this%sorted(s%nrows), source = .true.)
  allocate(this%keep_struct(s%nrows), source = .true.)

  ! load existing data from sparse matrix
  call this%load()

  ! not finished
  this%finished = .false.
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

subroutine SPBUILD_LOAD(this)
  !! Load existing values from sparse matrix into builder.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  ! local variables
  integer          :: i, j
  integer, pointer :: ia0(:), ia1(:), ja(:)
  TT,      pointer :: a(:)

  ! make sure reset was called after finish
  ASSERT(.not. this%finished)

  ! get csr representation of matrix (null() if empty)
  call this%sp%get_csr(ia0, ia1, ja, a)

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
    if (associated(ia0)) then
      do j = ia0(i), ia1(i)-1
        call this%cols(i)%push(ja(j))
        call this%vals(i)%push(a(j))
        call this%filled(i)%push(.true.)
      end do
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

  ASSERT(.not. this%finished)
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

  ! make sure reset was called after finish
  ASSERT(.not. this%finished)

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

  ! make sure reset was called after finish
  ASSERT(.not. this%finished)

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

  ! make sure reset was called after finish
  ASSERT(.not. this%finished)

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

subroutine SPBUILD_FINISH(this)
  !! Create sparse matrix.

  class(SPBUILD_TYPE), intent(inout) :: this
    !! Sparse matrix builder

  ! local variables
  integer              :: i, j
  integer, allocatable :: ia0(:), ia1(:)
  type(vector_int)     :: ja
  type(VECTOR_TYPE)    :: a

  ! sort column indices
  call this%sort()

  ! allocate csr
  allocate (ia0(this%sp%nrows))
  allocate (ia1(this%sp%nrows))
  call ja%init(n = 0, c = sum(this%cols%n))
  call a%init( n = 0, c = sum(this%cols%n))

  ! init csr
  do i = 1, this%sp%nrows
    ! row start pointer
    ia0(i) = ja%n + 1

    ! add values to csr (consider only filled elements)
    do j = 1, this%cols(i)%n
      if (.not. this%filled(i)%d(j)) cycle
      call ja%push(this%cols(i)%d(j))
      call a%push( this%vals(i)%d(j))
    end do

    ! row end pointer
    ia1(i) = ja%n + 1
  end do

  ! save csr in sparse matrix
  call this%sp%set_csr(ia0, ia1, ja%d(1:ja%n), a%d(1:a%n))

  ! set flag to indicate csr was built
  this%finished = .true.
end subroutine SPBUILD_FINISH

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

  ! not finished, enable rebuilding
  this%finished = .false.
end subroutine

subroutine BAND_INIT(this, ncols, ndiags, d0, symm)
  !! Initialize band matrix

  class(BAND_TYPE),  intent(out) :: this
    !! Band matrix
  integer,           intent(in)  :: ncols
    !! Number of columns = Number of rows
  integer,           intent(in)  :: ndiags
    !! Total number of diagonals. Must be positive uneven integer.
  TT, optional,      intent(in)  :: d0(-(ndiags-1)/2:,:)
    !! Initial data.
  logical, optional, intent(in)  :: symm
    !! Optional symmetry flag (default: false)

  ! init base
  call this%MATRIX_INIT("", ncols, ncols = ncols, symm = symm)

  ! set number of diagonals
  ASSERT(ndiags > 0)
  this%ndiags = ndiags

  ! min and max rows
  this%i0 = -(ndiags - 1) / 2
  this%i1 = +(ndiags - 1) / 2
  ASSERT(this%i1 < ncols)

  ! allocate data
  if (present(d0)) then
    allocate (this%d(this%i0:this%i1,1:ncols), source = d0)
  else
    allocate (this%d(this%i0:this%i1,1:ncols))
    this%d = 0
  end if
end subroutine

subroutine BAND_DESTRUCT(this)
  !! Destruct band matrix by deallocating its memory.

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix

  if (allocated(this%d   )) deallocate(this%d   )
  if (allocated(this%f   )) deallocate(this%f   )
  if (allocated(this%ipiv)) deallocate(this%ipiv)
end subroutine

subroutine BAND_RESET(this, only_factorization)
  !! Reset band matrix to zero (optionally only the factorization).

  class(BAND_TYPE),  intent(inout) :: this
    !! Band matrix
  logical, optional, intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  ! local variables
  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  if (.not. only_fact) then
    this%d = 0
  end if

  ! reset factorization
  if (allocated(this%f   )) deallocate(this%f   )
  if (allocated(this%ipiv)) deallocate(this%ipiv)
end subroutine

subroutine BAND_SCALE(this, fact)
  !! Scale band matrix by a factor.

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix
  TT,               intent(in)    :: fact
    !! Scaling factor

  ASSERT(.not. allocated(this%f))
  ASSERT(.not. allocated(this%ipiv))

  this%d = this%d * fact
end subroutine

subroutine BAND_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply band matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(BAND_TYPE),           intent(in)    :: this
    !! Band matrix
  TT,                         intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                         intent(inout) :: y(:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  if (this%ndiags == 1) then
    if (present(fact_y)) then
      y = this%d(:,0) * x + fact_y * y
    else
      y = this%d(:,0) * x
    end if
  else
    call gbmv(this%d(this%i0:this%i1,:), x, y, beta = fact_y, trans = trans)
  end if
end subroutine

subroutine BAND_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply band matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(BAND_TYPE),           intent(in)    :: this
    !! Band matrix
  TT,                         intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                         intent(inout) :: y(:,:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer :: i

  if (this%ndiags == 1) then
    if (present(fact_y)) then
      do i = 1, size(x, 2)
        y(:,i) = this%d(:,0) * x(:,i) + fact_y * y(:,i)
      end do
    else
      do i = 1, size(x, 2)
        y(:,i) = this%d(:,0) * x(:,i)
      end do
    end if
  else
    do i = 1, size(x, 2)
      call gbmv(this%d(this%i0:this%i1,:), x(:,i), y(:,i), beta = fact_y, trans = trans)
    end do
  end if
end subroutine BAND_MUL_MAT

subroutine BAND_FACTORIZE(this)
  !! Factorize band matrix. Call this before solving!

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix

  ! local variables
  integer :: info

  if (this%ndiags == 1) then ! diagonal
    if (.not. allocated(this%f)) allocate(this%f(0:0,this%ncols))
    this%f(0,:) = 1 / this%d(0,:)
  elseif (this%ndiags == 3) then ! tridiagonal
    if (.not. allocated(this%f)) allocate(this%f(-2:+1,this%ncols))
    if (.not. allocated(this%ipiv)) allocate(this%ipiv(this%ncols))

    associate (d => this%d, f => this%f, n => this%ncols)
      f = 0
      f(-1,1:n-1) = d(-1,2:n  )
      f( 0,1:n  ) = d( 0,1:n  )
      f(+1,1:n-1) = d(+1,1:n-1)

      ! factorization
      call gttrf(f(+1,1:n-1), f(0,1:n), f(-1,1:n-1), f(-2,3:n), ipiv = this%ipiv, info = info)
      if (info /= 0) call program_error("Tridiagonal Factorization error in LAPACK", code = info)
    end associate
  else ! general band matrix
    if (.not. allocated(this%f)) allocate(this%f((3 * this%ndiags - 1) / 2, this%ncols))
    if (.not. allocated(this%ipiv)) allocate(this%ipiv(this%nrows))

    ! copy data to f
    this%f = this%d(this%i0:this%i1,:)

    ! factorization
    call gbtrf(this%f, ipiv = this%ipiv, info = info)
    if (info /= 0) call program_error("Bandmatrix Factorization error in LAPACK", code = info)
  end if
end subroutine

subroutine BAND_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first!

  class(BAND_TYPE),           intent(in)  :: this
    !! Band matrix.
  TT,                         intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                         intent(out) :: x(:)
    !! Store result in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  TT :: rhs_(size(rhs),1), x_(size(x),1)

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

subroutine BAND_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(BAND_TYPE),           intent(in)  :: this
    !! Band matrix
  TT,                         intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                         intent(out) :: x(:,:)
    !! Store results in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  integer :: i, info
  logical :: iter_refine_

  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine

  if (this%ndiags == 1) then ! diagonal
    do i = 1, size(x, 2)
      x(:,i) = this%f(0,:) * rhs(:,i)
    end do
  elseif (this%ndiags == 3) then ! tridiagonal
    associate (d => this%d, f => this%f, n => this%ncols)
      x = rhs
      call gttrs(f(+1,1:n-1), f(0,1:n), f(-1,1:n-1), f(-2,3:n), x, this%ipiv, trans = trans, info = info)
      if (info /= 0) call program_error("Tridiagonal Solution error in LAPACK", code = info)

      if (iter_refine_) then
        call gtrfs(d(+1,1:n-1), d(0,1:n), d(-1,2:n), f(+1,1:n-1), f(0,1:n), f(-1,2:n), f(-2,3:n), this%ipiv, rhs, x, trans = trans, info = info)
        if (info /= 0) call program_error("Tridiagonal Iterative Refinement error in LAPACK", code = info)
      end if
    end associate
  else ! general band matrix
    x = rhs
    call gbtrs(this%f, x, this%ipiv, trans = trans, info = info)
    if (info /= 0) call program_error("Band Matrix Solution error in LAPACK", code = info)

    if (iter_refine_) then
      call gbrfs(this%d(this%i0:this%i1,:), this%f, this%ipiv, rhs, x, trans = trans, info = info)
      if (info /= 0) call program_error("Band Matrix Iterative Refinement error in LAPACK", code = info)
    end if
  end if
end subroutine

subroutine BAND_TO_DENSE(this, d, i0, j0)
  !! Insert band matrix into a dense matrix.

  class(BAND_TYPE),  intent(in)    :: this
    !! Band matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, i, j, k

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements into dense matrix
  do j = 1, this%ncols
    do i = this%i0, this%i1
      k = j + i ! actual row
      if ((k < 1) .or. (k > this%nrows)) cycle

      d%d(i0_+k-1,j0_+j-1) = this%d(i,j)
    end do
  end do
end subroutine

subroutine BAND_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert band matrix into a sparse matrix.

  class(BAND_TYPE),   intent(in)    :: this
    !! Band matrix (unchanged)
  type(SPBUILD_TYPE), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional,  intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional,  intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional,  intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional,  intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, i, j
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values
  do i = 1, this%nrows
    do j = i + this%i0, i + this%i1
      if ((j < 1) .or. (j > this%ncols)) cycle

      ! check structure
      if (present(struct)) then
        if (.not. struct(i,j)) cycle
      elseif (drop_zeros_) then
        if (this%d(i-j,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+j-1, this%d(i-j,j), search = .false.)
    end do
  end do
end subroutine

#ifdef TCMPLX
subroutine BAND_TO_REAL(this, b)
  !! Convert band complex matrix to real by ignoring imaginary part.

  class(band_cmplx), intent(in)  :: this
    !! Band complex matrix
  class(band_real),  intent(out) :: b
    !! Output real(this)

  call b%init(this%nrows, this%ndiags, d0 = real(this%d))
end subroutine
#else
subroutine BAND_TO_CMPLX(this, b)
  !! Convert band real matrix to complex.

  class(band_real), intent(in)  :: this
    !! Band real matrix
  class(band_cmplx),  intent(out) :: b
    !! Output cmplx(this)

  call b%init(this%nrows, this%ndiags, d0 = cmplx(this%d))
end subroutine
#endif

subroutine BAND_ADD_BAND(this, b, fact)
  !! Add two band matrices, store result in first band matrix.

  class(BAND_TYPE), intent(inout) :: this
    !! First band matrix (updated)
  type(BAND_TYPE),  intent(in)    :: b
    !! Second band matrix (unchanged)
  TT, optional,     intent(in)    :: fact
    !! Second band matrix scaling factor (default: 1)

  ! local variables
  integer         :: i
  TT              :: fact_
  TT, allocatable :: tmp(:,:)

  ASSERT(.not. allocated(this%f))
  ASSERT(.not. allocated(this%ipiv))

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  if (b%ndiags > this%ndiags) then
    allocate(tmp(b%i0:b%i1,1:this%ncols))
    tmp = 0
    do i = this%i0, this%i1
      tmp(i,:) = this%d(i,:)
    end do
    tmp = tmp + fact_ * b%d
    call move_alloc(tmp, this%d)

    this%ndiags = b%ndiags
    this%i0 = b%i0
    this%i1 = b%i1
  else
    this%d(b%i0:b%i1,1:b%ncols) = this%d(b%i0:b%i1,1:b%ncols) + fact_ * b%d
  end if
end subroutine

subroutine BAND_ADD_BAND3(this, b2, b3, fact1, fact2)
  !! Add two band matrices, store result in third band matrix.

  class(BAND_TYPE), intent(in)  :: this
    !! Band first matrix (unchanged)
  type(BAND_TYPE),  intent(in)  :: b2
    !! Band second matrix (unchanged)
  type(BAND_TYPE),  intent(out) :: b3
    !! Output band third matrix
  TT, optional,     intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,     intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ! local variables
  TT :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  if (this%ndiags >= b2%ndiags) then
    call b3%init(this%nrows, this%ndiags, d0 = fact1_ * this%d)
    b3%d(b2%i0:b2%i1,1:b2%ncols) = b3%d(b2%i0:b2%i1,1:b2%ncols) + fact2_ * b2%d
  else
    call b3%init(b2%nrows, b2%ndiags, d0 = fact2_ * b2%d)
    b3%d(this%i0:this%i1,1:this%ncols) = b3%d(this%i0:this%i1,1:this%ncols) + fact1_ * this%d
  end if
end subroutine

subroutine BAND_OUTPUT(this, file, funit, fmt)
  !! Write band matrix to file or stream.

  class(BAND_TYPE),           intent(in) :: this
    !! Band matrix
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

  ! output values
  do i = 1, this%nrows
    do j = i + this%i0, i + this%i1
      if ((j < 1) .or. (j > this%ncols)) cycle

      ! write row and column
      write(funit_, "(1I0,1A,1I0,1A)", advance = "no") i, " ", j, " "

      ! write value
      write(funit_, fmt_) this%d(i-j,j)
    end do
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

subroutine BLOCK_INIT(this, row_dim, col_dim, diag, tridiag, symm)
  !! Initialize Block matrix.

  class(BLOCK_TYPE), intent(out) :: this
    !! Block matrix
  integer,           intent(in)  :: row_dim(:)
    !! Number of Rows for each block
  integer, optional, intent(in)  :: col_dim(:)
    !! Number of Columns for each block
  logical, optional, intent(in)  :: diag
    !! diagonal block matrix flag (default: false)
  logical, optional, intent(in)  :: tridiag
    !! tridiagonal block matrix flag (default: false)
  logical, optional, intent(in)  :: symm
    !! Optional symmetry flag (default: false)

  ! local variables
  integer, allocatable :: col_dim_(:)
  integer              :: ib, i1, jb, j1

  ! init base
  if (present(col_dim)) then
    call this%MATRIX_INIT("", sum(row_dim), ncols = sum(col_dim), symm = symm)
    allocate (col_dim_(size(col_dim)), source = col_dim)
  else
    call this%MATRIX_INIT("", sum(row_dim), symm = symm)
    allocate (col_dim_(size(row_dim)), source = row_dim)
  end if

  ! number of blocks
  this%nbrows = size(row_dim)
  this%nbcols = size(col_dim_)

  ! set row indices
  allocate (this%i0(this%nbrows))
  allocate (this%i1(this%nbrows))
  i1 = 0
  do ib = 1, this%nbrows
    this%i0(ib) = i1 + 1
    this%i1(ib) = i1 + row_dim(ib)
    i1 = this%i1(ib)
  end do

  ! set col indices
  allocate (this%j0(this%ncols))
  allocate (this%j1(this%ncols))
  j1 = 0
  do jb = 1, this%ncols
    this%j0(jb) = j1 + 1
    this%j1(jb) = j1 + row_dim(jb)
    j1 = this%j1(jb)
  end do

  ! allocate data
  allocate (this%b(this%nbrows,this%nbcols))

  ! diag, tridiag
  this%diag = .false.
  if (present(diag)) this%diag = diag
  this%tridiag = .false.
  if (present(tridiag)) this%tridiag = tridiag

  if (this%diag .and. this%tridiag) call program_error("Block matrix can not be both diagonal and tridiagonal!")
  if (this%diag .or. this%tridiag) then
    if (this%nbrows /= this%nbcols) call program_error("(Tri)Diagonal Block matrix must have same number of blocks per row and column!")
  end if

  ! not factorized
  this%factorized = .false.
end subroutine

subroutine BLOCK_DESTRUCT(this)
  !! Destruct block matrix.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix

  ! local variables
  integer :: ib, jb

  if (allocated(this%i0)) deallocate(this%i0)
  if (allocated(this%i1)) deallocate(this%i1)
  if (allocated(this%j0)) deallocate(this%j0)
  if (allocated(this%j1)) deallocate(this%j1)

  if (allocated(this%b)) then
    do ib = 1, this%nbrows
      do jb = 1, this%nbcols
        if (.not. allocated(this%b(ib,jb)%p)) cycle

        call this%b(ib,jb)%p%destruct()
      end do
    end do
    deallocate(this%b)
  end if
end subroutine

subroutine BLOCK_RESET(this, only_factorization)
  !! Reset block matrix to zero (optionally only factorization).

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix
  logical, optional, intent(in)    :: only_factorization
    !! reset only factorization data (default: false)

  ! local variables
  integer :: ib, jb

  do ib = 1, this%nbrows
    do jb = 1, this%nbcols
      if (.not. allocated(this%b(ib,jb)%p)) cycle
      call this%b(ib,jb)%p%reset(only_factorization = only_factorization)
    end do
  end do
end subroutine

subroutine BLOCK_SCALE(this, fact)
  !! Scale block matrix by factor.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix
  TT,                intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: ib, jb

  ASSERT(.not. this%factorized)

  do ib = 1, this%nbrows
    do jb = 1, this%nbcols
      if (.not. allocated(this%b(ib,jb)%p)) cycle
      call this%b(ib,jb)%p%scale(fact)
    end do
  end do
end subroutine

subroutine BLOCK_MUL_VEC(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(BLOCK_TYPE),          intent(in)    :: this
    !! Block matrix
  TT,                         intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                         intent(inout) :: y(:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  character(len=1) :: trans_
  TT               :: one
  integer          :: ib, jb

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = 'N'
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ .eq. 'N') then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. allocated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_vec(x(j0(ib):j1(ib)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. allocated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_vec(x(i0(jb):i1(jb)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

subroutine BLOCK_MUL_MAT(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(BLOCK_TYPE),          intent(in)    :: this
    !! Block matrix
  TT,                         intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                         intent(inout) :: y(:,:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  character(len=1) :: trans_
  TT               :: one
  integer          :: ib, jb

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = 'N'
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ .eq. 'N') then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. allocated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_mat(x(j0(ib):j1(ib),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. allocated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_mat(x(i0(jb):i1(jb),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

subroutine BLOCK_FACTORIZE(this)
  !! Factorize diagonal Block Matrix.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix

  ! local variables
  integer :: ib

  ASSERT(.not. this%factorized)

  if (this%diag) then
    do ib = 1, this%nbrows
      if (.not. allocated(this%b(ib,ib)%p)) call program_error("Can not factorize diagonal block matrix with empty block on main diag!")
      call this%b(ib,ib)%p%factorize()
    end do
  else
    call program_error("Factorization for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if

  this%factorized = .true.
end subroutine

subroutine BLOCK_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve block diagonal system: this * x = rhs. Call this%factorize() first.

  class(BLOCK_TYPE),          intent(in)  :: this
    !! Block matrix
  TT,                         intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                         intent(out) :: x(:)
    !! Store result in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement

  ! local variables
  integer :: ib

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_vec(rhs(this%i0(ib):this%i1(ib)), x(this%j0(ib):this%j1(ib)), trans = trans, iter_refine = iter_refine)
    end do
  else
    call program_error("Solving for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if
end subroutine

subroutine BLOCK_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(BLOCK_TYPE),          intent(in)  :: this
    !! Block matrix
  TT,                         intent(in)  :: rhs(:,:)
    !! right hand side of equation system
  TT,                         intent(out) :: x(:,:)
    !! store result in x
  character(len=1), optional, intent(in)  :: trans
    !! transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! perform iterative refinement

  ! local variables
  integer :: ib

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_mat(rhs(this%i0(ib):this%i1(ib),:), x(this%j0(ib):this%j1(ib),:), trans = trans, iter_refine = iter_refine)
    end do
  else
    call program_error("Solving for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if
end subroutine

subroutine BLOCK_TO_DENSE(this, d, i0, j0)
  !! Insert block matrix into a dense matrix.

  class(BLOCK_TYPE), intent(in)    :: this
    !! Block matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, ib, jb

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements blockwise
  associate (bi0 => this%i0, bj0 => this%j0)
    do ib = 1, this%nbrows
      if (this%diag) then
        if (.not. allocated(this%b(ib,ib)%p)) cycle
        call this%b(ib,ib)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1)
      elseif (this%tridiag) then
        do jb = ib - 1, ib + 1
          if ((jb < 1) .or. (jb > this%nbcols)) cycle
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          call this%b(ib,jb)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1)
        end do
      else
        do jb = 1, this%nbcols
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          call this%b(ib,jb)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1)
        end do
      end if
    end do
  end associate
end subroutine

subroutine BLOCK_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert block matrix into a sparse matrix.

  class(BLOCK_TYPE),  intent(in)    :: this
    !! Block matrix (unchanged)
  type(SPBUILD_TYPE), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional, intent(in)     :: i0
    !! Optional: row offset in sparse matrix
  integer, optional, intent(in)     :: j0
    !! Optional: col offset in sparse matrix
  logical, optional, intent(in)     :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)     :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, ib, jb
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements blockwise
  associate (bi0 => this%i0, bi1 => this%i1, bj0 => this%j0, bj1 => this%j1)
    do ib = 1, this%nbrows
      if (this%diag) then
        if (.not. allocated(this%b(ib,ib)%p)) cycle
        if (present(struct)) then
          call this%b(ib,ib)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1, struct = struct(bi0(ib):bi1(ib), bj0(ib):bj1(ib)))
        else
          call this%b(ib,ib)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1, drop_zeros = drop_zeros)
        end if
      elseif (this%tridiag) then
        do jb = ib - 1, ib + 1
          if ((jb < 1) .or. (jb > this%nbcols)) cycle
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          if (present(struct)) then
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        end do
      else
        do jb = 1, this%nbcols
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          if (present(struct)) then
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        end do
      end if
    end do
  end associate
end subroutine

subroutine BLOCK_GET_GLOBAL_IDX(this, ib, jb, il, jl, ig, jg)
  !! Get global index from block and local index.

  class(BLOCK_TYPE), intent(in)  :: this
    !! Block matrix
  integer,           intent(in)  :: ib
    !! Block row
  integer,           intent(in)  :: jb
    !! Block col
  integer,           intent(in)  :: il
    !! Local row
  integer,           intent(in)  :: jl
    !! Local col
  integer,           intent(out) :: ig
    !! Output global row
  integer,           intent(out) :: jg
    !! Output global col

  ig = this%i0(ib) + il - 1
  jg = this%j0(jb) + jl - 1
end subroutine

subroutine BLOCK_GET_LOCAL_IDX(this, ig, jg, ib, jb, il, jl)
  !! Get block and local index from global index.

  class(BLOCK_TYPE), intent(in)  :: this
    !! Block matrix
  integer,           intent(in)  :: ig
    !! Global row
  integer,           intent(in)  :: jg
    !! Global col
  integer,           intent(out) :: ib
    !! Output block row
  integer,           intent(out) :: jb
    !! Output block col
  integer,           intent(out) :: il
    !! Output local row
  integer,           intent(out) :: jl
    !! Output local col

  ! find block row
  do ib = 1, this%nbrows
    if (this%i1(ib) >= ig) exit
  end do
  il = ig + 1 - this%i0(ib)

  ! find block col
  do jb = 1, this%ncols
    if (this%j1(jb) >= jg) exit
  end do
  jl = jg + 1 - this%j0(jb)
end subroutine

subroutine BLOCK_GET_DENSE(this, ib, jb, p)
  !! Return pointer to dense block.

  class(BLOCK_TYPE), target, intent(in)  :: this
    !! Block matrix
  integer,                   intent(in)  :: ib
    !! Block row
  integer,                   intent(in)  :: jb
    !! Block col
  type(DENSE_TYPE), pointer, intent(out) :: p
    !! Output pointer to dense block or null() if not dense

  nullify(p)
  if (allocated(this%b(ib,jb)%p)) then
    select type (b => this%b(ib,jb)%p)
    class is (DENSE_TYPE)
      p => b
    end select
  end if
end subroutine

subroutine BLOCK_GET_SPARSE(this, ib, jb, p)
  !! Return pointer to sparse block.

  class(BLOCK_TYPE), target,  intent(in)  :: this
    !! Block matrix
  integer,                    intent(in)  :: ib
    !! Block row
  integer,                    intent(in)  :: jb
    !! Block col
  type(SPARSE_TYPE), pointer, intent(out) :: p
    !! Output pointer to sparse block or null() if not dense

  nullify(p)
  if (allocated(this%b(ib,jb)%p)) then
    select type (b => this%b(ib,jb)%p)
    class is (SPARSE_TYPE)
      p => b
    end select
  end if
end subroutine

subroutine BLOCK_GET_BAND(this, ib, jb, p)
  !! Return pointer to band block.

  class(BLOCK_TYPE), target, intent(in)  :: this
    !! Block matrix
  integer,                   intent(in)  :: ib
    !! Block row
  integer,                   intent(in)  :: jb
    !! Block col
  type(BAND_TYPE), pointer,  intent(out) :: p
    !! Output pointer to band block or null() if not dense

  nullify(p)
  if (allocated(this%b(ib,jb)%p)) then
    select type (b => this%b(ib,jb)%p)
    class is (BAND_TYPE)
      p => b
    end select
  end if
end subroutine

subroutine BLOCK_GET_BLOCK(this, ib, jb, p)
  !! Return pointer to block block.

  class(BLOCK_TYPE), target, intent(in)  :: this
    !! Block matrix
  integer,                   intent(in)  :: ib
    !! Block row
  integer,                   intent(in)  :: jb
    !! Block col
  type(BLOCK_TYPE), pointer, intent(out) :: p
    !! Output pointer to block block or null() if not dense

  nullify(p)
  if (allocated(this%b(ib,jb)%p)) then
    select type (b => this%b(ib,jb)%p)
    class is (BLOCK_TYPE)
      p => b
    end select
  end if
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

#undef VECTOR_TYPE

#undef MATRIX_TYPE
#undef MATRIX_ALC_TYPE
#undef MATRIX_INIT
#undef MATRIX_DESTRUCT
#undef MATRIX_RESET
#undef MATRIX_SCALE
#undef MATRIX_MUL_VEC
#undef MATRIX_MUL_MAT
#undef MATRIX_FACTORIZE
#undef MATRIX_SOLVE_VEC
#undef MATRIX_SOLVE_MAT
#undef MATRIX_TO_DENSE
#undef MATRIX_TO_SPARSE

#undef DENSE_TYPE
#undef DENSE_INIT
#undef DENSE_DESTRUCT
#undef DENSE_RESET
#undef DENSE_SCALE
#undef DENSE_MUL_VEC
#undef DENSE_MUL_MAT
#undef DENSE_FACTORIZE
#undef DENSE_SOLVE_VEC
#undef DENSE_SOLVE_MAT
#undef DENSE_TO_DENSE
#undef DENSE_TO_SPARSE
#undef DENSE_TO_REAL
#undef DENSE_TO_CMPLX
#undef DENSE_ADD_DENSE
#undef DENSE_ADD_DENSE3
#undef DENSE_ADD_SPARSE
#undef DENSE_ADD_SPARSE3
#undef DENSE_ADD_BAND
#undef DENSE_ADD_BAND3
#undef DENSE_MUL_DENSE
#undef DENSE_EIG
#undef DENSE_OUTPUT

#undef SPARSE_TYPE
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
#undef SPARSE_ADD_SPARSE
#undef SPARSE_ADD_SPARSE3
#undef SPARSE_ADD_BAND
#undef SPARSE_ADD_BAND3
#undef SPARSE_MUL_SPARSE
#undef SPARSE_GET_CSR
#undef SPARSE_SET_CSR
#undef SPARSE_DESTRUCT_CSR
#undef SPARSE_SET_MV_HINT
#undef SPARSE_SET_MM_HINT
#undef SPARSE_OUTPUT

#undef SPBUILD_TYPE
#undef SPBUILD_INIT
#undef SPBUILD_DESTRUCT
#undef SPBUILD_LOAD
#undef SPBUILD_SET
#undef SPBUILD_SET_ROW
#undef SPBUILD_ADD
#undef SPBUILD_SORT
#undef SPBUILD_FINISH
#undef SPBUILD_RESET

#undef BAND_TYPE
#undef BAND_INIT
#undef BAND_DESTRUCT
#undef BAND_RESET
#undef BAND_SCALE
#undef BAND_MUL_VEC
#undef BAND_MUL_MAT
#undef BAND_FACTORIZE
#undef BAND_SOLVE_VEC
#undef BAND_SOLVE_MAT
#undef BAND_TO_DENSE
#undef BAND_TO_SPARSE
#undef BAND_TO_REAL
#undef BAND_TO_CMPLX
#undef BAND_ADD_BAND
#undef BAND_ADD_BAND3
#undef BAND_OUTPUT

#undef BLOCK_TYPE
#undef BLOCK_INIT
#undef BLOCK_DESTRUCT
#undef BLOCK_RESET
#undef BLOCK_SCALE
#undef BLOCK_MUL_VEC
#undef BLOCK_MUL_MAT
#undef BLOCK_FACTORIZE
#undef BLOCK_SOLVE_VEC
#undef BLOCK_SOLVE_MAT
#undef BLOCK_TO_DENSE
#undef BLOCK_TO_SPARSE
#undef BLOCK_GET_GLOBAL_IDX
#undef BLOCK_GET_LOCAL_IDX
#undef BLOCK_SET_BLOCK_TYPE
#undef BLOCK_GET_DENSE
#undef BLOCK_GET_SPARSE
#undef BLOCK_GET_BAND
#undef BLOCK_GET_BLOCK