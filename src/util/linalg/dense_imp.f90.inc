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
#define BLOCK_TYPE           CONCAT(block,T)

#define DENSE_INIT_BY_DATA   CONCAT3(dense,T,init_by_data)
#define DENSE_INIT_BY_SIZE   CONCAT3(dense,T,init_by_size)
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
#define DENSE_TRANSPOSE      CONCAT3(dense,T,transpose)
#define DENSE_TRANSPOSE2     CONCAT3(dense,T,transpose2)
#define DENSE_ADD_DENSE      CONCAT3(dense,T,add_dense)
#define DENSE_ADD_DENSE3     CONCAT3(dense,T,add_dense3)
#define DENSE_ADD_SPARSE     CONCAT3(dense,T,add_sparse)
#define DENSE_ADD_SPARSE3    CONCAT3(dense,T,add_sparse3)
#define DENSE_ADD_BAND       CONCAT3(dense,T,add_band)
#define DENSE_ADD_BAND3      CONCAT3(dense,T,add_band3)
#define DENSE_MUL_DENSE      CONCAT3(dense,T,mul_dense)
#define DENSE_EIG            CONCAT3(dense,T,eig)
#define DENSE_OUTPUT         CONCAT3(dense,T,output)
#define DENSE_EYE            CONCAT(eye,T)

subroutine DENSE_INIT_BY_DATA(this, d)
  !! Initialize dense matrix.

  class(DENSE_TYPE), intent(out) :: this
    !! Dense matrix
  TT,                intent(in)  :: d(:,:)
    !! Provide initial data.

  ! init base
  call this%init("", size(d, 1), ncols=size(d, 2))

  this%d = d
end subroutine

subroutine DENSE_INIT_BY_SIZE(this, nrows, ncols)
  !! Initialize dense matrix.

  class(DENSE_TYPE), intent(out) :: this
    !! Dense matrix
  integer,           intent(in)  :: nrows
    !! Number of rows (Default: d0rows)
  integer, optional, intent(in)  :: ncols
    !! Optional: Number of columns. (Default: ncols = nrows)

  integer :: ncols_

  ncols_ = nrows
  if (present(ncols)) ncols_ = ncols

  ! init base
  call this%init("", nrows, ncols=ncols_)

  ! allocate data
  allocate(this%d(nrows,ncols_))
  this%d = 0
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
    !! If present, scale existing data in y and add result of multiplication
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
  if (present(struct)) then
    ASSERT(size(this%d, 1) == size(struct, 1) .and. size(this%d, 2) == size(struct, 2))
  end if

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

  call d%init(real(this%d))
end subroutine
#else
subroutine DENSE_TO_CMPLX(this, d)
  !! Convert dense real matrix to complex with zero imaginary part.

  class(dense_real),  intent(in)  :: this
    !! Dense real matrix
  class(dense_cmplx), intent(out) :: d
    !! Output cmplx(this)

  call d%init(cmplx(this%d))
end subroutine
#endif

subroutine DENSE_TRANSPOSE(this)
  !! Transpose dense matrix

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix

  ! local variables
  TT, allocatable :: d(:,:)

  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  ! transpose matrix
  allocate (d(this%ncols,this%nrows), source = transpose(this%d))
  call move_alloc(d, this%d)
  this%nrows = size(this%d,1)
  this%ncols = size(this%d,2)
end subroutine

subroutine DENSE_TRANSPOSE2(this, d2)
  !! Transpose matrix and store result in second matrix

  class(DENSE_TYPE), intent(in) :: this
    !! First dense matrix (unchanged)
  type(DENSE_TYPE), intent(out) :: d2
    !! Output transposed matrix

  call d2%init(transpose(this%d))
end subroutine

subroutine DENSE_ADD_DENSE(this, d, fact)
  !! Add two dense matrices, store the result in first matrix.
  !!
  !! operation:  M1 <- M1 + fact*M2

  class(DENSE_TYPE), intent(inout) :: this
    !! First dense matrix M1 (updated)
  type(DENSE_TYPE),  intent(in)    :: d
    !! Second dense matrix M2 (unchanged)
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
  !!
  !! operation:  d3 <- fact1 d1 + fact2 d2

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix d1 (unchanged)
  type(DENSE_TYPE),  intent(in)  :: d2
    !! Dense second matrix d2 (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix d3
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
  call d3%init(fact1_ * this%d + fact2_ * d2%d)
end subroutine

subroutine DENSE_ADD_SPARSE(this, s, fact)
  !! Add dense and sparse matrix, store result in dense matrix.
  !!
  !! operation:  D <- D + fact S

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix D (updated)
  type(SPARSE_TYPE), intent(in)    :: s
    !! Sparse matrix S (unchanged)
  TT, optional,      intent(in)    :: fact
    !! Sparse matrix scaling factor (default: 1)

  ! local variables
  integer :: i, j
  TT      :: fact_

  ASSERT(this%nrows == s%nrows)
  ASSERT(this%ncols == s%ncols)
  ASSERT(.not. allocated(this%lu))
  ASSERT(.not. allocated(this%ipiv))

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! add contents of sparse matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      this%d(i,s%ja(j)) = this%d(i,s%ja(j)) + fact_ * s%a(j)
    end do
  end do
end subroutine

subroutine DENSE_ADD_SPARSE3(this, s2, d3, fact1, fact2)
  !! Add dense and sparse matrix, store result in another dense matrix.
  !!
  !! operation:  d3 <- fact1 d1 + fact2 s2

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix d1 (unchanged)
  type(SPARSE_TYPE), intent(in)  :: s2
    !! Sparse second matrix s2 (unchanged)
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
  call d3%init(fact1_ * this%d)

  ! add sparse matrix to d3
  call d3%add_sparse(s2, fact = fact2)
end subroutine

subroutine DENSE_ADD_BAND(this, b, fact)
  !! Add dense and band matrix, store result in dense matrix.
  !!
  !! operation:  d <- d + fact b

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix d (updated)
  type(BAND_TYPE),   intent(in)    :: b
    !! Band matrix b (unchanged)
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
  !!
  !! operation:  d3 <- fact1 d1 + fact2 b2

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix d1 (unchanged)
  type(BAND_TYPE),   intent(in)  :: b2
    !! Band second matrix b2 (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix d3
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
  call d3%init(fact1_ * this%d)

  ! add band matrix to d3
  call d3%add_band(b2, fact = fact2)
end subroutine

subroutine DENSE_MUL_DENSE(this, d2, d3)
  !! Multiply two dense matrices, store result in a third dense matrix.

  class(DENSE_TYPE), intent(in)  :: this
    !! Dense first matrix d1 (unchanged)
  type(DENSE_TYPE),  intent(in)  :: d2
    !! Dense second matrix d2 (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix d3

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

function DENSE_EYE(nrows, ncols) result(eye)
  !! dense identity matrix

  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(DENSE_TYPE)                :: eye

  ! local variables
  integer         :: i

  call eye%init(nrows, ncols)
  do i = 1, min(eye%nrows, eye%ncols)
    eye%d(i,i) = 1
  end do
end function

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
#undef BLOCK_TYPE

#undef DENSE_INIT_BY_DATA
#undef DENSE_INIT_BY_SIZE
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
#undef DENSE_TRANSPOSE
#undef DENSE_TRANSPOSE2
#undef DENSE_ADD_DENSE
#undef DENSE_ADD_DENSE3
#undef DENSE_ADD_SPARSE
#undef DENSE_ADD_SPARSE3
#undef DENSE_ADD_BAND
#undef DENSE_ADD_BAND3
#undef DENSE_MUL_DENSE
#undef DENSE_EIG
#undef DENSE_OUTPUT
#undef DENSE_EYE