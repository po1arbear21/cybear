#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define BAND_TYPE            CONCAT(band,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define BLOCK_TYPE           CONCAT(block,T)

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
#define BAND_EYE             CONCAT3(band,T,eye)

subroutine BAND_INIT(this, ncols, ndiags, d0)
  !! Initialize band matrix

  class(BAND_TYPE),  intent(out) :: this
    !! Band matrix
  integer,           intent(in)  :: ncols
    !! Number of columns = Number of rows
  integer,           intent(in)  :: ndiags
    !! Total number of diagonals. Must be positive uneven integer.
  TT, optional,      intent(in)  :: d0(-(ndiags-1)/2:,:)
    !! Initial data.

  ! init base
  call this%init("", ncols, ncols = ncols)

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
end subroutine

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

function BAND_EYE(nrows) result(eye)
  !! band identity matrix

  integer, intent(in)  :: nrows
    !! Number of rows
  type(BAND_TYPE)      :: eye

  call eye%init(nrows, 1)
  eye%d = 1
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
#undef BAND_TYPE
#undef SPARSE_TYPE
#undef SPBUILD_TYPE
#undef BLOCK_TYPE

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
#undef BAND_EYE