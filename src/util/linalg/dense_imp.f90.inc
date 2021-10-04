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
#define DENSE_TRANSPOSE      CONCAT3(dense,T,transpose)
#define DENSE_TRANSPOSE2     CONCAT3(dense,T,transpose2)
#define DENSE_MUL_DENSE      CONCAT3(dense,T,mul_dense)
#define DENSE_EIG            CONCAT3(dense,T,eig)
#define DENSE_OUTPUT         CONCAT3(dense,T,output)
#define DENSE_EYE            CONCAT3(dense,eye,T)

module subroutine DENSE_INIT_BY_DATA(this, d)
  !! Initialize dense matrix.

  class(DENSE_TYPE), intent(out) :: this
    !! Dense matrix
  TT,                intent(in)  :: d(:,:)
    !! Provide initial data.

  ! init base
  call this%init("", size(d, 1), ncols=size(d, 2))

  this%d = d
end subroutine

module subroutine DENSE_INIT_BY_SIZE(this, nrows, ncols)
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
  allocate (this%d(nrows,ncols_))
  this%d = 0
end subroutine

module subroutine DENSE_DESTRUCT(this)
  !! Destruct dense matrix; i.e. deallocate its members.

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix

  if (allocated(this%d   )) deallocate (this%d   )
  if (allocated(this%lu  )) deallocate (this%lu  )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine DENSE_RESET(this, only_factorization)
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
  if (allocated(this%lu  )) deallocate (this%lu  )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine DENSE_SCALE(this, fact)
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

module subroutine DENSE_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(DENSE_TYPE),   intent(in)    :: this
    !! Dense matrix
  TT,                  intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                  intent(inout) :: y(:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  call gemv(this%d, x, y, beta = fact_y, trans = trans)
end subroutine

module subroutine DENSE_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(DENSE_TYPE),   intent(in)    :: this
    !! Dense matrix
  TT,                  intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                  intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! Scale existing data in y and add result to it (default: 0)
  character, optional, intent(in)    :: trans
    !! Transposed flag for this matrix, must be 'N', 'C' or 'T' (default: 'N')

  call gemm(this%d, x, y, transa = trans, beta = fact_y)
end subroutine

module subroutine DENSE_FACTORIZE(this)
  !! Calculate LU decomposition.
  !! Call this subroutine before solving!

  class(DENSE_TYPE), intent(inout) :: this
    !! Dense matrix

  ! local variables
  integer :: info

  ! allocate memory for factorization
  if (.not. allocated(this%lu  )) allocate (this%lu(  this%nrows,this%ncols))
  if (.not. allocated(this%ipiv)) allocate (this%ipiv(this%nrows))

  ! factorization with lapack
  this%lu = this%d
  call getrf(this%lu, ipiv = this%ipiv, info = info)
  if (info /= 0) call program_error("Dense Factorization error in LAPACK", code = info)

  this%factorized = .true.
end subroutine

module subroutine DENSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for a single right-hand side.
  !! Solves: this * x = rhs, where x and rhs are raw vectors. Call this%factorize() first!

  class(DENSE_TYPE),   intent(in)  :: this
    !! Dense matrix
  TT,                  intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                  intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  TT :: rhs_(size(rhs),1), x_(size(x),1)

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

module subroutine DENSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for multiple right-hand sides.
  !! Solves: this * x = rhs, where x and rhs are raw matrices. Call this%factorize() first!

  class(DENSE_TYPE),   intent(in)  :: this
    !! Dense matrix
  TT,                  intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                  intent(out) :: x(:,:)
    !! Store results in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  integer :: info
  logical :: iter_refine_

  ASSERT(this%factorized)

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

module subroutine DENSE_TRANSPOSE(this)
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

module subroutine DENSE_TRANSPOSE2(this, d2)
  !! Transpose matrix and store result in second matrix

  class(DENSE_TYPE), intent(in) :: this
    !! First dense matrix (unchanged)
  type(DENSE_TYPE), intent(out) :: d2
    !! Output transposed matrix

  call d2%init(transpose(this%d))
end subroutine

module subroutine DENSE_MUL_DENSE(this, d2, d3)
  !! Multiply two dense matrices, store result in a third dense matrix.
  !!  Operation: d3 <- this * d2

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

module subroutine DENSE_EIG(this, e, R, L, sort)
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
  if (sort_) then
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

module subroutine DENSE_OUTPUT(this, file, funit, fmt)
  !! Write dense matrix to file or stream.

  class(DENSE_TYPE),      intent(in) :: this
    !! Dense matrix
  character(*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,      optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(*), optional, intent(in) :: fmt
    !! Value format string

  ! local variables
  integer                   :: i, j, funit_, status
  character(:), allocatable :: fmt_

  allocate (character(0) :: fmt_) ! remove gfortran warning
  if (present(fmt)) then
    fmt_ = fmt
  else
#ifdef TCMPLX
    fmt_ = '(A,ES24.16,SP,ES24.16,"i")'
#else
    fmt_ = '(A,ES24.16)'
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

module function DENSE_EYE(nrows, ncols) result(eye)
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
#undef DENSE_TRANSPOSE
#undef DENSE_TRANSPOSE2
#undef DENSE_MUL_DENSE
#undef DENSE_EIG
#undef DENSE_OUTPUT
#undef DENSE_EYE
