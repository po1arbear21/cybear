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
#define BAND_ADD_ELEM        CONCAT3(band,T,add_elem)
#define BAND_GET_ELEM        CONCAT3(band,T,get_elem)
#define BAND_SET_ELEM        CONCAT3(band,T,set_elem)
#define BAND_OUTPUT          CONCAT3(band,T,output)
#define BAND_EYE             CONCAT3(band,eye,T)

subroutine BAND_INIT(this, nrows, nlower, nupper, d0)
  !! Initialize band matrix

  class(BAND_TYPE),  intent(out) :: this
    !! Band matrix
  integer,           intent(in)  :: nrows
    !! Number of columns = Number of rows
  integer,           intent(in)  :: nlower
    !! Number of lower subdiagonals. >= 0
  integer, optional, intent(in)  :: nupper
    !! Number of upper superdiagonals. >= 0
    !! (default: nlower)
  TT, optional,      intent(in)  :: d0(:,:)
    !! Initial data.

  integer :: nupper_

  ! check input data
  ASSERT(nlower >= 0)
  nupper_ = nlower
  if (present(nupper)) nupper_ = nupper
  ASSERT(nupper_ >= 0)
  ASSERT(nlower+nupper_+1 <= nrows)

  ! init base
  call this%init("", nrows, ncols=nrows)

  this%nupper = nupper_
  this%nlower = nlower

  ! allocate data
  if (present(d0)) then
    ASSERT(size(d0, 1) == nlower+nupper_+1)
    ASSERT(size(d0, 2) == nrows)
    allocate(this%d(-nupper_:nlower,nrows), source=d0)
  else
    allocate(this%d(-nupper_:nlower,nrows))
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
  !! Operation: y <- this * x + fact_y * y

  class(BAND_TYPE),    intent(in)    :: this
    !! Band matrix
  TT,                  intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                  intent(inout) :: y(:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  if (this%nupper + this%nlower == 0) then    ! diagonal
    if (present(fact_y)) then
      y = this%d(0,:) * x + fact_y * y
    else
      y = this%d(0,:) * x
    end if
  else
    call gbmv(this%d, x, y, kl=this%nlower, beta=fact_y, trans=trans)
  end if
end subroutine

subroutine BAND_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply band matrix by a raw matrix.
  !! Operation: y <- this * x + fact_y * y

  class(BAND_TYPE),    intent(in)    :: this
    !! Band matrix
  TT,                  intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                  intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  integer :: i

  if (this%nupper + this%nlower == 0) then      ! diagonal
    if (present(fact_y)) then
      do i = 1, size(x, 2)
        y(:,i) = this%d(0,:) * x(:,i) + fact_y * y(:,i)
      end do
    else
      do i = 1, size(x, 2)
        y(:,i) = this%d(0,:) * x(:,i)
      end do
    end if
  else
    do i = 1, size(x, 2)
      call gbmv(this%d, x(:,i), y(:,i), kl=this%nlower, beta=fact_y, trans=trans)
    end do
  end if
end subroutine

subroutine BAND_FACTORIZE(this)
  !! Factorize band matrix. Call this before solving!

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix

  integer :: info

  if      (this%nupper + this%nlower == 0) then                   ! diagonal
    ! do nothing for diagonal matrix. no factorization needed.
  else if (this%nupper == 1 .and. this%nlower == 1) then          ! tridiagonal
    if (.not. allocated(this%f))    allocate(this%f(-2:1,this%ncols))
    if (.not. allocated(this%ipiv)) allocate(this%ipiv(this%ncols))

    this%f(-1:1,:) = this%d

    ! factorization
    associate(n   => this%ncols,                &
              du2 => this%f(-2,3:this%ncols  ), &
              du  => this%f(-1,2:this%ncols  ), &
              d   => this%f( 0,1:this%ncols  ), &
              dl  => this%f( 1,1:this%ncols-1)  )

      call gttrf(dl, d, du, du2, ipiv=this%ipiv, info=info)
      if (info /= 0) call program_error("Tridiagonal Factorization error in LAPACK", code = info)
    end associate
  else                                                            ! general band matrix
    associate(kl => this%nlower, ku => this%nupper)
      if (.not. allocated(this%f))    allocate(this%f(-(kl+ku):kl, this%ncols))
      if (.not. allocated(this%ipiv)) allocate(this%ipiv(this%nrows))

      this%f(-ku:kl,:) = this%d

      ! factorization
      call gbtrf(this%f, kl=kl, ipiv=this%ipiv, info=info)
      if (info /= 0) call program_error("Bandmatrix Factorization error in LAPACK", code = info)
    end associate
  end if
end subroutine

subroutine BAND_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first!

  class(BAND_TYPE),    intent(in)  :: this
    !! Band matrix.
  TT,                  intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                  intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  TT :: rhs_(size(rhs),1), x_(size(x),1)

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

subroutine BAND_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(BAND_TYPE),    intent(in)  :: this
    !! Band matrix
  TT,                  intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                  intent(out) :: x(:,:)
    !! Store results in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  integer :: i, info
  logical :: iter_refine_

  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine

  if      (this%nupper + this%nlower == 0) then                   ! diagonal
    do i = 1, size(x, 2)
      x(:,i) = rhs(:,i) / this%d(0,:)
    end do
  else if (this%nupper == 1 .and. this%nlower == 1) then          ! tridiagonal
    ! called factorize before?
    ASSERT(allocated(this%f) .and. allocated(this%ipiv))

    associate(n   => this%ncols,                &
              du  => this%d(-1,2:this%ncols  ), &
              d   => this%d( 0,1:this%ncols  ), &
              dl  => this%d( 1,1:this%ncols-1), &
              du2 => this%f(-2,3:this%ncols  ), &
              duf => this%f(-1,2:this%ncols  ), &
              df  => this%f( 0,1:this%ncols  ), &
              dlf => this%f( 1,1:this%ncols-1)  )

      x = rhs
      call gttrs(dlf, df, duf, du2, x, this%ipiv, trans=trans, info=info)
      if (info /= 0) call program_error("Tridiagonal Solution error in LAPACK", code = info)

      if (iter_refine_) then
        call gtrfs(dl, d, du, dlf, df, duf, du2, this%ipiv, rhs, x, trans=trans, info=info)
        if (info /= 0) call program_error("Tridiagonal Iterative Refinement error in LAPACK", code = info)
      end if
    end associate
  else                                                            ! general band matrix
    ! called factorize before?
    ASSERT(allocated(this%f) .and. allocated(this%ipiv))

    x = rhs
    call gbtrs(this%f, x, this%ipiv, kl=this%nlower, trans=trans, info=info)
    if (info /= 0) call program_error("Band Matrix Solution error in LAPACK", code = info)

    if (iter_refine_) then
      call gbrfs(this%d, this%f, this%ipiv, rhs, x, kl=this%nlower, trans=trans, info=info)
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

  integer :: i0_, j0_, i, j, k

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements into dense matrix
  do j = 1, this%ncols
    do i = -this%nupper, this%nlower
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
  integer :: i0_, j0_, i, j, k
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
  do j = 1, this%ncols
    do i = -this%nupper, this%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > this%nrows)) cycle

      ! check structure
      if      (present(struct)) then
        if (.not. struct(k,j)) cycle
      else if (drop_zeros_) then
        if (this%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+k-1, j0_+j-1, this%d(i,j), search=.false.)
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

  call b%init(this%nrows, this%nlower, nupper=this%nupper, d0=real(this%d))
end subroutine
#else
subroutine BAND_TO_CMPLX(this, b)
  !! Convert band real matrix to complex.

  class(band_real),  intent(in)  :: this
    !! Band real matrix
  class(band_cmplx), intent(out) :: b
    !! Output cmplx(this)

  call b%init(this%nrows, this%nlower, nupper=this%nupper, d0=cmplx(this%d))
end subroutine
#endif

subroutine BAND_ADD_BAND(this, b, fact)
  !! Add two band matrices.
  !! Operation: this <- this + fact * b

  class(BAND_TYPE), intent(inout) :: this
    !! First band matrix (updated)
  type(BAND_TYPE),  intent(in)    :: b
    !! Second band matrix (unchanged)
  TT, optional,     intent(in)    :: fact
    !! Second band matrix scaling factor (default: 1)

  TT              :: fact_
  TT, allocatable :: tmp(:,:)

  ASSERT(.not. allocated(this%f))
  ASSERT(.not. allocated(this%ipiv))
  ASSERT(this%nrows == b%nrows .and. this%ncols == b%ncols)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  if (this%nupper >= b%nupper .and. this%nlower >= b%nlower) then
    this%d(-b%nupper:b%nlower,1:this%ncols) = this%d(-b%nupper:b%nlower,1:this%ncols) + fact_ * b%d
  else
    allocate(tmp(-max(this%nupper,b%nupper):max(this%nlower,b%nlower),1:this%ncols))
    tmp = 0
    tmp(-this%nupper:this%nlower,1:this%ncols) = this%d
    tmp(-b%nupper:b%nlower,1:this%ncols) = this%d(-b%nupper:b%nlower,1:this%ncols) + fact_*b%d
    call move_alloc(tmp, this%d)
    this%nlower = max(this%nlower,b%nlower)
    this%nupper = max(this%nupper,b%nupper)
  end if
end subroutine

subroutine BAND_ADD_BAND3(this, b2, b3, fact1, fact2)
  !! Add two band matrices, store result in third band matrix.
  !! Operation: b3 <- fact1*this + fact2*b2

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

  ! case selection: hopefully no tmp array needs to be allocated
  if (this%nupper >= b2%nupper .and. this%nlower >= b2%nlower) then
    call b3%init(this%nrows, this%nlower, nupper=this%nupper, d0=fact1_*this%d)     ! b3 <- fact1*this
    call b3%add_band(b2, fact=fact2_)                                               ! b3 += fact2*b2
  else
    call b3%init(b2%nrows, b2%nlower, nupper=b2%nupper, d0=fact2_*b2%d)             ! b3 <- fact2*b2
    call b3%add_band(this, fact=fact1_)                                             ! b3 += fact1*this
  end if
end subroutine

subroutine BAND_OUTPUT(this, file, funit, fmt)
  !! Write band matrix to file or stream.

  class(BAND_TYPE),       intent(in) :: this
    !! Band matrix
  character(*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,      optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(*), optional, intent(in) :: fmt
    !! Value format string

  ! local variables
  integer                   :: i, j, k, status, funit_
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

  ! output values
  do j = 1, this%ncols
    do i = -this%nupper, this%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > this%nrows)) cycle

      ! write row and column
      write(funit_, "(I0,A,I0,A)", advance="no") k, " ", j, " "

      ! write value
      write(funit_, fmt_) this%d(i,j)
    end do
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

subroutine BAND_ADD_ELEM(this, i, j, x)
  class(BAND_TYPE), intent(inout) :: this
  integer,          intent(in)    :: i
    !! row index of full matrix to represent
  integer,          intent(in)    :: j
    !! col index of full matrix to represent
  TT,               intent(in)    :: x
    !! matrix element to add

  this%d(i-j,j) = this%get_elem(i,j) + x
end subroutine

function BAND_GET_ELEM(this, i, j) result(x)
  class(BAND_TYPE), intent(in) :: this
  integer,          intent(in) :: i
    !! row index of full matrix to represent
  integer,          intent(in) :: j
    !! col index of full matrix to represent
  TT                           :: x
    !! matrix element

  ASSERT(-this%nupper <= i-j .and. i-j <= this%nlower)
  ASSERT(           1 <= j   .and.   j <= this%ncols )

  x = this%d(i-j,j)
end function

subroutine BAND_SET_ELEM(this, i, j, x)
  class(BAND_TYPE), intent(inout) :: this
  integer,          intent(in)    :: i
    !! row index of full matrix to represent
  integer,          intent(in)    :: j
    !! col index of full matrix to represent
  TT,               intent(in)    :: x
    !! matrix element to set

  ASSERT(-this%nupper <= i-j .and. i-j <= this%nlower)
  ASSERT(           1 <= j   .and.   j <= this%ncols )

  this%d(i-j,j) = x
end subroutine

function BAND_EYE(nrows) result(eye)
  !! band identity matrix

  integer, intent(in) :: nrows
    !! Number of rows
  type(BAND_TYPE)     :: eye

  call eye%init(nrows, 0)
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
#undef SPARSE_TYPE
#undef SPBUILD_TYPE
#undef BAND_TYPE
#undef HESSENBERG_TYPE
#undef TRIANG_TYPE
#undef BLOCK_TYPE

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
#undef BAND_ADD_ELEM
#undef BAND_GET_ELEM
#undef BAND_SET_ELEM
#undef BAND_OUTPUT
#undef BAND_EYE
