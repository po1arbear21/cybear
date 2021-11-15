#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_PTR_TYPE        CONCAT(band_ptr,T)
#define MATRIX_TYPE          CONCAT(matrix,T)

#define BAND_INIT            CONCAT3(band,T,init)
#define BAND_DESTRUCT        CONCAT3(band,T,destruct)
#define BAND_RESET           CONCAT3(band,T,reset)
#define BAND_RESET_ROW       CONCAT3(band,T,reset_row)
#define BAND_SCALE           CONCAT3(band,T,scale)
#define BAND_MUL_VEC         CONCAT3(band,T,mul_vec)
#define BAND_MUL_MAT         CONCAT3(band,T,mul_mat)
#define BAND_FACTORIZE       CONCAT3(band,T,factorize)
#define BAND_SOLVE_VEC       CONCAT3(band,T,solve_vec)
#define BAND_SOLVE_MAT       CONCAT3(band,T,solve_mat)
#define BAND_ADD_ELEM        CONCAT3(band,T,add_elem)
#define BAND_GET_ELEM        CONCAT3(band,T,get_elem)
#define BAND_SET_ELEM        CONCAT3(band,T,set_elem)
#define BAND_SET_DIAG_ARR    CONCAT3(band,T,set_diag_arr)
#define BAND_SET_DIAG_VAL    CONCAT3(band,T,set_diag_val)
#define BAND_OUTPUT          CONCAT3(band,T,output)
#define BAND_EYE             CONCAT3(band,eye,T)

module subroutine BAND_INIT(this, nrows, nlower, nupper, d0)
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
  TT,      optional, intent(in)  :: d0(:,:)
    !! Initial data.

  integer :: nupper_

  ! check input data
  ASSERT(nlower >= 0)
  nupper_ = nlower
  if (present(nupper)) nupper_ = nupper
  ASSERT(nupper_ >= 0)
  ASSERT(nlower+nupper_+1 <= 2*nrows-1)

  ! init base
  call this%init("", nrows, ncols=nrows)

  this%nupper = nupper_
  this%nlower = nlower

  ! allocate data
  if (present(d0)) then
    ASSERT(size(d0, 1) == nlower+nupper_+1)
    ASSERT(size(d0, 2) == nrows)
    allocate (this%d(-nupper_:nlower,nrows), source=d0)
  else
    allocate (this%d(-nupper_:nlower,nrows))
    this%d = 0
  end if
end subroutine

module subroutine BAND_DESTRUCT(this)
  !! Destruct band matrix by deallocating its memory.

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix

  if (allocated(this%d   )) deallocate (this%d   )
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)

  this%factorized = .false.
end subroutine

module subroutine BAND_RESET(this, only_factorization)
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
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)

  this%factorized = .false.
end subroutine

module subroutine BAND_RESET_ROW(this, i)
  !! Reset i-th row to zero and its factorization.

  class(BAND_TYPE),  intent(inout) :: this
    !! Band matrix
  integer,           intent(in)    :: i
    !! i-th row

  integer :: j, j0, j1
#ifdef TCMPLX
  complex :: zero
#else
  real    :: zero
#endif

  zero = 0

  j0 = max(1,          i-this%nlower)
  j1 = min(this%ncols, i+this%nupper)

  do j = j0, j1
    call this%set_elem(i, j, zero)
  end do

  ! reset factorization
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)

  this%factorized = .false.
end subroutine

module subroutine BAND_SCALE(this, fact)
  !! Scale band matrix by a factor.

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix
  TT,               intent(in)    :: fact
    !! Scaling factor

  ASSERT(.not. allocated(this%f   ))
  ASSERT(.not. allocated(this%ipiv))

  this%d = this%d * fact
end subroutine

module subroutine BAND_MUL_VEC(this, x, y, fact_y, trans)
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

module subroutine BAND_MUL_MAT(this, x, y, fact_y, trans)
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

module subroutine BAND_FACTORIZE(this)
  !! Factorize band matrix. Call this before solving!

  class(BAND_TYPE), intent(inout) :: this
    !! Band matrix

  integer :: info

  if      (this%nupper + this%nlower == 0) then                   ! diagonal
    ! do nothing for diagonal matrix. no factorization needed.
  else if (this%nupper == 1 .and. this%nlower == 1) then          ! tridiagonal
    if (.not. allocated(this%f))    allocate (this%f(-2:1,this%ncols))
    if (.not. allocated(this%ipiv)) allocate (this%ipiv(this%ncols))

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
      if (.not. allocated(this%f))    allocate (this%f(-(kl+ku):kl, this%ncols))
      if (.not. allocated(this%ipiv)) allocate (this%ipiv(this%nrows))

      this%f(-ku:kl,:) = this%d

      ! factorization
      call gbtrf(this%f, kl=kl, ipiv=this%ipiv, info=info)
      if (info /= 0) call program_error("Bandmatrix Factorization error in LAPACK", code = info)
    end associate
  end if

  this%factorized = .true.
end subroutine

module subroutine BAND_SOLVE_VEC(this, rhs, x, trans, iter_refine)
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

module subroutine BAND_SOLVE_MAT(this, rhs, x, trans, iter_refine)
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

  else if ((this%nupper == 1) .and. (this%nlower == 1)) then      ! tridiagonal
    ASSERT(this%factorized)

    associate(n   => this%ncols,                &
      &       du  => this%d(-1,2:this%ncols  ), &
      &       d   => this%d( 0,1:this%ncols  ), &
      &       dl  => this%d( 1,1:this%ncols-1), &
      &       du2 => this%f(-2,3:this%ncols  ), &
      &       duf => this%f(-1,2:this%ncols  ), &
      &       df  => this%f( 0,1:this%ncols  ), &
      &       dlf => this%f( 1,1:this%ncols-1)  )

      x = rhs
      call gttrs(dlf, df, duf, du2, x, this%ipiv, trans=trans, info=info)
      if (info /= 0) call program_error("Tridiagonal Solution error in LAPACK", code = info)

      if (iter_refine_) then
        call gtrfs(dl, d, du, dlf, df, duf, du2, this%ipiv, rhs, x, trans=trans, info=info)
        if (info /= 0) call program_error("Tridiagonal Iterative Refinement error in LAPACK", code = info)
      end if
    end associate

  else                                                            ! general band matrix
    ASSERT(this%factorized)

    x = rhs
    call gbtrs(this%f, x, this%ipiv, kl=this%nlower, trans=trans, info=info)
    if (info /= 0) call program_error("Band Matrix Solution error in LAPACK", code = info)

    if (iter_refine_) then
      call gbrfs(this%d, this%f, this%ipiv, rhs, x, kl=this%nlower, trans=trans, info=info)
      if (info /= 0) call program_error("Band Matrix Iterative Refinement error in LAPACK", code = info)
    end if
  end if
end subroutine

module subroutine BAND_OUTPUT(this, file, funit, fmt)
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

module subroutine BAND_ADD_ELEM(this, i, j, x)
  class(BAND_TYPE), intent(inout) :: this
  integer,          intent(in)    :: i
    !! row index of full matrix to represent
  integer,          intent(in)    :: j
    !! col index of full matrix to represent
  TT,               intent(in)    :: x
    !! matrix element to add

  ASSERT((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  ASSERT((           1 <= j  ) .and. (  j <= this%ncols ))

  this%d(i-j,j) = this%get_elem(i,j) + x
end subroutine

module function BAND_GET_ELEM(this, i, j) result(x)
  class(BAND_TYPE), intent(in) :: this
  integer,          intent(in) :: i
    !! row index of full matrix to represent
  integer,          intent(in) :: j
    !! col index of full matrix to represent
  TT                           :: x
    !! matrix element

  ASSERT((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  ASSERT((           1 <= j  ) .and. (  j <= this%ncols ))

  x = this%d(i-j,j)
end function

module subroutine BAND_SET_ELEM(this, i, j, x)
  class(BAND_TYPE), intent(inout) :: this
  integer,          intent(in)    :: i
    !! row index of full matrix to represent
  integer,          intent(in)    :: j
    !! col index of full matrix to represent
  TT,               intent(in)    :: x
    !! matrix element to set

  ASSERT((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  ASSERT((           1 <= j  ) .and. (  j <= this%ncols ))

  this%d(i-j,j) = x
end subroutine

module subroutine BAND_SET_DIAG_VAL(this, val, k)
  !! inserts same value on the k-th diagonal.

  class(BAND_TYPE), intent(inout)        :: this
  TT,               intent(in)           :: val
    !! value
  integer,          intent(in), optional :: k
    !! k-th diagonal.
    !! k>0 lower diagonal, k<0 upper diagonal.
    !! ATTENTION: opposite to numbering in matlab!
    !! default: 0 == main diagonal

  integer :: k_, i

  k_ = 0
  if (present(k)) k_ = k

  ! create array of length: n-|k|
  call this%set_diag([(val, i=1,this%nrows-abs(k_))], k=k_)
end subroutine

module subroutine BAND_SET_DIAG_ARR(this, vals, k)
  !! inserts values on the k-th diagonal.

  class(BAND_TYPE), intent(inout)        :: this
  TT,               intent(in)           :: vals(:)
    !! values
  integer,          intent(in), optional :: k
    !! k-th diagonal.
    !! k>0 lower diagonal, k<0 upper diagonal.
    !! ATTENTION: opposite to numbering in matlab!
    !! default: 0 == main diagonal

  integer :: k_, i1, i2, n

  k_ = 0
  if (present(k)) k_ = k

  ASSERT((-k_ <= this%nupper) .and. (k_ <= this%nlower))

  n  = this%nrows
  i1 = merge(1, 1-k_, k_>=0)
  i2 = merge(n, n-k_, k_<=0)

  this%d(k_,i1:i2) = vals
end subroutine

module function BAND_EYE(nrows) result(eye)
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

#undef BAND_TYPE
#undef BAND_PTR_TYPE
#undef MATRIX_TYPE

#undef BAND_INIT
#undef BAND_DESTRUCT
#undef BAND_RESET
#undef BAND_RESET_ROW
#undef BAND_SCALE
#undef BAND_MUL_VEC
#undef BAND_MUL_MAT
#undef BAND_FACTORIZE
#undef BAND_SOLVE_VEC
#undef BAND_SOLVE_MAT
#undef BAND_ADD_ELEM
#undef BAND_GET_ELEM
#undef BAND_SET_ELEM
#undef BAND_SET_DIAG_ARR
#undef BAND_SET_DIAG_VAL
#undef BAND_OUTPUT
#undef BAND_EYE
