module subroutine band_{}T{}_init(this, nrows, nlower, nupper, d0)
  !! Initialize band matrix

  class(band_{}T),      intent(out) :: this
    !! Band matrix
  integer,              intent(in)  :: nrows
    !! Number of columns = Number of rows
  integer,              intent(in)  :: nlower
    !! Number of lower subdiagonals. >= 0
  integer,    optional, intent(in)  :: nupper
    !! Number of upper superdiagonals. >= 0
    !! (default: nlower)
  m4_type(T), optional, intent(in)  :: d0(:,:)
    !! Initial data.

  integer :: nupper_

  ! check input data
  m4_assert(nlower >= 0)
  nupper_ = nlower
  if (present(nupper)) nupper_ = nupper
  m4_assert(nupper_ >= 0)
  m4_assert(nlower+nupper_+1 <= 2*nrows-1)

  ! init base
  call this%init("", nrows, ncols=nrows)

  this%nupper = nupper_
  this%nlower = nlower

  ! allocate data
  if (present(d0)) then
    m4_assert(size(d0, 1) == nlower+nupper_+1)
    m4_assert(size(d0, 2) == nrows)
    allocate (this%d(-nupper_:nlower,nrows), source=d0)
  else
    allocate (this%d(-nupper_:nlower,nrows))
    this%d = 0
  end if

  this%solver = default_bsolver
end subroutine

module subroutine band_{}T{}_destruct(this)
  !! Destruct band matrix by deallocating its memory.

  class(band_{}T), intent(inout) :: this
    !! Band matrix

  if (allocated(this%d   )) deallocate (this%d   )
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)

  this%factorized = .false.
end subroutine

module subroutine band_{}T{}_reset(this, only_factorization)
  !! Reset band matrix to zero (optionally only the factorization).

  class(band_{}T),  intent(inout) :: this
    !! Band matrix
  logical, optional, intent(in)   :: only_factorization
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

module subroutine band_{}T{}_reset_row(this, i)
  !! Reset i-th row to zero and its factorization.

  class(band_{}T),  intent(inout) :: this
    !! Band matrix
  integer,          intent(in)    :: i
    !! i-th row

  integer :: j, j0, j1
  m4_ifelse(T,cmplx,{complex},{real}) :: zero

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

module subroutine band_{}T{}_scale(this, fact)
  !! Scale band matrix by a factor.

  class(band_{}T), intent(inout) :: this
    !! Band matrix
  m4_type(T),      intent(in)    :: fact
    !! Scaling factor

  m4_assert(.not. allocated(this%f   ))
  m4_assert(.not. allocated(this%ipiv))

  this%d = this%d * fact
end subroutine

module subroutine band_{}T{}_mul_vec(this, x, y, fact_y, trans)
  !! Multiply band matrix by a raw vector.
  !! Operation: y <- this * x + fact_y * y

  class(band_{}T),      intent(in)    :: this
    !! Band matrix
  m4_type(T),           intent(in)    :: x(:)
    !! Raw vector to multiply by
  m4_type(T),           intent(inout) :: y(:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
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

module subroutine band_{}T{}_mul_mat(this, x, y, fact_y, trans)
  !! Multiply band matrix by a raw matrix.
  !! Operation: y <- this * x + fact_y * y

  class(band_{}T),      intent(in)    :: this
    !! Band matrix
  m4_type(T),           intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  m4_type(T),           intent(inout) :: y(:,:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
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

module subroutine band_{}T{}_factorize(this)
  !! Factorize band matrix. Call this before solving!

  class(band_{}T), intent(inout) :: this
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
      select case (this%solver)
      case (BSOLVER_LAPACK)
        call gbtrf(this%f, kl=kl, ipiv=this%ipiv, info=info)
        if (info /= 0) call program_error("Bandmatrix Factorization error in LAPACK", code = info)
      m4_divert(m4_ifdef({m4_spike},0,-1))
      case (BSOLVER_SPIKE)
        call this%sp%init(m4_ifelse(T,cmplx,.true.,.false.), this%nrows, max(this%nlower, this%nupper))
        call spike_factorize(this%sp, this%nrows, this%nlower, this%nupper, this%f, this%ipiv, info = info)
        if (info /= 0) call program_error("Bandmatrix Factorization error in SPIKE", code = info)
      m4_divert(0)
      end select
    end associate
  end if

  this%factorized = .true.
end subroutine

module subroutine band_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first!

  class(band_{}T),     intent(in)  :: this
    !! Band matrix.
  m4_type(T),          intent(in)  :: rhs(:)
    !! Right hand side of equation system
  m4_type(T),          intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  m4_type(T) :: rhs_(size(rhs),1), x_(size(x),1)

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

module subroutine band_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(band_{}T),     intent(in)  :: this
    !! Band matrix
  m4_type(T),          intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  m4_type(T),          intent(out) :: x(:,:)
    !! Store results in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true, disabled for SPIKE)

  integer :: i, info
  logical :: iter_refine_

  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine

  if (this%nupper + this%nlower == 0) then ! diagonal
    do i = 1, size(x, 2)
      x(:,i) = rhs(:,i) / this%d(0,:)
    end do
  else if ((this%nupper == 1) .and. (this%nlower == 1)) then ! tridiagonal
    m4_assert(this%factorized)

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
  else ! general band matrix
    m4_assert(this%factorized)

    x = rhs
    select case (this%solver)
    case (BSOLVER_LAPACK)
      call gbtrs(this%f, x, this%ipiv, kl=this%nlower, trans=trans, info=info)
      if (info /= 0) call program_error("Band Matrix Solution error in LAPACK", code = info)
      if (iter_refine_) then
        call gbrfs(this%d, this%f, this%ipiv, rhs, x, kl=this%nlower, trans=trans, info=info)
        if (info /= 0) call program_error("Band Matrix Iterative Refinement error in LAPACK", code = info)
      end if
    m4_divert(m4_ifdef({m4_spike},0,-1))
    case (BSOLVER_SPIKE)
      call spike_solve(this%sp, this%nrows, this%nlower, this%nupper, this%f, this%ipiv, x)
    m4_divert(0)
    end select
  end if
end subroutine

module subroutine band_{}T{}_output(this, file, funit, fmt)
  !! Write band matrix to file or stream.

  class(band_{}T),        intent(in) :: this
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
    fmt_ = m4_ifelse(T,cmplx,{'(ES25.16E3,SP,ES25.16E3,"i")'},{'(ES25.16E3)'})
  end if

  ! open file if needed
  m4_assert(.not. (present(file) .and. present(funit)))
  m4_assert(present(file) .or. present(funit))
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

module subroutine band_{}T{}_add_elem(this, i, j, x)
  class(band_{}T), intent(inout) :: this
  integer,         intent(in)    :: i
    !! row index of full matrix to represent
  integer,         intent(in)    :: j
    !! col index of full matrix to represent
  m4_type(T),      intent(in)    :: x
    !! matrix element to add

  m4_assert((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  m4_assert((           1 <= j  ) .and. (  j <= this%ncols ))

  this%d(i-j,j) = this%get_elem(i,j) + x
end subroutine

module function band_{}T{}_get_elem(this, i, j) result(x)
  class(band_{}T), intent(in) :: this
  integer,         intent(in) :: i
    !! row index of full matrix to represent
  integer,         intent(in) :: j
    !! col index of full matrix to represent
  m4_type(T)                  :: x
    !! matrix element

  m4_assert((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  m4_assert((           1 <= j  ) .and. (  j <= this%ncols ))

  x = this%d(i-j,j)
end function

module subroutine band_{}T{}_set_elem(this, i, j, x)
  class(band_{}T), intent(inout) :: this
  integer,         intent(in)    :: i
    !! row index of full matrix to represent
  integer,         intent(in)    :: j
    !! col index of full matrix to represent
  m4_type(T),      intent(in)    :: x
    !! matrix element to set

  m4_assert((-this%nupper <= i-j) .and. (i-j <= this%nlower))
  m4_assert((           1 <= j  ) .and. (  j <= this%ncols ))

  this%d(i-j,j) = x
end subroutine

module subroutine band_{}T{}_set_diag_val(this, val, k)
  !! inserts same value on the k-th diagonal.

  class(band_{}T),   intent(inout) :: this
  m4_type(T),        intent(in)    :: val
    !! value
  integer, optional, intent(in)    :: k
    !! k-th diagonal.
    !! k>0 lower diagonal, k<0 upper diagonal.
    !! Am4_type(T)ENTION: opposite to numbering in matlab!
    !! default: 0 == main diagonal

  integer :: k_, i

  k_ = 0
  if (present(k)) k_ = k

  ! create array of length: n-|k|
  call this%set_diag([(val, i=1,this%nrows-abs(k_))], k=k_)
end subroutine

module subroutine band_{}T{}_set_diag_arr(this, vals, k)
  !! inserts values on the k-th diagonal.

  class(band_{}T),   intent(inout) :: this
  m4_type(T),        intent(in)    :: vals(:)
    !! values
  integer, optional, intent(in)    :: k
    !! k-th diagonal.
    !! k>0 lower diagonal, k<0 upper diagonal.
    !! Am4_type(T)ENTION: opposite to numbering in matlab!
    !! default: 0 == main diagonal

  integer :: k_, i1, i2, n

  k_ = 0
  if (present(k)) k_ = k

  m4_assert((-k_ <= this%nupper) .and. (k_ <= this%nlower))

  n  = this%nrows
  i1 = merge(1, 1-k_, k_>=0)
  i2 = merge(n, n-k_, k_<=0)

  this%d(k_,i1:i2) = vals
end subroutine

module function band_eye_{}T{}(nrows) result(eye)
  !! band identity matrix

  integer, intent(in) :: nrows
    !! Number of rows
  type(band_{}T)      :: eye

  call eye%init(nrows, 0)
  eye%d = 1
end function
