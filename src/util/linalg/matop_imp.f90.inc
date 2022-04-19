subroutine matop_{}T{}_init(this, dum, nrows, ncols)
  !! initialize matrix operation

  class(matop_{}T),  intent(out) :: this
  character(*),      intent(in)  :: dum
    !! Dummy, unused
  integer,           intent(in)  :: nrows
    !! number of output values
  integer, optional, intent(in)  :: ncols
    !! optional number of input values (default: ncols = nrows)

  m4_ignore(dum)

  this%nrows = nrows
  this%ncols = nrows
  if (present(ncols)) this%ncols = ncols
end subroutine

subroutine matop_{}T{}_exec1(this, x, y)
  !! virtual exec1 procedure. just a placeholder s.t. all extending matops have to redefine it.

  class(matop_{}T),   intent(in)  :: this
  m4_type(T), target, intent(in)  :: x(:)
    !! input vector
  m4_type(T), target, intent(out) :: y(:)
    !! output vector

  m4_ignore(this)
  m4_ignore(x)
  y = 0 ! to remove gfortran warning for intent(out) variable
  call program_error("virtual parent method. shouldnt be called.")
end subroutine

subroutine matop_{}T{}_exec2(this, x, y)
  !! virtual exec2 procedure. just a placeholder s.t. all extending matops have to redefine it.

  class(matop_{}T),   intent(in)  :: this
  m4_type(T), target, intent(in)  :: x(:,:)
    !! input vector
  m4_type(T), target, intent(out) :: y(:,:)
    !! output vector

  m4_ignore(this)
  m4_ignore(x)
  y = 0 ! to remove gfortran warning for intent(out) variable
  call program_error("virtual parent method. shouldnt be called.")
end subroutine

subroutine single_matop_{}T{}_init(this, mat, inv)
  !! initialize single matrix operation

  class(SINGLE_matop_{}T),   intent(out) :: this
  class(matrix_{}T), target, intent(in)  :: mat
    !! matrix
  logical, optional,         intent(in)  :: inv
    !! optional inverse flag (default: false)

  ! init base
  call this%init("", mat%nrows, ncols = mat%ncols)

  ! set members
  this%mat => mat
  this%inv = .false.
  if (present(inv)) this%inv = inv
end subroutine

subroutine single_matop_{}T{}_exec1(this, x, y)
  !! apply single matrix operation to vector

  class(SINGLE_matop_{}T), intent(in)  :: this
  m4_type(T), target,      intent(in)  :: x(:)
    !! input vector
  m4_type(T), target,      intent(out) :: y(:)
    !! output vector

  ! apply matrix
  if (this%inv) then
    call this%mat%solve_vec(x, y)
  else
    call this%mat%mul_vec(x, y)
  end if
end subroutine

subroutine single_matop_{}T{}_exec2(this, x, y)
  !! apply single matrix operation to matrix

  class(SINGLE_matop_{}T), intent(in)  :: this
  m4_type(T), target,      intent(in)  :: x(:,:)
    !! input matrix
  m4_type(T), target,      intent(out) :: y(:,:)
    !! output matrix

  ! apply matrix
  if (this%inv) then
    call this%mat%solve_mat(x, y)
  else
    call this%mat%mul_mat(x, y)
  end if
end subroutine

subroutine chain_matop_{}T{}_init(this, nmat, nrows, ncols)
  !! initialize chain matrix operation

  class(CHAIN_matop_{}T), intent(out) :: this
  integer,                intent(in)  :: nmat
    !! number of matrices chained together in this operation
  integer,                intent(in)  :: nrows
    !! number of output values
  integer, optional,      intent(in)  :: ncols
    !! optional number of input values (default: ncols = nrows)

  ! init base
  call this%init("", nrows, ncols)

  ! allocate matrix pointers and inverse flags
  allocate (this%mat(nmat))
  allocate (this%inv(nmat), source = .false.)
end subroutine

subroutine chain_matop_{}T{}_set_mat(this, imat, mat, inv)
  !! set matrix pointer and inverse flag

  class(CHAIN_matop_{}T),    intent(inout) :: this
  integer,                   intent(in)    :: imat
    !! matrix index (mat(1) is applied first, mat(nmat) last)
  class(matrix_{}T), target, intent(in)    :: mat
    !! applied matrix
  logical, optional,         intent(in)    :: inv
    !! inverse flag (default: false)

  this%mat(imat)%p => mat
  this%inv(imat) = .false.
  if (present(inv)) this%inv(imat) = inv
end subroutine

subroutine chain_matop_{}T{}_exec1(this, x, y)
  !! apply chain matrix operation to vector

  class(CHAIN_matop_{}T), intent(in)  :: this
  m4_type(T), target,     intent(in)  :: x(:)
    !! input vector
  m4_type(T), target,     intent(out) :: y(:)
    !! output vector

  ! local variables
  integer                 :: i, ntmp
  m4_type(T), allocatable :: tmp1(:), tmp2(:)

  ! get temporary space size
  ntmp = 0
  do i = 1, size(this%mat)
    if (this%mat(i)%p%nrows > ntmp) then
      ntmp = this%mat(i)%p%nrows
    end if
    if (this%mat(i)%p%ncols > ntmp) then
      ntmp = this%mat(i)%p%ncols
    end if
  end do

  ! allocate and init temporary space
  allocate (tmp1(ntmp), tmp2(ntmp))
  tmp1(1:size(x)) = x

  ! apply matrices
  do i = 1, size(this%mat)
    if (this%inv(i)) then
      call this%mat(i)%p%solve_vec(tmp1(1:this%mat(i)%p%nrows), tmp2(1:this%mat(i)%p%ncols))
      tmp1(1:this%mat(i)%p%ncols) = tmp2(1:this%mat(i)%p%ncols)
    else
      call this%mat(i)%p%mul_vec(tmp1(1:this%mat(i)%p%ncols), tmp2(1:this%mat(i)%p%nrows))
      tmp1(1:this%mat(i)%p%nrows) = tmp2(1:this%mat(i)%p%nrows)
    end if
  end do

  ! output
  y(1:size(y)) = tmp1(1:size(y))
end subroutine

subroutine chain_matop_{}T{}_exec2(this, x, y)
  !! apply chain matrix operation to matrix

  class(CHAIN_matop_{}T), intent(in)  :: this
  m4_type(T), target,     intent(in)  :: x(:,:)
    !! input matrix
  m4_type(T), target,     intent(out) :: y(:,:)
    !! output matrix

  ! local variables
  integer                 :: i, ntmp, nrhs
  m4_type(T), allocatable :: tmp1(:,:), tmp2(:,:)

  ! get temporary space size
  ntmp = 0
  do i = 1, size(this%mat)
    if (this%mat(i)%p%nrows > ntmp) then
      ntmp = this%mat(i)%p%nrows
    end if
    if (this%mat(i)%p%ncols > ntmp) then
      ntmp = this%mat(i)%p%ncols
    end if
  end do
  nrhs = size(x, 2)

  ! allocate and init temporary space
  allocate (tmp1(ntmp,nrhs), tmp2(ntmp,nrhs))
  tmp1(1:size(x,1),:) = x

  ! apply matrices
  do i = 1, size(this%mat)
    if (this%inv(i)) then
      call this%mat(i)%p%solve_mat(tmp1(1:this%mat(i)%p%nrows,:), tmp2(1:this%mat(i)%p%ncols,:))
      tmp1(1:this%mat(i)%p%ncols,:) = tmp2(1:this%mat(i)%p%ncols,:)
    else
      call this%mat(i)%p%mul_mat(tmp1(1:this%mat(i)%p%ncols,:), tmp2(1:this%mat(i)%p%nrows,:))
      tmp1(1:this%mat(i)%p%nrows,:) = tmp2(1:this%mat(i)%p%nrows,:)
    end if
  end do

  ! output
  y(1:size(y),:) = tmp1(1:size(y),:)
end subroutine
