#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_PTR_TYPE      CONCAT(matrix_ptr,T)

#define MATOP_TYPE           CONCAT( matop,T)
#define MATOP_INIT           CONCAT3(matop,T,init)
#define MATOP_EXEC1          CONCAT3(matop,T,exec1)
#define MATOP_EXEC2          CONCAT3(matop,T,exec2)

#define SINGLE_MATOP_TYPE    CONCAT(single_matop,T)
#define SINGLE_MATOP_INIT    CONCAT3(single_matop,T,init)
#define SINGLE_MATOP_EXEC1   CONCAT3(single_matop,T,exec1)
#define SINGLE_MATOP_EXEC2   CONCAT3(single_matop,T,exec2)

#define CHAIN_MATOP_TYPE     CONCAT(chain_matop,T)
#define CHAIN_MATOP_INIT     CONCAT3(chain_matop,T,init)
#define CHAIN_MATOP_SET_MAT  CONCAT3(chain_matop,T,set_mat)
#define CHAIN_MATOP_EXEC1    CONCAT3(chain_matop,T,exec1)
#define CHAIN_MATOP_EXEC2    CONCAT3(chain_matop,T,exec2)

subroutine MATOP_INIT(this, dum, nrows, ncols)
  !! initialize matrix operation

  class(MATOP_TYPE), intent(out) :: this
  character(*),      intent(in)  :: dum
    !! Dummy, unused
  integer,           intent(in)  :: nrows
    !! number of output values
  integer, optional, intent(in)  :: ncols
    !! optional number of input values (default: ncols = nrows)

  IGNORE(dum)

  this%nrows = nrows
  this%ncols = nrows
  if (present(ncols)) this%ncols = ncols
end subroutine

subroutine SINGLE_MATOP_INIT(this, mat, inv)
  !! initialize single matrix operation

  class(SINGLE_MATOP_TYPE),   intent(out) :: this
  class(MATRIX_TYPE), target, intent(in)  :: mat
    !! matrix
  logical, optional,          intent(in)  :: inv
    !! optional inverse flag (default: false)

  ! init base
  call this%init("", mat%nrows, ncols = mat%ncols)

  ! set members
  this%mat => mat
  this%inv = .false.
  if (present(inv)) this%inv = inv
end subroutine

subroutine SINGLE_MATOP_EXEC1(this, x, y)
  !! apply single matrix operation to vector

  class(SINGLE_MATOP_TYPE), intent(inout) :: this
  TT,                       intent(in)    :: x(:)
    !! input vector
  TT,                       intent(out)   :: y(:)
    !! output vector

  ! apply matrix
  if (this%inv) then
    call this%mat%solve_vec(x, y)
  else
    call this%mat%mul_vec(x, y)
  end if
end subroutine

subroutine SINGLE_MATOP_EXEC2(this, x, y)
  !! apply single matrix operation to matrix

  class(SINGLE_MATOP_TYPE), intent(inout) :: this
  TT,                       intent(in)    :: x(:,:)
    !! input matrix
  TT,                       intent(out)   :: y(:,:)
    !! output matrix

  ! apply matrix
  if (this%inv) then
    call this%mat%solve_mat(x, y)
  else
    call this%mat%mul_mat(x, y)
  end if
end subroutine

subroutine CHAIN_MATOP_INIT(this, nmat, nrows, ncols)
  !! initialize chain matrix operation

  class(CHAIN_MATOP_TYPE), intent(out) :: this
  integer,                 intent(in)  :: nmat
    !! number of matrices chained together in this operation
  integer,                 intent(in)  :: nrows
    !! number of output values
  integer, optional,       intent(in)  :: ncols
    !! optional number of input values (default: ncols = nrows)

  ! init base
  call this%init("", nrows, ncols)

  ! allocate matrix pointers and inverse flags
  allocate (this%mat(nmat))
  allocate (this%inv(nmat), source = .false.)
end subroutine

subroutine CHAIN_MATOP_SET_MAT(this, imat, mat, inv)
  !! set matrix pointer and inverse flag

  class(CHAIN_MATOP_TYPE),    intent(inout) :: this
  integer,                    intent(in)    :: imat
    !! matrix index (mat(1) is applied first, mat(nmat) last)
  class(MATRIX_TYPE), target, intent(in)    :: mat
    !! applied matrix
  logical, optional,          intent(in)    :: inv
    !! inverse flag (default: false)

  this%mat(imat)%p => mat
  this%inv(imat) = .false.
  if (present(inv)) this%inv(imat) = inv
end subroutine

subroutine CHAIN_MATOP_EXEC1(this, x, y)
  !! apply chain matrix operation to vector

  class(CHAIN_MATOP_TYPE), intent(inout) :: this
  TT,                      intent(in)    :: x(:)
    !! input vector
  TT,                      intent(out)   :: y(:)
    !! output vector

  ! local variables
  integer         :: i, ntmp
  TT, allocatable :: tmp1(:), tmp2(:)

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

subroutine CHAIN_MATOP_EXEC2(this, x, y)
  !! apply chain matrix operation to matrix

  class(CHAIN_MATOP_TYPE), intent(inout) :: this
  TT,                      intent(in)    :: x(:,:)
    !! input matrix
  TT,                      intent(out)   :: y(:,:)
    !! output matrix

  ! local variables
  integer         :: i, ntmp, nrhs
  TT, allocatable :: tmp1(:,:), tmp2(:,:)

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

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef MATRIX_TYPE
#undef MATRIX_PTR_TYPE

#undef MATOP_TYPE
#undef MATOP_INIT
#undef MATOP_EXEC1
#undef MATOP_EXEC2

#undef SINGLE_MATOP_TYPE
#undef SINGLE_MATOP_INIT
#undef SINGLE_MATOP_EXEC1
#undef SINGLE_MATOP_EXEC2

#undef CHAIN_MATOP_TYPE
#undef CHAIN_MATOP_INIT
#undef CHAIN_MATOP_SET_MAT
#undef CHAIN_MATOP_EXEC1
#undef CHAIN_MATOP_EXEC2
