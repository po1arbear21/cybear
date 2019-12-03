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

#define HESSENBERG_INIT      CONCAT3(hessenberg,T,init)
#define HESSENBERG_DESTRUCT  CONCAT3(hessenberg,T,destruct)
#define HESSENBERG_RESET     CONCAT3(hessenberg,T,reset)
#define HESSENBERG_SCALE     CONCAT3(hessenberg,T,scale)
#define HESSENBERG_MUL_VEC   CONCAT3(hessenberg,T,mul_vec)
#define HESSENBERG_MUL_MAT   CONCAT3(hessenberg,T,mul_mat)
#define HESSENBERG_FACTORIZE CONCAT3(hessenberg,T,factorize)
#define HESSENBERG_SOLVE_VEC CONCAT3(hessenberg,T,solve_vec)
#define HESSENBERG_SOLVE_MAT CONCAT3(hessenberg,T,solve_mat)
#define HESSENBERG_TO_DENSE  CONCAT3(hessenberg,T,to_dense)
#define HESSENBERG_TO_SPARSE CONCAT3(hessenberg,T,to_sparse)
#ifdef TCMPLX
#define HESSENBERG_TO_REAL   CONCAT3(hessenberg,T,to_real)
#else
#define HESSENBERG_TO_CMPLX  CONCAT3(hessenberg,T,to_cmplx)
#endif

subroutine HESSENBERG_INIT(this, nrows, upper)
  !! Initialize hessenberg matrix

  class(HESSENBERG_TYPE), intent(out) :: this
    !! Hessenberg matrix
  integer,                intent(in)  :: nrows
    !! Number of rows
  logical,                intent(in)  :: upper
    !! Upper or lower hessenberg matrix

  ! init base
  call this%init("", nrows)

  ! set upper/lower flag
  this%upper = upper

  ! allocate memory
  allocate (this%d(nrows,nrows))
  this%d = 0
end subroutine

subroutine HESSENBERG_DESTRUCT(this)
  !! Destruct hessenberg matrix by deallocating its memory.

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix

  if (allocated(this%d)) deallocate(this%d)
end subroutine

subroutine HESSENBERG_RESET(this, only_factorization)
  !! Reset hessenberg matrix to zero (optionally only the factorization).

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix
  logical, optional,      intent(in)    :: only_factorization
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

subroutine HESSENBERG_SCALE(this, fact)
  !! Scale hessenberg matrix by a factor.

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix
  TT,                     intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: i, j

  ASSERT(.not. allocated(this%f))
  ASSERT(.not. allocated(this%ipiv))

  this%d = this%d * fact

  if (this%upper) then
    do i = 1, this%nrows
      j = max(i-1,1)
      this%d(i,j:this%nrows) = this%d(i,j:this%nrows) * fact
    end do
  else
    do i = 1, this%nrows
      j = min(i+1,this%nrows)
      this%d(i,1:j) = this%d(i,1:j) * fact
    end do
  end if
end subroutine

subroutine HESSENBERG_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply hessenberg matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(HESSENBERG_TYPE),     intent(in)    :: this
    !! Hessenberg matrix
  TT,                         intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                         intent(inout) :: y(:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer          :: i
  character(len=1) :: uplo, trans_
  TT, allocatable  :: tmp(:)

  uplo = 'L'
  if (this%upper) uplo = 'U'
  trans_ = 'N'
  if (present(trans)) trans_ = trans

#ifdef TCMPLX
  ASSERT((trans_ == 'N') .or. (trans_ == 'T') .or. (trans_ == 'C'))
#else
  ASSERT((trans_ == 'N') .or. (trans_ == 'T'))
#endif

  ! triangular part
  allocate (tmp(size(x)), source = x)
  call trmv(this%d, tmp, uplo = uplo, trans = trans_)

  ! off-diagonal
  if (trans_ == 'N') then
    if (this%upper) then
      do i = 2, this%nrows
        tmp(i) = tmp(i) + this%d(i,i-1) * x(i-1)
      end do
    else
      do i = 1, this%nrows - 1
        tmp(i) = tmp(i) + this%d(i,i+1) * x(i+1)
      end do
    end if
  elseif (trans_ == 'T') then
    if (this%upper) then
      do i = 1, this%nrows - 1
        tmp(i) = tmp(i) + this%d(i+1,i) * x(i+1)
      end do
    else
      do i = 2, this%nrows
        tmp(i) = tmp(i) + this%d(i-1,i) * x(i-1)
      end do
    end if
#ifdef TCMPLX
  else ! trans_ == 'C'
    if (this%upper) then
      do i = 1, this%nrows - 1
        tmp(i) = tmp(i) + conjg(this%d(i+1,i)) * x(i+1)
      end do
    else
      do i = 2, this%nrows
        tmp(i) = tmp(i) + conjg(this%d(i-1,i)) * x(i-1)
      end do
    end if
#endif
  end if

  if (present(fact_y)) then
    y = fact_y * y + tmp
  else
    y = tmp
  end if
end subroutine

subroutine HESSENBERG_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply hessenberg matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(HESSENBERG_TYPE),     intent(in)    :: this
    !! Hessenberg matrix
  TT,                         intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                         intent(inout) :: y(:,:)
    !! Store result in y
  TT,               optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character(len=1), optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer          :: i
  character(len=1) :: uplo, trans_
  TT, allocatable  :: tmp(:,:)

  uplo = 'L'
  if (this%upper) uplo = 'U'
  trans_ = 'N'
  if (present(trans)) trans_ = trans

#ifdef TCMPLX
  ASSERT((trans_ == 'N') .or. (trans_ == 'T') .or. (trans_ == 'C'))
#else
  ASSERT((trans_ == 'N') .or. (trans_ == 'T'))
#endif

  ! triangular part
  allocate (tmp(size(x,1),size(x,2)), source = x)
  call trmm(this%d, tmp, uplo = uplo, transa = trans_)

  ! off-diagonal
  if (trans_ == 'N') then
    if (this%upper) then
      do i = 2, this%nrows
        tmp(i,:) = tmp(i,:) + this%d(i,i-1) * x(i-1,:)
      end do
    else
      do i = 1, this%nrows - 1
        tmp(i,:) = tmp(i,:) + this%d(i,i+1) * x(i+1,:)
      end do
    end if
  elseif (trans_ == 'T') then
    if (this%upper) then
      do i = 1, this%nrows - 1
        tmp(i,:) = tmp(i,:) + this%d(i+1,i) * x(i+1,:)
      end do
    else
      do i = 2, this%nrows
        tmp(i,:) = tmp(i,:) + this%d(i-1,i) * x(i-1,:)
      end do
    end if
#ifdef TCMPLX
  else ! trans_ == 'C'
    if (this%upper) then
      do i = 1, this%nrows - 1
        tmp(i,:) = tmp(i,:) + conjg(this%d(i+1,i)) * x(i+1,:)
      end do
    else
      do i = 2, this%nrows
        tmp(i,:) = tmp(i,:) + conjg(this%d(i-1,i)) * x(i-1,:)
      end do
    end if
#endif
  end if

  if (present(fact_y)) then
    y = fact_y * y + tmp
  else
    y = tmp
  end if
end subroutine

subroutine HESSENBERG_FACTORIZE(this)
  !! Factorize hessenberg matrix

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix

  ! local variables
  integer :: i, n, itmp

  ! system size
  n = this%nrows

  ! allocate factorization if necessary
  if (.not. allocated(this%f)) allocate(this%f(n,n))
  if (.not. allocated(this%ipiv)) allocate(this%ipiv(n))

  ! init permutation vector
  do i = 1, n
    this%ipiv(i) = i
  end do

  ! init factorization by original matrix
  this%f = this%d

  if (this%upper) then ! LU decomposition
    do i = 1, n - 1
      ! partial pivoting
      if (abs(this%f(i,i)) < abs(this%f(i+1,i))) then
        ! swap (i)-th and (i+1)-th row
        call swap(this%f(i,i:n), this%f(i+1,i:n))

        ! update permutation vector
        itmp           = this%ipiv(i  )
        this%ipiv(i  ) = this%ipiv(i+1)
        this%ipiv(i+1) = itmp
      end if

      ! update L, U
      this%f(i+1,i    ) =                     this%f(i+1,i) / this%f(i,i    )
      this%f(i+1,i+1:n) = this%f(i+1,i+1:n) - this%f(i+1,i) * this%f(i,i+1:n)
    end do
  else ! UL decomposition
    do i = n, 2, -1
      ! partial pivoting
      if (abs(this%f(i,i)) < abs(this%f(i-1,i))) then
        ! swap (i)-th and (i-1)-th row
        call swap(this%f(i,1:i), this%f(i-1,1:i))

        ! update permutation vector
        itmp           = this%ipiv(i  )
        this%ipiv(i  ) = this%ipiv(i-1)
        this%ipiv(i-1) = itmp
      end if

      ! update U, L
      this%f(i-1,i    ) =                     this%f(i-1,i) / this%f(i,i    )
      this%f(i-1,1:i-1) = this%f(i-1,1:i-1) - this%f(i-1,i) * this%f(i,1:i-1)
    end do
  end if
end subroutine

subroutine HESSENBERG_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first!

  class(HESSENBERG_TYPE),     intent(in)  :: this
    !! Hessenberg matrix
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

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed hessenberg matrix not implemented!")
  end if

  ! solve with matrix routine
  rhs_(:,1) = rhs
  call this%solve_mat(rhs_, x_, trans = trans, iter_refine = iter_refine)
  x = x_(:,1)
end subroutine

subroutine HESSENBERG_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(HESSENBERG_TYPE),     intent(in)  :: this
    !! Hessenberg matrix
  TT,                         intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                         intent(out) :: x(:,:)
    !! Store results in x
  character(len=1), optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,          optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  integer :: i, j, n, nrhs
  logical :: iter_refine_

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed hessenberg matrix not implemented!")
  end if

  ! system size
  n = this%nrows

  ! number of right-hand sides
  nrhs = size(rhs, 2)

  ! init solution vector
  x = rhs(this%ipiv,:)

  if (this%upper) then
    ! solve for L
    do i = 2, n
      x(i,:) = x(i,:) - this%f(i,i-1) * x(i-1,:)
    end do

    ! solve for U
    call trsm(this%f, x, uplo = 'U')
  else
    ! solve for U
    do i = n-1, 1, -1
      x(i,:) = x(i,:) - this%f(i,i+1) * x(i+1,:)
    end do

    ! solve for L
    call trsm(this%f, x, uplo = 'L')
  end if

  ! iterative refinement
  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine
  if (iter_refine_) then
#ifdef TCMPLX
    block
      real,    allocatable :: a(:), b(:), rr(:,:), ri(:,:)
      complex, allocatable :: r(:,:), d(:,:)

      ! allocate temporary memory
      allocate (a(0:2*n+1), b(0:2*n+1), rr(n,nrhs), ri(n,nrhs))
      allocate (r(n,nrhs), d(n,nrhs))

      ! calculate residual with high precision
      a(0) = -1d0
      a(1) = 0d0
      if (this%upper) then
        do j = 1, nrhs
          do i = 1, n
            ! rr = ar * br + ai * (- bi)
            a(2:2*min(n-i+2,n)  :2) = real(this%d(i,max(i-1,1):n))
            a(3:2*min(n-i+2,n)+1:2) = imag(this%d(i,max(i-1,1):n))
            b(0) = real(rhs(i,j))
            b(1) = -imag(rhs(i,j))
            b(2:2*min(n-i+2,n)  :2) = real(x(max(i-1,1):n,j))
            b(3:2*min(n-i+2,n)+1:2) = -imag(x(max(i-1,1):n,j))
            rr(i,j) = hp_dot(a(0:2*min(n-i+2,n)), b(0:2*min(n-i+2,n)))

            ! ri = ar * bi + ai * br
            b(0) = imag(rhs(i,j))
            b(1) = real(rhs(i,j))
            b(2:2*min(n-i+2,n)  :2) = imag(x(max(i-1,1):n,j))
            b(3:2*min(n-i+2,n)+1:2) = real(x(max(i-1,1):n,j))
            ri(i,j) = hp_dot(a(0:2*min(n-i+2,n)), b(0:2*min(n-i+2,n)))

            ! r = rr + 1i * ri
            r(i,j) = cmplx(rr(i,j), ri(i,j))
          end do
        end do
      else
        do j = 1, nrhs
          do i = 1, n
            ! rr = ar * br - ai * bi
            a(2:2*min(i+1,n)  :2) = real(this%d(i,1:min(i+1,n)))
            a(3:2*min(i+1,n)+1:2) = imag(this%d(i,1:min(i+1,n)))
            b(0) = real(rhs(i,j))
            b(1) = -imag(rhs(i,j))
            b(2:2*min(i+1,n)  :2) = real(x(1:min(i+1,n),j))
            b(3:2*min(i+1,n)+1:2) = -imag(x(1:min(i+1,n),j))
            rr(i,j) = hp_dot(a(0:2*min(i+1,n)), b(0:2*min(i+1,n)))

            ! ri = ar * bi + ai * br
            b(0) = imag(rhs(i,j))
            b(1) = real(rhs(i,j))
            b(2:2*min(i+1,n)  :2) = imag(x(1:min(i+1,n),j))
            b(3:2*min(i+1,n)+1:2) = real(x(1:min(i+1,n),j))
            rr(i,j) = hp_dot(a(0:2*min(i+1,n)), b(0:2*min(i+1,n)))

            ! r = rr + 1i * ri
            r(i,j) = cmplx(rr(i,j), ri(i,j))
          end do
        end do
      end if

      ! solve system again for residual as right-hand side
      call this%solve_mat(r, d, iter_refine = .false.)

      ! update solution
      x = x - d
    end block
#else
    block
      real, allocatable :: a(:), b(:), r(:,:), d(:,:)

      ! allocate temporary memory
      allocate (a(0:n), b(0:n), r(n,nrhs), d(n,nrhs))

      ! calculate residual with high precision
      a(0) = -1d0
      if (this%upper) then
        do j = 1, nrhs
          do i = 1, n
            a(1:min(n-i+2,n)) = this%d(i,max(i-1,1):n)
            b(0) = rhs(i,j)
            b(1:min(n-i+2,n)) = x(max(i-1,1):n,j)
            r(i,j) = hp_dot(a(0:min(n-i+2,n)), b(0:min(n-i+2,n)))
          end do
        end do
      else
        do j = 1, nrhs
          do i = 1, n
            a(1:min(i+1,n)) = this%d(i,1:min(i+1,n))
            b(0) = rhs(i,j)
            b(1:min(i+1,n)) = x(1:min(i+1,n),j)
            r(i,j) = hp_dot(a(0:min(i+1,n)), b(0:min(i+1,n)))
          end do
        end do
      end if

      ! solve system again for residual as right-hand side
      call this%solve_mat(r, d, iter_refine = .false.)

      ! update solution
      x = x - d
    end block
#endif
  end if
end subroutine

subroutine HESSENBERG_TO_DENSE(this, d, i0, j0)
  !! Insert hessenberg matrix into a dense matrix.

  class(HESSENBERG_TYPE), intent(in)    :: this
    !! Hessenberg matrix (unchanged)
  type(DENSE_TYPE),       intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional,      intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional,      intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, i, n

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! system size
  n = this%nrows

  ! insert elements into dense matrix
  if (this%upper) then
    d%d(i0_,j0_:j0_+n-1) = this%d(1,1:n)
    do i = 2, n
      d%d(i0_+i-1,j0_+i-2:j0_+n-1) = this%d(i,i-1:n)
    end do
  else
    do i = 1, n - 1
      d%d(i0_+i-1,j0_:j0_+i) = this%d(i,1:i+1)
    end do
    d%d(i0_+n-1,j0_:j0_+n-1) = this%d(n,1:n)
  end if
end subroutine

subroutine HESSENBERG_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert hessenberg matrix into a sparse matrix.

  class(HESSENBERG_TYPE), intent(in)    :: this
    !! Hessenberg matrix (unchanged)
  type(SPBUILD_TYPE),     intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional,      intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional,      intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional,      intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional,      intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, i, j, n
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! system size
  n = this%nrows

  ! insert values
  if (this%upper) then
    do i = 1, n
      do j = max(i-1,1), n
        ! check structure
        if (present(struct)) then
          if (.not. struct(i,j)) cycle
        elseif (drop_zeros_) then
          if (this%d(i,j) == 0) cycle
        end if

        ! insert value without searching (append to rows)
        call sb%set(i0_+i-1,j0_+j-1, this%d(i,j), search = .false.)
      end do
    end do
  else
    do i = 1, n
      do j = 1, min(i+1,n)
        ! check structure
        if (present(struct)) then
          if (.not. struct(i,j)) cycle
        elseif (drop_zeros_) then
          if (this%d(i,j) == 0) cycle
        end if

        ! insert value without searching (append to rows)
        call sb%set(i0_+i-1,j0_+j-1, this%d(i,j), search = .false.)
      end do
    end do
  end if
end subroutine

#ifdef TCMPLX
subroutine HESSENBERG_TO_REAL(this, h)
  !! Convert hessenberg complex matrix to real by ignoring imaginary part.

  class(hessenberg_cmplx), intent(in)  :: this
    !! Hessenberg complex matrix
  class(hessenberg_real),  intent(out) :: h
    !! Output real(this)

  ! local variables
  integer :: i, j

  call h%init(this%nrows, this%upper)

  if (this%upper) then
    do i = 1, this%nrows
      j = max(i-1,1)
      h%d(i,j:this%nrows) = real(this%d(i,j:this%nrows))
    end do
  else
    do i = 1, this%nrows
      j = min(i+1,this%nrows)
      h%d(i,1:j) = real(this%d(i,1:j))
    end do
  end if
end subroutine
#else
subroutine HESSENBERG_TO_CMPLX(this, h)
  !! Convert hessenberg real matrix to complex.

  class(hessenberg_real),  intent(in)  :: this
    !! Hessenberg real matrix
  class(hessenberg_cmplx), intent(out) :: h
    !! Output cmplx(this)

  ! local variables
  integer :: i, j

  call h%init(this%nrows, this%upper)

  if (this%upper) then
    do i = 1, this%nrows
      j = max(i-1,1)
      h%d(i,j:this%nrows) = cmplx(this%d(i,j:this%nrows))
    end do
  else
    do i = 1, this%nrows
      j = min(i+1,this%nrows)
      h%d(i,1:j) = cmplx(this%d(i,1:j))
    end do
  end if
end subroutine
#endif

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

#undef HESSENBERG_INIT
#undef HESSENBERG_DESTRUCT
#undef HESSENBERG_RESET
#undef HESSENBERG_SCALE
#undef HESSENBERG_MUL_VEC
#undef HESSENBERG_MUL_MAT
#undef HESSENBERG_FACTORIZE
#undef HESSENBERG_SOLVE_VEC
#undef HESSENBERG_SOLVE_MAT
#undef HESSENBERG_TO_DENSE
#undef HESSENBERG_TO_SPARSE
#undef HESSENBERG_TO_REAL
#undef HESSENBERG_TO_CMPLX