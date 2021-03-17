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

#define HESSENBERG_INIT      CONCAT3(hessenberg,T,init)
#define HESSENBERG_DESTRUCT  CONCAT3(hessenberg,T,destruct)
#define HESSENBERG_RESET     CONCAT3(hessenberg,T,reset)
#define HESSENBERG_SCALE     CONCAT3(hessenberg,T,scale)
#define HESSENBERG_MUL_VEC   CONCAT3(hessenberg,T,mul_vec)
#define HESSENBERG_MUL_MAT   CONCAT3(hessenberg,T,mul_mat)
#define HESSENBERG_FACTORIZE CONCAT3(hessenberg,T,factorize)
#define HESSENBERG_SOLVE_VEC CONCAT3(hessenberg,T,solve_vec)
#define HESSENBERG_SOLVE_MAT CONCAT3(hessenberg,T,solve_mat)
#define HESSENBERG_IS_TRIANG CONCAT3(hessenberg,T,is_triang)

module subroutine HESSENBERG_INIT(this, nrows, upper)
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

module subroutine HESSENBERG_DESTRUCT(this)
  !! Destruct hessenberg matrix by deallocating its memory.

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix

  if (allocated(this%d   )) deallocate (this%d   )
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine HESSENBERG_RESET(this, only_factorization)
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
  if (allocated(this%f   )) deallocate (this%f   )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine HESSENBERG_SCALE(this, fact)
  !! Scale hessenberg matrix by a factor.

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix
  TT,                     intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: i, j

  ASSERT(.not. allocated(this%f))
  ASSERT(.not. allocated(this%ipiv))

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

module subroutine HESSENBERG_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply hessenberg matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(HESSENBERG_TYPE), intent(in)    :: this
    !! Hessenberg matrix
  TT,                     intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                     intent(inout) :: y(:)
    !! Store result in y
  TT,        optional,    intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional,    intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer         :: i
  character       :: uplo, trans_
  TT, allocatable :: tmp(:)

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

module subroutine HESSENBERG_MUL_MAT(this, x, y, fact_y, trans)
  !! Multiply hessenberg matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(HESSENBERG_TYPE), intent(in)    :: this
    !! Hessenberg matrix
  TT,                     intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                     intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional,    intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional,    intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  integer         :: i
  character       :: uplo, trans_
  TT, allocatable :: tmp(:,:)

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

module subroutine HESSENBERG_FACTORIZE(this)
  !! Factorize hessenberg matrix

  class(HESSENBERG_TYPE), intent(inout) :: this
    !! Hessenberg matrix

  ! local variables
  integer :: i, n, itmp

  ! system size
  n = this%nrows

  ! allocate factorization if necessary
  if (.not. allocated(this%f   )) allocate (this%f(   n,n))
  if (.not. allocated(this%ipiv)) allocate (this%ipiv(n  ))

  ! init permutation vector
  this%ipiv = [(i, i=1,n)]

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

  this%factorized = .true.
end subroutine

module subroutine HESSENBERG_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first!

  class(HESSENBERG_TYPE), intent(in)  :: this
    !! Hessenberg matrix
  TT,                     intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                     intent(out) :: x(:)
    !! Store result in x
  character, optional,    intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional,    intent(in)  :: iter_refine
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

module subroutine HESSENBERG_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(HESSENBERG_TYPE), intent(in)  :: this
    !! Hessenberg matrix
  TT,                     intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                     intent(out) :: x(:,:)
    !! Store results in x
  character, optional,    intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional,    intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  ! local variables
  integer :: i, j, n, nrhs
  logical :: iter_refine_

  ASSERT(this%factorized)

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
      a(0) = -1.0
      a(1) = 0.0
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
      a(0) = -1.0
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

logical module function HESSENBERG_IS_TRIANG(this, abs_tol)
  !! Checks if Hessenberg matrix is a triangular matrix.
  !! Only checks if offdiagonal part across Hessenerg's triangular part is zero.

  class(HESSENBERG_TYPE), intent(in) :: this
    !! Hessenberg matrix
  real, optional,         intent(in) :: abs_tol
    !! Absolute value's threshold (default: 0.0).

  integer :: i
  real    :: abs_tol_

  abs_tol_ = 0.0
  if (present(abs_tol)) abs_tol_ = abs_tol

  if (this%upper) then
    HESSENBERG_IS_TRIANG = all([(abs(this%d(i,i-1)) <= abs_tol_, i=2, this%nrows)])
  else
    HESSENBERG_IS_TRIANG = all([(abs(this%d(i,i+1)) <= abs_tol_, i=1, this%nrows-1)])
  end if
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

#undef HESSENBERG_INIT
#undef HESSENBERG_DESTRUCT
#undef HESSENBERG_RESET
#undef HESSENBERG_SCALE
#undef HESSENBERG_MUL_VEC
#undef HESSENBERG_MUL_MAT
#undef HESSENBERG_FACTORIZE
#undef HESSENBERG_SOLVE_VEC
#undef HESSENBERG_SOLVE_MAT
#undef HESSENBERG_IS_TRIANG
