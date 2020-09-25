#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DENSE_TYPE           CONCAT(dense,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define TRIANG_TYPE          CONCAT(triang,T)

#define TRIANG_INIT          CONCAT3(triang,T,init)
#define TRIANG_DESTRUCT      CONCAT3(triang,T,destruct)
#define TRIANG_RESET         CONCAT3(triang,T,reset)
#define TRIANG_SCALE         CONCAT3(triang,T,scale)
#define TRIANG_MUL_VEC       CONCAT3(triang,T,mul_vec)
#define TRIANG_MUL_MAT_ADD   CONCAT3(triang,T,mul_mat_add)
#define TRIANG_MUL_MAT_SIDES CONCAT3(triang,T,mul_mat_sides)
#define TRIANG_FACTORIZE     CONCAT3(triang,T,factorize)
#define TRIANG_SOLVE_VEC     CONCAT3(triang,T,solve_vec)
#define TRIANG_SOLVE_MAT     CONCAT3(triang,T,solve_mat)
#define TRIANG_TO_DENSE      CONCAT3(triang,T,to_dense)
#define TRIANG_TO_SPARSE     CONCAT3(triang,T,to_sparse)
#define TRIANG_TO_REAL       CONCAT3(triang,T,to_real)
#define TRIANG_TO_CMPLX      CONCAT3(triang,T,to_cmplx)

subroutine TRIANG_INIT(this, n, upper)
  !! Initialize (square) triangular matrix

  class(TRIANG_TYPE), intent(out) :: this
    !! Triangular matrix
  integer,            intent(in)  :: n
    !! Number of rows/cols
  logical,            intent(in)  :: upper
    !! Upper or lower triangular matrix

  ! init base
  call this%init("", n)

  ! set upper/lower flag
  this%upper = upper

  ! allocate memory
  allocate (this%d(n,n))
  this%d = 0                ! inserts 0, regardless of type

  this%factorized = .true.
end subroutine

subroutine TRIANG_DESTRUCT(this)
  !! Destruct triangular matrix by deallocating its memory.

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix

  if (allocated(this%d)) deallocate (this%d)
end subroutine

subroutine TRIANG_RESET(this, only_factorization)
  !! Reset triangular matrix to zero

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix
  logical, optional,  intent(in)    :: only_factorization
    !! ignored

  IGNORE(only_factorization)

  this%d = 0
end subroutine

subroutine TRIANG_SCALE(this, fact)
  !! Scale triangular matrix by a factor.

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix
  TT,                 intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: i

  if (this%upper) then
    do i = 1, this%nrows
      this%d(i,i:this%nrows) = this%d(i,i:this%nrows) * fact
    end do
  else
    do i = 1, this%nrows
      this%d(i,1:i) = this%d(i,1:i) * fact
    end do
  end if
end subroutine

subroutine TRIANG_MUL_VEC(this, x, y, fact_y, trans)
  !! Multiply triangular matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(TRIANG_TYPE),  intent(in)    :: this
    !! triangular matrix
  TT,                  intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                  intent(inout) :: y(:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N') for Not/Transposed/Hermitian transp.

  character       :: uplo
  TT, allocatable :: y0(:)

  uplo = 'L'
  if (this%upper) uplo = 'U'

  ! BLAS doesnt offer fact_y feature...
  if (present(fact_y)) allocate (y0(size(y)), source = y)
  y = x
  call trmv(this%d, y, uplo=uplo, trans=trans)
  if (present(fact_y)) y = y + fact_y*y0
end subroutine

subroutine TRIANG_MUL_MAT_ADD(this, x, y, fact_y, trans)
  !! Multiply triangular matrix by a raw matrix.
  !! Performs operation: y <- this * x + fact_y * y
  !!
  !! Inefficient: does not make use of triangular structure.

  class(TRIANG_TYPE),  intent(in)    :: this
    !! Triangular matrix
  TT,                  intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                  intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! Scale existing data in y and add result to it (default: 0)
  character, optional, intent(in)    :: trans
    !! Transposed flag for this matrix, must be 'N', 'C' or 'T' (default: 'N')

  call gemm(this%d, x, y, transa=trans, beta=fact_y)
end subroutine

subroutine TRIANG_MUL_MAT_SIDES(this, y, fact, trans, right)
  !! Multiply triangular matrix by a raw matrix.
  !! Performs operation:
  !!
  !!        y <- { fact * op(x) * y    if right=false
  !!             { fact * y * op(x)    if right=true
  !!
  !!     op(x) = { x                   if trans=N
  !!             { x^T                 if trans=T
  !!             { x^H                 if trans=C

  class(TRIANG_TYPE),  intent(in)    :: this
    !! triangular matrix x
  TT,                  intent(inout) :: y(:,:)
    !! raw matrix y
  TT,        optional, intent(in)    :: fact
    !! Scale product (default: 1)
  character, optional, intent(in)    :: trans
    !! Transposed flag for this matrix, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)    :: right
    !! use x as right sided factor? (default: .false.)

  character :: uplo, side
  logical   :: right_

  right_ = .false.
  if (present(right)) right_ = right

  uplo = 'L'
  if (this%upper) uplo = 'U'
  side = 'L'
  if (right_) side = 'R'

  call trmm(this%d, y, side=side, uplo=uplo, transa=trans, alpha=fact)
end subroutine

subroutine TRIANG_FACTORIZE(this)
  !! Throws error: triangular matrix doesnt need to be factorized! One can solve with it.

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix

  ASSERT(.false.)   ! already triangular
  IGNORE(this)
end subroutine

subroutine TRIANG_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for a single right-hand side.
  !! Solves: this * x = rhs, where x and rhs are raw vectors.

  class(TRIANG_TYPE),  intent(in)  :: this
    !! Triangular matrix
  TT,                  intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                  intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Neglected: Perform iterative refinement. Error if .true. (default: .false.)

  character :: uplo

  IGNORE(iter_refine)
  ASSERT(.not. present(iter_refine))

  uplo = 'L'
  if (this%upper) uplo = 'U'

  x = rhs
  call trsv(this%d, x, uplo=uplo, trans=trans)
end subroutine

subroutine TRIANG_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for multiple right-hand sides.
  !! Solves: this * x = rhs, where x and rhs are raw matrices.

  class(TRIANG_TYPE),  intent(in)  :: this
    !! Triangular matrix
  TT,                  intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  TT,                  intent(out) :: x(:,:)
    !! Store results in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Neglected: Perform iterative refinement. Error if .true. (default: .false.)

  character :: uplo

  IGNORE(iter_refine)
  ASSERT(.not. present(iter_refine))

  uplo = 'L'
  if (this%upper) uplo = 'U'

  x = rhs
  call trsm(this%d, x, uplo=uplo, transa=trans)
end subroutine

subroutine TRIANG_TO_DENSE(this, d, i0, j0)
  !! Insert triangular matrix into a dense matrix.

  class(TRIANG_TYPE),  intent(in)    :: this
    !! triangular matrix (unchanged)
  type(DENSE_TYPE),    intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional,   intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional,   intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  integer :: i0_, j0_, n

  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  n = this%nrows
  d%d(i0_:i0_+n-1,j0_:j0_+n-1) = this%d
end subroutine

subroutine TRIANG_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert triangular matrix into a sparse matrix.

  class(TRIANG_TYPE), intent(in)    :: this
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
  integer :: i0_, j0_, i, j, jmin, jmax
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
    if (this%upper) then
      jmin = i
      jmax = this%nrows
    else
      jmin = 1
      jmax = i
    end if

    do j = jmin, jmax
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,j)) cycle
      elseif (drop_zeros_) then
        if (this%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+j-1, this%d(i,j), search=.false.)
    end do
  end do
end subroutine

#ifdef TCMPLX
subroutine TRIANG_TO_REAL(this, a)
  !! Convert triangular complex matrix to real by ignoring imaginary part.

  class(TRIANG_cmplx), intent(in)  :: this
    !! TRIANG complex matrix
  type(TRIANG_real),   intent(out) :: a
    !! Output real(this)

  call a%init(this%nrows, this%upper)
  a%d = real(this%d)
end subroutine
#else
subroutine TRIANG_TO_CMPLX(this, a)
  !! Convert triangular real matrix to complex.

  class(TRIANG_real), intent(in)  :: this
    !! Triangular real matrix
  type(TRIANG_cmplx), intent(out) :: a
    !! Output cmplx(this)

  call a%init(this%nrows, this%upper)
  a%d = this%d
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

#undef DENSE_TYPE
#undef SPBUILD_TYPE
#undef TRIANG_TYPE

#undef TRIANG_INIT
#undef TRIANG_DESTRUCT
#undef TRIANG_RESET
#undef TRIANG_SCALE
#undef TRIANG_MUL_VEC
#undef TRIANG_MUL_MAT_ADD
#undef TRIANG_MUL_MAT_SIDES
#undef TRIANG_FACTORIZE
#undef TRIANG_SOLVE_VEC
#undef TRIANG_SOLVE_MAT
#undef TRIANG_TO_DENSE
#undef TRIANG_TO_SPARSE
#undef TRIANG_TO_TRIANG
#undef TRIANG_TO_REAL
#undef TRIANG_TO_CMPLX
