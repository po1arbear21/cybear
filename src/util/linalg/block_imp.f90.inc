#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)
#define MATRIX_TYPE          CONCAT(matrix,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define BAND_TYPE            CONCAT(band,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define TRIANG_TYPE          CONCAT(triang,T)
#define BLOCK_TYPE           CONCAT(block,T)

#define BLOCK_INIT           CONCAT3(block,T,init)
#define BLOCK_DESTRUCT       CONCAT3(block,T,destruct)
#define BLOCK_RESET          CONCAT3(block,T,reset)
#define BLOCK_SCALE          CONCAT3(block,T,scale)
#define BLOCK_MUL_VEC        CONCAT3(block,T,mul_vec)
#define BLOCK_MUL_MAT        CONCAT3(block,T,mul_mat)
#define BLOCK_FACTORIZE      CONCAT3(block,T,factorize)
#define BLOCK_SOLVE_VEC      CONCAT3(block,T,solve_vec)
#define BLOCK_SOLVE_MAT      CONCAT3(block,T,solve_mat)
#define BLOCK_TO_DENSE       CONCAT3(block,T,to_dense)
#define BLOCK_TO_SPARSE      CONCAT3(block,T,to_sparse)
#define BLOCK_GET_GLOBAL_IDX CONCAT3(block,T,get_global_idx)
#define BLOCK_GET_LOCAL_IDX  CONCAT3(block,T,get_local_idx)
#define BLOCK_SET_BLOCK_TYPE CONCAT3(block,T,set_block_type)
#define BLOCK_GET_DENSE      CONCAT3(block,T,get_dense)
#define BLOCK_GET_SPARSE     CONCAT3(block,T,get_sparse)
#define BLOCK_GET_BAND       CONCAT3(block,T,get_band)
#define BLOCK_GET_BLOCK      CONCAT3(block,T,get_block)
#define BLOCK_SET            CONCAT3(block,T,set)

subroutine BLOCK_INIT(this, row_dim, col_dim, diag, tridiag)
  !! Initialize Block matrix.

  class(BLOCK_TYPE), intent(out) :: this
    !! Block matrix
  integer,           intent(in)  :: row_dim(:)
    !! Number of Rows for each block
  integer, optional, intent(in)  :: col_dim(:)
    !! Number of Columns for each block (default: row_dim)
  logical, optional, intent(in)  :: diag
    !! diagonal block matrix flag (default: false)
  logical, optional, intent(in)  :: tridiag
    !! tridiagonal block matrix flag (default: false)

  ! local variables
  integer, allocatable :: col_dim_(:)
  integer              :: ib, i1, jb, j1

  col_dim_ = row_dim
  if (present(col_dim)) col_dim_ = col_dim

  ! init base
  call this%init("", sum(row_dim), ncols=sum(col_dim_))

  ! number of blocks
  this%nbrows = size(row_dim)
  this%nbcols = size(col_dim_)

  ! set row indices
  allocate (this%i0(this%nbrows))
  allocate (this%i1(this%nbrows))
  i1 = 0
  do ib = 1, this%nbrows
    this%i0(ib) = i1 + 1
    this%i1(ib) = i1 + row_dim(ib)
    i1 = this%i1(ib)
  end do

  ! set col indices
  allocate (this%j0(this%nbcols))
  allocate (this%j1(this%nbcols))
  j1 = 0
  do jb = 1, this%nbcols
    this%j0(jb) = j1 + 1
    this%j1(jb) = j1 + col_dim_(jb)
    j1 = this%j1(jb)
  end do

  ! allocate data
  allocate (this%b(this%nbrows,this%nbcols))

  ! diag, tridiag
  this%diag = .false.
  if (present(diag)) this%diag = diag
  this%tridiag = .false.
  if (present(tridiag)) this%tridiag = tridiag

  if (this%diag .and. this%tridiag) call program_error("Block matrix cannot be both diagonal and tridiagonal!")
  if (this%diag .or.  this%tridiag) then
    if (this%nbrows /= this%nbcols) call program_error("(Tri)Diagonal Block matrix must have same number of blocks per row and column!")
  end if

  ! not factorized
  this%factorized = .false.
end subroutine

subroutine BLOCK_DESTRUCT(this)
  !! Destruct block matrix.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix

  ! local variables
  integer :: ib, jb

  if (allocated(this%i0)) deallocate (this%i0)
  if (allocated(this%i1)) deallocate (this%i1)
  if (allocated(this%j0)) deallocate (this%j0)
  if (allocated(this%j1)) deallocate (this%j1)

  if (allocated(this%b)) then
    do ib = 1, this%nbrows
      do jb = 1, this%nbcols
        if (.not. allocated(this%b(ib,jb)%p)) cycle

        call this%b(ib,jb)%p%destruct()
      end do
    end do
    deallocate (this%b)
  end if

  this%factorized = .false.
end subroutine

subroutine BLOCK_RESET(this, only_factorization)
  !! Reset block matrix to zero (optionally only factorization).

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix
  logical, optional, intent(in)    :: only_factorization
    !! reset only factorization data (default: false)

  ! local variables
  integer :: ib, jb

  do ib = 1, this%nbrows
    do jb = 1, this%nbcols
      if (.not. allocated(this%b(ib,jb)%p)) cycle
      call this%b(ib,jb)%p%reset(only_factorization = only_factorization)
    end do
  end do

  this%factorized = .false.
end subroutine

subroutine BLOCK_SCALE(this, fact)
  !! Scale block matrix by factor.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix
  TT,                intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: ib, jb

  ASSERT(.not. this%factorized)

  do ib = 1, this%nbrows
    do jb = 1, this%nbcols
      if (.not. allocated(this%b(ib,jb)%p)) cycle
      call this%b(ib,jb)%p%scale(fact)
    end do
  end do
end subroutine

subroutine BLOCK_MUL_VEC(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(BLOCK_TYPE),   intent(in)    :: this
    !! Block matrix
  TT,                  intent(in)    :: x(:)
    !! Raw vector to multiply by
  TT,                  intent(inout) :: y(:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  character :: trans_
  TT        :: one
  integer   :: ib, jb

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = 'N'
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ == 'N') then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. allocated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_vec(x(j0(ib):j1(ib)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. allocated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_vec(x(i0(jb):i1(jb)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

subroutine BLOCK_MUL_MAT(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(BLOCK_TYPE),   intent(in)    :: this
    !! Block matrix
  TT,                  intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  TT,                  intent(inout) :: y(:,:)
    !! Store result in y
  TT,        optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character, optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  ! local variables
  character :: trans_
  TT        :: one
  integer   :: ib, jb

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = 'N'
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ == 'N') then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. allocated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_mat(x(j0(ib):j1(ib),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. allocated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_mat(x(i0(jb):i1(jb),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. allocated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

subroutine BLOCK_FACTORIZE(this)
  !! Factorize diagonal Block Matrix.

  class(BLOCK_TYPE), intent(inout) :: this
    !! Block matrix

  ! local variables
  integer :: ib

  ASSERT(.not. this%factorized)

  if (this%diag) then
    do ib = 1, this%nbrows
      if (.not. allocated(this%b(ib,ib)%p)) call program_error("Can not factorize diagonal block matrix with empty block on main diag!")
      call this%b(ib,ib)%p%factorize()
    end do
  else
    call program_error("Factorization for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if

  this%factorized = .true.
end subroutine

subroutine BLOCK_SOLVE_VEC(this, rhs, x, trans, iter_refine)
  !! Solve block diagonal system: this * x = rhs. Call this%factorize() first.

  class(BLOCK_TYPE),   intent(in)  :: this
    !! Block matrix
  TT,                  intent(in)  :: rhs(:)
    !! Right hand side of equation system
  TT,                  intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement

  ! local variables
  integer :: ib

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_vec(rhs(this%i0(ib):this%i1(ib)), x(this%j0(ib):this%j1(ib)), trans = trans, iter_refine = iter_refine)
    end do
  else
    call program_error("Solving for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if
end subroutine

subroutine BLOCK_SOLVE_MAT(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(BLOCK_TYPE),   intent(in)  :: this
    !! Block matrix
  TT,                  intent(in)  :: rhs(:,:)
    !! right hand side of equation system
  TT,                  intent(out) :: x(:,:)
    !! store result in x
  character, optional, intent(in)  :: trans
    !! transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! perform iterative refinement

  ! local variables
  integer :: ib

  ASSERT(this%factorized)

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_mat(rhs(this%i0(ib):this%i1(ib),:), x(this%j0(ib):this%j1(ib),:), trans = trans, iter_refine = iter_refine)
    end do
  else
    call program_error("Solving for non-diagonal block matrix not implemented! Convert to dense or sparse instead!")
  end if
end subroutine

subroutine BLOCK_TO_DENSE(this, d, i0, j0)
  !! Insert block matrix into a dense matrix.

  class(BLOCK_TYPE), intent(in)    :: this
    !! Block matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, ib, jb

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements blockwise
  associate (bi0 => this%i0, bj0 => this%j0)
    do ib = 1, this%nbrows
      if (this%diag) then
        if (.not. allocated(this%b(ib,ib)%p)) cycle
        call this%b(ib,ib)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1)
      elseif (this%tridiag) then
        do jb = ib - 1, ib + 1
          if ((jb < 1) .or. (jb > this%nbcols)) cycle
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          call this%b(ib,jb)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1)
        end do
      else
        do jb = 1, this%nbcols
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          call this%b(ib,jb)%p%to_dense(d, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1)
        end do
      end if
    end do
  end associate
end subroutine

subroutine BLOCK_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
  !! Insert block matrix into a sparse matrix.

  class(BLOCK_TYPE),  intent(in)    :: this
    !! Block matrix (unchanged)
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
  integer :: i0_, j0_, ib, jb
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert elements blockwise
  associate (bi0 => this%i0, bi1 => this%i1, bj0 => this%j0, bj1 => this%j1)
    do ib = 1, this%nbrows
      if (this%diag) then
        if (.not. allocated(this%b(ib,ib)%p)) cycle
        if (present(struct)) then
          call this%b(ib,ib)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1, struct = struct(bi0(ib):bi1(ib), bj0(ib):bj1(ib)))
        else
          call this%b(ib,ib)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(ib)-1, drop_zeros = drop_zeros)
        end if
      elseif (this%tridiag) then
        do jb = ib - 1, ib + 1
          if ((jb < 1) .or. (jb > this%nbcols)) cycle
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          if (present(struct)) then
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        end do
      else
        do jb = 1, this%nbcols
          if (.not. allocated(this%b(ib,jb)%p)) cycle
          if (present(struct)) then
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call this%b(ib,jb)%p%to_sparse(sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        end do
      end if
    end do
  end associate
end subroutine

subroutine BLOCK_GET_GLOBAL_IDX(this, ib, jb, il, jl, ig, jg)
  !! Get global index from block and local index.

  class(BLOCK_TYPE), intent(in)  :: this
    !! Block matrix
  integer,           intent(in)  :: ib
    !! Block row
  integer,           intent(in)  :: jb
    !! Block col
  integer,           intent(in)  :: il
    !! Local row
  integer,           intent(in)  :: jl
    !! Local col
  integer,           intent(out) :: ig
    !! Output global row
  integer,           intent(out) :: jg
    !! Output global col

  ig = this%i0(ib) + il - 1
  jg = this%j0(jb) + jl - 1
end subroutine

subroutine BLOCK_GET_LOCAL_IDX(this, ig, jg, ib, jb, il, jl)
  !! Get block and local index from global index.

  class(BLOCK_TYPE), intent(in)  :: this
    !! Block matrix
  integer,           intent(in)  :: ig
    !! Global row
  integer,           intent(in)  :: jg
    !! Global col
  integer,           intent(out) :: ib
    !! Output block row
  integer,           intent(out) :: jb
    !! Output block col
  integer,           intent(out) :: il
    !! Output local row
  integer,           intent(out) :: jl
    !! Output local col

  ! find block row
  do ib = 1, this%nbrows
    if (this%i1(ib) >= ig) exit
  end do
  il = ig + 1 - this%i0(ib)

  ! find block col
  do jb = 1, this%ncols
    if (this%j1(jb) >= jg) exit
  end do
  jl = jg + 1 - this%j0(jb)
end subroutine

subroutine BLOCK_GET_DENSE(this, ib, jb, p)
  !! Return pointer to dense block.

  class(BLOCK_TYPE), target, intent(inout) :: this
    !! Block matrix
  integer,                   intent(in)    :: ib
    !! Block row
  integer,                   intent(in)    :: jb
    !! Block col
  type(DENSE_TYPE), pointer, intent(out)   :: p
    !! Output pointer to dense block

  ! allocate if block not already allocated
  if (.not. allocated(this%b(ib,jb)%p)) allocate (DENSE_TYPE :: this%b(ib,jb)%p)

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
    class is (DENSE_TYPE)
      p => b
    class default
      call program_error("Wrong matrix type given! (given: dense, block's type: something else)")
  end select
end subroutine

subroutine BLOCK_GET_SPARSE(this, ib, jb, p)
  !! Return pointer to sparse block.

  class(BLOCK_TYPE), target,  intent(inout) :: this
    !! Block matrix
  integer,                    intent(in)    :: ib
    !! Block row
  integer,                    intent(in)    :: jb
    !! Block col
  type(SPARSE_TYPE), pointer, intent(out)   :: p
    !! Output pointer to sparse block

  ! allocate if block not already allocated
  if (.not. allocated(this%b(ib,jb)%p)) allocate (SPARSE_TYPE :: this%b(ib,jb)%p)

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
    class is (SPARSE_TYPE)
      p => b
    class default
      call program_error("Wrong matrix type given! (given: SPARSE, block's type: something else)")
  end select
end subroutine

subroutine BLOCK_GET_BAND(this, ib, jb, p)
  !! Return pointer to BAND block.

  class(BLOCK_TYPE), target, intent(inout) :: this
    !! Block matrix
  integer,                   intent(in)    :: ib
    !! Block row
  integer,                   intent(in)    :: jb
    !! Block col
  type(BAND_TYPE), pointer,  intent(out)   :: p
    !! Output pointer to BAND block

  ! allocate if block not already allocated
  if (.not. allocated(this%b(ib,jb)%p)) allocate (BAND_TYPE :: this%b(ib,jb)%p)

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
    class is (BAND_TYPE)
      p => b
    class default
      call program_error("Wrong matrix type given! (given: BAND, block's type: something else)")
  end select
end subroutine

subroutine BLOCK_GET_BLOCK(this, ib, jb, p)
  !! Return pointer to BLOCK block.

  class(BLOCK_TYPE), target, intent(inout) :: this
    !! Block matrix
  integer,                   intent(in)    :: ib
    !! Block row
  integer,                   intent(in)    :: jb
    !! Block col
  type(BLOCK_TYPE), pointer, intent(out)   :: p
    !! Output pointer to BLOCK block

  ! allocate if block not already allocated
  if (.not. allocated(this%b(ib,jb)%p)) allocate (BLOCK_TYPE :: this%b(ib,jb)%p)

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
    class is (BLOCK_TYPE)
      p => b
    class default
      call program_error("Wrong matrix type given! (given: BLOCK, block's type: something else)")
  end select
end subroutine

subroutine BLOCK_SET(this, ib, jb, b, fact)
  !! Set one specific block.

  class(BLOCK_TYPE),  intent(inout) :: this
    !! Block matrix
  integer,            intent(in)    :: ib
    !! Block row
  integer,            intent(in)    :: jb
    !! Block col
  class(MATRIX_TYPE), intent(in)    :: b
  TT, optional,       intent(in)    :: fact
    !! Scaling factor (default: 1)

  ! check that ib,jb are in correct ranges (by block counts and by diag/tridiag flags)
  ASSERT((ib >= 1) .and. (ib <= this%nbrows))
  ASSERT((jb >= 1) .and. (jb <= this%nbcols))
  if (this%diag) then
    ASSERT(ib == jb)
  else if (this%tridiag) then
    ASSERT(abs(ib-jb) <= 1)
  end if

  ! destruct old block
  if (allocated(this%b(ib,jb)%p)) then
    call this%b(ib,jb)%p%destruct()
  end if

  ! set block
  this%b(ib,jb)%p = b

  ! scaling
  if (present(fact)) call this%b(ib,jb)%p%scale(fact)
end subroutine

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

#undef BLOCK_INIT
#undef BLOCK_DESTRUCT
#undef BLOCK_RESET
#undef BLOCK_SCALE
#undef BLOCK_MUL_VEC
#undef BLOCK_MUL_MAT
#undef BLOCK_FACTORIZE
#undef BLOCK_SOLVE_VEC
#undef BLOCK_SOLVE_MAT
#undef BLOCK_TO_DENSE
#undef BLOCK_TO_SPARSE
#undef BLOCK_GET_GLOBAL_IDX
#undef BLOCK_GET_LOCAL_IDX
#undef BLOCK_SET_BLOCK_TYPE
#undef BLOCK_GET_DENSE
#undef BLOCK_GET_SPARSE
#undef BLOCK_GET_BAND
#undef BLOCK_GET_BLOCK
#undef BLOCK_SET
