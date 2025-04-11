module subroutine block_{}T{}_init(this, row_dim, col_dim, diag, tridiag)
  !! Initialize Block matrix.

  class(block_{}T),  intent(out) :: this
    !! Block matrix
  integer,           intent(in)  :: row_dim(:)
    !! Number of Rows for each block
  integer, optional, intent(in)  :: col_dim(:)
    !! Number of Columns for each block (default: row_dim)
  logical, optional, intent(in)  :: diag
    !! diagonal block matrix flag (default: false)
  logical, optional, intent(in)  :: tridiag
    !! tridiagonal block matrix flag (default: false)

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
  allocate (this%b(  this%nbrows,this%nbcols))
  allocate (this%alc(this%nbrows,this%nbcols), source = .false.)

  ! diag, tridiag
  this%diag = .false.
  if (present(diag)) this%diag = diag
  this%tridiag = .false.
  if (present(tridiag)) this%tridiag = tridiag

  if (this%diag .and. this%tridiag) call program_error("Block matrix cannot be both diagonal and tridiagonal!")
  if (this%diag .or.  this%tridiag) then
    if (this%nbrows /= this%nbcols) call program_error("(Tri)Diagonal Block matrix must have same number of blocks per row and column!")
  end if
end subroutine

module subroutine block_{}T{}_destruct(this)
  !! Destruct block matrix.

  class(block_{}T), intent(inout) :: this
    !! Block matrix

  integer :: ib, jb

  if (allocated(this%i0)) deallocate (this%i0)
  if (allocated(this%i1)) deallocate (this%i1)
  if (allocated(this%j0)) deallocate (this%j0)
  if (allocated(this%j1)) deallocate (this%j1)

  if (allocated(this%b)) then
    do ib = 1, this%nbrows; do jb = 1, this%nbcols
      if (.not. this%alc(ib,jb)) cycle
      if (.not. associated(this%b(ib,jb)%p)) cycle

      call this%b(ib,jb)%p%destruct()
      deallocate (this%b(ib,jb)%p)
      this%alc(ib,jb) = .false.
    end do; end do
    deallocate (this%b)
  end if

  call this%sparse_conversion%destruct()
  call this%dense_conversion%destruct()

  this%factorized = .false.
end subroutine

module subroutine block_{}T{}_delete_block(this, ib, jb)
  !! Deletes block (ib, jb).
  !!
  !! Necessary in case one "get"s a block but doesnt insert any data (e.g. 0 matrix).
  !! Then delete that block. Otherwise an error is thrown.

  class(block_{}T), intent(inout) :: this
    !! Block matrix
  integer,          intent(in)    :: ib, jb
    !! row/col index.

  m4_assert(this%block_allowed(ib, jb))

  if (this%alc(ib,jb)) then
    if (associated(this%b(ib,jb)%p)) deallocate (this%b(ib,jb)%p)
    this%alc(ib,jb) = .false.
  else
    this%b(ib,jb)%p => null()
  end if

  this%factorized = .false.
end subroutine

module subroutine block_{}T{}_reset(this, only_factorization)
  !! Reset block matrix to zero (optionally only factorization).

  class(block_{}T),  intent(inout) :: this
    !! Block matrix
  logical, optional, intent(in)    :: only_factorization
    !! reset only factorization data (default: false)

  integer :: ib, jb

  do ib = 1, this%nbrows; do jb = 1, this%nbcols
    if (.not. associated(this%b(ib,jb)%p)) cycle
    call this%b(ib,jb)%p%reset(only_factorization = only_factorization)
  end do; end do

  call this%sparse_conversion%destruct()
  call this%dense_conversion%destruct()

  this%factorized = .false.
end subroutine

module subroutine block_{}T{}_scale(this, fact)
  !! Scale block matrix by factor.

  class(block_{}T), intent(inout) :: this
    !! Block matrix
  m4_type(T),       intent(in)    :: fact
    !! Scaling factor

  integer :: ib, jb

  m4_assert(.not. this%factorized)
  m4_assert(this%check_sizes())

  do ib = 1, this%nbrows; do jb = 1, this%nbcols
    if (.not. associated(this%b(ib,jb)%p)) cycle
    if (.not. this%alc(ib,jb)) call program_error("Trying to scale block matrix with non-owned blocks")
    call this%b(ib,jb)%p%scale(fact)
  end do; end do
end subroutine

module subroutine block_{}T{}_mul_vec(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(block_{}T),     intent(in)    :: this
    !! Block matrix
  m4_type(T),           intent(in)    :: x(:)
    !! Raw vector to multiply by
  m4_type(T),           intent(inout) :: y(:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
    !! Transposed flag, must be "N", "T" or "C" (default: "N")

  character  :: trans_
  m4_type(T) :: one
  integer    :: ib, jb

  m4_assert(this%check_sizes())

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = "N"
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ == "N") then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. associated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_vec(x(j0(ib):j1(ib)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(j0(jb):j1(jb)), y(i0(ib):i1(ib)), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. associated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_vec(x(i0(jb):i1(jb)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_vec(x(i0(ib):i1(ib)), y(j0(jb):j1(jb)), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

module subroutine block_{}T{}_mul_mat(this, x, y, fact_y, trans)
  !! Perform operation: y = this * x + fact_y * y

  class(block_{}T),     intent(in)    :: this
    !! Block matrix
  m4_type(T),           intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  m4_type(T),           intent(inout) :: y(:,:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
    !! Transposed flag, must be "N", "T" or "C" (default: "N")

  character  :: trans_
  m4_type(T) :: one
  integer    :: ib, jb

  m4_assert(this%check_sizes())

  ! optional arguments
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if
  trans_ = "N"
  if (present(trans)) trans_ = trans

  ! real or complex 1
  one = 1

  associate (i0 => this%i0, i1 => this%i1, j0 => this%j0, j1 => this%j1)
    if (trans_ == "N") then
      do ib = 1, this%nbrows
        if (this%diag) then
          if (.not. associated(this%b(ib,ib)%p)) cycle
          call this%b(ib,ib)%p%mul_mat(x(j0(ib):j1(ib),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do jb = ib - 1, ib + 1
            if ((jb < 1) .or. (jb > this%nbcols)) cycle
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        else
          do jb = 1, this%nbcols
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(j0(jb):j1(jb),:), y(i0(ib):i1(ib),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    else
      do jb = 1, this%nbcols
        if (this%diag) then
          if (.not. associated(this%b(jb,jb)%p)) cycle
          call this%b(jb,jb)%p%mul_mat(x(i0(jb):i1(jb),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
        elseif (this%tridiag) then
          do ib = jb - 1, jb + 1
            if ((ib < 1) .or. (ib > this%nbrows)) cycle
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        else
          do ib = 1, this%nbrows
            if (.not. associated(this%b(ib,jb)%p)) cycle
            call this%b(ib,jb)%p%mul_mat(x(i0(ib):i1(ib),:), y(j0(jb):j1(jb),:), fact_y = one, trans = trans)
          end do
        end if
      end do
    end if
  end associate
end subroutine

module subroutine block_{}T{}_factorize(this, solver)
  !! Factorize diagonal Block Matrix.

  class(block_{}T),  intent(inout) :: this
  integer, optional, intent(in)    :: solver

  integer :: ib, solver_

  call this%reset(only_factorization = .true.)

  solver_ = SPSOLVER_PARDISO
  if(present(solver)) solver_ = solver

  if (this%diag) then
    do ib = 1, this%nbrows
      if (.not. associated(this%b(ib,ib)%p)) call program_error("Can not factorize diagonal block matrix with empty block on main diag!")
      call this%b(ib,ib)%p%factorize()
    end do
  else
    select case(solver_)
    case(SPSOLVER_PARDISO)
      call matrix_convert(this, this%sparse_conversion)
      call this%sparse_conversion%factorize(SPSOLVER_PARDISO)
    m4_divert(m4_ifdef({m4_klu2},0,-1))
    case (SPSOLVER_KLU2)
      call matrix_convert(this, this%sparse_conversion)
      call this%sparse_conversion%factorize(SPSOLVER_KLU2)
    m4_divert(m4_ifdef({m4_mumps},0,-1))
    case (SPSOLVER_MUMPS)
      call matrix_convert(this, this%sparse_conversion)
      call this%sparse_conversion%factorize(SPSOLVER_MUMPS)
    m4_divert(m4_ifdef({m4_ilupack},0,-1))
    case (SPSOLVER_ILUPACK)
      call matrix_convert(this, this%sparse_conversion)
      call this%sparse_conversion%factorize(SPSOLVER_ILUPACK)
    m4_divert(0)
    case (BSOLVER_LAPACK)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(m4_ifdef({m4_spike},0,-1))
    case (BSOLVER_SPIKE)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(0)
    case (DSOLVER_LAPACK)
      call matrix_convert(this, this%dense_conversion)
      call this%dense_conversion%factorize(DSOLVER_LAPACK)
    case default
      call program_error("Invalid solver for block matrix selected")
    end select
  end if

  this%solver = solver_
  this%factorized = .true.
end subroutine

module subroutine block_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
  !! Solve block diagonal system: this * x = rhs. Call this%factorize() first.

  class(block_{}T),    intent(in)  :: this
    !! Block matrix
  m4_type(T),          intent(in)  :: rhs(:)
    !! Right hand side of equation system
  m4_type(T),          intent(out) :: x(:)
    !! Store result in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be "N", "C" or "T" (default: "N")
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement

  integer :: ib

  m4_assert(this%check_sizes())

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_vec(rhs(this%i0(ib):this%i1(ib)), x(this%j0(ib):this%j1(ib)), trans = trans, iter_refine = iter_refine)
    end do
  else
    select case(this%solver)
    case(SPSOLVER_PARDISO)
      call this%sparse_conversion%solve_vec(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_klu2},0,-1))
    case (SPSOLVER_KLU2)
      call this%sparse_conversion%solve_vec(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_mumps},0,-1))
    case (SPSOLVER_MUMPS)
      call this%sparse_conversion%solve_vec(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_ilupack},0,-1))
    case (SPSOLVER_ILUPACK)
      call this%sparse_conversion%solve_vec(rhs, x, trans, iter_refine)
    m4_divert(0)
    case (BSOLVER_LAPACK)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(m4_ifdef({m4_spike},0,-1))
    case (BSOLVER_SPIKE)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(0)
    case (DSOLVER_LAPACK)
      call this%dense_conversion%solve_vec(rhs, x, trans, iter_refine)
    end select
  end if
end subroutine

module subroutine block_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
  !! Solve system: this * x = rhs. Call this%factorize() first.

  class(block_{}T),    intent(in)  :: this
    !! Block matrix
  m4_type(T),          intent(in)  :: rhs(:,:)
    !! right hand side of equation system
  m4_type(T),          intent(out) :: x(:,:)
    !! store result in x
  character, optional, intent(in)  :: trans
    !! transposed flag, must be "N", "C" or "T" (default: "N")
  logical,   optional, intent(in)  :: iter_refine
    !! perform iterative refinement

  integer :: ib

  m4_assert(this%factorized)
  m4_assert(this%check_sizes())

  if (this%diag) then
    do ib = 1, this%nbrows
      call this%b(ib,ib)%p%solve_mat(rhs(this%i0(ib):this%i1(ib),:), x(this%j0(ib):this%j1(ib),:), trans = trans, iter_refine = iter_refine)
    end do
  else
    select case(this%solver)
    case(SPSOLVER_PARDISO)
      call this%sparse_conversion%solve_mat(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_klu2},0,-1))
    case (SPSOLVER_KLU2)
      call this%sparse_conversion%solve_mat(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_mumps},0,-1))
    case (SPSOLVER_MUMPS)
      call this%sparse_conversion%solve_mat(rhs, x, trans, iter_refine)
    m4_divert(m4_ifdef({m4_ilupack},0,-1))
    case (SPSOLVER_ILUPACK)
      call this%sparse_conversion%solve_mat(rhs, x, trans, iter_refine)
    m4_divert(0)
    case (BSOLVER_LAPACK)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(m4_ifdef({m4_spike},0,-1))
    case (BSOLVER_SPIKE)
      call program_error("Solving block matrices with band solver not supported yet because there is no conversion routine from block to band")
    m4_divert(0)
    case (DSOLVER_LAPACK)
      call this%dense_conversion%solve_mat(rhs, x, trans, iter_refine)
    end select
  end if
end subroutine

module subroutine block_{}T{}_solve_lsqr(this, rhs, damp, x, istop, atol, btol)
  !! Convert block matrix to sparse and solve Linear equation system this * x = rhs in the least squares sense.
  !! Works also for singular, over- and underdetermined systems. Only implemented for real matrices!
  !! Iterative solver algorithm, implementation taken from Github.
  class(block_{}T), intent(in)  :: this
    !! Block matrix
  m4_type(T),       intent(in)  :: rhs(:)
    !! Right hand side
  real,             intent(in)  :: damp
    !! Damping factor, (0.0 for simply solving this * x = rhs)
  m4_type(T),       intent(out) :: x(:)
    !! Solution vector
  integer,          intent(out) :: istop
    !! Exit code (see comments in lsqr.f90)
  real, optional,   intent(in)  :: atol
    !! Relative tolerance of the matrix elements
  real, optional,   intent(in)  :: btol
    !! Relative tolerance of the right hand side elements

  type(sparse_{}T) :: sparse

  ! check sizes
  m4_assert(size(x)   == this%ncols)
  m4_assert(size(rhs) == this%nrows)
  m4_assert(this%check_sizes())

  m4_ifelse(T, cmplx, {
    m4_ignore(rhs)
    m4_ignore(damp)
    m4_ignore(x)
    m4_ignore(istop)
    m4_ignore(atol)
    m4_ignore(btol)
    m4_ignore(sparse)
    call program_error("lsqr solver only supported for real matrices")}, {

    call matrix_convert(this, sparse)
    call sparse%solve_lsqr(rhs, damp, x, istop, atol, btol)
    call sparse%destruct()

    })
end subroutine

module subroutine block_{}T{}_get_global_idx(this, ib, jb, il, jl, ig, jg)
  !! Get global index from block and local index.

  class(block_{}T), intent(in)  :: this
    !! Block matrix
  integer,          intent(in)  :: ib
    !! Block row
  integer,          intent(in)  :: jb
    !! Block col
  integer,          intent(in)  :: il
    !! Local row
  integer,          intent(in)  :: jl
    !! Local col
  integer,          intent(out) :: ig
    !! Output global row
  integer,          intent(out) :: jg
    !! Output global col

  ig = this%i0(ib) + il - 1
  jg = this%j0(jb) + jl - 1
end subroutine

module subroutine block_{}T{}_get_local_idx(this, ig, jg, ib, jb, il, jl)
  !! Get block and local index from global index.

  class(block_{}T), intent(in)  :: this
    !! Block matrix
  integer,          intent(in)  :: ig
    !! Global row
  integer,          intent(in)  :: jg
    !! Global col
  integer,          intent(out) :: ib
    !! Output block row
  integer,          intent(out) :: jb
    !! Output block col
  integer,          intent(out) :: il
    !! Output local row
  integer,          intent(out) :: jl
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

module subroutine block_{}T{}_get_dense(this, ib, jb, p)
  !! Return pointer to dense block.

  class(block_{}T),         intent(inout) :: this
    !! Block matrix
  integer,                  intent(in)    :: ib
    !! Block row
  integer,                  intent(in)    :: jb
    !! Block col
  type(dense_{}T), pointer, intent(out)   :: p
    !! Output pointer to dense block

  m4_assert(this%block_allowed(ib, jb))

  ! allocate if block not already allocated
  if (.not. associated(this%b(ib,jb)%p)) then
    this%alc(ib,jb) = .true.
    allocate (dense_{}T :: this%b(ib,jb)%p)
  end if

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
  class is (dense_{}T)
    p => b
  class default
    call program_error("Wrong matrix type given! (given: dense; block type: something else)")
  end select
end subroutine

module subroutine block_{}T{}_get_sparse(this, ib, jb, p)
  !! Return pointer to sparse block.

  class(block_{}T),          intent(inout) :: this
    !! Block matrix
  integer,                   intent(in)    :: ib
    !! Block row
  integer,                   intent(in)    :: jb
    !! Block col
  type(sparse_{}T), pointer, intent(out)   :: p
    !! Output pointer to sparse block

  m4_assert(this%block_allowed(ib, jb))

  ! allocate if block not already allocated
  if (.not. associated(this%b(ib,jb)%p)) then
    this%alc(ib,jb) = .true.
    allocate (sparse_{}T :: this%b(ib,jb)%p)
  end if

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
  class is (sparse_{}T)
    p => b
  class default
    call program_error("Wrong matrix type given! (given: SPARSE; block type: something else)")
  end select
end subroutine

module subroutine block_{}T{}_get_band(this, ib, jb, p)
  !! Return pointer to BAND block.

  class(block_{}T),        intent(inout) :: this
    !! Block matrix
  integer,                 intent(in)    :: ib
    !! Block row
  integer,                 intent(in)    :: jb
    !! Block col
  type(band_{}T), pointer, intent(out)   :: p
    !! Output pointer to BAND block

  m4_assert(this%block_allowed(ib, jb))

  ! allocate if block not already allocated
  if (.not. associated(this%b(ib,jb)%p)) then
    this%alc(ib,jb) = .true.
    allocate (band_{}T :: this%b(ib,jb)%p)
  end if

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
  class is (band_{}T)
    p => b
  class default
    call program_error("Wrong matrix type given! (given: BAND; block type: something else)")
  end select
end subroutine

module subroutine block_{}T{}_get_block(this, ib, jb, p)
  !! Return pointer to BLOCK block.

  class(block_{}T),         intent(inout) :: this
    !! Block matrix
  integer,                  intent(in)    :: ib
    !! Block row
  integer,                  intent(in)    :: jb
    !! Block col
  type(block_{}T), pointer, intent(out)   :: p
    !! Output pointer to BLOCK block

  m4_assert(this%block_allowed(ib, jb))

  ! allocate if block not already allocated
  if (.not. associated(this%b(ib,jb)%p)) then
    this%alc(ib,jb) = .true.
    allocate (block_{}T :: this%b(ib,jb)%p)
  end if

  ! set intent out pointer to block
  select type (b => this%b(ib,jb)%p)
  class is (block_{}T)
    p => b
  class default
    call program_error("Wrong matrix type given! (given: BLOCK; block type: something else)")
  end select
end subroutine

module subroutine block_{}T{}_set(this, ib, jb, b, fact)
  !! Set one specific block (deep copy, shallow copy for non-owned subpointers)
  class(block_{}T),     intent(inout) :: this
    !! Block matrix
  integer,              intent(in)    :: ib
    !! Block row
  integer,              intent(in)    :: jb
    !! Block col
  class(matrix_{}T),    intent(in)    :: b
  m4_type(T), optional, intent(in)    :: fact
    !! Scaling factor (default: 1)

  m4_assert(this%block_allowed(ib, jb))

  ! destruct old block
  if (this%alc(ib,jb)) then
    if (associated(this%b(ib,jb)%p)) then
      call this%b(ib,jb)%p%destruct()
      deallocate (this%b(ib,jb)%p)
    end if
    this%alc(ib,jb) = .false.
  end if

  ! copy block
  this%alc(ib,jb) = .true.
  select type (b)
  class is (block_{}T)
    allocate (block_{}T :: this%b(ib,jb)%p)
    select type (this_b => this%b(ib,jb)%p)
    type is (block_{}T)
      call b%copy(this_b)
    end select
  class default
    allocate (this%b(ib,jb)%p, source = b)
  end select

  ! scaling
  if (present(fact)) call this%b(ib,jb)%p%scale(fact)
end subroutine

module subroutine block_{}T{}_set_ptr(this, ib, jb, b)
  !! Set one specific block by pointer
  class(block_{}T),          intent(inout) :: this
    !! Block matrix
  integer,                   intent(in)    :: ib
    !! Block row
  integer,                   intent(in)    :: jb
    !! Block col
  class(matrix_{}T), target, intent(in)    :: b

  m4_assert(this%block_allowed(ib, jb))

  ! destruct old block
  if (this%alc(ib,jb)) then
    if (associated(this%b(ib,jb)%p)) then
      call this%b(ib,jb)%p%destruct()
      deallocate (this%b(ib,jb)%p)
    end if
    this%alc(ib,jb) = .false.
  end if

  ! set pointer
  this%b(ib,jb)%p => b
end subroutine

module subroutine block_{}T{}_copy(this, copy)
  !! create deep copy of this block matrix (shallow copy for non-owned pointers)
  class(block_{}T), intent(in)  :: this
  type(block_{}T),  intent(out) :: copy
    !! output deep copy

  integer :: ib, jb

  ! initialize copy
  call copy%init(this%i1-this%i0+1, this%j1-this%j0+1, diag = this%diag, tridiag = this%tridiag)

  do ib = 1, this%nbrows; do jb = 1, this%nbcols
    if (this%alc(ib,jb)) then
      ! deep copy
      copy%alc(ib,jb) = .true.

      select type (p => this%b(ib,jb)%p)
      class is (block_{}T)
        allocate (block_{}T :: copy%b(ib,jb)%p)
        select type (copy_b => copy%b(ib,jb)%p)
        type is (block_{}T)
          call p%copy(copy_b)
        end select

      class default
        allocate (copy%b(ib,jb)%p, source = p)

      end select
    else
      ! shallow copy (pointer only)
      copy%b(ib,jb)%p => this%b(ib,jb)%p
    end if
  end do; end do
end subroutine

module subroutine block_{}T{}_copy_deep(this, copy)
  !! create deep copy of this block matrix (deep copy also for non-owned pointers)
  class(block_{}T), intent(in)  :: this
  type(block_{}T),  intent(out) :: copy
    !! output deep copy

  integer :: ib, jb

  ! initialize copy
  call copy%init(this%i1-this%i0+1, this%j1-this%j0+1, diag = this%diag, tridiag = this%tridiag)

  do ib = 1, this%nbrows; do jb = 1, this%nbcols
    if (.not. associated(this%b(ib,jb)%p)) cycle
    copy%alc(ib,jb) = .true.

    select type (p => this%b(ib,jb)%p)
    class is (block_{}T)
      allocate (block_{}T :: copy%b(ib,jb)%p)
      select type (copy_b => copy%b(ib,jb)%p)
      type is (block_{}T)
        call p%copy_deep(copy_b)
      end select
    class default
      allocate (copy%b(ib,jb)%p, source = p)
    end select
  end do; end do
end subroutine

module function block_{}T{}_block_allowed(this, ib, jb) result(tf)
  !! Checks if block given by row ib, col jb is allowed (bounds, diag, tridiag?)

  class(block_{}T), intent(in) :: this
    !! Block matrix
  integer,          intent(in) :: ib
    !! Block row
  integer,          intent(in) :: jb
    !! Block col
  logical                      :: tf
    !! Is block (ib,jb) allowed to be used?

  tf = ((ib >= 1) .and. (ib <= this%nbrows)) .and. ((jb >= 1) .and. (jb <= this%nbcols))

  if (this%diag) then
    tf = tf .and. (ib == jb)
  else if (this%tridiag) then
    tf = tf .and. (abs(ib-jb) <= 1)
  end if
end function

module function block_{}T{}_check_sizes(this) result(tf)
  !! Checks if the allocated matrix blocks have correct sizes.

  class(block_{}T), intent(in) :: this
  logical                      :: tf
    !! Have all allocated blocks the correct dimensions?

  integer :: ib, jb

  do ib = 1, this%nbrows; do jb = 1, this%nbcols
    if (.not. associated(this%b(ib,jb)%p)) cycle

    if (this%b(ib,jb)%p%nrows /= this%i1(ib) - this%i0(ib) + 1) then
      call program_error("block (ib,jb)=("//int2str(ib)//","//int2str(jb)//") has a wrong number of rows.")
    end if
    if (this%b(ib,jb)%p%ncols /= this%j1(jb) - this%j0(jb) + 1) then
      call program_error("block (ib,jb)=("//int2str(ib)//","//int2str(jb)//") has a wrong number of cols.")
    end if
  end do; end do

  tf = .true.
end function
