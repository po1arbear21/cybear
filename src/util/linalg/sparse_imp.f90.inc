module subroutine sparse_{}T{}_init(this, nrows, ncols)
  !! Initialize sparse matrix.

  class(sparse_{}T), intent(out) :: this
    !! Sparse matrix
  integer,           intent(in)  :: nrows
    !! Number of rows
  integer, optional, intent(in)  :: ncols
    !! Optional: Number of columns (default: ncols = nrows)

  ! init base
  call this%init("", nrows, ncols = ncols)
end subroutine

module subroutine sparse_{}T{}_destruct(this)
  !! Destruct sparse matrix by deallocating memory and deleting factorization.

  class(sparse_{}T), intent(inout) :: this
    !! Sparse matrix

  call this%reset()
end subroutine

module subroutine sparse_{}T{}_reset(this, only_factorization)
  !! Reset sparse matrix to empty state (optionally only the factorization).

  class(sparse_{}T), intent(inout) :: this
    !! Sparse matrix
  logical, optional, intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  ! deallocate CSR data
  if (.not. only_fact) then
    if (allocated(this%ia)) deallocate (this%ia)
    if (allocated(this%ja)) deallocate (this%ja)
    if (allocated(this%a )) deallocate (this%a )
  end if

  ! destruct factorization
  if (this%solver_handle > 0) then
    select case (this%solver)
    case (SPSOLVER_PARDISO)
      call destruct_pardiso_handle(this%solver_handle)
    m4_divert(m4_ifdef({m4_klu2},0,-1))
    case (SPSOLVER_KLU2)
      call destruct_klu2_handle_{}m4_ifelse(T,cmplx,{c},{r}){}(this%solver_handle)
    m4_divert(m4_ifdef({m4_mumps},0,-1))
    case (SPSOLVER_MUMPS)
      call destruct_mumps_handle_{}m4_ifelse(T,cmplx,{c},{r}){}(this%solver_handle)
    m4_divert(m4_ifdef({m4_ilupack},0,-1))
    case (SPSOLVER_ILUPACK)
      call destruct_ilupack_handle(this%solver_handle)
    m4_divert(0)
    end select
  end if

  this%factorized = .false.
end subroutine

module subroutine sparse_{}T{}_scale(this, fact)
  !! Scale sparse matrix by a factor.

  class(sparse_{}T), intent(inout) :: this
    !! Sparse matrix
  m4_type(T),        intent(in)    :: fact
    !! Scaling factor

  ! h == 0 is default value for unfactorized
  m4_assert(this%solver_handle == 0)

  ! scale values
  if (.not. this%is_empty()) this%a = fact * this%a
end subroutine

module subroutine sparse_{}T{}_mul_vec_slice(this, x, y, fact_y, i0, i1)
  !! Multiply sparse matrix by a raw vector where only some rows of y are computed.
  !! Performs operation: y = this * x + fact_y * y

  class(sparse_{}T),    intent(in)    :: this
    !! Sparse matrix
  m4_type(T),           intent(in)    :: x(:)
    !! Raw vector to multiply by
  m4_type(T),           intent(inout) :: y(:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  integer,    optional, intent(in)    :: i0
    !! first row index. default: 1
  integer,    optional, intent(in)    :: i1
    !! last row index. default: this%nrows

  integer             :: i0_, i1_, i
  integer(SPARSE_IDX) :: j

  i0_ = 1
  if (present(i0)) i0_ = i0
  i1_ = this%nrows
  if (present(i1)) i1_ = i1

  m4_assert(size(y) == i1_-i0_+1)
  m4_assert((i0_ >= 1) .and. (i0_ <= this%nrows))
  m4_assert((i1_ >= 1) .and. (i1_ <= this%nrows))
  m4_assert(i0_ <= i1_)

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! matrix vector multiplication
  !$omp parallel do default(none) private(i, j) shared(this, x, y, i0_, i1_) schedule(dynamic)
  do i = 1, i1_-i0_+1
    do j = this%ia(i+i0_-1), this%ia(i+i0_)-1
      y(i) = y(i) + this%a(j) * x(this%ja(j))
    end do
  end do
  !$omp end parallel do
end subroutine

module subroutine sparse_{}T{}_mul_vec(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(sparse_{}T),    intent(in)    :: this
    !! Sparse matrix
  m4_type(T),           intent(in)    :: x(:)
    !! Raw vector to multiply by
  m4_type(T),           intent(inout) :: y(:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
    !! Transpose flag, must be 'N', '{T}' or 'C' (default: 'N')
    !! 'N':   y =                     this  * x + fact_y * y
    !! '{T}':   y =           transpose(this) * x + fact_y * y
    !! 'C':   y = conjugate_transpose(this) * x + fact_y * y

  character           :: trans_
  integer             :: i
  integer(SPARSE_IDX) :: j

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! check sizes
  m4_assert(size(x) == this%ncols)
  m4_assert(size(y) == this%nrows)

  ! transpose flag
  trans_ = 'N'
  if (present(trans)) trans_ = trans
  m4_assert((trans_ == 'N') .or. (trans_ == '{{T}}')m4_ifelse(T,cmplx,{ .or. (trans_ == 'C')},))

  if (trans_ == 'N') then
    ! matrix vector multiplication
    !$omp parallel do default(none) private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        y(i) = y(i) + this%a(j) * x(this%ja(j))
      end do
    end do
    !$omp end parallel do

  else
    ! fixme y_tmp is unnecessary! use directly y whereever there is y_tmp now. (in reduction and assignments). yields catastrophic compiler error for ifort ...
    block
      m4_type(T), allocatable :: y_tmp(:)

      allocate (y_tmp(size(y)))
      y_tmp = 0

      m4_divert(m4_ifelse(T,cmplx,0,-1))
        if (trans_ == 'C') then
          ! matrix vector multiplication with complex conjugated transposed matrix
          !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
          do i = 1, this%nrows
            do j = this%ia(i), this%ia(i+1)-1
              y_tmp(this%ja(j)) = y_tmp(this%ja(j)) + conjg(this%a(j)) * x(i)
            end do
          end do
          !$omp end parallel do
        end if
      m4_divert(0)
      if (trans_ == '{T}') then
        ! matrix vector multiplication with transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j)) = y_tmp(this%ja(j)) + this%a(j) * x(i)
          end do
        end do
        !$omp end parallel do
      end if

      y = y + y_tmp
    end block
  end if
end subroutine

module subroutine sparse_{}T{}_mul_mat(this, x, y, fact_y, trans)
  !! Multiply sparse matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(sparse_{}T),    intent(in)    :: this
    !! Sparse matrix
  m4_type(T),           intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  m4_type(T),           intent(inout) :: y(:,:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
    !! Transpose flag for this matrix, must be 'N', '{T}' or 'C' (default: 'N')
    !! 'N':   y =                     this  * x + fact_y * y
    !! '{T}':   y =           transpose(this) * x + fact_y * y
    !! 'C':   y = conjugate_transpose(this) * x + fact_y * y

  character           :: trans_
  integer             :: i
  integer(SPARSE_IDX) :: j

  ! y scaling factor
  if (present(fact_y)) then
    y = fact_y * y
  else
    y = 0
  end if

  ! ignore x if matrix is empty
  if (this%is_empty()) return

  ! check sizes
  m4_assert(size(x, dim=1) == this%ncols)
  m4_assert(size(x, dim=2) == size(y, dim=2))
  m4_assert(size(y, dim=1) == this%nrows)

  ! transpose flag
  trans_ = 'N'
  if (present(trans)) trans_ = trans
  m4_assert((trans_ == 'N') .or. (trans_ == '{{T}}')m4_ifelse(T,cmplx,{ .or. (trans_ == 'C')},))

  if (trans_ == 'N') then
    ! matrix matrix multiplication
    !$omp parallel do default(none) private(i, j) shared(this, x, y) schedule(dynamic)
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        y(i,:) = y(i,:) + this%a(j) * x(this%ja(j),:)
      end do
    end do
    !$omp end parallel do

  else
    ! fixme y_tmp is unnecessary! use directly y whereever there is y_tmp now. (in reduction and assignments). yields catastrophic compiler error for ifort ...
    block
      m4_type(T), allocatable :: y_tmp(:,:)

      allocate (y_tmp(size(y, dim=1), size(y, dim=2)))
      y_tmp = 0

      m4_divert(m4_ifelse(T,cmplx,0,-1))
        if (trans_ == 'C') then
          ! matrix matrix multiplication with complex conjugated transposed matrix
          !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
          do i = 1, this%nrows
            do j = this%ia(i), this%ia(i+1)-1
              y_tmp(this%ja(j),:) = y_tmp(this%ja(j),:) + conjg(this%a(j)) * x(i,:)
            end do
          end do
          !$omp end parallel do
        end if
      m4_divert(0)
      if (trans_ == '{T}') then
        ! matrix matrix multiplication with transposed matrix
        !$omp parallel do default(none) private(i, j) shared(this, x) reduction(+: y_tmp) schedule(dynamic)
        do i = 1, this%nrows
          do j = this%ia(i), this%ia(i+1)-1
            y_tmp(this%ja(j),:) = y_tmp(this%ja(j),:) + this%a(j) * x(i,:)
          end do
        end do
        !$omp end parallel do
      end if
      y = y + y_tmp
    end block
  end if
end subroutine

module subroutine sparse_{}T{}_init_solver(this, solver)
  !! creates solver handle
  !!
  !! this routine exists s.t. one can init ilupack, then get its handler and change parameters.

  class(sparse_{}T), intent(inout) :: this
  integer,           intent(in)    :: solver
    !! solver used for factorization

  ! make sure matrix is unfactorized and not empty
  m4_assert(this%solver_handle == 0)
  m4_assert(.not. this%is_empty())

  select case (solver)
  case (SPSOLVER_PARDISO)
    this%solver_handle = create_pardiso_handle(this%nrows, m4_ifelse(T,cmplx,{.true.},{.false.}))
  m4_divert(m4_ifdef({m4_klu2},0,-1))
  case (SPSOLVER_KLU2)
    this%solver_handle = create_klu2_handle_{}m4_ifelse(T,cmplx,{c},{r}){}()
  m4_divert(m4_ifdef({m4_mumps},0,-1))
  case (SPSOLVER_MUMPS)
    this%solver_handle = create_mumps_handle_{}m4_ifelse(T,cmplx,{c},{r}){}()
  m4_divert(m4_ifdef({m4_ilupack},0,-1))
  case (SPSOLVER_ILUPACK)
    this%solver_handle = create_ilupack_handle(this%ia, this%ja, this%a)
  m4_divert(0)
  end select
end subroutine

module subroutine sparse_{}T{}_factorize(this, solver)
  !! Factorize Matrix using PARDISO, KLU2, MUMPS or ILUPACK.

  class(sparse_{}T), intent(inout) :: this
  integer, optional, intent(in)    :: solver
    !! solver used for factorization (default: SPSOLVER_PARDISO)

  integer :: solver_

  call this%reset(only_factorization = .true.)

  solver_ = SPSOLVER_PARDISO
  if(present(solver)) solver_ = solver

  call this%init_solver(solver_)

  select case (solver_)
  case (SPSOLVER_PARDISO)
    call pardiso_factorize(this%solver_handle, this%ia, this%ja, this%a)
  m4_divert(m4_ifdef({m4_klu2},0,-1))
  case (SPSOLVER_KLU2)
    call klu2_factorize(this%solver_handle, this%ia, this%ja, this%a)
  m4_divert(m4_ifdef({m4_mumps},0,-1))
  case (SPSOLVER_MUMPS)
    call mumps_factorize(this%solver_handle, this%ia, this%ja, this%a)
  m4_divert(m4_ifdef({m4_ilupack},0,-1))
  case (SPSOLVER_ILUPACK)
    call ilupack_factorize(this%solver_handle)
  m4_divert(0)
  case default
    call program_error("Trying to solve sparse matrix with non-sparse solver!")
  end select

  this%solver = solver_
  this%factorized = .true.
end subroutine

module subroutine sparse_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using the chosen solver.
  !!
  !! Call this%factorize() first!

  class(sparse_{}T),   intent(in)  :: this
    !! Sparse matrix
  m4_type(T),          intent(in)  :: rhs(:)
    !! Right hand side
  m4_type(T),          intent(out) :: x(:)
    !! Solution vector
  character, optional, intent(in)  :: trans
    !! Not implemented for '{T}' or 'C' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  ! matrix must be factorized before solving
  m4_assert(this%factorized)

  ! check sizes
  m4_assert(size(x)   == this%ncols)
  m4_assert(size(rhs) == this%nrows)

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed sparse matrix not implemented!")
  end if
  if (present(iter_refine)) then
    if (.not. iter_refine) call program_error("Solving with sparse matrix without iterative refinement not implemented!")
  end if

  ! solve system with pardiso, klu2, mumps or ilupack
  select case (this%solver)
  case (SPSOLVER_PARDISO)
    call pardiso_solve(this%solver_handle, this%ia, this%ja, this%a, rhs, x)
  m4_divert(m4_ifdef({m4_klu2},0,-1))
  case (SPSOLVER_KLU2)
    call klu2_solve(this%solver_handle, rhs, x)
  m4_divert(m4_ifdef({m4_mumps},0,-1))
  case (SPSOLVER_MUMPS)
    call mumps_solve(this%solver_handle, rhs, x)
  m4_divert(m4_ifdef({m4_ilupack},0,-1))
  case (SPSOLVER_ILUPACK)
    call ilupack_solve(this%solver_handle, rhs, x)
  m4_divert(0)
  end select
end subroutine

module subroutine sparse_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
  !! Solve Linear equation system this * x = rhs using the chosen solver.
  !!
  !! Call this%factorize() first.

  class(sparse_{}T),   intent(in)  :: this
    !! Sparse matrix
  m4_type(T),          intent(in)  :: rhs(:,:)
    !! Right hand side
  m4_type(T),          intent(out) :: x(:,:)
    !! Solution vector
  character, optional, intent(in)  :: trans
    !! Not implemented for '{T}' or 'C' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Not implemented for false (default: true)

  m4_type(T), allocatable :: rhs_(:), x_(:)

  ! matrix must be factorized before solving
  m4_assert(this%factorized)

  ! check sizes
  m4_assert(size(x,   dim=1) == this%ncols)
  m4_assert(size(x,   dim=2) == size(rhs, dim=2))
  m4_assert(size(rhs, dim=1) == this%nrows)

  ! optional arguments
  if (present(trans)) then
    if (trans /= 'N') call program_error("Solving with transposed sparse matrix not implemented!")
  end if
  if (present(iter_refine)) then
    if (.not. iter_refine) call program_error("Solving with sparse matrix without iterative refinement not implemented!")
  end if

  ! allocate space for temporary copy
  allocate (rhs_(size(rhs,1) * size(rhs,2)))
  allocate (  x_(size(  x,1) * size(  x,2)))

  ! solve with pardiso, klu2, mumps or ilupack, pass right-hand sides one after another in 1D array
  rhs_ = reshape(rhs, [size(rhs,1) * size(rhs,2)])
  select case (this%solver)
  case (SPSOLVER_PARDISO)
    call pardiso_solve(this%solver_handle, this%ia, this%ja, this%a, rhs_, x_)
  m4_divert(m4_ifdef({m4_klu2},0,-1))
  case (SPSOLVER_KLU2)
    call klu2_solve(this%solver_handle, rhs_, x_)
  m4_divert(m4_ifdef({m4_mumps},0,-1))
  case (SPSOLVER_MUMPS)
    call mumps_solve(this%solver_handle, rhs_, x_)
  m4_divert(m4_ifdef({m4_ilupack},0,-1))
  case (SPSOLVER_ILUPACK)
    call ilupack_solve(this%solver_handle, rhs_, x_)
  m4_divert(0)
  end select
  x = reshape(x_, [size(x,1), size(x,2)])
end subroutine

module subroutine sparse_{}T{}_solve_lsqr(this, rhs, damp, x, istop, atol, btol)
  !! Solve Linear equation system this * x = rhs in the least squares sense.
  !! Works also for singular, over- and underdetermined systems. Only implemented for real matrices!
  !! Iterative solver algorithm, implementation taken from Github.
  class(sparse_{}T), intent(in)  :: this
    !! Sparse matrix
  m4_type(T),        intent(in)  :: rhs(:)
    !! Right hand side
  real,              intent(in)  :: damp
    !! Damping factor, (0.0 for simply solving this * x = rhs)
  m4_type(T),        intent(out) :: x(:)
    !! Solution vector
  integer,           intent(out) :: istop
    !! Exit code (see comments in lsqr.f90)
  real, optional,    intent(in)  :: atol
    !! Relative tolerance of the matrix elements
  real, optional,    intent(in)  :: btol
    !! Relative tolerance of the right hand side elements

  integer, allocatable :: irow(:)
  type(lsqr_solver_ez) :: lsqr

  ! check sizes
  m4_assert(size(x)   == this%ncols)
  m4_assert(size(rhs) == this%nrows)

  m4_ifelse(T, cmplx, {
    m4_ignore(rhs)
    m4_ignore(damp)
    m4_ignore(x)
    m4_ignore(istop)
    m4_ignore(atol)
    m4_ignore(btol)
    m4_ignore(irow)
    m4_ignore(lsqr)
    call program_error("lsqr solver only supported for real matrices")}, {

    allocate (irow(this%nnz(only_nonzeros = .true.)))
    do i = 1, this%nrows
      irow(this%ia(i):this%ia(i+1)-1) = i
    end do
    call lsqr%initialize(this%nrows, this%ncols, this%a, irow, this%ja, atol = atol, btol = btol)
    call lsqr%solve(rhs, damp, x, istop)

    })
end subroutine

module subroutine sparse_{}T{}_transpose(this)
  !! Transpose sparse matrix.

  class(sparse_{}T), intent(inout) :: this
    !! Sparse matrix

  type(sparse_{}T) :: s2

  ! h == 0 is default value for unfactorized
  m4_assert(this%solver_handle == 0)

  ! transpose matrix
  call this%transpose(s2)

  ! copy transposed matrix
  this%nrows = s2%nrows
  this%ncols = s2%ncols
  if (.not. s2%is_empty()) then
    this%ia = s2%ia
    this%ja = s2%ja
    this%a  = s2%a
  end if
end subroutine

module subroutine sparse_{}T{}_transpose2(this, s2)
  !! Transpose sparse matrix, store result in second matrix.

  class(sparse_{}T), intent(in)  :: this
    !! First sparse matrix
  type(sparse_{}T),  intent(out) :: s2
    !! Output transposed sparse matrix

  integer                       :: i
  integer(SPARSE_IDX)           :: j, k
  type(vector_int), allocatable :: col(:)
  type(vector_{}T), allocatable :: val(:)

  ! init s2
  call s2%init(this%ncols, ncols = this%nrows)

  ! do not copy values of an empty matrix
  if (this%is_empty()) return

  ! allocate memory
  allocate (col(this%ncols))
  allocate (val(this%ncols))
  do i = 1, this%ncols
    call col(i)%init(n = 0, c = 8)
    call val(i)%init(n = 0, c = 8)
  end do

  ! transpose matrix
  do i = 1, this%nrows
    do j = this%ia(i), this%ia(i+1)-1
      call col(this%ja(j))%push(i)
      call val(this%ja(j))%push(this%a(j))
    end do
  end do

  ! allocate memory for csr
  allocate (s2%ia(s2%nrows + 1))
  allocate (s2%ja(size(this%ja,kind=SPARSE_IDX)))
  allocate (s2%a( size(this%a ,kind=SPARSE_IDX)))

  ! construct csr
  k        = 1
  s2%ia(1) = k
  do i = 1, size(col)
    do j = 1, col(i)%n
      s2%ja(k) = col(i)%d(j)
      s2%a( k) = val(i)%d(j)
      k        = k + 1
    end do
    s2%ia(i+1) = k
  end do
end subroutine

module subroutine sparse_{}T{}_input(this, file, funit)
  !! Read sparse matrix from file or stream (does not work for large matrices with more than 2^31 non-zero entries)

  class(sparse_{}T),     intent(out) :: this
    !! Sparse matrix
  character(*), optional, intent(in)  :: file
    !! File name (alternative: provide funit)
  integer,      optional, intent(in)  :: funit
    !! File unit (alternative: provide file name)

  integer                   :: i, row, col, nrows, ncols, status, funit_
  character(:), allocatable :: fmt_
  real                      :: rval{}m4_ifelse(T,cmplx,{, ival},)
  type(vector_int)          :: rows, cols
  type(vector_{}T)          :: vals
  type(spbuild_{}T)         :: sb

  allocate (character(0) :: fmt_) ! remove gfortran warning

  ! open file if needed
  m4_assert(.not. (present(file) .and. present(funit)))
  m4_assert(present(file) .or. present(funit))
  if (present(file)) then
    open (newunit = funit_, file = file, status = "old", action = "read", iostat = status)
    if (status /= 0) call program_error("File could not be opened for reading")
  else
    funit_ = funit
  end if

  ! read and save into vectors
  call rows%init(0, c = 1024)
  call cols%init(0, c = 1024)
  call vals%init(0, c = 1024)
  do
    read (funit_, *, iostat = status) row, col, rval{}m4_ifelse(T,cmplx,{, ival},)

    if (status < 0) then
      ! end of file reached
      close (funit_)
      exit
    elseif (status > 0) then
      ! IO error
      close(funit_)
      call program_error("IO Error")
    end if

    call rows%push(row)
    call cols%push(col)
    m4_ifelse(T,cmplx,{call vals%push(cmplx(rval, ival))},{call vals%push(rval)})
  end do

  ! close file if openend here
  if (present(file)) close(funit_)

  ! get dimensions
  nrows = maxval(rows%d(1:rows%n))
  ncols = maxval(cols%d(1:cols%n))

  ! create matrix
  call this%init(nrows, ncols = ncols)
  call sb%init(this)
  do i = 1, rows%n
    call sb%add(rows%d(i), cols%d(i), vals%d(i))
  end do
  call sb%save()
end subroutine

module subroutine sparse_{}T{}_output(this, file, funit, fmt)
  !! Write sparse matrix to file or stream

  class(sparse_{}T),      intent(in) :: this
    !! Sparse matrix
  character(*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,      optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(*), optional, intent(in) :: fmt
    !! Value format string

  integer                   :: i, status, funit_
  integer(SPARSE_IDX)       :: j
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

  ! output csr
  if (.not. this%is_empty()) then
    do i = 1, this%nrows
      do j = this%ia(i), this%ia(i+1)-1
        ! write row and column
        write (funit_, "(I0,A,I0,A)", advance = "no") i, " ", this%ja(j), " "

        ! write value
        write (funit_, fmt_) this%a(j)
      end do
    end do
  end if

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

module function sparse_{}T{}_nnz(this, only_nonzeros) result(nnz)
  !! computes number of non-zero elements

  class(sparse_{}T), intent(in) :: this
  logical, optional, intent(in) :: only_nonzeros
    !! assume that all entries of sparse matrix are non-zero. thus, is much faster.
    !! default: false
  integer(SPARSE_IDX)           :: nnz

  logical :: only_nonzeros_

  if (this%is_empty()) then
    nnz = 0
  else
    only_nonzeros_ = .false.
    if (present(only_nonzeros)) only_nonzeros_ = only_nonzeros

    if (only_nonzeros_) then
      nnz = this%ia(this%nrows+1)-1
    else
      nnz = count(abs(this%a) > 0.0, kind=SPARSE_IDX)
    end if
  end if
end function

module function sparse_{}T{}_is_empty(this) result(empty)
  !! check if sparse matrix is empty (CSR not allocated)

  class(sparse_{}T), intent(in) :: this
  logical                       :: empty
    !! return true if CSR is not allocated, false otherwise

  if ((this%nrows == 0) .or. &
    & (this%ncols == 0) .or. &
    & (.not. allocated(this%ia)) .or. &
    & (.not. allocated(this%ja)) .or. &
    & (.not. allocated(this%a))) then
    empty = .true.
  else
    empty = .false.
  end if
end function

module function sparse_eye_{}T{}(nrows, ncols) result(eye)
  !! sparse identity matrix

  integer,           intent(in)  :: nrows
    !! Number of rows
  integer, optional, intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(sparse_{}T)               :: eye
    !! return sparse identity matrix

  integer :: i, n

  call eye%init(nrows, ncols = ncols)
  n = min(eye%nrows, eye%ncols)

  ! allocate csr memory
  allocate (eye%ia(eye%nrows + 1))
  allocate (eye%ja(n))
  allocate (eye%a( n))

  eye%ia(1) = 1
  do i = 1, n
    eye%ja(i  ) = i
    eye%a( i  ) = 1
    eye%ia(i+1) = int(i+1, kind = SPARSE_IDX)
  end do
  do i = n + 1, eye%nrows
    eye%ia(i+1) = eye%ia(n+1)
  end do
end function

module function sparse_zero_{}T{}(nrows, ncols) result(zero)
  !! create empty sparse matrix

  integer,           intent(in) :: nrows
    !! Number of rows
  integer, optional, intent(in) :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(sparse_{}T)              :: zero
    !! return empty sparse matrix

  call zero%init(nrows, ncols = ncols)

  ! allocate csr memory
  allocate (zero%ia(nrows+1), source = int(1, kind = SPARSE_IDX))
  allocate (zero%ja(0))
  allocate (zero%a( 0))
end function

module subroutine spbuild_{}T{}_init(this, s, n_per_row, assume_unique)
  !! Initialize sparse matrix builder.

  class(spbuild_{}T),       intent(out) :: this
    !! Sparse matrix builder
  type(sparse_{}T), target, intent(in)  :: s
    !! Sparse matrix this builder is tied to
  integer, optional,        intent(in)  :: n_per_row
    !! Optional: provide hint for number of elements per row (default: 8)
  logical, optional,        intent(in)  :: assume_unique
    !! Optional: user ensures no duplicate entries (default: .false.)

  integer :: i, c

  ! number of elements per row
  c = 8
  if (present(n_per_row)) c = n_per_row

  ! assume unique
  this%assume_unique = .false.
  if (present(assume_unique)) this%assume_unique = assume_unique

  ! save pointer to matrix
  this%sp => s

  ! allocate data
  allocate (this%cols(  s%nrows))
  allocate (this%vals(  s%nrows))
  allocate (this%filled(s%nrows))
  do i = 1, s%nrows
    call this%cols(i)%init(  n = 0, c = c)
    call this%vals(i)%init(  n = 0, c = c)
    call this%filled(i)%init(n = 0, c = c)
  end do

  ! allocate row flags
  allocate (this%sorted(s%nrows),      source = .true.)
  allocate (this%keep_struct(s%nrows), source = .true.)

  ! load existing data from sparse matrix
  call this%load()
end subroutine

module subroutine spbuild_{}T{}_destruct(this)
  !! Deallocate memory.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder

  nullify (this%sp)
  if (allocated(this%cols       )) deallocate (this%cols       )
  if (allocated(this%vals       )) deallocate (this%vals       )
  if (allocated(this%filled     )) deallocate (this%filled     )
  if (allocated(this%sorted     )) deallocate (this%sorted     )
  if (allocated(this%keep_struct)) deallocate (this%keep_struct)
end subroutine

module subroutine spbuild_{}T{}_reset(this)
  !! Delete all values of each row

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder

  integer :: i

  do i = 1, this%sp%nrows
    call this%reset_row(i)
  end do
end subroutine

module subroutine spbuild_{}T{}_reset_row(this, row)
  !! Delete all values of a row

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder
  integer,            intent(in)    :: row
    !! row index

  if (this%keep_struct(row)) then
    ! reset data while keeping structure
    this%vals(  row)%d(1:this%vals(  row)%n) = 0
    this%filled(row)%d(1:this%filled(row)%n) = .false.

  else
    ! reset data without keeping structure
    call this%cols(  row)%resize(0)
    call this%vals(  row)%resize(0)
    call this%filled(row)%resize(0)
    this%sorted(row)      = .true.
    this%keep_struct(row) = .true.
  end if
end subroutine

module subroutine spbuild_{}T{}_set(this, row, col, val, fact, search)
  !! Insert/Update value.

  class(spbuild_{}T),   intent(inout) :: this
    !! Sparse matrix builder
  integer,              intent(in)    :: row
    !! Row index
  integer,              intent(in)    :: col
    !! Column index
  m4_type(T),           intent(in)    :: val
    !! New value
  m4_type(T), optional, intent(in)    :: fact
    !! Scale existing element before adding new value (default: 0)
  logical,    optional, intent(in)    :: search
    !! Search for existing value instead of appending to end of row (default: true)

  m4_type(T) :: fact_
  logical    :: search_
  integer    :: i, l, r, m

  m4_assert(row > 0)
  m4_assert(row <= this%sp%nrows)
  m4_assert(col > 0)
  m4_assert(col <= this%sp%ncols)

  ! optional args
  search_ = .true.
  if (present(search)) search_ = search
  fact_ = 0
  if (present(fact)) fact_ = fact

  ! get index i for new value
  if (search_) then
    ! decide if binary or linear search
    if (this%sorted(row)) then
      ! binary search in sorted row
      i = -1
      l = 1
      r = this%cols(row)%n
      do while (l <= r)
        m = (l + r) / 2
        if (this%cols(row)%d(m) < col) then
          l = m + 1
        elseif (this%cols(row)%d(m) > col) then
          r = m - 1
        else
          i = m
          exit
        end if
      end do

      ! not found, append to end
      if (i == -1) i = this%cols(row)%n + 1
    else
      ! linear search
      do i = 1, this%cols(row)%n
        if (this%cols(row)%d(i) == col) exit
      end do
    end if
  else
    ! append element without searching
    i = this%cols(row)%n + 1

    ! do not keep structure when resetting
    this%keep_struct(row) = .false.
  end if

  ! Insert/Update value
  if (i > this%cols(row)%n) then ! append to end
    ! update sorted flag
    if (this%sorted(row) .and. (this%cols(row)%n > 0)) then
      if (col <= this%cols(row)%back()) this%sorted(row) = .false.
    end if

    ! make sure that appending without searching happens in order
    ! user can disable this check if he ensures no duplicates
    if (.not. this%assume_unique) then
      m4_assert(search_ .or. this%sorted(row))
    end if

    ! append new value
    call this%cols(row)%push(col)
    call this%vals(row)%push(val)
    call this%filled(row)%push(.true.)
  else
    ! update existing value
    if (this%filled(row)%d(i)) then
      this%vals(row)%d(i) = fact_ * this%vals(row)%d(i) + val
    else
      this%vals(row)%d(i) = val
      this%filled(row)%d(i) = .true.
    end if
  end if
end subroutine

module subroutine spbuild_{}T{}_set_row(this, row, vals, j0, j1)
  !! Insert continuous block of values.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder
  integer,            intent(in)    :: row
    !! Row index
  m4_type(T),         intent(in)    :: vals(:)
    !! Values to insert (j1 - j0 + 1 elements)
  integer, optional,  intent(in)    :: j0
    !! Start column index (default: 1)
  integer, optional,  intent(in)    :: j1
    !! End column index (default: this%sp%ncols)

  integer :: j, j0_, j1_

  ! optional arguments
  j0_ = 1
  if (present(j0)) j0_ = j0
  j1_ = this%sp%ncols
  if (present(j1)) j1_ = j1
  m4_assert(j0_ > 0)
  m4_assert(j0_ <= this%sp%ncols)
  m4_assert(j1_ > 0)
  m4_assert(j1_ <= this%sp%ncols)

  ! setting whole row -> reset previous data of that row.
  ! thus, making sure row is sorted.
  if ((1 == j0_) .and. (this%sp%ncols == j1_)) then
    call this%cols(  row)%reset
    call this%vals(  row)%reset
    call this%filled(row)%reset
    this%sorted(row)      = .true.
    this%keep_struct(row) = .true.
  end if

  ! insert values
  do j = j0_, j1_
    call this%set(row, j, vals(j-j0_+1), search = .false.)
  end do
end subroutine

module subroutine spbuild_{}T{}_add(this, row, col, val)
  !! Insert new element or add to existing value.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder
  integer,            intent(in)    :: row
    !! Row index
  integer,            intent(in)    :: col
    !! Column index
  m4_type(T),         intent(in)    :: val
    !! New value

  m4_type(T) :: fact

  ! set scaling factor to 1 both for real and complex values
  fact = 1

  ! add value
  call this%set(row, col, val, fact = fact)
end subroutine

module subroutine spbuild_{}T{}_sort(this)
  !! Sort column indices.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder

  integer              :: i, n
  integer, allocatable :: perm(:)

  ! sort rows in parallel
  !$omp parallel default(none) private(i, n, perm) shared(this)

  ! thread local memory: allocate enough memory for permutation array
  allocate (perm(maxval(this%cols%n)))

  !$omp do schedule(dynamic)
  do i = 1, this%sp%nrows
    ! do nothing if row is already sorted
    if (this%sorted(i)) cycle
    this%sorted(i) = .true.

    ! number of elements in row
    n = this%cols(i)%n

    ! sort column indices and save permutation
    call qsort(this%cols(i)%d(1:n), perm = perm(1:n))

    ! apply permutation to values
    this%vals(  i)%d(1:n) = this%vals(  i)%d(perm(1:n))
    this%filled(i)%d(1:n) = this%filled(i)%d(perm(1:n))
  end do
  !$omp end do

  ! deallocate thread local memory (otherwise: memory leak!!!)
  deallocate (perm)

  !$omp end parallel
end subroutine

module subroutine spbuild_{}T{}_load(this)
  !! Load existing values from sparse matrix into builder.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder

  integer             :: i
  integer(SPARSE_IDX) :: j

  ! load all rows
  do i = 1, this%sp%nrows
    ! delete values
    call this%cols(  i)%resize(0)
    call this%vals(  i)%resize(0)
    call this%filled(i)%resize(0)

    ! reset row flags
    this%sorted(i)      = .true.
    this%keep_struct(i) = .true.

    ! add values from csr (if not empty)
    if (allocated(this%sp%ia)) then
      do j = this%sp%ia(i), this%sp%ia(i+1)-1
        call this%cols(  i)%push(this%sp%ja(j))
        call this%vals(  i)%push(this%sp%a( j))
        call this%filled(i)%push(.true.       )
      end do
    end if
  end do
end subroutine

module subroutine spbuild_{}T{}_save(this)
  !! Save values in sparse matrix.

  class(spbuild_{}T), intent(inout) :: this
    !! Sparse matrix builder

  integer             :: i, j
  integer(SPARSE_IDX) :: k

  ! reset sparse matrix. needed in case it has old values
  call this%sp%reset()

  ! sort column indices
  call this%sort()

  ! get total number of elements
  k = 0
  do i = 1, this%sp%nrows
    do j = 1, this%cols(i)%n
      if (.not. this%filled(i)%d(j)) cycle
      k = k + 1
    end do
  end do

  ! allocate memory
  allocate (this%sp%ia(this%sp%nrows + 1))
  allocate (this%sp%ja(k))
  allocate (this%sp%a( k))

  ! init csr
  this%sp%ia(1) = 1
  k = 0
  do i = 1, this%sp%nrows
    do j = 1, this%cols(i)%n
      if (.not. this%filled(i)%d(j)) cycle
      k             = k + 1
      this%sp%ja(k) = this%cols(i)%d(j)
      this%sp%a( k) = this%vals(i)%d(j)
    end do

    ! next row pointer
    this%sp%ia(i+1) = k + 1
  end do
end subroutine
