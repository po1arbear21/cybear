module subroutine dense_{}T{}_init_by_data(this, d)
  !! Initialize dense matrix.

  class(dense_{}T), intent(out) :: this
    !! Dense matrix
  m4_type(T),       intent(in)  :: d(:,:)
    !! Provide initial data.

  ! init base
  call this%init("", size(d, 1), ncols=size(d, 2))

  this%d = d
end subroutine

module subroutine dense_{}T{}_init_by_size(this, nrows, ncols)
  !! Initialize dense matrix.

  class(dense_{}T),  intent(out) :: this
    !! Dense matrix
  integer,           intent(in)  :: nrows
    !! Number of rows (Default: d0rows)
  integer, optional, intent(in)  :: ncols
    !! Optional: Number of columns. (Default: ncols = nrows)

  integer :: ncols_

  ncols_ = nrows
  if (present(ncols)) ncols_ = ncols

  ! init base
  call this%init("", nrows, ncols=ncols_)

  ! allocate data
  allocate (this%d(nrows,ncols_))
  this%d = 0
end subroutine

module subroutine dense_{}T{}_destruct(this)
  !! Destruct dense matrix; i.e. deallocate its members.

  class(dense_{}T), intent(inout) :: this
    !! Dense matrix

  if (allocated(this%d   )) deallocate (this%d   )
  if (allocated(this%lu  )) deallocate (this%lu  )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine dense_{}T{}_reset(this, only_factorization)
  !! Reset dense matrix to zero (optionally only the factorization).

  class(dense_{}T),  intent(inout) :: this
    !! Dense matrix
  logical, optional, intent(in)    :: only_factorization
    !! Reset only factorization data (default: false)

  logical :: only_fact

  only_fact = .false.
  if (present(only_factorization)) only_fact = only_factorization

  ! reset data
  if (.not. only_fact) then
    this%d = 0
  end if

  ! reset factorization
  if (allocated(this%lu  )) deallocate (this%lu  )
  if (allocated(this%ipiv)) deallocate (this%ipiv)
  this%factorized = .false.
end subroutine

module subroutine dense_{}T{}_scale(this, fact)
  !! Scale dense matrix by a factor.

  class(dense_{}T), intent(inout) :: this
    !! Dense matrix
  m4_type(T),       intent(in)    :: fact
    !! Scaling factor

  m4_assert(.not. allocated(this%lu))
  m4_assert(.not. allocated(this%ipiv))

  ! scale data
  this%d = fact * this%d
end subroutine

module subroutine dense_{}T{}_mul_vec(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw vector.
  !! Performs operation: y = this * x + fact_y * y

  class(dense_{}T),     intent(in)    :: this
    !! Dense matrix
  m4_type(T),           intent(in)    :: x(:)
    !! Raw vector to multiply by
  m4_type(T),           intent(inout) :: y(:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! If present, scale existing data in y and add result of multiplication
  character,  optional, intent(in)    :: trans
    !! Transposed flag, must be 'N', 'T' or 'C' (default: 'N')

  call gemv(this%d, x, y, beta = fact_y, trans = trans)
end subroutine

module subroutine dense_{}T{}_mul_mat(this, x, y, fact_y, trans)
  !! Multiply dense matrix by a raw matrix.
  !! Performs operation: y = this * x + fact_y * y

  class(dense_{}T),     intent(in)    :: this
    !! Dense matrix
  m4_type(T),           intent(in)    :: x(:,:)
    !! Raw matrix to multiply by
  m4_type(T),           intent(inout) :: y(:,:)
    !! Store result in y
  m4_type(T), optional, intent(in)    :: fact_y
    !! Scale existing data in y and add result to it (default: 0)
  character,  optional, intent(in)    :: trans
    !! Transposed flag for this matrix, must be 'N', 'C' or 'T' (default: 'N')

  call gemm(this%d, x, y, transa = trans, beta = fact_y)
end subroutine

module subroutine dense_{}T{}_factorize(this, solver)
  !! Calculate LU decomposition.
  !! Call this subroutine before solving!

  class(dense_{}T),  intent(inout) :: this
  integer, optional, intent(in)    :: solver

  integer :: info, solver_

  call this%reset(only_factorization = .true.)

  solver_ = DSOLVER_LAPACK
  if(present(solver)) solver_ = solver

  select case(solver_)
  case(DSOLVER_LAPACK)
    ! allocate memory for factorization
    if (.not. allocated(this%lu  )) allocate (this%lu(  this%nrows,this%ncols))
    if (.not. allocated(this%ipiv)) allocate (this%ipiv(this%nrows))

    ! factorization with lapack
    this%lu = this%d
    call getrf(this%lu, ipiv = this%ipiv, info = info)
    if (info /= 0) call program_error("Dense Factorization error in LAPACK", code = info)
  case default
    call program_error("Trying to solve dense matrix with non-dense solver!")
  end select

  this%solver = solver_
  this%factorized = .true.
end subroutine

module subroutine dense_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for a single right-hand side.
  !! Solves: this * x = rhs, where x and rhs are raw vectors. Call this%factorize() first!

  class(dense_{}T),    intent(in)  :: this
    !! Dense matrix
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

module subroutine dense_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
  !! Solve linear equation system for multiple right-hand sides.
  !! Solves: this * x = rhs, where x and rhs are raw matrices. Call this%factorize() first!

  class(dense_{}T),    intent(in)  :: this
    !! Dense matrix
  m4_type(T),          intent(in)  :: rhs(:,:)
    !! Right hand sides of equation system
  m4_type(T),          intent(out) :: x(:,:)
    !! Store results in x
  character, optional, intent(in)  :: trans
    !! Transposed flag, must be 'N', 'C' or 'T' (default: 'N')
  logical,   optional, intent(in)  :: iter_refine
    !! Perform iterative refinement (default: true)

  integer :: info
  logical :: iter_refine_

  m4_assert(this%factorized)

  ! optional arguments
  iter_refine_ = .true.
  if (present(iter_refine)) iter_refine_ = iter_refine

  ! solve with lapack
  x = rhs
  call getrs(this%lu, this%ipiv, x, trans = trans, info = info)
  if (info /= 0) call program_error("Dense solving error in LAPACK", code = info)

  ! iterative refinement
  if (iter_refine_) then
    call gerfs(this%d, this%lu, this%ipiv, rhs, x, trans = trans, info = info)
    if (info /= 0) call program_error("Dense iterative refinement error in LAPACK", code = info)
  end if
end subroutine

module subroutine dense_{}T{}_transpose(this)
  !! Transpose dense matrix

  class(dense_{}T), intent(inout) :: this
    !! Dense matrix

  m4_type(T), allocatable :: d(:,:)

  m4_assert(.not. allocated(this%lu))
  m4_assert(.not. allocated(this%ipiv))

  ! transpose matrix
  allocate (d(this%ncols,this%nrows), source = transpose(this%d))
  call move_alloc(d, this%d)
  this%nrows = size(this%d,1)
  this%ncols = size(this%d,2)
end subroutine

module subroutine dense_{}T{}_transpose2(this, d2)
  !! Transpose matrix and store result in second matrix

  class(dense_{}T), intent(in)  :: this
    !! First dense matrix (unchanged)
  type(dense_{}T),  intent(out) :: d2
    !! Output transposed matrix

  call d2%init(transpose(this%d))
end subroutine

module subroutine dense_{}T{}_eig(this, e, R, L, sort)
  !! Perform a dense eigenvalue decomposition.

  class(dense_{}T),  intent(in)  :: this
    !! Dense matrix
  complex,           intent(out) :: e(:)
    !! Output eigenvalues (complex in general)
  complex, optional, intent(out) :: R(:,:)
    !! Output right eigenvectors in columns (complex in general)
  complex, optional, intent(out) :: L(:,:)
    !! Output left  eigenvectors in columns (complex in general)
  logical, optional, intent(in)  :: sort
    !! Sort by magnitude of eigenvalues (default: false)

  integer :: info
  logical :: sort_
  m4_divert(m4_ifelse(T,cmplx,0,-1))
    complex, allocatable :: a(:,:)

    ! allocate memory
    allocate (a(this%nrows,this%nrows), source = this%d)

    ! get eigenvalue decomposition
    call geev(a, e, vl = L, vr = R, info = info)
    if (info /= 0) call program_error("error in complex LAPACK geev", code = info)
  m4_divert(m4_ifelse(T,real,0,-1))
    integer           :: i
    real, allocatable :: ar(:,:), wr(:), wi(:), vl(:,:), vr(:,:)

    ! allocate memory
    allocate (ar(this%nrows,this%nrows), source = this%d)
    allocate (wr(this%nrows))
    allocate (wi(this%nrows))

    if (present(L) .and. present(R)) then
      allocate (vl(this%nrows,this%nrows))
      allocate (vr(this%nrows,this%nrows))
      call geev(ar, wr, wi, vl = vl, vr = vr, info = info)
    elseif (present(L)) then
      allocate (vl(this%nrows,this%nrows))
      call geev(ar, wr, wi, vl = vl, info = info)
    elseif (present(R)) then
      allocate (vr(this%nrows,this%nrows))
      call geev(ar, wr, wi, vr = vr, info = info)
    else
      call geev(ar, wr, wi, info = info)
    end if
    if (info /= 0) call program_error("error in real LAPACK geev", code = info)

    ! construct eigenvalues and eigenvectors
    i = 1
    do while (i <= this%nrows)
      if (wi(i) == 0) then
        e(i) = wr(i)
        if (present(L)) L(:,i) = vl(:,i)
        if (present(R)) R(:,i) = vr(:,i)
        i = i + 1
      else
        e(i  ) = cmplx(wr(i  ), wi(i  ))
        e(i+1) = cmplx(wr(i+1), wi(i+1))
        if (present(L)) then
          L(:,i  ) = cmplx(vl(:,i), + vl(:,i+1))
          L(:,i+1) = cmplx(vl(:,i), - vl(:,i+1))
        end if
        if (present(R)) then
          R(:,i  ) = cmplx(vr(:,i), + vr(:,i+1))
          R(:,i+1) = cmplx(vr(:,i), - vr(:,i+1))
        end if
        i = i + 2
      end if
    end do
  m4_divert(0)

  ! sort
  sort_ = .false.
  if (present(sort)) sort_ = sort
  if (sort_) then
    block
      real,    allocatable :: abs_e(:)
      integer, allocatable :: perm(:)

      allocate (abs_e(this%nrows), source = abs(e))
      allocate (perm(this%nrows))

      ! sort eigenvalues
      call qsort(abs_e, perm = perm)

      e = e(perm)
      if (present(L)) L = L(:,perm)
      if (present(R)) R = R(:,perm)
    end block
  end if
end subroutine

module subroutine dense_{}T{}_output(this, file, funit, fmt)
  !! Write dense matrix to file or stream.

  class(dense_{}T),       intent(in) :: this
    !! Dense matrix
  character(*), optional, intent(in) :: file
    !! Filename (alternative: provide funit)
  integer,      optional, intent(in) :: funit
    !! File unit (alternative: provide file name)
  character(*), optional, intent(in) :: fmt
    !! Value format string

  integer                   :: i, j, funit_, status
  character(:), allocatable :: fmt_

  allocate (character(0) :: fmt_) ! remove gfortran warning
  if (present(fmt)) then
    fmt_ = fmt
  else
    fmt_ = m4_ifelse(T,cmplx,{'(A,ES25.16E3,SP,ES25.16E3,"i")'},{'(A,ES25.16E3)'})
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

  do i = 1, this%nrows
    do j = 1, this%ncols
      if (j > 1) then
        write(funit_, fmt_, advance = "no") ", ", this%d(i,j)
      else
        write(funit_, fmt_, advance = "no") "", this%d(i,j)
      end if
    end do

    ! next line
    write(funit_,*)
  end do

  ! close file if openend here
  if (present(file)) close(funit_)
end subroutine

module function dense_eye_{}T{}(nrows, ncols) result(eye)
  !! dense identity matrix

  integer,           intent(in) :: nrows
    !! Number of rows
  integer, optional, intent(in) :: ncols
    !! Optional number of columns (default: ncols = nrows)
  type(dense_{}T)               :: eye

  integer :: i

  call eye%init(nrows, ncols)
  do i = 1, min(eye%nrows, eye%ncols)
    eye%d(i,i) = 1
  end do
end function
