module subroutine add_band2_{}T{}(b1, b2, fact)
  !! Add two band matrices.
  !! Operation: b2 <- fact * b1 + b2

  type(band_{}T),       intent(in)    :: b1
    !! First band matrix (unchanged)
  type(band_{}T),       intent(inout) :: b2
    !! Second band matrix (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! First band matrix' scaling factor (default: 1)

  m4_type(T)              :: fact_
  m4_type(T), allocatable :: tmp(:,:)

  m4_assert(.not. allocated(b1%f))
  m4_assert(.not. allocated(b2%ipiv))
  m4_assert(b2%nrows == b1%nrows .and. b2%ncols == b1%ncols)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  if ((b2%nupper >= b1%nupper) .and. (b2%nlower >= b1%nlower)) then
    b2%d(-b1%nupper:b1%nlower,1:b2%ncols) = b2%d(-b1%nupper:b1%nlower,1:b2%ncols) + fact_ * b1%d
  else
    allocate (tmp(-max(b2%nupper,b1%nupper):max(b2%nlower,b1%nlower),1:b2%ncols))
    tmp = 0
    tmp(-b2%nupper:b2%nlower,1:b2%ncols) = b2%d
    tmp(-b1%nupper:b1%nlower,1:b2%ncols) = b2%d(-b1%nupper:b1%nlower,1:b2%ncols) + fact_*b1%d
    call move_alloc(tmp, b2%d)
    b2%nlower = max(b2%nlower, b1%nlower)
    b2%nupper = max(b2%nupper, b1%nupper)
  end if
end subroutine

module subroutine add_band3_{}T{}(b1, b2, b3, fact1, fact2)
  !! Add two band matrices, store result in third band matrix.
  !! Operation: b3 <- fact1*b1 + fact2*b2

  type(band_{}T),       intent(in)  :: b1
    !! Band first matrix (unchanged)
  type(band_{}T),       intent(in)  :: b2
    !! Band second matrix (unchanged)
  type(band_{}T),       intent(out) :: b3
    !! Output band third matrix
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_type(T) :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! case selection: hopefully no tmp array needs to be allocated
  if ((b1%nupper >= b2%nupper) .and. (b1%nlower >= b2%nlower)) then
    call b3%init(b1%nrows, b1%nlower, nupper=b1%nupper, d0=fact1_*b1%d)     ! b3 <- fact1*b1
    call matrix_add(b2, b3, fact=fact2_)                                           ! b3 += fact2*b2
  else
    call b3%init(b2%nrows, b2%nlower, nupper=b2%nupper, d0=fact2_*b2%d)     ! b3 <- fact2*b2
    call matrix_add(b1, b3, fact=fact1_)                                           ! b3 += fact1*b1
  end if
end subroutine

module subroutine add_band_dense_{}T{}(b, d, fact)
  !! Add dense and band matrix, store result in dense matrix.
  !!
  !! operation:  d <- fact * b + d

  type(band_{}T),       intent(in)    :: b
    !! Band matrix b (unchanged)
  type(dense_{}T),      intent(inout) :: d
    !! Dense matrix d (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  integer    :: i, j, k
  m4_type(T) :: fact_

  m4_assert(d%nrows == b%nrows)
  m4_assert(d%ncols == b%ncols)
  m4_assert(.not. allocated(d%lu  ))
  m4_assert(.not. allocated(d%ipiv))

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! add contents of band matrix
  do j = 1, b%ncols
    do i = -b%nupper, b%nlower
      ! get actual column
      k = j + i
      if ((k < 1) .or. (k > b%nrows)) cycle

      d%d(k,j) = d%d(k,j) + fact_ * b%d(i,j)
    end do
  end do
end subroutine

module subroutine add_band_dense2_{}T{}(b1, d2, d3, fact1, fact2)
  !! Add dense and band matrix, store result in a third dense matrix.
  !!
  !! operation:  d3 <- fact1 * b1 + fact2 * d2

  type(band_{}T),       intent(in)  :: b1
    !! Band first matrix b1 (unchanged)
  type(dense_{}T),      intent(in)  :: d2
    !! Dense second matrix d1 (unchanged)
  type(dense_{}T),      intent(out) :: d3
    !! Output dense third matrix d3
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_type(T) :: fact2_

  m4_assert(d2%nrows == b1%nrows)
  m4_assert(d2%ncols == b1%ncols)

  ! optional arguments
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! d3 <- fact2 * d2
  call d3%init(fact2_ * d2%d)

  ! d3 += fact1 * b1
  call matrix_add(b1, d3, fact = fact1)
end subroutine

module subroutine add_band_sparse_{}T{}(b, s, fact)
  !! Add sparse and band matrix, store result in sparse matrix.
  !!
  !! operation: s <- fact * b + s

  type(band_{}T),       intent(in)    :: b
    !! Band matrix (unchanged)
  type(sparse_{}T),     intent(inout) :: s
    !! Sparse matrix (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  integer                          :: i, j2
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, k, b_nrows, b_nlower, b_nupper
  integer(SPARSE_IDX), allocatable :: ia(:)
  m4_type(T)                       :: fact_
  m4_type(T),          allocatable :: a(:)

  m4_assert(s%nrows == b%nrows)
  m4_assert(s%ncols == b%ncols)

  ! h == 0 is default value for unfactorized
  m4_assert(s%solver_handle == 0)

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! if s is empty, set to zero matrix first
  if (s%is_empty()) then
    allocate (s%ia(s%nrows+1), source = int(1, kind = SPARSE_IDX))
    allocate (s%ja(0))
    allocate (s%a( 0))
  end if

  ! allocate memory
  allocate (ia(s%nrows + 1))
  b_nrows  = int(b%nrows ,kind=SPARSE_IDX)
  b_nlower = int(b%nlower,kind=SPARSE_IDX)
  b_nupper = int(b%nupper,kind=SPARSE_IDX)
  k = b_nrows + b_nlower * b_nrows - b_nlower * (b_nlower + 1) / 2 + b_nupper * b_nrows - b_nupper * (b_nupper + 1) / 2
  k = k + size(s%ja,kind=SPARSE_IDX)
  allocate (ja(k))
  allocate ( a(k))

  ! add matrices
  k     = 1
  ia(1) = k
  do i = 1, s%nrows
    j1 = s%ia(i)
    j2 = max(i - b%nlower, 1)
    do while ((j1 < s%ia(i+1)) .and. (j2 <= min(i + b%nupper, b%ncols)))
      if (s%ja(j1) < j2) then
        ja(k) = s%ja(j1)
        a( k) = s%a( j1)
        k     = k  + 1
        j1    = j1 + 1
      elseif (s%ja(j1) > j2) then
        ja(k) = j2
        a( k) = fact_ * b%d(i-j2,j2)
        k     = k  + 1
        j2    = j2 + 1
      else
        ja(k) = j2
        a( k) = s%a(j1) + fact_ * b%d(i-j2,j2)
        k     = k  + 1
        j1    = j1 + 1
        j2    = j2 + 1
      end if
    end do
    do while (j1 < s%ia(i+1))
      ja(k) = s%ja(j1)
      a( k) = s%a( j1)
      k     = k  + 1
      j1    = j1 + 1
    end do
    do while (j2 <= min(i + b%nupper, b%ncols))
      ja(k) = j2
      a( k) = fact_ * b%d(i-j2,j2)
      k     = k  + 1
      j2    = j2 + 1
    end do
    ia(i+1) = k
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, s%ia)
  deallocate (s%ja)
  deallocate (s%a )
  allocate (s%ja(k-1), source = ja(1:k-1))
  allocate (s%a( k-1), source =  a(1:k-1))
end subroutine

module subroutine add_band_sparse2_{}T{}(b1, s2, s3, fact1, fact2)
  !! Add sparse and band matrix, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * b1 + fact2 * s2

  type(band_{}T),       intent(in)  :: b1
    !! Band first matrix (unchanged)
  type(sparse_{}T),     intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(sparse_{}T),     intent(out) :: s3
    !! Output sparse third matrix
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_assert(s2%nrows == b1%nrows)
  m4_assert(s2%ncols == b1%ncols)

  ! s3 <- fact2 * s2
  call s3%init(s2%nrows, ncols = s2%ncols)
  if (.not. s2%is_empty()) then
    s3%ia = s2%ia
    s3%ja = s2%ja
    if (present(fact2)) then
      s3%a = fact2 * s2%a
    else
      s3%a = s2%a
    end if
  end if

  ! s3 += fact1 * b1
  call matrix_add(b1, s3, fact = fact1)
end subroutine

module subroutine add_dense2_{}T{}(d1, d2, fact)
  !! Add two dense matrices, store the result in second matrix.
  !!
  !! operation:  d2 <- fact * d1 + d2

  type(dense_{}T),      intent(in)    :: d1
    !! First dense matrix d1 (unchanged)
  type(dense_{}T),      intent(inout) :: d2
    !! Second dense matrix d2 (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! Second dense matrix scaling factor (default: 1)

  m4_assert(.not. allocated(d2%lu))
  m4_assert(.not. allocated(d2%ipiv))

  if (present(fact)) then
    d2%d = fact * d1%d + d2%d
  else
    d2%d =        d1%d + d2%d
  end if
end subroutine

module subroutine add_dense3_{}T{}(d1, d2, d3, fact1, fact2)
  !! Add two dense matrices, store the result in third dense matrix.
  !!
  !! operation:  d3 <- fact1 * d1 + fact2 * d2

  type(dense_{}T),      intent(in)  :: d1
    !! Dense first matrix d1 (unchanged)
  type(dense_{}T),      intent(in)  :: d2
    !! Dense second matrix d2 (unchanged)
  type(dense_{}T),      intent(out) :: d3
    !! Output dense third matrix d3
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_type(T) :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! init output matrix
  call d3%init(fact1_ * d1%d + fact2_ * d2%d)
end subroutine

module subroutine add_sparse_dense_{}T{}(s, d, fact)
  !! Add dense and sparse matrix, store result in dense matrix.
  !!
  !! operation:  D <- fact * S + D

  type(sparse_{}T),     intent(in)    :: s
    !! Sparse matrix S (unchanged)
  type(dense_{}T),      intent(inout) :: d
    !! Dense matrix D (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! Sparse matrix scaling factor (default: 1)

  integer             :: i
  integer(SPARSE_IDX) :: j
  m4_type(T)          :: fact_

  m4_assert(d%nrows == s%nrows)
  m4_assert(d%ncols == s%ncols)
  m4_assert(.not. allocated(d%lu))
  m4_assert(.not. allocated(d%ipiv))

  ! do nothing if s is empty
  if (s%is_empty()) return

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! add contents of sparse matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      d%d(i,s%ja(j)) = d%d(i,s%ja(j)) + fact_ * s%a(j)
    end do
  end do
end subroutine

module subroutine add_sparse_dense2_{}T{}(s1, d2, d3, fact1, fact2)
  !! Add dense and sparse matrix, store result in another dense matrix.
  !!
  !! operation:  d3 <- fact1 * s1 + fact2 * d2

  type(sparse_{}T),     intent(in)  :: s1
    !! Sparse first matrix s1 (unchanged)
  type(dense_{}T),      intent(in)  :: d2
    !! Dense second matrix d1 (unchanged)
  type(dense_{}T),      intent(out) :: d3
    !! Output dense third matrix
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_type(T) :: fact2_

  m4_assert(d2%nrows == s1%nrows)
  m4_assert(d2%ncols == s1%ncols)

  ! optional argument
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! d3 <- fact2 * d2
  call d3%init(fact2_ * d2%d)

  ! d3 += fact1 * s1
  if (.not. s1%is_empty()) call matrix_add(s1, d3, fact = fact1)
end subroutine

module subroutine add_sparse2_{}T{}(s1, s2, fact)
  !! Add two sparse matrices, store result in first matrix.
  !!
  !! operation: s2 <- fact * s1 + s2

  type(sparse_{}T),     intent(in)    :: s1
    !! First sparse matrix (unchanged)
  type(sparse_{}T),     intent(inout) :: s2
    !! Second sparse matrix (updated)
  m4_type(T), optional, intent(in)    :: fact
    !! Second sparse matrix scaling factor (default: 1)

  integer                          :: i
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, j2, k
  integer(SPARSE_IDX), allocatable :: ia(:)
  m4_type(T)                       :: fact_
  m4_type(T),          allocatable :: a(:)

  m4_assert(s2%nrows == s1%nrows)
  m4_assert(s2%ncols == s1%ncols)

  ! h == 0 is default value for unfactorized
  m4_assert(s2%solver_handle == 0)

  ! do nothing if s1 is empty
  if (s1%is_empty()) return

  ! optional arguments
  fact_ = 1
  if (present(fact)) fact_ = fact

  ! simply copy and scale s1 if s2 matrix is empty
  if (s2%is_empty()) then
    s2%ia = s1%ia
    s2%ja = s1%ja
    s2%a  = s1%a * fact_
    return
  end if

  ! allocate memory
  allocate (ia(s2%nrows + 1))
  allocate (ja(size(s2%ja,kind=SPARSE_IDX) + size(s1%ja,kind=SPARSE_IDX)))
  allocate ( a(size(s2%ja,kind=SPARSE_IDX) + size(s1%ja,kind=SPARSE_IDX)))

  ! add matrices
  k     = 1
  ia(1) = k
  do i = 1, s2%nrows
    j1 = s2%ia(i)
    j2 = s1%ia(i)
    do while ((j1 < s2%ia(i+1)) .and. (j2 < s1%ia(i+1)))
      if (s2%ja(j1) < s1%ja(j2)) then
        ja(k) = s2%ja(j1)
        a( k) = s2%a( j1)
        k     = k  + 1
        j1    = j1 + 1
      else if (s2%ja(j1) > s1%ja(j2)) then
        ja(k) = s1%ja(j2)
        a( k) = fact_ * s1%a(j2)
        k     = k  + 1
        j2    = j2 + 1
      else
        ja(k) = s2%ja(j1)
        a( k) = s2%a( j1) + fact_ * s1%a(j2)
        k     = k  + 1
        j1    = j1 + 1
        j2    = j2 + 1
      end if
    end do
    do while (j1 < s2%ia(i+1))
      ja(k) = s2%ja(j1)
      a( k) = s2%a( j1)
      k     = k  + 1
      j1    = j1 + 1
    end do
    do while (j2 < s1%ia(i+1))
      ja(k) = s1%ja(j2)
      a( k) = fact_ * s1%a(j2)
      k     = k  + 1
      j2    = j2 + 1
    end do
    ia(i+1) = k
  end do

  ! delete old csr data and replace it
  call move_alloc(ia, s2%ia)
  deallocate (s2%ja)
  deallocate (s2%a )
  allocate (s2%ja(k-1), source = ja(1:k-1))
  allocate (s2%a (k-1), source =  a(1:k-1))
end subroutine

module subroutine add_sparse3_{}T{}(s1, s2, s3, fact1, fact2)
  !! Add two sparse matrices, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * s1 + fact2 * s2

  type(sparse_{}T),     intent(in)  :: s1
    !! Sparse first matrix (unchanged)
  type(sparse_{}T),     intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(sparse_{}T),     intent(out) :: s3
    !! Output sparse third matrix
  m4_type(T), optional, intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  m4_type(T), optional, intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  m4_assert(s1%nrows == s2%nrows)
  m4_assert(s1%ncols == s2%ncols)

  ! init result matrix
  call s3%init(s1%nrows, ncols = s1%ncols)

  ! s3 <- fact1 * s1
  if (.not. s1%is_empty()) then
    s3%ia = s1%ia
    s3%ja = s1%ja
    if (present(fact1)) then
      s3%a = fact1 * s1%a
    else
      s3%a = s1%a
    end if
  end if

  ! s3 += fact2 * s2
  call matrix_add(s2, s3, fact = fact2)
end subroutine

module subroutine mul_dense_dense_{}T{}(d1, d2, result)
  !! multiply two dense matrices: result <- d1 * d2
  type(dense_{}T), intent(in)  :: d1
    !! First dense matrix (unchanged)
  type(dense_{}T), intent(in)  :: d2
    !! Second dense matrix (unchanged)
  type(dense_{}T), intent(out) :: result
    !! Resulting dense matrix

  ! init result matrix
  call result%init(d1%nrows, ncols = d2%ncols)

  ! matrix multiplication
  call gemm(d1%d, d2%d, result%d)
end subroutine

module subroutine mul_dense_sparse_{}T{}(d1, s2, result)
  !! multiply dense matrix by sparse matrix: result <- d1 * s2
  type(dense_{}T),  intent(in)  :: d1
    !! First matrix (dense, unchanged)
  type(sparse_{}T), intent(in)  :: s2
    !! Second matrix (sparse, unchanged)
  type(dense_{}T),  intent(out) :: result
    !! Resulting dense matrix

  type(dense_{}T) :: tmp

  ! convert sparse matrix to dense
  call tmp%init(s2%nrows, ncols = s2%ncols)
  call matrix_convert(s2, tmp)

  ! perform dense-dense multiplication
  call matrix_mul(d1, tmp, result)
end subroutine

module subroutine mul_sparse_dense_{}T{}(s1, d2, result)
  !! multiply sparse matrix by dense matrix: result <- s1 * d2
  type(sparse_{}T), intent(in)  :: s1
    !! First matrix (sparse, unchanged)
  type(dense_{}T),  intent(in)  :: d2
    !! Second matrix (dense, unchanged)
  type(dense_{}T),  intent(out) :: result
    !! Resulting dense matrix

  ! use sparse-raw matrix multiplication routine
  call result%init(s1%nrows, ncols = d2%ncols)
  call s1%mul_mat(d2%d, result%d)
end subroutine

module subroutine mul_sparse_sparse_{}T{}(s1, s2, result)
  !! multiply two sparse matrices: result <- s1 * s2
  type(sparse_{}T), intent(in)  :: s1
    !! First sparse matrix (unchanged)
  type(sparse_{}T), intent(in)  :: s2
    !! Second sparse matrix (unchanged)
  type(sparse_{}T), intent(out) :: result
    !! Resulting sparse matrix

  integer                          :: i, nthreads, ithread
  integer(SPARSE_IDX)              :: j, k, icol0, icol1, ncol
  integer(SPARSE_IDX), allocatable :: ops(:)
  integer,             allocatable :: i1(:), col(:)
  logical,             allocatable :: fil(:)
  m4_type(T)                       :: next_val
  m4_type(T),          allocatable :: val(:), val1(:)

  ! check dimensions
  m4_assert(s1%ncols == s2%nrows)

  ! init result
  call result%init(s1%nrows, ncols = s2%ncols)
  if (s1%is_empty() .or. s2%is_empty()) then
    ! return empty matrix if either s1 or s2 is empty
    return
  end if
  allocate (result%ia(result%nrows + 1), source = int(0, kind = SPARSE_IDX))

  ! multiplication operation counters for each result row
  allocate (ops(0:s1%nrows), source = int(0, kind = SPARSE_IDX))

  ! work in parallel
  !$omp parallel default(none) &
  !$omp private(i, j, k, ithread, col, icol0, icol1, ncol, val, val1, next_val, fil) &
  !$omp shared(s1, s2, result, sort_, nthreads, i1, ops)

    !$omp single
    ! get actual number of threads
    nthreads = omp_get_num_threads()

    ! allocate upper work chunk index bounds
    allocate (i1(0:nthreads), source = 0)
    !$omp end single nowait

    ! get current thread id
    ithread = omp_get_thread_num() + 1

    ! count number of multiplications per row
    !$omp do schedule(dynamic)
    do i = 1, s1%nrows
      do j = s1%ia(i), s1%ia(i+1)-1
        ops(i) = ops(i) + s2%ia(s1%ja(j)+1) - s2%ia(s1%ja(j))
      end do
    end do
    !$omp end do

    ! prefix sum (parallelization here probably not worth it)
    !$omp single
    do i = 2, s1%nrows
      ops(i) = ops(i) + ops(i-1)
    end do
    !$omp end single

    ! give every thread roughly the same number of ops
    if (ithread < nthreads) then
      ! equidistant index
      j = nint(ops(s1%nrows) * (real(ithread) / real(nthreads)), kind = SPARSE_IDX)

      ! find row index
      i1(ithread) = bin_search(ops, j) - 1
    else
      i1(ithread) = s1%nrows
    end if
    !$omp barrier

    ! get upper limit on memory required
    j = ops(i1(ithread)) - ops(i1(ithread-1))
    j = max(0, int(min(int(j, kind=int64), int(i1(ithread) - i1(ithread-1), kind=int64) * int(s2%ncols, kind=int64)), kind=int32))

    ! allocate thread local memory
    allocate (col(j))
    allocate (val(j))
    allocate (val1(s2%ncols))
    allocate (fil( s2%ncols), source = .false.)

    ! reset column counters
    icol1 = 0
    ncol  = 0

    ! calculate result rows
    do i = i1(ithread-1) + 1, i1(ithread)
      ! loop over columns of first matrix
      do j = s1%ia(i), s1%ia(i+1)-1
        ! loop over columns of second matrix where the row is equal to column of the first matrix
        do k = s2%ia(s1%ja(j)), s2%ia(s1%ja(j)+1)-1
          ! calculate new value
          next_val = s1%a(j) * s2%a(k)

          ! check whether element already exists
          if (fil(s2%ja(k))) then
            ! update existing value
            val1(s2%ja(k)) = val1(s2%ja(k)) + next_val
          else
            ! set new value and fill flag
            val1(s2%ja(k)) = next_val
            fil( s2%ja(k)) = .true.

            ! save column index
            ncol = ncol + 1
            col(ncol) = s2%ja(k)
          end if
        end do
      end do

      ! update column pointers
      icol0 = icol1 + 1
      icol1 = ncol

      ! sort column indices
      call qsort(col(icol0:icol1))

      ! save column pointers (at this point: indices to thread-local memory)
      result%ia(i+1) = icol1 + 1

      ! save values
      val(icol0:icol1) = val1(col(icol0:icol1))

      ! reset fill flag (where it was set)
      if (i < i1(ithread)) fil(col(icol0:icol1)) = .false.
    end do
    !$omp barrier

    !$omp single
      ! update column pointers (change into indices to global memory)
      result%ia(1) = 1
      do i = 1, nthreads
        result%ia(i1(i-1)+2:i1(i)+1) = result%ia(i1(i-1)+2:i1(i)+1) + (result%ia(i1(i-1)+1) - 1)
      end do

      ! allocate global arrays
      allocate (result%ja(result%ia(result%nrows+1)-1))
      allocate (result%a (result%ia(result%nrows+1)-1))
    !$omp end single

    ! copy thread-local data to global arrays
    result%ja(result%ia(i1(ithread-1)+1):result%ia(i1(ithread)+1)-1) = col(1:ncol)
    result%a( result%ia(i1(ithread-1)+1):result%ia(i1(ithread)+1)-1) = val(1:ncol)

    ! deallocate thread-local memory (avoid memory leak due to openmp bug)
    deallocate (val1, fil, col, val)
  !$omp end parallel
end subroutine

! fixme add testing
module subroutine diag_arr_band_{}T{}(arr_in, b_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  m4_type(T),     intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(band_{}T), intent(out) :: b_out
    !! dense output matrix

  call b_out%init(size(arr_in), 0, d0=reshape(arr_in, [1, size(arr_in)]))
end subroutine

module subroutine diag_arr_dense_{}T{}(arr_in, d_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  m4_type(T),      intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(dense_{}T), intent(out) :: d_out
    !! dense output matrix

  integer :: i

  call d_out%init(size(arr_in))
  do i = 1, size(arr_in)
    d_out%d(i,i) = arr_in(i)
  end do
end subroutine

module subroutine diag_arr_sparse_{}T{}(arr_in, s_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  m4_type(T),       intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(sparse_{}T), intent(out) :: s_out
    !! sparse output matrix

  integer :: i, n

  n = size(arr_in)
  call s_out%init(n)

  s_out%a  = arr_in
  s_out%ia = [(int(i, kind = SPARSE_IDX), i = 1, n+1)]
  s_out%ja = [(    i,                     i = 1, n  )]
end subroutine

! fixme add testing
module subroutine diag_band_arr_{}T{}(b_in, arr_out)
  !! extract main diagonal of a band matrix

  type(band_{}T), intent(in)  :: b_in
    !! band input matrix
  m4_type(T),     intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  arr_out = b_in%d(0,:)
end subroutine

module subroutine diag_dense_arr_{}T{}(d_in, arr_out)
  !! extract main diagonal of a dense matrix

  type(dense_{}T), intent(in)  :: d_in
    !! dense input matrix
  m4_type(T),      intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  integer :: i, n

  n = min(d_in%nrows, d_in%ncols)
  m4_assert(n == size(arr_out))

  arr_out = [(d_in%d(i,i), i = 1, n)]
end subroutine

module subroutine diag_sparse_arr_{}T{}(s_in, arr_out)
  !! extract main diagonal of a sparse matrix

  type(sparse_{}T), intent(in)  :: s_in
    !! sparse input matrix
  m4_type(T),       intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  integer             :: i, n
  integer(SPARSE_IDX) :: j

  n = min(s_in%nrows, s_in%ncols)
  m4_assert(n == size(arr_out))

  ! empty matrix => zero diagonal
  if (s_in%is_empty()) then
    arr_out = 0
    return
  end if

  !$omp parallel do   &
  !$omp default(none) &
  !$omp private(i,j)  &
  !$omp shared(n,arr_out,s_in)
  do i = 1, n
    ! get value closest to main diagonal
    j = s_in%ia(i) + (bin_search(s_in%ja(s_in%ia(i):s_in%ia(i+1)-1), i) - 1)

    ! output value only, if found exactly
    if (s_in%ja(j) == i) then
      arr_out(i) = s_in%a(j)
    else
      arr_out(i) = 0
    end if
  end do
  !$omp end parallel do
end subroutine

! fixme needs testing for complex case
module subroutine approx_sparse_{}T{}(thres, s1, s2)
  !! take row-wise approximation of matrix:
  !!
  !!      s2_ij = / 0,        if |s1_ij| <  thres * |s1_ii|
  !!              \ s1_ij,    else

  real,             intent(in)          :: thres
    !! relative threshold
  type(sparse_{}T), intent(in),  target :: s1
    !! sparse input matrix
  type(sparse_{}T), intent(out), target :: s2
    !! sparse output matrix

  integer           :: i, j, j0
  type(spbuild_{}T) :: sb

  ! load data from s1
  call sb%init(s1)

  ! apply threshold to each row
  do i = 1, s1%nrows
    ! find diagonal index j0, i.e. sb%cols(i)%d(j0) == i
    j0 = BIN_SEARCH(sb%cols(i)%d(1:sb%cols(i)%n), i)

    ! did we find diagonal index? if not keep all of rows elements
    if (sb%cols(i)%d(j0) /= i) cycle

    ! apply threshold
    do j = 1, sb%vals(i)%n
      if (abs(sb%vals(i)%d(j)) < thres * abs(sb%vals(i)%d(j0))) sb%filled(i)%d(j) = .false.
    end do
  end do

  ! save values in sb into s2
  call s2%init(s1%nrows, ncols=s1%ncols)
  sb%sp => s2
  call sb%save()
end subroutine
