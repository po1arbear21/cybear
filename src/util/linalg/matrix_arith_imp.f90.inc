#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

#define BAND_TYPE            CONCAT(band,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define SPARSE_TYPE          CONCAT(sparse,T)

#define ADD_BAND2            CONCAT(add_band2,T)
#define ADD_BAND3            CONCAT(add_band3,T)
#define ADD_BAND_DENSE       CONCAT(add_band_dense,T)
#define ADD_BAND_DENSE2      CONCAT(add_band_dense2,T)
#define ADD_BAND_SPARSE      CONCAT(add_band_sparse,T)
#define ADD_BAND_SPARSE2     CONCAT(add_band_sparse2,T)
#define ADD_DENSE2           CONCAT(add_dense2,T)
#define ADD_DENSE3           CONCAT(add_dense3,T)
#define ADD_SPARSE_DENSE     CONCAT(add_sparse_dense,T)
#define ADD_SPARSE_DENSE2    CONCAT(add_sparse_dense2,T)
#define ADD_SPARSE2          CONCAT(add_sparse2,T)
#define ADD_SPARSE3          CONCAT(add_sparse3,T)

#define DIAG_ARR_BAND        CONCAT(diag_arr_band,T)
#define DIAG_ARR_DENSE       CONCAT(diag_arr_dense,T)
#define DIAG_ARR_SPARSE      CONCAT(diag_arr_sparse,T)
#define DIAG_BAND_ARR        CONCAT(diag_band_arr,T)
#define DIAG_DENSE_ARR       CONCAT(diag_dense_arr,T)
#define DIAG_SPARSE_ARR      CONCAT(diag_sparse_arr,T)

! fixme needs testing
module subroutine ADD_BAND2(b1, b2, fact)
  !! Add two band matrices.
  !! Operation: b2 <- fact * b1 + b2

  type(BAND_TYPE), intent(in)    :: b1
    !! First band matrix (unchanged)
  type(BAND_TYPE), intent(inout) :: b2
    !! Second band matrix (updated)
  TT, optional,    intent(in)    :: fact
    !! First band matrix' scaling factor (default: 1)

  TT              :: fact_
  TT, allocatable :: tmp(:,:)

  ASSERT(.not. allocated(b1%f))
  ASSERT(.not. allocated(b2%ipiv))
  ASSERT(b2%nrows == b1%nrows .and. b2%ncols == b1%ncols)

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

! fixme needs testing
module subroutine ADD_BAND3(b1, b2, b3, fact1, fact2)
  !! Add two band matrices, store result in third band matrix.
  !! Operation: b3 <- fact1*b1 + fact2*b2

  type(BAND_TYPE), intent(in)  :: b1
    !! Band first matrix (unchanged)
  type(BAND_TYPE), intent(in)  :: b2
    !! Band second matrix (unchanged)
  type(BAND_TYPE), intent(out) :: b3
    !! Output band third matrix
  TT, optional,    intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,    intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  TT :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! case selection: hopefully no tmp array needs to be allocated
  if ((b1%nupper >= b2%nupper) .and. (b1%nlower >= b2%nlower)) then
    call b3%init(b1%nrows, b1%nlower, nupper=b1%nupper, d0=fact1_*b1%d)     ! b3 <- fact1*b1
    call add(b2, b3, fact=fact2_)                                           ! b3 += fact2*b2
  else
    call b3%init(b2%nrows, b2%nlower, nupper=b2%nupper, d0=fact2_*b2%d)     ! b3 <- fact2*b2
    call add(b1, b3, fact=fact1_)                                           ! b3 += fact1*b1
  end if
end subroutine

module subroutine ADD_BAND_DENSE(b, d, fact)
  !! Add dense and band matrix, store result in dense matrix.
  !!
  !! operation:  d <- fact * b + d

  type(BAND_TYPE),  intent(in)    :: b
    !! Band matrix b (unchanged)
  type(DENSE_TYPE), intent(inout) :: d
    !! Dense matrix d (updated)
  TT, optional,     intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  integer :: i, j, k
  TT      :: fact_

  ASSERT(d%nrows == b%nrows)
  ASSERT(d%ncols == b%ncols)
  ASSERT(.not. allocated(d%lu  ))
  ASSERT(.not. allocated(d%ipiv))

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

module subroutine ADD_BAND_DENSE2(b1, d2, d3, fact1, fact2)
  !! Add dense and band matrix, store result in a third dense matrix.
  !!
  !! operation:  d3 <- fact1 * b1 + fact2 * d2

  type(BAND_TYPE),  intent(in)  :: b1
    !! Band first matrix b1 (unchanged)
  type(DENSE_TYPE), intent(in)  :: d2
    !! Dense second matrix d1 (unchanged)
  type(DENSE_TYPE), intent(out) :: d3
    !! Output dense third matrix d3
  TT, optional,     intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,     intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  TT :: fact2_

  ASSERT(d2%nrows == b1%nrows)
  ASSERT(d2%ncols == b1%ncols)

  ! optional arguments
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! d3 <- fact2 * d2
  call d3%init(fact2_ * d2%d)

  ! d3 += fact1 * b1
  call add(b1, d3, fact = fact1)
end subroutine

module subroutine ADD_BAND_SPARSE(b, s, fact)
  !! Add sparse and band matrix, store result in sparse matrix.
  !!
  !! operation: s <- fact * b + s

  type(BAND_TYPE),   intent(in)    :: b
    !! Band matrix (unchanged)
  type(SPARSE_TYPE), intent(inout) :: s
    !! Sparse matrix (updated)
  TT, optional,      intent(in)    :: fact
    !! Band matrix scaling factor (default: 1)

  integer                          :: i, j2
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, k, b_nrows, b_nlower, b_nupper
  integer(SPARSE_IDX), allocatable :: ia(:)
  TT                               :: fact_
  TT,                  allocatable :: a(:)

  ASSERT(s%nrows == b%nrows)
  ASSERT(s%ncols == b%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(s%solver_handle == 0)

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

module subroutine ADD_BAND_SPARSE2(b1, s2, s3, fact1, fact2)
  !! Add sparse and band matrix, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * b1 + fact2 * s2

  type(BAND_TYPE),   intent(in)  :: b1
    !! Band first matrix (unchanged)
  type(SPARSE_TYPE), intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE), intent(out) :: s3
    !! Output sparse third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ASSERT(s2%nrows == b1%nrows)
  ASSERT(s2%ncols == b1%ncols)

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
  call add(b1, s3, fact = fact1)
end subroutine

module subroutine ADD_DENSE2(d1, d2, fact)
  !! Add two dense matrices, store the result in second matrix.
  !!
  !! operation:  d2 <- fact * d1 + d2

  type(DENSE_TYPE), intent(in)    :: d1
    !! First dense matrix d1 (unchanged)
  type(DENSE_TYPE), intent(inout) :: d2
    !! Second dense matrix d2 (updated)
  TT, optional,     intent(in)    :: fact
    !! Second dense matrix scaling factor (default: 1)

  ASSERT(.not. allocated(d2%lu))
  ASSERT(.not. allocated(d2%ipiv))

  if (present(fact)) then
    d2%d = fact * d1%d + d2%d
  else
    d2%d =        d1%d + d2%d
  end if
end subroutine

module subroutine ADD_DENSE3(d1, d2, d3, fact1, fact2)
  !! Add two dense matrices, store the result in third dense matrix.
  !!
  !! operation:  d3 <- fact1 * d1 + fact2 * d2

  type(DENSE_TYPE), intent(in)  :: d1
    !! Dense first matrix d1 (unchanged)
  type(DENSE_TYPE), intent(in)  :: d2
    !! Dense second matrix d2 (unchanged)
  type(DENSE_TYPE), intent(out) :: d3
    !! Output dense third matrix d3
  TT, optional,     intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,     intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  TT :: fact1_, fact2_

  ! optional arguments
  fact1_ = 1
  if (present(fact1)) fact1_ = fact1
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! init output matrix
  call d3%init(fact1_ * d1%d + fact2_ * d2%d)
end subroutine

module subroutine ADD_SPARSE_DENSE(s, d, fact)
  !! Add dense and sparse matrix, store result in dense matrix.
  !!
  !! operation:  D <- fact * S + D

  type(SPARSE_TYPE), intent(in)    :: s
    !! Sparse matrix S (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Dense matrix D (updated)
  TT, optional,      intent(in)    :: fact
    !! Sparse matrix scaling factor (default: 1)

  integer             :: i
  integer(SPARSE_IDX) :: j
  TT                  :: fact_

  ASSERT(d%nrows == s%nrows)
  ASSERT(d%ncols == s%ncols)
  ASSERT(.not. allocated(d%lu))
  ASSERT(.not. allocated(d%ipiv))

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

module subroutine ADD_SPARSE_DENSE2(s1, d2, d3, fact1, fact2)
  !! Add dense and sparse matrix, store result in another dense matrix.
  !!
  !! operation:  d3 <- fact1 * s1 + fact2 * d2

  type(SPARSE_TYPE), intent(in)  :: s1
    !! Sparse first matrix s1 (unchanged)
  type(DENSE_TYPE),  intent(in)  :: d2
    !! Dense second matrix d1 (unchanged)
  type(DENSE_TYPE),  intent(out) :: d3
    !! Output dense third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  TT :: fact2_

  ASSERT(d2%nrows == s1%nrows)
  ASSERT(d2%ncols == s1%ncols)

  ! optional argument
  fact2_ = 1
  if (present(fact2)) fact2_ = fact2

  ! d3 <- fact2 * d2
  call d3%init(fact2_ * d2%d)

  ! d3 += fact1 * s1
  if (.not. s1%is_empty()) call add(s1, d3, fact = fact1)
end subroutine

module subroutine ADD_SPARSE2(s1, s2, fact)
  !! Add two sparse matrices, store result in first matrix.
  !!
  !! operation: s2 <- fact * s1 + s2

  type(SPARSE_TYPE), intent(in)    :: s1
    !! First sparse matrix (unchanged)
  type(SPARSE_TYPE), intent(inout) :: s2
    !! Second sparse matrix (updated)
  TT, optional,      intent(in)    :: fact
    !! Second sparse matrix scaling factor (default: 1)

  integer                          :: i
  integer,             allocatable :: ja(:)
  integer(SPARSE_IDX)              :: j1, j2, k
  integer(SPARSE_IDX), allocatable :: ia(:)
  TT                               :: fact_
  TT,                  allocatable :: a(:)

  ASSERT(s2%nrows == s1%nrows)
  ASSERT(s2%ncols == s1%ncols)

  ! h == 0 is default value for unfactorized
  ASSERT(s2%solver_handle == 0)

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

module subroutine ADD_SPARSE3(s1, s2, s3, fact1, fact2)
  !! Add two sparse matrices, store result in sparse third matrix.
  !!
  !! operation: s3 <- fact1 * s1 + fact2 * s2

  type(SPARSE_TYPE), intent(in)  :: s1
    !! Sparse first matrix (unchanged)
  type(SPARSE_TYPE), intent(in)  :: s2
    !! Sparse second matrix (unchanged)
  type(SPARSE_TYPE), intent(out) :: s3
    !! Output sparse third matrix
  TT, optional,      intent(in)  :: fact1
    !! First matrix scaling factor (default: 1)
  TT, optional,      intent(in)  :: fact2
    !! Second matrix scaling factor (default: 1)

  ASSERT(s1%nrows == s2%nrows)
  ASSERT(s1%ncols == s2%ncols)

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
  call add(s2, s3, fact = fact2)
end subroutine

! fixme add testing
module subroutine DIAG_ARR_BAND(arr_in, b_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  TT,              intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(BAND_TYPE), intent(out) :: b_out
    !! dense output matrix

  call b_out%init(size(arr_in), 0, d0=reshape(arr_in, [1, size(arr_in)]))
end subroutine

module subroutine DIAG_ARR_DENSE(arr_in, d_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  TT,               intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(DENSE_TYPE), intent(out) :: d_out
    !! dense output matrix

  integer :: i

  call d_out%init(size(arr_in))
  do i = 1, size(arr_in)
    d_out%d(i,i) = arr_in(i)
  end do
end subroutine

module subroutine DIAG_ARR_SPARSE(arr_in, s_out)
  !! creates a quadratic diagonal matrix given the diagonals as an array.

  TT,                intent(in)  :: arr_in(:)
    !! diagonal elements of matrix
  type(SPARSE_TYPE), intent(out) :: s_out
    !! sparse output matrix

  integer :: i, n

  n = size(arr_in)
  call s_out%init(n)

  s_out%a  = arr_in
  s_out%ia = [(int(i, kind = SPARSE_IDX), i = 1, n+1)]
  s_out%ja = [(    i,                     i = 1, n  )]
end subroutine

! fixme add testing
module subroutine DIAG_BAND_ARR(b_in, arr_out)
  !! extract main diagonal of a band matrix

  type(BAND_TYPE), intent(in)  :: b_in
    !! band input matrix
  TT,              intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  arr_out = b_in%d(0,:)
end subroutine

module subroutine DIAG_DENSE_ARR(d_in, arr_out)
  !! extract main diagonal of a dense matrix

  type(DENSE_TYPE), intent(in)  :: d_in
    !! dense input matrix
  TT,               intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  integer :: i, n

  n = min(d_in%nrows, d_in%ncols)
  ASSERT(n == size(arr_out))

  arr_out = [(d_in%d(i,i), i = 1, n)]
end subroutine

module subroutine DIAG_SPARSE_ARR(s_in, arr_out)
  !! extract main diagonal of a sparse matrix

  type(SPARSE_TYPE), intent(in)  :: s_in
    !! sparse input matrix
  TT,                intent(out) :: arr_out(:)
    !! diagonal elements of matrix

  integer             :: i, n
  integer(SPARSE_IDX) :: j

  n = min(s_in%nrows, s_in%ncols)
  ASSERT(n == size(arr_out))

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

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef BAND_TYPE
#undef DENSE_TYPE
#undef SPARSE_TYPE

#undef ADD_BAND2
#undef ADD_BAND3
#undef ADD_BAND_DENSE
#undef ADD_BAND_DENSE2
#undef ADD_BAND_SPARSE
#undef ADD_BAND_SPARSE2
#undef ADD_DENSE2
#undef ADD_DENSE3
#undef ADD_SPARSE_DENSE
#undef ADD_SPARSE_DENSE2
#undef ADD_SPARSE2
#undef ADD_SPARSE3

#undef DIAG_ARR_BAND
#undef DIAG_ARR_DENSE
#undef DIAG_ARR_SPARSE
#undef DIAG_BAND_ARR
#undef DIAG_DENSE_ARR
#undef DIAG_SPARSE_ARR
