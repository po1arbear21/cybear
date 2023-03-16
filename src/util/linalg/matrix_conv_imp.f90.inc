module subroutine sparse_to_band_{}T{}(s, b, i0, j0)
  !! Insert sparse matrix into a band matrix.
  type(sparse_{}T),  intent(in)    :: s
    !! Sparse matrix (unchanged)
  type(band_{}T),    intent(inout) :: b
    !! Band matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_, nlower, nupper
  integer(SPARSE_IDX) :: j

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! initialize
  if (.not. allocated(b%d)) then
    m4_assert((.not. present(i0)) .and. (.not. present(j0)))
    m4_assert(s%nrows == s%ncols)

    ! get nlower, nupper
    nlower = 0
    nupper = 0
    do i = 1, s%nrows
      do j = s%ia(i), s%ia(i+1)-1
        nupper = max(nupper, s%ja(j) - i)
        nlower = max(nlower, i - s%ja(j))
      end do
    end do

    ! initialize band matrix
    call b%init(s%nrows, nlower, nupper = nupper)
  end if

  ! insert values into band matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      call b%set_elem(i0_+i-1, j0_+s%ja(j)-1, s%a(j))
    end do
  end do
end subroutine

module subroutine band_to_dense_{}T{}(b, d, i0, j0)
  !! Insert band matrix into a dense matrix.

  type(band_{}T),    intent(in)    :: b
    !! Band matrix (unchanged)
  type(dense_{}T),   intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  integer :: i0_, j0_, i, j, k

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! initialize if necessary
  if (.not. allocated(d%d)) then
    m4_assert((.not. present(i0)) .and. (.not. present(j0)))
    call d%init(b%nrows)
  end if

  ! insert elements into dense matrix
  do j = 1, b%ncols
    do i = -b%nupper, b%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > b%nrows)) cycle

      d%d(i0_+k-1,j0_+j-1) = b%d(i,j)
    end do
  end do
end subroutine

module recursive subroutine block_to_dense_{}T{}(b, d, i0, j0)
  !! Insert block matrix into a dense matrix.

  class(block_{}T),  intent(in)    :: b
    !! Block matrix (unchanged)
  type(dense_{}T),   intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, ib, jb, jb0(2)

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  m4_assert(b%check_sizes())

  ! initialize if necessary
  if (.not. allocated(d%d)) then
    m4_assert((.not. present(i0)) .and. (.not. present(j0)))
    call d%init(b%nrows, ncols = b%ncols)
  end if

  ! insert elements blockwise
  associate (bi0 => b%i0, bj0 => b%j0)
    do ib = 1, b%nbrows
      ! determine range for column blocks
      if (b%diag) then
        jb0 = ib
      else if (b%tridiag) then
        jb0 = [ib-1, ib+1]
      else
        jb0 = [1, b%nbcols]
      end if

      ! loop column blocks
      do jb = jb0(1), jb0(2)
        if ((jb < 1) .or. (jb > b%nbcols)) cycle
        if (.not. associated(b%b(ib,jb)%p)) cycle

        select type (m => b%b(ib,jb)%p)
          class is (block_{}T)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (band_{}T)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (dense_{}T)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (sparse_{}T)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
        end select
      end do
    end do
  end associate
end subroutine

module subroutine dense_to_dense_{}T{}(d1, d2, i0, j0)
  !! Insert dense matrix into another dense matrix.

  type(dense_{}T),   intent(in)    :: d1
    !! First dense matrix (unchanged)
  type(dense_{}T),   intent(inout) :: d2
    !! Second dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in second matrix (default 1)

  ! local variables
  integer :: i0_, j0_, i1, j1

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! initialize if necessary
  if (.not. allocated(d2%d)) then
    m4_assert((.not. present(i0)) .and. (.not. present(j0)))
    call d2%init(d1%nrows, ncols = d1%ncols)
  end if

  ! end indices
  i1 = i0_ + d1%nrows - 1
  j1 = j0_ + d1%ncols - 1

  d2%d(i0_:i1, j0_:j1) = d1%d
end subroutine

module subroutine sparse_to_dense_{}T{}(s, d, i0, j0)
  !! Insert sparse matrix into a dense matrix.

  type(sparse_{}T),  intent(in)    :: s
    !! Sparse matrix (unchanged)
  type(dense_{}T),   intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j

  m4_assert(.not. s%is_empty())
  m4_assert(allocated(d%d))

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! initialize if necessary
  if (.not. allocated(d%d)) then
    m4_assert((.not. present(i0)) .and. (.not. present(j0)))
    call d%init(s%nrows, ncols = s%ncols)
  end if

  ! insert values into dense matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      d%d(i0_+i-1,j0_+s%ja(j)-1) = s%a(j)
    end do
  end do
end subroutine

module subroutine block_to_sparse_{}T{}(b, s, struct, drop_zeros)
  !! Create sparse matrix from block matrix.
  class(block_{}T),  intent(in)  :: b
    !! Block matrix (unchanged)
  type(sparse_{}T),  intent(out) :: s
    !! Output: Sparse matrix
  logical, optional, intent(in)  :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)  :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  type(spbuild_{}T) :: sb

  call s%init(b%i1(b%nbrows), ncols=b%j1(b%nbcols))
  call sb%init(s)
  call matrix_convert(b, sb, struct=struct, drop_zeros=drop_zeros)
  call sb%save()
end subroutine

module subroutine band_to_spbuild_{}T{}(b, sb, i0, j0, struct, drop_zeros)
  !! Insert band matrix into a sparse builder.
  type(band_{}T),    intent(in)    :: b
    !! Band matrix (unchanged)
  type(spbuild_{}T), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional, intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, i, j, k
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  m4_assert(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  m4_assert(associated(sb%sp))

  ! insert values
  do j = 1, b%ncols
    do i = -b%nupper, b%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > b%nrows)) cycle

      ! check structure
      if (present(struct)) then
        if (.not. struct(k,j)) cycle
      else if (drop_zeros_) then
        if (b%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+k-1, j0_+j-1, b%d(i,j), search=.false.)
    end do
  end do
end subroutine

module subroutine block_to_spbuild_{}T{}(b, sb, i0, j0, struct, drop_zeros)
  !! Insert block matrix into a sparse builder.
  class(block_{}T),  intent(in)    :: b
    !! Block matrix (unchanged)
  type(spbuild_{}T), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional, intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, ib, jb, jb0(2)
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  m4_assert(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  m4_assert(b%check_sizes())
  m4_assert(associated(sb%sp))

  ! insert elements blockwise
  associate (bi0 => b%i0, bi1 => b%i1, bj0 => b%j0, bj1 => b%j1)
    do ib = 1, b%nbrows
      if      (b%diag) then
        jb0 = ib
      else if (b%tridiag) then
        jb0(1) = max(ib-1, 1        )
        jb0(2) = min(ib+1, b%nbcols)
      else
        jb0 = [1, b%nbcols]
      end if

      do jb = jb0(1), jb0(2)
        if (.not. associated(b%b(ib,jb)%p)) cycle

        select type (m => b%b(ib,jb)%p)
        class is (block_{}T)
          if (present(struct)) then
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        type is (band_{}T)
          if (present(struct)) then
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        type is (dense_{}T)
          if (present(struct)) then
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        type is (sparse_{}T)
          if (present(struct)) then
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
          else
            call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
          end if
        end select
      end do
    end do
  end associate
end subroutine

module subroutine dense_to_spbuild_{}T{}(d, sb, i0, j0, struct, drop_zeros)
  !! Insert dense matrix into a sparse builder.
  type(dense_{}T),   intent(in)    :: d
    !! Dense matrix (unchanged)
  type(spbuild_{}T), intent(inout) :: sb
    !! Sparse matrix represented by sparse matrix builder (updated)
  integer, optional, intent(in)    :: i0
    !! Optional: row offset in sparse matrix
  integer, optional, intent(in)    :: j0
    !! Optional: col offset in sparse matrix
  logical, optional, intent(in)    :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)    :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  ! local variables
  integer :: i0_, j0_, i, j
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  m4_assert(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  m4_assert(i0_ > 0)
  m4_assert(i0_ + d%nrows -1 <= sb%sp%nrows)
  j0_ = 1
  if (present(j0)) j0_ = j0
  m4_assert(j0_ > 0)
  m4_assert(j0_ + d%ncols -1 <= sb%sp%ncols)
  if (present(struct)) then
    m4_assert(size(d%d, 1) == size(struct, 1) .and. size(d%d, 2) == size(struct, 2))
  end if

  m4_assert(associated(sb%sp))

  ! insert values
  do i = 1, d%nrows
    do j = 1, d%ncols
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,j)) cycle
      elseif (drop_zeros_) then
        if (d%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+j-1, d%d(i,j), search = .false.)
    end do
  end do
end subroutine

module subroutine sparse_to_spbuild_{}T{}(s, sb, i0, j0, struct, drop_zeros)
  !! Insert sparse matrix into another sparse builder.
  type(sparse_{}T),  intent(in)    :: s
    !! First sparse matrix (unchanged)
  type(spbuild_{}T), intent(inout) :: sb
    !! Second sparse matrix represented by sparse matrix builder (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in sparse matrix
  integer, optional, intent(in)    :: j0
    !! optional: col offset in sparse matrix
  logical, optional, intent(in)    :: struct(:,:)
    !! optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)    :: drop_zeros
    !! optional: determine sparsity structure by checking for zero entries

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j
  logical             :: drop_zeros_

  ! do nothing if empty
  if (s%is_empty()) return

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  m4_assert(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  m4_assert(associated(sb%sp))

  ! insert values into second sparse matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      ! check structure
      if (present(struct)) then
        if (.not. struct(i,s%ja(j))) cycle
      else if (drop_zeros_) then
        if (s%a(j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+i-1, j0_+s%ja(j)-1, s%a(j), search = .false.)
    end do
  end do
end subroutine

m4_divert(m4_ifelse(T,cmplx,0,-1))

  module subroutine to_real_band(bc_in, br_out)
    !! Convert band complex matrix to real by ignoring imaginary part.
    type(band_cmplx), intent(in)  :: bc_in
      !! Band complex matrix
    type(band_real),  intent(out) :: br_out
      !! Output real(bc_in)

    call br_out%init(bc_in%nrows, bc_in%nlower, nupper=bc_in%nupper, d0=real(bc_in%d))
  end subroutine

  module subroutine to_real_dense(dc_in, dr_out)
    !! Convert dense complex matrix to real by ignoring imaginary part.
    type(dense_cmplx), intent(in)  :: dc_in
      !! Dense complex matrix
    type(dense_real),  intent(out) :: dr_out
      !! Output real(dc_in)

    call dr_out%init(real(dc_in%d))
  end subroutine

  module subroutine to_real_sparse(sc_in, sr_out)
    !! Convert sparse complex matrix to real by ignoring imaginary part
    type(sparse_cmplx), intent(in)  :: sc_in
      !! Sparse complex matrix
    type(sparse_real),  intent(out) :: sr_out
      !! Output real(sc_in)

    ! init sr_out
    call sr_out%init(sc_in%nrows, ncols = sc_in%ncols)

    ! copy csr
    if (.not. sc_in%is_empty()) then
      allocate (sr_out%ia(size(sc_in%ia,kind=SPARSE_IDX)), source = sc_in%ia     )
      allocate (sr_out%ja(size(sc_in%ja,kind=SPARSE_IDX)), source = sc_in%ja     )
      allocate (sr_out%a( size(sc_in%a ,kind=SPARSE_IDX)), source = real(sc_in%a))
    end if
  end subroutine

m4_divert(m4_ifelse(T,real,0,-1))

  module subroutine to_cmplx_band(br_in, bc_out)
    !! Convert band real matrix to complex.
    type(band_real),  intent(in)  :: br_in
      !! Band real matrix
    type(band_cmplx), intent(out) :: bc_out
      !! Output cmplx(br_in)

    call bc_out%init(br_in%nrows, br_in%nlower, nupper=br_in%nupper, d0=cmplx(br_in%d))
  end subroutine

  module subroutine to_cmplx_dense(dr_in, dc_out)
    !! Convert dense real matrix to complex with zero imaginary part.
    type(dense_real),  intent(in)  :: dr_in
      !! Dense real matrix
    type(dense_cmplx), intent(out) :: dc_out
      !! Output cmplx(dr_in)

    call dc_out%init(cmplx(dr_in%d))
  end subroutine

  module subroutine to_cmplx_sparse(sr_in, sc_out)
    !! Convert sparse real matrix to complex
    type(sparse_real),  intent(in)  :: sr_in
      !! Sparse real matrix
    type(sparse_cmplx), intent(out) :: sc_out
      !! Output cmplx(sr_in)

    ! init sc_out
    call sc_out%init(sr_in%nrows, ncols = sr_in%ncols)

    ! copy csr
    if (.not. sr_in%is_empty()) then
      allocate (sc_out%ia(size(sr_in%ia,kind=SPARSE_IDX)), source = sr_in%ia      )
      allocate (sc_out%ja(size(sr_in%ja,kind=SPARSE_IDX)), source = sr_in%ja      )
      allocate (sc_out%a( size(sr_in%a ,kind=SPARSE_IDX)), source = cmplx(sr_in%a))
    end if
  end subroutine

m4_divert(0)
