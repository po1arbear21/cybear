#define PASTE(X)          X
#define PASTE2(X)         PASTE(X)_
#define CONCATHELP(X, Y)  PASTE2(X)Y
#define CONCAT(X, Y)      CONCATHELP(X,Y)

#define BAND_TYPE         CONCAT(band,T)
#define BAND_PTR_TYPE     CONCAT(band_ptr,T)
#define BLOCK_TYPE        CONCAT(block,T)
#define BLOCK_PTR_TYPE    CONCAT(block_ptr,T)
#define DENSE_TYPE        CONCAT(dense,T)
#define DENSE_PTR_TYPE    CONCAT(dense_ptr,T)
#define MATRIX_TYPE       CONCAT(matrix,T)
#define MATRIX_PTR_TYPE   CONCAT(matrix_ptr,T)
#define SPARSE_TYPE       CONCAT(sparse,T)
#define SPARSE_PTR_TYPE   CONCAT(sparse_ptr,T)
#define SPBUILD_TYPE      CONCAT(spbuild,T)

#define SPARSE_TO_BAND    CONCAT(sparse_to_band,T)

#define BAND_TO_DENSE     CONCAT(band_to_dense,T)
#define BLOCK_TO_DENSE    CONCAT(block_to_dense,T)
#define DENSE_TO_DENSE    CONCAT(dense_to_dense,T)
#define SPARSE_TO_DENSE   CONCAT(sparse_to_dense,T)

#define BLOCK_TO_SPARSE   CONCAT(block_to_sparse,T)

#define BAND_TO_SPBUILD   CONCAT(band_to_spbuild,T)
#define BLOCK_TO_SPBUILD  CONCAT(block_to_spbuild,T)
#define DENSE_TO_SPBUILD  CONCAT(dense_to_spbuild,T)
#define SPARSE_TO_SPBUILD CONCAT(sparse_to_spbuild,T)

module subroutine SPARSE_TO_BAND(s, b, i0, j0)
  !! Insert sparse matrix into a band matrix.
  type(SPARSE_TYPE), intent(in)    :: s
    !! Sparse matrix (unchanged)
  type(BAND_TYPE),   intent(inout) :: b
    !! Band matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j

  ASSERT(.not. s%is_empty())
  ASSERT(allocated(b%d))

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values into band matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      call b%set_elem(i0_+i-1, j0_+s%ja(j)-1, s%a(j))
    end do
  end do
end subroutine

module subroutine BAND_TO_DENSE(b, d, i0, j0)
  !! Insert band matrix into a dense matrix.

  type(BAND_TYPE),   intent(in)    :: b
    !! Band matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
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

  ! insert elements into dense matrix
  do j = 1, b%ncols
    do i = -b%nupper, b%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > b%nrows)) cycle

      d%d(i0_+k-1,j0_+j-1) = b%d(i,j)
    end do
  end do
end subroutine

module recursive subroutine BLOCK_TO_DENSE(b, d, i0, j0)
  !! Insert block matrix into a dense matrix.

  class(BLOCK_TYPE), intent(in)    :: b
    !! Block matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
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

  ASSERT(b%check_sizes())

  ! insert elements blockwise
  associate (bi0 => b%i0, bj0 => b%j0)
    do ib = 1, b%nbrows
      ! determine range for column blocks
      if      (b%diag   ) then
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
          class is (BLOCK_TYPE)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (BAND_TYPE)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (DENSE_TYPE)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
          type is (SPARSE_TYPE)
            call matrix_convert(m, d, i0=i0_+bi0(ib)-1, j0=j0_+bj0(jb)-1)
        end select
      end do
    end do
  end associate
end subroutine

module subroutine DENSE_TO_DENSE(d1, d2, i0, j0)
  !! Insert dense matrix into another dense matrix.

  type(DENSE_TYPE),  intent(in)    :: d1
    !! First dense matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d2
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

  ! end indices
  i1 = i0_ + d1%nrows - 1
  j1 = j0_ + d1%ncols - 1

  d2%d(i0_:i1, j0_:j1) = d1%d
end subroutine

module subroutine SPARSE_TO_DENSE(s, d, i0, j0)
  !! Insert sparse matrix into a dense matrix.

  type(SPARSE_TYPE), intent(in)    :: s
    !! Sparse matrix (unchanged)
  type(DENSE_TYPE),  intent(inout) :: d
    !! Dense matrix (updated)
  integer, optional, intent(in)    :: i0
    !! optional: row offset in second matrix (default 1)
  integer, optional, intent(in)    :: j0
    !! optional: col offset in second matrix (default 1)

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j

  ASSERT(.not. s%is_empty())
  ASSERT(allocated(d%d))

  ! optional arguments
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ! insert values into dense matrix
  do i = 1, s%nrows
    do j = s%ia(i), s%ia(i+1)-1
      d%d(i0_+i-1,j0_+s%ja(j)-1) = s%a(j)
    end do
  end do
end subroutine

module subroutine BLOCK_TO_SPARSE(b, s, struct, drop_zeros)
  !! Create sparse matrix from block matrix.
  class(BLOCK_TYPE), intent(in)  :: b
    !! Block matrix (unchanged)
  type(SPARSE_TYPE), intent(out) :: s
    !! Output: Sparse matrix
  logical, optional, intent(in)  :: struct(:,:)
    !! Optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional, intent(in)  :: drop_zeros
    !! Optional: determine sparsity structure by checking for zero entries

  type(SPBUILD_TYPE) :: sb

  call s%init(b%i1(b%nbrows), ncols=b%j1(b%nbcols))
  call sb%init(s)
  call BLOCK_TO_SPBUILD(b, sb, struct=struct, drop_zeros=drop_zeros)
  call sb%save()
end subroutine

module subroutine BAND_TO_SPBUILD(b, sb, i0, j0, struct, drop_zeros)
  !! Insert band matrix into a sparse builder.
  type(BAND_TYPE),    intent(in)    :: b
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
  integer :: i0_, j0_, i, j, k
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
  do j = 1, b%ncols
    do i = -b%nupper, b%nlower
      k = j + i ! actual row
      if ((k < 1) .or. (k > b%nrows)) cycle

      ! check structure
      if      (present(struct)) then
        if (.not. struct(k,j)) cycle
      else if (drop_zeros_) then
        if (b%d(i,j) == 0) cycle
      end if

      ! insert value without searching (append to rows)
      call sb%set(i0_+k-1, j0_+j-1, b%d(i,j), search=.false.)
    end do
  end do
end subroutine

module subroutine BLOCK_TO_SPBUILD(b, sb, i0, j0, struct, drop_zeros)
  !! Insert block matrix into a sparse builder.
  class(BLOCK_TYPE),  intent(in)    :: b
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
  integer :: i0_, j0_, ib, jb, jb0(2)
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

  ASSERT(b%check_sizes())

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
          class is (BLOCK_TYPE)
            if (present(struct)) then
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
            else
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
            end if
          type is (BAND_TYPE)
            if (present(struct)) then
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
            else
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
            end if
          type is (DENSE_TYPE)
            if (present(struct)) then
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, struct = struct(bi0(ib):bi1(ib), bj0(jb):bj1(jb)))
            else
              call matrix_convert(m, sb, i0 = i0_+bi0(ib)-1, j0 = j0_+bj0(jb)-1, drop_zeros = drop_zeros)
            end if
          type is (SPARSE_TYPE)
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

module subroutine DENSE_TO_SPBUILD(d, sb, i0, j0, struct, drop_zeros)
  !! Insert dense matrix into a sparse builder.
  type(DENSE_TYPE),   intent(in)    :: d
    !! Dense matrix (unchanged)
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
  integer :: i0_, j0_, i, j
  logical :: drop_zeros_

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  ASSERT(i0_ > 0)
  ASSERT(i0_ + d%nrows -1 <= sb%sp%nrows)
  j0_ = 1
  if (present(j0)) j0_ = j0
  ASSERT(j0_ > 0)
  ASSERT(j0_ + d%ncols -1 <= sb%sp%ncols)
  if (present(struct)) then
    ASSERT(size(d%d, 1) == size(struct, 1) .and. size(d%d, 2) == size(struct, 2))
  end if

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

module subroutine SPARSE_TO_SPBUILD(s, sb, i0, j0, struct, drop_zeros)
  !! Insert sparse matrix into another sparse builder.
  type(SPARSE_TYPE),  intent(in)    :: s
    !! First sparse matrix (unchanged)
  type(SPBUILD_TYPE), intent(inout) :: sb
    !! Second sparse matrix represented by sparse matrix builder (updated)
  integer, optional,  intent(in)    :: i0
    !! optional: row offset in sparse matrix
  integer, optional,  intent(in)    :: j0
    !! optional: col offset in sparse matrix
  logical, optional,  intent(in)    :: struct(:,:)
    !! optional: use sparsity structure (true means an element exists; false means it does not)
  logical, optional,  intent(in)    :: drop_zeros
    !! optional: determine sparsity structure by checking for zero entries

  integer             :: i, i0_, j0_
  integer(SPARSE_IDX) :: j
  logical             :: drop_zeros_

  ! do nothing if empty
  if (s%is_empty()) return

  ! optional arguments
  drop_zeros_ = .false.
  if (present(drop_zeros)) drop_zeros_ = drop_zeros
  ASSERT(.not. (present(struct) .and. drop_zeros_))
  i0_ = 1
  if (present(i0)) i0_ = i0
  j0_ = 1
  if (present(j0)) j0_ = j0

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

#ifdef TCMPLX
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
#endif

#ifndef TCMPLX
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
#endif

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef BAND_TYPE
#undef BAND_PTR_TYPE
#undef BLOCK_TYPE
#undef BLOCK_PTR_TYPE
#undef DENSE_TYPE
#undef DENSE_PTR_TYPE
#undef MATRIX_TYPE
#undef MATRIX_PTR_TYPE
#undef SPARSE_TYPE
#undef SPARSE_PTR_TYPE
#undef SPBUILD_TYPE

#undef SPARSE_TO_BAND

#undef BAND_TO_DENSE
#undef BLOCK_TO_DENSE
#undef DENSE_TO_DENSE
#undef SPARSE_TO_DENSE

#undef BLOCK_TO_SPARSE

#undef BAND_TO_SPBUILD
#undef BLOCK_TO_SPBUILD
#undef DENSE_TO_SPBUILD
#undef SPARSE_TO_SPBUILD
