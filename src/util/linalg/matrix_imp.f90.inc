subroutine matrix_{}T{}_init(this, dum, nrows, ncols)
  !! Initialize matrix.

  class(matrix_{}T), intent(out) :: this
    !! Matrix
  character(*),      intent(in)  :: dum
    !! Dummy, unused
  integer,           intent(in)  :: nrows
    !! Number of rows
  integer, optional, intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)

  m4_ignore(dum)

  this%nrows = nrows
  this%ncols = nrows
  if (present(ncols)) this%ncols = ncols

  m4_assert(this%nrows >= 0)
  m4_assert(this%ncols >= 0)

  this%factorized = .false.
end subroutine
