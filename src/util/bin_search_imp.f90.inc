function bin_search_{}T{}(x, x1, mode) result(i)
  !! find index of nearest value in sorted array
  m4_type(T),        intent(in) :: x(:)
    !! sorted array
  m4_type(T),        intent(in) :: x1
    !! value to find
  integer, optional, intent(in) :: mode
    !! binary search mode (default: BS_NEAR)
  integer                       :: i
    !! return array index

  integer :: i0, i1, mode_

  ! mode
  mode_ = BS_NEAR
  if (present(mode)) mode_ = mode
  m4_assert(any(mode_ == [BS_NEAR, BS_LESS, BS_GREAT]))

  ! starting interval is the whole array
  i0 = 1
  i1 = size(x)

  ! test if x1 is outside of interval [x(1), x(end)]
  if (x1 <= x(i0)) then
    i = i0
    return
  end if
  if (x1 >= x(i1)) then
    i = i1
    return
  end if

  ! binary search
  do while (i1 > (i0 + 1))
    i = (i1 + i0) / 2
    if      (x(i) < x1) then
      i0 = i
    else if (x(i) > x1) then
      i1 = i
    else
      return
    end if
  end do

  select case (mode_)
    case (BS_NEAR)
      ! pick index of value that is closer
      if ((2 * x1) < (x(i0) + x(i1))) then
        i = i0
      else
        i = i1
      end if

    case (BS_LESS)
      ! pick smaller index
      i = i0

    case (BS_GREAT)
      ! pick larger index
      i = i1
  end select
end function

m4_undefine({T})
