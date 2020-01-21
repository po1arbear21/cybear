module util_m
  use error_m
  use iso_c_binding
  implicit none

  ! binary search modes
  integer, parameter :: BS_NEAR  = 1
  integer, parameter :: BS_LESS  = 2
  integer, parameter :: BS_GREAT = 3

  interface hash_int
    module procedure :: hash_int32, hash_int64
  end interface

contains

  pure function cstrlen(cstr) result(len)
    !! get length of c string
    character(len=1), intent(in) :: cstr(*)
    integer                      :: len

    ! go through cstr until null character is found
    len = 1
    do while (cstr(len) .ne. c_null_char)
      len = len + 1
    end do

    ! discard trailing null character
    len = len - 1
  end function

  pure function c2fstring(cstr) result(fstr)
    !! convert c string to fortran string
    character(len=1), intent(in) :: cstr(*)
    character(len=cstrlen(cstr)) :: fstr

    fstr = transfer(cstr(1:len(fstr)), fstr)
  end function

  pure function f2cstring(fstr) result(cstr)
    !! convert fortran string to c string
    character(len=*), intent(in) :: fstr
    character(len=1)             :: cstr(len_trim(fstr) + 1)

    ! local variables
    integer :: i

    do i = 1, size(cstr) - 1
      cstr(i) = fstr(i:i)
    end do
    cstr(size(cstr)) = c_null_char
  end function

  function int2str(i) result(str)
    !! convert integer to string
    integer,      intent(in)             :: i
    character(:),            allocatable :: str

    character(24) :: tmp

    ! write to temporary string
    write (tmp, "(I24)") i

    ! adjustl: leading blanks cut and appended at end
    ! trim:    remove trailing blanks
    str = trim(adjustl(tmp))
  end function

  function select_int(flags, ints) result(t)
    !! select an integer according to the flag that is set (no flags set => 0, multiple set => - 1)
    logical, intent(in) :: flags(:)
    integer, intent(in) :: ints(:)
    integer             :: t

    ! local variables
    integer             :: i, j

    t = 0
    j = 0
    do i = 1, size(flags)
      if (flags(i)) then
        t = ints(i)
        j = j + 1
      end if
    end do

    ! ambiguous
    if (j .gt. 1) then
      t = -1
    end if
  end function

  function hash_int32(i) result(h)
    !! 32-bit integer hash function (nullprogram.com/blog/2018/07/31/)
    integer(kind=4), intent(in) :: i
    integer(kind=4)             :: h

    h = i + z'7f4a7c15'
    h = ieor(h, shiftr(h, 17))
    h = h * z'ed5ad4bb'
    h = ieor(h, shiftr(h, 11))
    h = h * z'ac4c1b51'
    h = ieor(h, shiftr(h, 15))
    h = h * z'31848bab'
    h = ieor(h, shiftr(h, 14))
  end function

  function hash_int64(i) result(h)
    !! 64-bit integer hash function (splitmix64)
    integer(kind=8), intent(in) :: i
    integer(kind=8)             :: h

    h = i + z'9e3779b97f4a7c15'
    h = ieor(h, shiftr(h, 30))
    h = h * z'bf58476d1ce4e5b9'
    h = ieor(h, shiftr(h, 27))
    h = h * z'94d049bb133111eb'
    h = ieor(h, shiftr(h, 31))
  end function

  function bin_search(x, x1, mode) result(i)
    !! find index of nearest value in sorted array
    real,              intent(in) :: x(:)
      !! sorted array
    real,              intent(in) :: x1
      !! value to find
    integer, optional, intent(in) :: mode
      !! binary search mode (default: BS_NEAR)
    integer                       :: i
      !! return array index

    ! local variables
    integer :: i0, i1, mode_

    ! mode
    mode_ = BS_NEAR
    if (present(mode)) mode_ = mode

    ! starting interval is the whole array
    i0 = 1
    i1 = size(x)

    ! test if x1 is outside of interval [x(1), x(end)]
    if (x1 .lt. x(i0)) then
      i = i0
      return
    end if
    if (x1 .gt. x(i1)) then
      i = i1
      return
    end if

    ! binary search
    do while (i1 .gt. (i0 + 1))
      i = (i1 + i0) / 2
      if (x(i) .lt. x1) then
        i0 = i
      elseif (x(i) .gt. x1) then
        i1 = i
      else
        return
      end if
    end do

    select case (mode_)
      case (BS_NEAR)
        ! pick index of value that is closer
        if ((2 * x1) .lt. (x(i0) + x(i1))) then
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

      case default
        i = 0
        write(*,*) mode
        call program_error("bin_search: unknown search mode!")

    end select
  end function

  function linspace(x0, x1, nx) result(x)
    !! create array of linear spaced values
    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    ! local variables
    integer :: i
    real    :: dx

    ! spacing between values
    dx = (x1 - x0) / (nx - 1)

    x(1 ) = x0
    do i = 2, nx - 1
      x(i) = x0 + (i - 1) * dx
    end do
    x(nx) = x1
  end function

  function logspace(x0, x1, nx) result(x)
    !! create array of logarithmic spaced values
    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    ! local variables
    integer :: i
    real    :: e0, e1, de

    e0 = log(x0)
    e1 = log(x1)
    de = (e1 - e0) / (nx - 1)

    x(1 ) = x0
    do i = 2, nx - 1
      x(i) = exp(e0 + (i - 1) * de)
    end do
    x(nx) = x1
  end function

end module
