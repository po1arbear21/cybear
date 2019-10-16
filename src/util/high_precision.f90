module high_precision_m
  implicit none

  type hp_real
    !! represents high precision value x + c

    real :: x
      !! principal value
    real :: c
      !! correction value
  end type

  interface operator (+)
    module procedure :: hp_real_add_hh
    module procedure :: hp_real_add_hd
    module procedure :: hp_real_add_dh
  end interface

  interface operator (-)
    module procedure :: hp_real_neg
    module procedure :: hp_real_sub_hh
    module procedure :: hp_real_sub_hd
    module procedure :: hp_real_sub_dh
  end interface

contains

  elemental function real_to_hp(r) result(h)
    !! convert real to high precision real

    real, intent(in) :: r
    type(hp_real)    :: h

    ! principal is r; correction zero
    h%x = r
    h%c = 0
  end function

  elemental function hp_to_real(h) result(r)
    !! convert high precision real to real

    type(hp_real), intent(in) :: h
    real                      :: r

    ! add principal and correction
    r = h%x + h%c
  end function

  elemental function hp_real_add_hh(h1, h2) result(h3)
    !! add two high precision reals

    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h1 + h2%x
    h3%c = h3%c + h2%c
  end function

  elemental function hp_real_add_hd(h1, r2) result(h3)
    !! add real to high precision real

    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    ! local variables
    real :: z

    ! add principals
    h3%x = h1%x + r2

    ! update correction
    z    = h3%x - h1%x
    h3%c = h1%c + (h1%x - (h3%x - z)) + (r2 - z)
  end function

  elemental function hp_real_add_dh(r1, h2) result(h3)
    !! add high precision real to real

    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h2 + r1
  end function

  elemental function hp_real_neg(h1) result(h2)
    !! negate high precision real

    type(hp_real), intent(in) :: h1
    type(hp_real)             :: h2

    ! negate principal and correction
    h2%x = - h1%x
    h2%c = - h1%c
  end function

  elemental function hp_real_sub_hh(h1, h2) result(h3)
    !! subtract two high precision reals

    type(hp_real), intent(in) :: h1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = h1 + (-h2)
  end function

  elemental function hp_real_sub_hd(h1, r2) result(h3)
    !! subtract real from high precision real

    type(hp_real), intent(in) :: h1
    real,          intent(in) :: r2
    type(hp_real)             :: h3

    h3 = h1 + (-r2)
  end function

  elemental function hp_real_sub_dh(r1, h2) result(h3)
    !! subtract high precision real from real

    real,          intent(in) :: r1
    type(hp_real), intent(in) :: h2
    type(hp_real)             :: h3

    h3 = r1 + (-h2)
  end function

  function hp_sum(a) result(s)
    !! compute sum with high precision

    real, intent(in) :: a(:)
      !! summands
    real             :: s
      !! result

    ! local variabls
    integer       :: i
    type(hp_real) :: t

    ! empty array
    if (size(a) < 1) then
      s = 0
      return
    end if

    ! init high precision real with first value
    t = real_to_hp(a(1))

    ! add rest of values to high precision real one at a time
    do i = 2, size(a)
      t = t + a(i)
    end do

    ! convert result back to real
    s = hp_to_real(t)
  end function

end module high_precision_m
