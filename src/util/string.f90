module string_m

  implicit none

  type string
    character(:), allocatable :: s
  end type

  interface operator(<)
    module procedure :: string_lt_string
  end interface

  interface operator(<=)
    module procedure :: string_le_string
  end interface

  interface operator(>)
    module procedure :: string_gt_string
  end interface

  interface operator(>=)
    module procedure :: string_ge_string
  end interface

  interface operator(==)
    module procedure :: string_eq_string
  end interface

  interface operator(//)
    module procedure :: string_cat_char
    module procedure :: char_cat_string
    module procedure :: string_cat_string
  end interface

contains

  function new_string(c) result(s)
    !! create new string from character array

    character(*), intent(in) :: c
    type(string)             :: s

    s%s = c
  end function

  function string_lt_string(s1, s2) result(r)
    !! return true if s1 < s2 in collating sense

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    logical                  :: r

    r = s1%s < s2%s
  end function

  function string_le_string(s1, s2) result(r)
    !! return true if s1 <= s2 in collating sense

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    logical                  :: r

    r = s1%s <= s2%s
  end function

  function string_gt_string(s1, s2) result(r)
    !! return true if s1 > s2 in collating sense

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    logical                  :: r

    r = s1%s > s2%s
  end function

  function string_ge_string(s1, s2) result(r)
    !! return true if s1 >= s2 in collating sense

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    logical                  :: r

    r = s1%s >= s2%s
  end function

  function string_eq_string(s1, s2) result(r)
    !! return true if s1 == s2 in collating sense

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    logical                  :: r

    r = s1%s == s2%s
  end function

  function string_cat_char(s, c) result(r)
    !! Concatenate two character sequences.
    !! In this version the left-hand side character sequences is by a string.

    type(string), intent(in) :: s
    character(*), intent(in) :: c
    type(string)             :: r

    r%s = s%s // c
  end function

  function char_cat_string(c, s) result(r)
    !! Concatenate two character sequences.
    !! In this version the right-hand side character sequences is by a string.

    character(*), intent(in) :: c
    type(string), intent(in) :: s
    type(string)             :: r

    r%s = c // s%s
  end function

  function string_cat_string(s1, s2) result(r)
    !! Concatenate two character sequences.
    !! In this version both character sequences are by a string.

    type(string), intent(in) :: s1
    type(string), intent(in) :: s2
    type(string)             :: r

    r%s = s1%s // s2%s
  end function

end module
