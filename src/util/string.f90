module string_m
  implicit none

  type string
    character(len=:), allocatable :: s
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

contains

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

end module