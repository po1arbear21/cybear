module test_string_m

  use string_m
  use test_case_m, only: test_case

  implicit none

  private
  public test_string

contains

  subroutine test_string()
    type(string)    :: str(4)
    type(test_case) :: tc

    call tc%init("string")

    str(1)%s = "Maxwell"
    str(2)%s = "Naxwell"
    str(3)%s = "maxwell"
    str(4)%s = " Maxwell"

    ! test1: less
    call tc%assert(str(1) <  str(2), "<")
    call tc%assert(str(1) <  str(3), "<")

    ! test2: less equal
    call tc%assert(str(1) <= str(1), "<=")
    call tc%assert(str(1) <= str(2), "<=")

    ! test3: greater
    call tc%assert(str(2) >  str(1), ">")
    call tc%assert(str(1) >  str(4), ">")

    ! test4: greater equal
    call tc%assert(str(2) >= str(1), ">=")
    call tc%assert(str(2) >= str(2), ">=")

    ! test5: equal
    call tc%assert(      (str(1) == str(1)), "==")
    call tc%assert(.not. (str(1) == str(2)), "==")
    call tc%assert(.not. (str(1) == str(3)), "==")
    call tc%assert(.not. (str(1) == str(4)), "==")

    call tc%finish()
  end subroutine

end module
