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

    ! test6: concat str <- str//chars
    str(1)%s = "tree"               ! lhs
    str(2)   = str(1) // "window"   ! result
    str(3)%s = "treewindow"         ! expected
    call tc%assert_eq(str(3), str(2), "concat str//chars")

    ! test7: concat str <- chars//str
    str(1)%s = "display"            ! rhs
    str(2)   = "mouse" // str(1)    ! result
    str(3)%s = "mousedisplay"       ! expected
    call tc%assert_eq(str(3), str(2), "concat chars//str")

    ! test8: concat str <- str//str
    str(1)%s = "pen"                ! lhs
    str(2)%s = "paper"              ! rhs
    str(3)   = str(1) // str(2)     ! result
    str(4)%s = "penpaper"           ! expected
    call tc%assert_eq(str(4), str(3), "concat str//str")

    call tc%finish()
  end subroutine

end module
