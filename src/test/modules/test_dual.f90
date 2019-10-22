module test_dual_m
  use test_case_m
  use dual_m
  implicit none

contains

  subroutine test_dual()
    type(test_case) :: tc

    print "(1A)", "test_dual"
    call tc%init("dual")



    call tc%finish()
  end subroutine test_dual

end module test_dual_m
