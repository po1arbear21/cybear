module test_high_precision_m
  use test_case_m
  use high_precision_m
  implicit none

contains

  subroutine test_high_precision()
    type(test_case) :: tc

    print "(1A)", "test_high_precision"
    call tc%init("high_precision")



    call tc%finish()
  end subroutine test_high_precision

end module test_high_precision_m
