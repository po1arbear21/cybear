module test_normalization_m
  use test_case_m
  use normalization_m
  implicit none

contains

  subroutine test_normalization()
    type(test_case) :: tc

    print "(1A)", "test_normalization"
    call tc%init("normalization")



    call tc%finish()
  end subroutine test_normalization

end module test_normalization_m
