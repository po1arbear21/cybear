module test_matrix_m
  use test_case_m
  use matrix_m
  implicit none

contains

  subroutine test_matrix()
    type(test_case) :: tc

    print "(1A)", "test_matrix"
    call tc%init("matrix")



    call tc%finish()
  end subroutine test_matrix

end module test_matrix_m
