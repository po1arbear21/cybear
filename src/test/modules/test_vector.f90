module test_vector_m
  use test_case_m
  use vector_m
  implicit none

contains

  subroutine test_vector()
    type(test_case) :: tc

    print "(1A)", "test_vector"
    call tc%init("vector")



    call tc%finish()
  end subroutine test_vector

end module test_vector_m
