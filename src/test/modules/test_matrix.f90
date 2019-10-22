module test_matrix_m
  use test_case_m
  use matrix_m
  implicit none

  private
  public test_matrix

  interface
    module subroutine test_dense()
    end subroutine

    module subroutine test_sparse()
    end subroutine

    module subroutine test_band()
    end subroutine

    module subroutine test_block()
    end subroutine
  end interface

contains

  subroutine test_matrix()
    call test_dense()
    call test_sparse()
  end subroutine test_matrix

end module test_matrix_m
