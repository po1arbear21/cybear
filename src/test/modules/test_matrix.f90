module test_matrix_m

  use example_matrices_m
  use matrix_m
  use test_case_m

  implicit none

  private
  public test_matrix

  interface
    module subroutine test_arith()
    end subroutine

    module subroutine test_band()
    end subroutine

    module subroutine test_block()
    end subroutine

    module subroutine test_conv()
    end subroutine

    module subroutine test_dense()
    end subroutine

    module subroutine test_sparse()
    end subroutine
  end interface

contains

  subroutine test_matrix()
    call test_arith()
    call test_band()
    call test_block()
    call test_conv()
    call test_dense()
    call test_sparse()
  end subroutine

end module
