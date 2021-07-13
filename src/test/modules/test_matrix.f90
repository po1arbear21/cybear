module test_matrix_m

  use example_matrices_m, only: example_matrix3, example_matrix4, example_matrix5, example_matrix6, matrix1, matrix2
  use matrix_m
  use test_case_m,        only: test_case

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

    module subroutine test_triang()
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
    call test_triang()
  end subroutine

end module
