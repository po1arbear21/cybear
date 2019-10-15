module test_matrix_m
  use test_case_m
  use matrix_m
  implicit none

contains

  subroutine test_matrix()
    call test_dense()
    call test_sparse()
    call test_band()
    call test_block()
  end subroutine test_matrix

  subroutine test_dense()
    type(test_case) :: tc

    print "(1A)", "test_dense"
    call tc%init("dense")



    call tc%finish()
  end subroutine test_dense

  subroutine test_sparse()
    type(test_case) :: tc

    print "(1A)", "test_sparse"
    call tc%init("sparse")



    call tc%finish()
  end subroutine test_sparse

  subroutine test_band()
    type(test_case) :: tc

    print "(1A)", "test_band"
    call tc%init("band")



    call tc%finish()
  end subroutine test_band

  subroutine test_block()
    type(test_case) :: tc

    print "(1A)", "test_block"
    call tc%init("block")



    call tc%finish()
  end subroutine test_block

end module test_matrix_m
