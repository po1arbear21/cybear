module test_feast_m
  use feast_m
  use test_case_m

  implicit none

  private
  public test_feast

contains

  subroutine test_feast()
    type(feast_option) :: opt
    type(test_case) :: tc

    print "(A)", "test_feast"
    call tc%init("feast")


    ! call feast_solve(opt)

    call tc%finish()
  end subroutine

end module
