m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mumps},0,-1))

module test_klu2_m

  use klu2_m,       only: create_klu2_handle, destruct_klu2_handle, klu2_factorize, klu2_solve
  use sparse_idx_m, only: SPARSE_IDX
  use test_case_m,  only: test_case

  implicit none

  private
  public test_klu2

contains

  subroutine test_klu2()
    type(test_case)     :: tc
    integer(SPARSE_IDX) :: ia(6)
    integer             :: ja(12), h
    real(kind=16)       :: a(12), b(5), x(5), x_exp(5)

    call tc%init("klu2")

    ia = [1, 3, 6, 9, 10, 13]
    ja = [1, 2, 1, 3, 5, 2, 3, 4, 3, 2, 3, 5]
    a  = [2, 3, 3, 4, 6, -1, -3, 2, 1, 4, 2, 1]
    b  = [8, 45, -3, 3, 19]
    x_exp = [1, 2, 3, 4, 5]

    h = create_klu2_handle()
    call klu2_factorize(h, ia, ja, a)
    call klu2_solve(h, b, x)
    call destruct_klu2_handle(h)

    print *, x

    call tc%assert_eq(x_exp, x, real(1e-31, kind=16), real(1e-31, kind=16), "klu2 solve")

    call tc%finish()
  end subroutine

end module

m4_divert(0)
