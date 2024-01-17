m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_spike},0,-1))

module test_spike_m

  use matrix_m,     only: band_real, band_cmplx, BSOLVER_LAPACK, BSOLVER_SPIKE
  use random_m,     only: random
  use test_case_m,  only: test_case

  implicit none

  private
  public test_spike

contains

  subroutine test_spike()
    integer, parameter :: n = 1000, kl = 11, ku = 7, nrhs = 4
    type(test_case)    :: tc
    type(random)       :: rnd

    call tc%init("spike")

    call rnd%init(int([1234567, 89101112], kind = 8), int([1, 2], kind = 8))

    call test_spike_real()
    call test_spike_cmplx()

    call tc%finish()

  contains

    subroutine test_spike_real()
      integer           :: k
      real, allocatable :: tmp(:), b(:,:), x_lapack(:,:), x_spike(:,:)
      type(band_real)   :: A

      call A%init(n, kl, ku)
      allocate (b(n,nrhs), x_lapack(n,nrhs), x_spike(n,nrhs), tmp(n))

      do k = -ku, 0
        call rnd%next_reals(tmp(1:(n+k)))
        A%d(k,(1-k):n) = tmp(1:(n+k))
      end do
      do k = 1, kl
        call rnd%next_reals(tmp(1:(n-k)))
        A%d(k,1:(n-k)) = tmp(1:(n-k))
      end do
      do k = 1, nrhs
        call rnd%next_reals(b(:,k))
      end do

      A%solver = BSOLVER_LAPACK
      call A%factorize()
      call A%solve_mat(b, x_lapack)

      A%solver = BSOLVER_SPIKE
      call A%factorize()
      call A%solve_mat(b, x_spike)

      call tc%assert_eq(x_lapack, x_spike, 1e-10, 1e-16, "real")
    end subroutine

    subroutine test_spike_cmplx()
      integer              :: k
      real,    allocatable :: tmp(:)
      complex, allocatable :: b(:,:), x_lapack(:,:), x_spike(:,:)
      type(band_cmplx)     :: A

      call A%init(n, kl, ku)
      allocate (b(n,nrhs), x_lapack(n,nrhs), x_spike(n,nrhs), tmp(2 * n))
      do k = -ku, 0
        call rnd%next_reals(tmp(1:2*(n+k)))
        A%d(k,(1-k):n) = cmplx(tmp(1:(n+k)), tmp((n+k+1):2*(n+k)))
      end do
      do k = 1, kl
        call rnd%next_reals(tmp(1:2*(n-k)))
        A%d(k,1:n-k) = cmplx(tmp(1:(n-k)), tmp((n-k+1):2*(n-k)))
      end do
      do k = 1, nrhs
        call rnd%next_reals(tmp(1:2*n))
        b(:,k) = cmplx(tmp(1:n), tmp((n+1):2*n))
      end do

      A%solver = BSOLVER_LAPACK
      call A%factorize()
      call A%solve_mat(b, x_lapack)

      A%solver = BSOLVER_SPIKE
      call A%factorize()
      call A%solve_mat(b, x_spike)

      call tc%assert_eq(x_lapack, x_spike, 1e-10, 1e-16, "cmplx")
    end subroutine

  end subroutine

end module

m4_divert(0)
