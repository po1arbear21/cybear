#include "../macro.f90.inc"

module arnoldi_m
  use blas95
  use error_m
  use matop_m
  use matrix_m
  implicit none

  private
  public :: arnoldi

contains

  subroutine arnoldi(A, b, H, Q, brkd, tol, m0)
    !! Computes a basis of the (n + 1)-Krylov subspace of A: the space spanned by {b, Ab, ..., A^n b}.

    class(matop_real) , intent(in)    :: A
      !! A: m  m array
    real              , intent(in)    :: b(:)
      !! size: m
    type(dense_real)  , intent(inout) :: Q
      !! Q: m x (n + 1) array, the columns are an orthonormal basis of the Krylov subspace.
    type(dense_real)  , intent(inout) :: H      ! fixme change type: dense to hessenberg
      !! h: (n + 1) x n array, A on basis Q. It is upper Hessenberg.
    integer           , intent(out)   :: brkd
      !! breakdown happend at k=brkd. brkd=0 means no breakdown happend.
    real,    optional , intent(in)    :: tol
      !! tolerance for breakdown check
      !!      H%d(k+1,k) < k*tol*nrm2(H%d(1:k+1,k))
      !! default: 1e-12
    integer, optional , intent(in)    :: m0
      !! extend arnoldi decomposition of dimension m0 to m.
      !! Q(:,j) for j=1..m0 is set.
      !! default: 0 (<=> no column of Q is set)

    integer           :: m, m0_, n, k, j
    real              :: tol_
    real, allocatable :: q_tmp(:), v(:)

    m0_ = 0
    if (present(m0)) m0_ = m0
    tol_ = 1e-12
    if (present(tol)) tol_ = tol

    brkd = 0

    m = A%nrows
    n = H%ncols

    ASSERT(size(b) == m)
    ASSERT(Q%nrows == m .and. Q%ncols == n+1)
    ASSERT(H%nrows == n+1)

    call Q%init(m, n+1)
    call H%init(n+1, n)

    q_tmp = b / nrm2(b)
    Q%d(1,:) = q_tmp

    allocate (v(m))
    do k = m0_+1, n
      Q%d(k,:) = q_tmp
      call A%exec(q_tmp, v)                         ! v <- A * qtmp

      do j = 1, k+1
        H%d(j,k) = dot(Q%d(j,:), v)
        call axpy(Q%d(j,:), v, a=-H%d(j,k))         ! v <- v - h[j, k] * Q[:, j]
      end do

      H%d(k+1,k) = nrm2(v)

      if (H%d(k+1,k) < k*tol_*nrm2(H%d(1:k+1,k))) then
        brkd = k
        exit
      end if

      q_tmp = v / H%d(k+1,k)
      Q%d(k+1,:) = q_tmp
    end do
    Q%d = transpose(Q%d)
  end subroutine

end module
