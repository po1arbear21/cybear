#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DENSE_TYPE           CONCAT(dense,T)

#define EXPM                 CONCAT3(expm,M,T)

subroutine EXPM(A, E)
  !! expm by scaling and squaring. factors out the trace of A -> should be better for scaling.

  type(DENSE_TYPE), intent(in)  :: A
    !! input matrix
  type(DENSE_TYPE), intent(out) :: E
    !! output: E=exp(A)

  integer, parameter :: q = 6
  integer            :: n, s, i, k
  logical            :: pm
  real               :: nrm, fact, tr
  type(DENSE_TYPE)   :: A_, Ak, Ak_, B, C, E_

  ! matrix dimension
  n = A%nrows
  ASSERT(n == A%ncols)

  ! trace
  tr = 0.0
  do i = 1, n
    tr = tr + A%d(i,i)
  end do

  ! A <- A - tr(A)/n*I
  call A_%init(A%d)
  do i = 1, n
    A_%d(i,i) = A%d(i,i) - tr/n
  end do

  ! compute the L-infinity norm
  nrm = 0
  do i = 1, n
    nrm = max(nrm, sum(abs(A_%d(i,1:n))))
  end do

  ! determine scaling factor
  s = max(0, int(log(nrm) / log(2.0)) + 2)

  ! scale matrix by 1 / 2^s
  call A_%scale(1.0/2.0 ** s)

  ! 0th+1st term of series:
  !   B <- I - 0.5 As
  !   C <- I + 0.5 As
  fact = 0.5
  call B%init(A_%d)
  call B%scale(-fact)
  call C%init(A_%d)
  call C%scale( fact)
  do i = 1, n
    B%d(i,i) = 1.0 + B%d(i,i)
    C%d(i,i) = 1.0 + C%d(i,i)
  end do

  ! 2nd to qth term of series
  pm = .true.
  call Ak%init(A_%d)
  do k = 2, q
    fact = (fact * (q - k + 1)) / (k * (2 * q - k + 1))

    ! Ak = As^k
    ! Ak = matmul(As, Ak)
    call A_%mul_dense(Ak, Ak_)
    Ak%d = Ak_%d

    ! B = sum_{k=1}^{q} (-1)^k * fact_k * As^k
    if (pm) then
      call B%add_dense(Ak, fact= fact)
    else
      call B%add_dense(Ak, fact=-fact)
    end if
    pm = .not. pm

    ! C = sum_{k=1}^{q} fact_k * As^k
    call C%add_dense(Ak, fact=fact)
  end do

  ! E = inv(B) * C
  call E%init(n)
  call B%factorize
  call B%solve_mat(C%d, E%d)

  ! exp(A) = exp(A/(2^s)) ^ (2^s)
  call E_%init(n)
  do k = 1, s
    ! E = matmul(E, E)
    call E%mul_mat(E%d, E_%d)
    E%d = E_%d
  end do

  ! exp(A) = exp(tr(A)/n) exp(A-tr(A)/n*I)
  call E%scale(exp(tr/n))
end subroutine

#undef T
#undef TT
#undef M
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3EXPM
#undef CONCAT3

#undef DENSE_TYPE

#undef EXPM
