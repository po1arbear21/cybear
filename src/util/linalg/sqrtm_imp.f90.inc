#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define MATRIX_TYPE          CONCAT(M,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define TRIANG_TYPE          CONCAT(triang,T)

#define SQRTM_MATRIX_TYPE    CONCAT3(sqrtm,M,T)

subroutine SQRTM_MATRIX_TYPE(M, sqrtM)
  !! Computes the matrix square root of a matrix M.
  !!
  !!  sqrtM = sqrtm(M) as in MATLAB
  !!
  !!  Uses Schur decomposition of M: M = Q V Q^H
  !!
  !!  compute square root of V
  !!    Eq.sys.: U^2 = V
  !!    Solution:        U_ii^2 = V_ii
  !!                 and U_ii U_ij + U_ij U_jj = V_ij - sum_{k=i+1}^{j-1} U_ik U_kj

  type(MATRIX_TYPE), intent(in)  :: M
    !! General matrix M
  type(DENSE_TYPE),  intent(out) :: sqrtM
    !! Square root of M: sqrtM = sqrtm(M)

  integer               :: n, j, i, k
  TT                    :: s
  type(DENSE_TYPE)      :: Q, tmp
  type(HESSENBERG_TYPE) :: V_quasi
  type(TRIANG_TYPE)     :: V, U

#ifdef TCMPLX
  ASSERT(.false.)                                                           ! fixme not implemented for complex matrices ....
#endif

  ! Schur decomposition of M: M = Q V Q^H
  n = M%nrows
  call schur(M, Q, V_quasi)

  ! is quasi-triangular V_quasi a triangular V ? aka only real EVals?
  ASSERT(V_quasi%is_triang())
  call V_quasi%to_triang(V)

  ! Compute square root U of upper triangular V, a column at a time.
  ASSERT(V%upper)
  call U%init(n, .true.)    ! upper triangular
  do j = 1, n
    U%d(j,j) = sqrt(V%d(j,j))                                               ! fixme some assertion that sqrt is well defined? type == real ??
    do i = j-1, 1, -1
      s = 0.0
      do k = i+1, j-1
        s = s + U%d(i,k)*U%d(k,j)
      end do
      U%d(i,j) = (V%d(i,j) - s) / (U%d(i,i) + U%d(j,j))
    end do
  end do

  ! build sqrtM = Q U Q^H
  call tmp%init(n)
  call sqrtM%init(n)

  call U%mul_mat(transpose(Q%d), tmp%d)       ! tmp   <- U Q^T              ! fixme: transpose should be hermitian transpose. only correct for real types
  call Q%mul_mat(tmp%d, sqrtM%d)              ! sqrtM <- Q tmp = Q U Q^T
end subroutine

#undef T
#undef TT
#undef M
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef MATRIX_TYPE
#undef DENSE_TYPE
#undef HESSENBERG_TYPE
#undef TRIANG_TYPE

#undef SQRTM_MATRIX_TYPE
