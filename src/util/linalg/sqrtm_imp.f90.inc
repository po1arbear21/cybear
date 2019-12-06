#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

#define DENSE_TYPE           CONCAT(dense,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)

#define SQRTM_DENSE          CONCAT(SQRTM_dense,T)
#define SQRTM_HESSENBERG     CONCAT(SQRTM_hessenberg,T)

! subroutine SQRTM_DENSE(A, V, U)
!   !! Computes SQRTM decomposition for dense matrices
!   !!
!   !!  A = V U V^H

!   type(DENSE_TYPE), intent(in)  :: A
!   type(DENSE_TYPE), intent(out) :: V
!     !! unitary matrix: V^-1 = V^H. for real: V^-1 = V^T.
!   type(DENSE_TYPE), intent(out) :: U     ! fixme change to triangular
!     !! upper tridiag matrix

!   integer               :: n, info
!   type(HESSENBERG_TYPE) :: H
!   TT, allocatable       :: tau(:)
!   real, allocatable     :: scale_(:)

!   ASSERT(A%nrows == A%ncols)

!   n = A%nrows
!   call H%init(n, upper=.true.)
!   H%d = A%d

!   ! balance A: A2 = PAP^T
!   allocate(scale_(n))
!   call gebal(H%d, scale=scale_, job='P', info=info)   ! only permutation allowed
!   ASSERT(info==0)

!   ! reduce A to upper hessenberg form: A = V H V^H
!   allocate(tau(n-1))
!   call gehrd(H%d, tau=tau, info=info)
!   ASSERT(info==0)

!   ! form V explicitly
!   call V%init(n)
!   V%d = H%d
! #ifdef TCMPLX
!   call unghr(V%d, tau=tau, info=info)
! #else
!   call orghr(V%d, tau=tau, info=info)
! #endif
!   ASSERT(info==0)

!   ! SQRTM decomp of Hessenberg: H = Z*U*Z^H   =>  A = (QZ)*U*(QZ)^H
!   call SQRTM(H, V, U, nonHessenberg=.true.)

!   ! transform EVecs
!   call gebak(V%d, scale_, job='P', side='R', info=info)
!   ASSERT(info==0)
! end subroutine



subroutine SQRTM_HESSENBERG(H, sqrtH)
  !! Computes the matrix square root of a Hessenberg matrix H.
  !!
  !!  sqrtH = sqrtm(H)
  !!
  !!  H Schur decomposition: H = Q U Q^H

  type(HESSENBERG_TYPE), intent(in)  :: H
    !! upper triangular matrix     ! fixme implement for lower triangular
  type(HESSENBERG_TYPE), intent(out) :: sqrtH
    !! upper triangular. sqrtH = sqrtm(H)

  type(TRIANG_TYPE)     :: U
  type(HESSENBERG_TYPE) :: U_quasi
  type(DENSE_TYPE)      :: Q
  integer :: n

  ASSERT(H%upper)

  n = H%nrows
  call U_quasi%init(n, upper=.true.)

  call schur(H, Q, U_quasi)

  ! is quasi-triangular U_quasi a triangular U ? aka only real EVals?
  ASSERT(U_quasi%is_triangular())
  call U%init(U_quasi)

  ...

end subroutine

subroutine SQRTM_TRIANG(A, sqrtA)
  !! Computes the matrix square root of a triangular matrix A.
  !!
  !!  sqrtA = sqrtm(A)
  !!
  !! Eq.sys.:
  !!    U^2 = T
  !!
  !!    =>    U_ii^2 = T_ii
  !!      and U_ii U_ij + U_ij U_jj = T_ij - sum_{k=i+1}^{j-1} U_ik U_kj

  type(TRIANG_TYPE), intent(in)  :: A
    !! upper triangular matrix     ! fixme implement for lower triangular
  type(TRIANG_TYPE), intent(out) :: sqrtA
    !! upper triangular. sqrtH = sqrtm(H)

  integer :: n, info, i, j, k
  TT      :: s

  ASSERT(A%upper)     ! fixme implement for lower triangular

  n = A%nrows

  call sqrtA%init(n, upper=A%upper)

  do j = 1, n
    ASSERT(sqrt(A%d(j,j)) > 0.0 )
    sqrtA%d(j,j) = sqrt(A%d(j,j))
    do i = j-1, 1, -1
      s = 0.0
      do k = i+1, j-1
        s = s + A%d(i,k)*A%d(k,j)
      end do
      sqrtA%d(i,j) = (A%d(i,j) - s) / (sqrtA%d(i,i) + sqrtA%d(j,j))
    end do
  end do
end subroutine

! subroutine SQRTM_DENSE(H, sqrtH)
!   !! Computes the matrix square root of a dense matrix H.
!   !!
!   !!  sqrtH = sqrtm(H)
!   !!
!   !! Eq.sys.:
!   !!    U^2 = T
!   !!
!   !!    =>    U_ii^2 = T_ii
!   !!      and U_ii U_ij + U_ij U_jj = T_ij - sum_{k=i+1}^{j-1} U_ik U_kj

!   type(DENSE_TYPE), intent(in)  :: H
!     !! upper Hessenberg     ! fixme implement for lower Hessenberg
!   type(DENSE_TYPE),      intent(out) :: sqrtH
!     !! upper triangular. sqrtH = sqrtm(H)

!   integer :: n, info, j

!   ASSERT(H%upper)     ! fixme implement for lower Hessenberg

!   n = H%nrows



!   ! is Hessenberg triangular?
!   do j = 1, n-1
!     ASSERT(abs(H%d(j+1,j)) < 5*epsilon(1.0))
!   end do

!   call sqrtH%init(n)
!   sqrtH%d = 0.0
!   do j = 1, n
!     sqrtH%d(j,j) = sqrt(H%d(j,j))
!     do i = j-1, 1, -1
!       s = 0.0
!       do k = i+1, j-1
!         s = s + H%d(i,k)*H%d(k,j)
!       end do
!       sqrtH%d(i,j) = (H%d(i,j) - s) / (sqrtH%d(i,i) + sqrtH%d(j,j))
!     end do
!   end do
! end subroutine

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef DENSE_TYPE
#undef HESSENBERG_TYPE

#undef SQRTM_DENSE
#undef SQRTM_HESSENBERG
