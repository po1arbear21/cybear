!===============================================================================
! Copyright 2005-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!      F95 interface for LAPACK routines
!*******************************************************************************
! This file was generated automatically!
!*******************************************************************************

PURE SUBROUTINE DTPRFB_F95(T,V,A,B,DIRECT,STOREV,SIDE,TRANS)
    ! Fortran77 call:
    ! DTPRFB( SIDE,TRANS,DIRECT,STOREV,M,N,K,L,V,LDV,T,LDT,A,LDA,B,LDB,WORK,LDWORK )
    ! SIDE='l','r' ; default: 'l'
    ! TRANS='n','t' ; default: 'n'
    ! DIRECT='f','b' ; default: 'f'
    ! STOREV='c','r' ; default: 'c'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TPRFB, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIRECT
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: STOREV
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: T(:,:)
    REAL(WP), INTENT(IN) :: V(:,:)
    REAL(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(INOUT) :: B(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TPRFB'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_DIRECT
    CHARACTER(LEN=1) :: O_STOREV
    CHARACTER(LEN=1) :: O_SIDE
    CHARACTER(LEN=1) :: O_TRANS
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDT
    INTEGER :: L
    INTEGER :: LDA
    INTEGER :: LDV
    INTEGER :: K
    INTEGER :: LDB
    INTEGER :: M
    INTEGER :: LDWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    INTEGER :: L_NN
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:,:)
    ! <<< Intrinsic functions >>>
    INTRINSIC INT, MAX, PRESENT, REAL, SIZE, SQRT
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(DIRECT)) THEN
        O_DIRECT = DIRECT
    ELSE
        O_DIRECT = 'F'
    ENDIF
    IF(PRESENT(SIDE)) THEN
        O_SIDE = SIDE
    ELSE
        O_SIDE = 'L'
    ENDIF
    IF(PRESENT(STOREV)) THEN
        O_STOREV = STOREV
    ELSE
        O_STOREV = 'C'
    ENDIF
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    L_NN = SIZE(T)
    IF(L_NN <= 0) THEN
        K = L_NN
    ELSE
        ! Packed matrix "T(K*(K+1)/2)", so: K=(-1+8*(SIZE(T)))/2
        K = INT((-1+SQRT(1+8*REAL(L_NN,WP)))*0.5)
    ENDIF
    L = SIZE(V,2)
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDT = MAX(1,SIZE(T,1))
    LDV = MAX(1,SIZE(V,1))
    M = SIZE(B,1)
    N = SIZE(B,2)
    LDWORK = MAX(K,M)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(WORK(LDWORK,N), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TPRFB( O_SIDE,O_TRANS,O_DIRECT,O_STOREV,M,N,K,L,V,LDV,T,LDT,A,LDA,B,LDB,WORK,LDWORK )
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
END SUBROUTINE DTPRFB_F95
