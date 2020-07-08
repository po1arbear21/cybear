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

PURE SUBROUTINE DTPMQRT_F95(V,T,A,B,K,NB,TRANS,SIDE,INFO)
    ! Fortran77 call:
    ! DTPMQRT( SIDE,TRANS,M,N,K,L,NB,V,LDV,T,LDT,A,LDA,B,LDB,WORK,INFO )
    ! SIDE='l','r' ; default: 'l'
    ! TRANS='n','c' ; default: 'n'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TPMQRT, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN) :: K
    INTEGER, INTENT(IN) :: NB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: V(:,:)
    REAL(WP), INTENT(IN) :: T(:,:)
    REAL(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(INOUT) :: B(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=6), PARAMETER :: SRNAME = 'TPMQRT'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_TRANS
    CHARACTER(LEN=1) :: O_SIDE
    INTEGER :: O_INFO
    INTEGER :: LDA
    INTEGER :: LDV
    INTEGER :: LDB
    INTEGER :: LDT
    INTEGER :: N
    INTEGER :: L
    INTEGER :: M
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(SIDE)) THEN
        O_SIDE = SIDE
    ELSE
        O_SIDE = 'L'
    ENDIF
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    L = SIZE(V,2)
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDT = MAX(1,SIZE(T,1))
    LDV = MAX(1,SIZE(V,1))
    M = SIZE(B,1)
    N = SIZE(B,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(WORK(MAX(M,N)*NB), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TPMQRT( SIDE,O_TRANS,M,N,K,L,NB,V,LDV,T,LDT,A,LDA,B,LDB,WORK,INFO )
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DTPMQRT_F95
