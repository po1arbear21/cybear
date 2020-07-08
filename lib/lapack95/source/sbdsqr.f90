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

PURE SUBROUTINE SBDSQR_F95(D,E,VT,U,C,UPLO,INFO)
    ! Fortran77 call:
    ! SBDSQR(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,LDC,WORK,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_BDSQR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: D(:)
    REAL(WP), INTENT(INOUT) :: E(:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: VT(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: U(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: C(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=6), PARAMETER :: SRNAME = 'RBDSQR'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: NCVT
    INTEGER :: NRU
    INTEGER :: NCC
    INTEGER :: LDVT
    INTEGER :: LDU
    INTEGER :: LDC
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_VT(:,:)
    REAL(WP), POINTER :: O_U(:,:)
    REAL(WP), POINTER :: O_C(:,:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    N = SIZE(D)
    IF(PRESENT(C)) THEN
        NCC = SIZE(C,2)
    ELSE
        NCC = 0
    ENDIF
    IF(PRESENT(VT)) THEN
        NCVT = SIZE(VT,2)
    ELSE
        NCVT = 0
    ENDIF
    IF(PRESENT(U)) THEN
        NRU = SIZE(U,1)
    ELSE
        NRU = 0
    ENDIF
    LDC = MAX(1,N)
    LDU = MAX(1,NRU)
    LDVT = MAX(1,N)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(C)) THEN
        O_C => C
    ELSE
        O_C => L_A2_REAL
    ENDIF
    IF(PRESENT(U)) THEN
        O_U => U
    ELSE
        O_U => L_A2_REAL
    ENDIF
    IF(PRESENT(VT)) THEN
        O_VT => VT
    ELSE
        O_VT => L_A2_REAL
    ENDIF
    ALLOCATE(WORK(4*N), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_BDSQR(O_UPLO,N,NCVT,NRU,NCC,D,E,O_VT,LDVT,O_U,LDU,O_C, &
     &                                                  LDC,WORK,O_INFO)
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
END SUBROUTINE SBDSQR_F95
