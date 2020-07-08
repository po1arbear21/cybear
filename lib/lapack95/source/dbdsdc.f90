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

PURE SUBROUTINE DBDSDC_F95(D,E,U,VT,Q,IQ,UPLO,INFO)
    ! Fortran77 call:
    ! DBDSDC(UPLO,COMPQ,N,D,E,U,LDU,VT,LDVT,Q,IQ,WORK,IWORK,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_BDSDC, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: D(:)
    REAL(WP), INTENT(INOUT) :: E(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: VT(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Q(:)
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IQ(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'BDSDC'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: COMPQ
    INTEGER :: N
    INTEGER :: LDU
    INTEGER :: LDVT
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_U(:,:)
    REAL(WP), POINTER :: O_VT(:,:)
    REAL(WP), POINTER :: O_Q(:)
    INTEGER, POINTER :: O_IQ(:)
    REAL(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    INTEGER, TARGET :: L_A1_INTE(1)
    REAL(WP), TARGET :: L_A1_REAL(1)
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
    IF((PRESENT(U).OR.PRESENT(VT)).AND.(PRESENT(Q).OR.PRESENT(IQ))) THEN
        O_INFO=-1001; GOTO 1001
    ELSEIF((PRESENT(U).AND.PRESENT(VT))) THEN
        COMPQ = 'I'
    ELSEIF((PRESENT(Q).AND.PRESENT(IQ))) THEN
        COMPQ = 'P'
    ELSEIF(PRESENT(U).OR.PRESENT(VT).OR.PRESENT(Q).OR.PRESENT(IQ)) THEN
        O_INFO=-1001; GOTO 1001
    ELSE
        COMPQ = 'N'
    ENDIF
    IF(PRESENT(U)) THEN
        LDU = MAX(1,SIZE(U,1))
    ELSE
        LDU = 1
    ENDIF
    IF(PRESENT(VT)) THEN
        LDVT = MAX(1,SIZE(VT,1))
    ELSE
        LDVT = 1
    ENDIF
    N = SIZE(D)
    IF((COMPQ.EQ.'I'.OR.COMPQ.EQ.'i')) THEN
        LWORK = MAX(1,3*N*N+4*N)
    ELSEIF((COMPQ.EQ.'P'.OR.COMPQ.EQ.'p')) THEN
        LWORK = MAX(1,6*N)
    ELSE
        LWORK = MAX(1,4*N)
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(IQ)) THEN
        O_IQ => IQ
    ELSE
        O_IQ => L_A1_INTE
    ENDIF
    IF(PRESENT(Q)) THEN
        O_Q => Q
    ELSE
        O_Q => L_A1_REAL
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
    ALLOCATE(IWORK(8*N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_BDSDC(O_UPLO,COMPQ,N,D,E,O_U,LDU,O_VT,LDVT,O_Q,O_IQ,   &
     &                                                WORK,IWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DBDSDC_F95
