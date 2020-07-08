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

PURE SUBROUTINE STPCON_F95(AP,RCOND,UPLO,DIAG,NORM,INFO)
    ! Fortran77 call:
    ! STPCON(NORM,UPLO,DIAG,N,AP,RCOND,WORK,IWORK,INFO)
    ! NORM='1','O','I'; default: '1'
    ! UPLO='U','L'; default: 'U'
    ! DIAG='N','U'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TPCON, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(OUT) :: RCOND
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: AP(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TPCON'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    CHARACTER(LEN=1) :: O_DIAG
    CHARACTER(LEN=1) :: O_NORM
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    INTEGER :: L_NN
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC INT, PRESENT, REAL, SIZE, SQRT
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(DIAG)) THEN
        O_DIAG = DIAG
    ELSE
        O_DIAG = 'N'
    ENDIF
    IF(PRESENT(NORM)) THEN
        O_NORM = NORM
    ELSE
        O_NORM = '1'
    ENDIF
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    L_NN = SIZE(AP)
    IF(L_NN <= 0) THEN
        N = L_NN
    ELSE
        ! Packed matrix "AP(N*(N+1)/2)", so: N=(-1+8*(SIZE(AP)))/2
        N = INT((-1+SQRT(1+8*REAL(L_NN,WP)))*0.5)
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(IWORK(N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(3*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TPCON(O_NORM,O_UPLO,O_DIAG,N,AP,RCOND,WORK,IWORK,      &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE STPCON_F95
