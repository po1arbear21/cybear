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

PURE SUBROUTINE DOPMTR_F95(AP,TAU,C,SIDE,UPLO,TRANS,INFO)
    ! Fortran77 call:
    ! DOPMTR(SIDE,UPLO,TRANS,M,N,AP,TAU,C,LDC,WORK,INFO)
    ! SIDE='L','R'; default: 'L'
    ! UPLO='U','L'; default: 'U'
    ! TRANS='N','T'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_OPMTR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: AP(:)
    REAL(WP), INTENT(IN) :: TAU(:)
    REAL(WP), INTENT(INOUT) :: C(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'OPMTR'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_SIDE
    CHARACTER(LEN=1) :: O_UPLO
    CHARACTER(LEN=1) :: O_TRANS
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDC
    INTEGER :: LWORK
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
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    LDC = MAX(1,SIZE(C,1))
    M = SIZE(C,1)
    N = SIZE(C,2)
    IF((O_SIDE.EQ.'L'.OR.O_SIDE.EQ.'l')) THEN
        LWORK = MAX(1,N)
    ELSE
        LWORK = MAX(1,M)
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_OPMTR(O_SIDE,O_UPLO,O_TRANS,M,N,AP,TAU,C,LDC,WORK,     &
     &                                                           O_INFO)
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
END SUBROUTINE DOPMTR_F95
