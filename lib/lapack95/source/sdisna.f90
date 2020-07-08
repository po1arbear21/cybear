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

PURE SUBROUTINE SDISNA_F95(D,SEP,JOB,MINMN,INFO)
    ! Fortran77 call:
    ! SDISNA(JOB,M,N,D,SEP,INFO)
    ! JOB='E','L','R'; default: 'E'
    ! MINMN='M','N'; default: 'M';  Superfluous if JOB='E'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_DISNA, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: MINMN
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: D(:)
    REAL(WP), INTENT(OUT) :: SEP(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'DISNA'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_JOB
    CHARACTER(LEN=1) :: O_MINMN
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    ! <<< Intrinsic functions >>>
    INTRINSIC PRESENT
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(JOB)) THEN
        O_JOB = JOB
    ELSE
        O_JOB = 'E'
    ENDIF
    IF(PRESENT(MINMN)) THEN
        O_MINMN = MINMN
    ELSE
        O_MINMN = 'M'
    ENDIF
    IF((O_MINMN.EQ.'M'.OR.O_MINMN.EQ.'m')) THEN
        M = SIZE(D)
    ELSE
        M = SIZE(D)+1
    ENDIF
    IF((O_MINMN.EQ.'M'.OR.O_MINMN.EQ.'m')) THEN
        N = SIZE(D)+1
    ELSE
        N = SIZE(D)
    ENDIF
    ! <<< Call lapack77 routine >>>
    CALL F77_DISNA(O_JOB,M,N,D,SEP,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SDISNA_F95
