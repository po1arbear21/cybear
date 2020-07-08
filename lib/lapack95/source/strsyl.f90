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

PURE SUBROUTINE STRSYL_F95(A,B,C,SCALE,TRANA,TRANB,ISGN,INFO)
    ! Fortran77 call:
    ! STRSYL(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,SCALE,INFO)
    ! TRANA='N','C','T'; default: 'N'
    ! TRANB='N','C','T'; default: 'N'
    ! ISGN=+1,-1; default: +1
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TRSYL, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(OUT) :: SCALE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANA
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANB
    INTEGER, INTENT(IN), OPTIONAL :: ISGN
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: A(:,:)
    REAL(WP), INTENT(IN) :: B(:,:)
    REAL(WP), INTENT(INOUT) :: C(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TRSYL'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_TRANA
    CHARACTER(LEN=1) :: O_TRANB
    INTEGER :: O_ISGN
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDC
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ISGN)) THEN
        O_ISGN = ISGN
    ELSE
        O_ISGN = 1
    ENDIF
    IF(PRESENT(TRANA)) THEN
        O_TRANA = TRANA
    ELSE
        O_TRANA = 'N'
    ENDIF
    IF(PRESENT(TRANB)) THEN
        O_TRANB = TRANB
    ELSE
        O_TRANB = 'N'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDC = MAX(1,SIZE(C,1))
    M = SIZE(A,2)
    N = SIZE(B,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_TRSYL(O_TRANA,O_TRANB,O_ISGN,M,N,A,LDA,B,LDB,C,LDC,SCALE,  &
     &                                                           O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE STRSYL_F95
