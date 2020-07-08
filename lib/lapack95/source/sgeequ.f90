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

PURE SUBROUTINE SGEEQU_F95(A,R,C,ROWCND,COLCND,AMAX,INFO)
    ! Fortran77 call:
    ! SGEEQU(M,N,A,LDA,R,C,ROWCND,COLCND,AMAX,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEEQU, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(OUT), OPTIONAL :: ROWCND
    REAL(WP), INTENT(OUT), OPTIONAL :: COLCND
    REAL(WP), INTENT(OUT), OPTIONAL :: AMAX
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: A(:,:)
    REAL(WP), INTENT(OUT) :: R(:)
    REAL(WP), INTENT(OUT) :: C(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GEEQU'
    ! <<< Local scalars >>>
    REAL(WP) :: O_ROWCND
    REAL(WP) :: O_COLCND
    REAL(WP) :: O_AMAX
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDA = MAX(1,SIZE(A,1))
    M = SIZE(A,1)
    N = SIZE(A,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_GEEQU(M,N,A,LDA,R,C,O_ROWCND,O_COLCND,O_AMAX,O_INFO)
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(AMAX)) THEN
        AMAX = O_AMAX
    ENDIF
    IF(PRESENT(COLCND)) THEN
        COLCND = O_COLCND
    ENDIF
    IF(PRESENT(ROWCND)) THEN
        ROWCND = O_ROWCND
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SGEEQU_F95
