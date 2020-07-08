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

PURE SUBROUTINE SSYCONV_F95(A,IPIV,E,UPLO,WAY,INFO)
    ! Fortran77 call:
    ! SSYCONV(UPLO,WAY,N,A,LDA,IPIV,E,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_SYCONV, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: WAY
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: A(:,:)
    INTEGER, INTENT(IN) :: IPIV(:)
    REAL(WP), INTENT(OUT) :: E(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=6), PARAMETER :: SRNAME = 'SYCONV'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    CHARACTER(LEN=1) :: O_WAY
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    IF(PRESENT(WAY)) THEN
        O_WAY = WAY
    ELSE
        O_WAY = 'C'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    N = SIZE(A,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_SYCONV(O_UPLO,O_WAY,N,A,LDA,IPIV,E,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SSYCONV_F95
