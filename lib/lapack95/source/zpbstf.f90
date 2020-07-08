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

PURE SUBROUTINE ZPBSTF_F95(BB,UPLO,INFO)
    ! Fortran77 call:
    ! ZPBSTF(UPLO,N,KB,BB,LDBB,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_PBSTF, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: BB(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'PBSTF'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: KB
    INTEGER :: LDBB
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    KB = SIZE(BB,1)-1
    LDBB = MAX(1,SIZE(BB,1))
    N = SIZE(BB,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_PBSTF(O_UPLO,N,KB,BB,LDBB,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZPBSTF_F95
