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

PURE SUBROUTINE CPBTRS1_F95(AB,B,UPLO,INFO)
    ! Fortran77 call:
    ! CPBTRS(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK1, ONLY: F77_PBTRS
    USE F77_LAPACK, ONLY: F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: AB(:,:)
    COMPLEX(WP), INTENT(INOUT) :: B(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'PBTRS'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: KD
    INTEGER :: NRHS
    INTEGER :: LDAB
    INTEGER :: LDB
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    KD = SIZE(AB,1)-1
    LDAB = MAX(1,SIZE(AB,1))
    LDB = MAX(1,SIZE(B,1))
    N = SIZE(AB,2)
    NRHS = 1
    ! <<< Call lapack77 routine >>>
    CALL F77_PBTRS(O_UPLO,N,KD,NRHS,AB,LDAB,B,LDB,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CPBTRS1_F95
