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

PURE SUBROUTINE CSYSWAPR_F95(A,I1,I2,UPLO)
    ! Fortran77 call:
    ! CSYSWAPR(UPLO,N,A,I1,I2)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_SYSWAPR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN) :: I1
    INTEGER, INTENT(IN) :: I2
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'SYSWAPR'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: N
    ! <<< Intrinsic functions >>>
    INTRINSIC PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    N = SIZE(A,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_SYSWAPR(O_UPLO,N,A,I1,I2)
END SUBROUTINE CSYSWAPR_F95
