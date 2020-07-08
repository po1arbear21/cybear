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

PURE SUBROUTINE DLAPMR_F95(X,K,FORWRD)
    ! Fortran77 call:
    ! DLAPMR(FORWRD,M,N,X,LDX,K)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_LAPMR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    LOGICAL, INTENT(IN), OPTIONAL :: FORWRD
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: X(:,:)
    INTEGER, INTENT(INOUT) :: K(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'LAPMR'
    ! <<< Local scalars >>>
    LOGICAL :: O_FORWRD
    INTEGER :: LDX
    INTEGER :: N
    INTEGER :: M
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(FORWRD)) THEN
        O_FORWRD = FORWRD
    ELSE
        O_FORWRD = .TRUE.
    ENDIF
    LDX = MAX(1,SIZE(X,1))
    M = SIZE(X,1)
    N = SIZE(X,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_LAPMR(O_FORWRD,M,N,X,LDX,K)
END SUBROUTINE DLAPMR_F95
