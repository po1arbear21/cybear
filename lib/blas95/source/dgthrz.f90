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
!      F95 interface for BLAS routines
!*******************************************************************************
! This file was generated automatically!
!*******************************************************************************

PURE SUBROUTINE DGTHRZ_F95(X,INDX,Y)
    ! Fortran77 call:
    ! DGTHRZ(NZ,Y,X,INDX)
    ! <<< Use statements >>>
    USE F77_BLAS, ONLY: F77_GTHRZ
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Array arguments >>>
    REAL(WP), INTENT(OUT) :: X(:)
    INTEGER, INTENT(IN) :: INDX(:)
    REAL(WP), INTENT(INOUT) :: Y(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GTHRZ'
    ! <<< Local scalars >>>
    INTEGER :: NZ
    ! <<< Intrinsic functions >>>
    INTRINSIC SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    NZ = SIZE(X)
    ! <<< Call blas77 routine >>>
    CALL F77_GTHRZ(NZ,Y,X,INDX)
END SUBROUTINE DGTHRZ_F95
