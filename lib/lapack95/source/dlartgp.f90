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

PURE SUBROUTINE DLARTGP_F95(F,G,CS,SN,R)
    ! Fortran77 call:
    ! DLARTGP( F,G,CS,SN,R )
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_LARTGP, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(IN) :: F
    REAL(WP), INTENT(IN) :: G
    REAL(WP), INTENT(OUT) :: CS
    REAL(WP), INTENT(OUT) :: SN
    REAL(WP), INTENT(OUT) :: R
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=6), PARAMETER :: SRNAME = 'LARTGP'
    ! <<< Local scalars >>>
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    ! <<< Call lapack77 routine >>>
    CALL F77_LARTGP( F,G,CS,SN,R )
END SUBROUTINE DLARTGP_F95
