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

PURE FUNCTION ICAMAX_F95(X)
    ! Fortran77 call:
    ! ICAMAX(N,X,INCX)
    ! <<< Use statements >>>
    USE F77_BLAS, ONLY: F77_IAMAX
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    INTEGER :: ICAMAX_F95
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: X(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'IAMAX'
    ! <<< Local scalars >>>
    INTEGER :: INCX
    INTEGER :: N
    ! <<< Intrinsic functions >>>
    INTRINSIC SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    INCX = 1
    N = SIZE(X)
    ! <<< Call blas77 routine >>>
    ICAMAX_F95 = F77_IAMAX(N,X,INCX)
END FUNCTION ICAMAX_F95
