!===============================================================================
! Copyright 2015-2020 Intel Corporation.
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

PURE FUNCTION SCABS1_F95(Z)
    ! Fortran77 call:
    ! SCABS(Z)
    ! <<< Use statements >>>
    USE F77_BLAS, ONLY: F77_CABS1
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    REAL(WP) :: SCABS1_F95
    ! <<< Arguments >>>
    COMPLEX(WP), INTENT(IN) :: Z
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'CABS1'
    ! <<< Executable statements >>>
    ! <<< Call blas77 routine >>>
    SCABS1_F95 = F77_CABS1(Z)
END FUNCTION SCABS1_F95
