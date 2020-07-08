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

PURE SUBROUTINE CGEQRT3_F95(A,T,INFO)
    ! Fortran77 call:
    ! CGEQRT3( M,N,A,LDA,T,LDT,INFO )
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEQRT3, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(OUT) :: T(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=6), PARAMETER :: SRNAME = 'GEQRT3'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: LDT
    INTEGER :: M
    INTEGER :: LDA
    INTEGER :: N
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDA = MAX(1,SIZE(A,1))
    LDT = MAX(1,SIZE(T,1))
    M = SIZE(A,1)
    N = SIZE(A,2)
    ! <<< Call lapack77 routine >>>
    CALL F77_GEQRT3( M,N,A,LDA,T,LDT,INFO )
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CGEQRT3_F95
