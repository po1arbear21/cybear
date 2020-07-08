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

PURE SUBROUTINE ZGTTRS1_F95(DL,D,DU,DU2,B,IPIV,TRANS,INFO)
    ! Fortran77 call:
    ! ZGTTRS(TRANS,N,NRHS,DL,D,DU,DU2,IPIV,B,LDB,INFO)
    ! TRANS='N','C','T'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK1, ONLY: F77_GTTRS
    USE F77_LAPACK, ONLY: F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: DL(:)
    COMPLEX(WP), INTENT(IN) :: D(:)
    COMPLEX(WP), INTENT(IN) :: DU(:)
    COMPLEX(WP), INTENT(IN) :: DU2(:)
    COMPLEX(WP), INTENT(INOUT) :: B(:)
    INTEGER, INTENT(IN) :: IPIV(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GTTRS'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_TRANS
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: NRHS
    INTEGER :: LDB
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    LDB = MAX(1,SIZE(B,1))
    N = SIZE(D)
    NRHS = 1
    ! <<< Call lapack77 routine >>>
    CALL F77_GTTRS(O_TRANS,N,NRHS,DL,D,DU,DU2,IPIV,B,LDB,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZGTTRS1_F95
