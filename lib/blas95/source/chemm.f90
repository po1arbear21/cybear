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

PURE SUBROUTINE CHEMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
    ! Fortran77 call:
    ! CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    ! SIDE='L','R'; default: 'L'
    ! UPLO='U','L'; default: 'U'
    ! Default ALPHA=1
    ! Default BETA=0
    ! <<< Use statements >>>
    USE F77_BLAS, ONLY: F77_HEMM
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
    COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: A(:,:)
    COMPLEX(WP), INTENT(IN) :: B(:,:)
    COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=4), PARAMETER :: SRNAME = 'HEMM'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_SIDE
    CHARACTER(LEN=1) :: O_UPLO
    COMPLEX(WP) :: O_ALPHA
    COMPLEX(WP) :: O_BETA
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDC
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ALPHA)) THEN
        O_ALPHA = ALPHA
    ELSE
        O_ALPHA = 1
    ENDIF
    IF(PRESENT(BETA)) THEN
        O_BETA = BETA
    ELSE
        O_BETA = 0
    ENDIF
    IF(PRESENT(SIDE)) THEN
        O_SIDE = SIDE
    ELSE
        O_SIDE = 'L'
    ENDIF
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDC = MAX(1,SIZE(C,1))
    M = SIZE(C,1)
    N = SIZE(C,2)
    ! <<< Call blas77 routine >>>
    CALL F77_HEMM(O_SIDE,O_UPLO,M,N,O_ALPHA,A,LDA,B,LDB,O_BETA,C,LDC)
END SUBROUTINE CHEMM_F95
