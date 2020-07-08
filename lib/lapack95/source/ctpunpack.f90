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

PURE SUBROUTINE MKL_CTPUNPACK_F95(AP,I,J,ROWS,COLS,A,UPLO,TRANS,INFO)
    ! Fortran77 call:
    ! MKL_CTPUNPACK(UPLO,TRANS,N,AP,I,J,ROWS,COLS,A,LDA,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_MKL_TPUNPACK, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(IN) :: J
    INTEGER, INTENT(IN) :: ROWS
    INTEGER, INTENT(IN) :: COLS
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL:: UPLO
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL:: TRANS
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: AP(:,:)
    COMPLEX(WP), INTENT(OUT) :: A(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=12), PARAMETER :: SRNAME = 'MKL_CTPUNPACK'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    CHARACTER(LEN=1) :: O_TRANS
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDA = MAX(1,SIZE(A,1))
    N = SIZE(A,2)
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    ! <<< Call lapack77 routine >>>
    CALL F77_MKL_TPUNPACK(O_UPLO,O_TRANS,N,AP,I,J,ROWS,COLS,A,LDA,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE MKL_CTPUNPACK_F95
