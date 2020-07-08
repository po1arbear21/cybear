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

PURE SUBROUTINE CPPEQU_F95(AP,S,SCOND,AMAX,UPLO,INFO)
    ! Fortran77 call:
    ! CPPEQU(UPLO,N,AP,S,SCOND,AMAX,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_PPEQU, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(OUT), OPTIONAL :: SCOND
    REAL(WP), INTENT(OUT), OPTIONAL :: AMAX
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: AP(:)
    REAL(WP), INTENT(OUT) :: S(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'PPEQU'
    ! <<< Local scalars >>>
    REAL(WP) :: O_SCOND
    REAL(WP) :: O_AMAX
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: L_NN
    ! <<< Intrinsic functions >>>
    INTRINSIC INT, PRESENT, REAL, SIZE, SQRT
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    L_NN = SIZE(AP)
    IF(L_NN <= 0) THEN
        N = L_NN
    ELSE
        ! Packed matrix "AP(N*(N+1)/2)", so: N=(-1+8*(SIZE(AP)))/2
        N = INT((-1+SQRT(1+8*REAL(L_NN,WP)))*0.5)
    ENDIF
    ! <<< Call lapack77 routine >>>
    CALL F77_PPEQU(O_UPLO,N,AP,S,O_SCOND,O_AMAX,O_INFO)
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(AMAX)) THEN
        AMAX = O_AMAX
    ENDIF
    IF(PRESENT(SCOND)) THEN
        SCOND = O_SCOND
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CPPEQU_F95
