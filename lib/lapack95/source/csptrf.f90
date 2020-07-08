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

PURE SUBROUTINE CSPTRF_F95(AP,UPLO,IPIV,INFO)
    ! Fortran77 call:
    ! CSPTRF(UPLO,N,AP,IPIV,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_SPTRF, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: AP(:)
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'SPTRF'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    INTEGER :: L_NN
    ! <<< Local arrays >>>
    INTEGER, POINTER :: O_IPIV(:)
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
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(IPIV)) THEN
        O_IPIV => IPIV
    ELSE
        ALLOCATE(O_IPIV(N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_SPTRF(O_UPLO,N,AP,O_IPIV,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(IPIV)) THEN
        DEALLOCATE(O_IPIV, STAT=L_STAT_DEALLOC)
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CSPTRF_F95
