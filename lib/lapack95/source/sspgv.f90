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

PURE SUBROUTINE SSPGV_F95(AP,BP,W,ITYPE,UPLO,Z,INFO)
    ! Fortran77 call:
    ! SSPGV(ITYPE,JOBZ,UPLO,N,AP,BP,W,Z,LDZ,WORK,INFO)
    ! ITYPE=1,2,3; default: 1
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_SPGV, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ITYPE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: AP(:)
    REAL(WP), INTENT(INOUT) :: BP(:)
    REAL(WP), INTENT(OUT) :: W(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=4), PARAMETER :: SRNAME = 'SPGV'
    ! <<< Local scalars >>>
    INTEGER :: O_ITYPE
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: JOBZ
    INTEGER :: N
    INTEGER :: LDZ
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    INTEGER :: L_NN
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_Z(:,:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC INT, MAX, PRESENT, REAL, SIZE, SQRT
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ITYPE)) THEN
        O_ITYPE = ITYPE
    ELSE
        O_ITYPE = 1
    ENDIF
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    IF(PRESENT(Z)) THEN
        JOBZ = 'V'
    ELSE
        JOBZ = 'N'
    ENDIF
    IF(PRESENT(Z)) THEN
        LDZ = MAX(1,SIZE(Z,1))
    ELSE
        LDZ = 1
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
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_REAL
    ENDIF
    ALLOCATE(WORK(3*N), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_SPGV(O_ITYPE,JOBZ,O_UPLO,N,AP,BP,W,O_Z,LDZ,WORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SSPGV_F95
