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

PURE SUBROUTINE SPTSVX1_F95(D,E,B,X,DF,EF,FACT,FERR,BERR,RCOND,INFO)
    ! Fortran77 call:
    ! SPTSVX(FACT,N,NRHS,D,E,DF,EF,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,
    !   INFO)
    ! FACT='N','F'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK1, ONLY: F77_PTSVX
    USE F77_LAPACK, ONLY: F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT
    REAL(WP), INTENT(OUT), OPTIONAL :: FERR
    REAL(WP), INTENT(OUT), OPTIONAL :: BERR
    REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: D(:)
    REAL(WP), INTENT(IN) :: E(:)
    REAL(WP), INTENT(IN) :: B(:)
    REAL(WP), INTENT(OUT) :: X(:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: DF(:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: EF(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'PTSVX'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_FACT
    REAL(WP) :: O_RCOND
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: NRHS
    INTEGER :: LDB
    INTEGER :: LDX
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_DF(:)
    REAL(WP), POINTER :: O_EF(:)
    REAL(WP), POINTER :: WORK(:)
    REAL(WP) :: R_BERR(1)
    REAL(WP) :: R_FERR(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(FACT)) THEN
        O_FACT = FACT
    ELSE
        O_FACT = 'N'
    ENDIF
    LDB = MAX(1,SIZE(B,1))
    LDX = MAX(1,SIZE(X,1))
    N = SIZE(D)
    NRHS = 1
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(DF)) THEN
        O_DF => DF
    ELSE
        ALLOCATE(O_DF(N), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(EF)) THEN
            O_EF => EF
        ELSE
            ALLOCATE(O_EF(N-1), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(2*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_PTSVX(O_FACT,N,NRHS,D,E,O_DF,O_EF,B,LDB,X,LDX,O_RCOND, &
     &                                        R_FERR,R_BERR,WORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalars for reduced arrays >>>
    IF(PRESENT(BERR)) THEN
        BERR = R_BERR(1)
    ENDIF
    IF(PRESENT(FERR)) THEN
        FERR = R_FERR(1)
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(RCOND)) THEN
        RCOND = O_RCOND
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(DF)) THEN
        DEALLOCATE(O_DF, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(EF)) THEN
        DEALLOCATE(O_EF, STAT=L_STAT_DEALLOC)
    ENDIF
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SPTSVX1_F95
