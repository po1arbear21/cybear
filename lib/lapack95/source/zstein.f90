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

PURE SUBROUTINE ZSTEIN_F95(D,E,W,IBLOCK,ISPLIT,Z,IFAILV,INFO)
    ! Fortran77 call:
    ! ZSTEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,IFAILV,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_STEIN, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: D(:)
    REAL(WP), INTENT(IN) :: E(:)
    REAL(WP), INTENT(IN) :: W(:)
    INTEGER, INTENT(IN) :: IBLOCK(:)
    INTEGER, INTENT(IN) :: ISPLIT(:)
    COMPLEX(WP), INTENT(OUT) :: Z(:,:)
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAILV(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'STEIN'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: M
    INTEGER :: LDZ
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    INTEGER, POINTER :: O_IFAILV(:)
    REAL(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDZ = MAX(1,SIZE(Z,1))
    M = SIZE(Z,2)
    N = SIZE(D)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(IFAILV)) THEN
        O_IFAILV => IFAILV
    ELSE
        ALLOCATE(O_IFAILV(M), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(IWORK(N), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(5*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_STEIN(N,D,E,M,W,IBLOCK,ISPLIT,Z,LDZ,WORK,IWORK,        &
     &                                                  O_IFAILV,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(IFAILV)) THEN
        DEALLOCATE(O_IFAILV, STAT=L_STAT_DEALLOC)
    ENDIF
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZSTEIN_F95
