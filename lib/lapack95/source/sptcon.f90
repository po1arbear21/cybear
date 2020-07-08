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

PURE SUBROUTINE SPTCON_F95(D,E,ANORM,RCOND,INFO)
    ! Fortran77 call:
    ! SPTCON(N,D,E,ANORM,RCOND,WORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_PTCON, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(IN) :: ANORM
    REAL(WP), INTENT(OUT) :: RCOND
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: D(:)
    REAL(WP), INTENT(IN) :: E(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'PTCON'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    N = SIZE(D)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(WORK(N), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_PTCON(N,D,E,ANORM,RCOND,WORK,O_INFO)
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
END SUBROUTINE SPTCON_F95
