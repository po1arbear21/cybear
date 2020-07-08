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

PURE SUBROUTINE SSTEVD_F95(D,E,Z,INFO)
    ! Fortran77 call:
    ! SSTEVD(JOBZ,N,D,E,Z,LDZ,WORK,LWORK,IWORK,LIWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_STEVD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: D(:)
    REAL(WP), INTENT(INOUT) :: E(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'STEVD'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: JOBZ
    INTEGER :: N
    INTEGER :: LDZ
    INTEGER :: LWORK
    INTEGER :: LIWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_Z(:,:)
    REAL(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    INTEGER :: S_IWORK(1)
    REAL(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
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
    N = SIZE(D)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_REAL
    ENDIF
    ! <<< Request work array(s) size >>>
    LIWORK = -1
    LWORK = -1
    CALL F77_STEVD(JOBZ,N,D,E,O_Z,LDZ,S_WORK,LWORK,S_IWORK,LIWORK,      &
     &                                                           O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LIWORK = S_IWORK(1)
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(IWORK(LIWORK), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_STEVD(JOBZ,N,D,E,O_Z,LDZ,WORK,LWORK,IWORK,LIWORK,      &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SSTEVD_F95
