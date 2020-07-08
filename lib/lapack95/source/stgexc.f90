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

PURE SUBROUTINE STGEXC_F95(A,B,IFST,ILST,Z,Q,INFO)
    ! Fortran77 call:
    ! STGEXC(WANTQ,WANTZ,N,A,LDA,B,LDB,Q,LDQ,Z,LDZ,IFST,ILST,WORK,LWORK,
    !   INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TGEXC, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(INOUT), OPTIONAL :: IFST
    INTEGER, INTENT(INOUT), OPTIONAL :: ILST
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(INOUT) :: B(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: Z(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: Q(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TGEXC'
    ! <<< Local scalars >>>
    INTEGER :: O_IFST
    INTEGER :: O_ILST
    INTEGER :: O_INFO
    LOGICAL :: WANTQ
    LOGICAL :: WANTZ
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDQ
    INTEGER :: LDZ
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_Z(:,:)
    REAL(WP), POINTER :: O_Q(:,:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    REAL(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(IFST)) THEN
        O_IFST = IFST
    ELSE
        O_IFST = 1
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    IF(PRESENT(Q)) THEN
        LDQ = MAX(1,SIZE(Q,1))
    ELSE
        LDQ = 1
    ENDIF
    IF(PRESENT(Z)) THEN
        LDZ = MAX(1,SIZE(Z,1))
    ELSE
        LDZ = 1
    ENDIF
    N = SIZE(A,2)
    IF(PRESENT(Q)) THEN
        WANTQ = .TRUE.
    ELSE
        WANTQ = .FALSE.
    ENDIF
    IF(PRESENT(Z)) THEN
        WANTZ = .TRUE.
    ELSE
        WANTZ = .FALSE.
    ENDIF
    IF(PRESENT(ILST)) THEN
        O_ILST = ILST
    ELSE
        O_ILST = N
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(Q)) THEN
        O_Q => Q
    ELSE
        O_Q => L_A2_REAL
    ENDIF
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_REAL
    ENDIF
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_TGEXC(WANTQ,WANTZ,N,A,LDA,B,LDB,O_Q,LDQ,O_Z,LDZ,O_IFST,    &
     &                                       O_ILST,S_WORK,LWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TGEXC(WANTQ,WANTZ,N,A,LDA,B,LDB,O_Q,LDQ,O_Z,LDZ,O_IFST,&
     &                                         O_ILST,WORK,LWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(IFST)) THEN
        IFST = O_IFST
    ENDIF
    IF(PRESENT(ILST)) THEN
        ILST = O_ILST
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE STGEXC_F95
