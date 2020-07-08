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

PURE SUBROUTINE CUNBDB_F95(X11,X12,X21,X22,THETA,PHI,TAUP1,TAUP2,TAUQ1, &
     &                                           TAUQ2,TRANS,SIGNS,INFO)
    ! Fortran77 call:
    ! CUNBDB(TRANS,SIGNS,M,P,Q,X11,LDX11,X12,LDX12,X21,LDX21,X22,LDX22,
    !   THETA,PHI,TAUP1,TAUP2,TAUQ1,TAUQ2,WORK,LWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_UNBDB, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIGNS
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: X11(:,:)
    COMPLEX(WP), INTENT(INOUT) :: X12(:,:)
    COMPLEX(WP), INTENT(INOUT) :: X21(:,:)
    COMPLEX(WP), INTENT(INOUT) :: X22(:,:)
    REAL(WP), INTENT(OUT) :: THETA(:)
    REAL(WP), INTENT(OUT) :: PHI(:)
    COMPLEX(WP), INTENT(OUT) :: TAUP1(:)
    COMPLEX(WP), INTENT(OUT) :: TAUP2(:)
    COMPLEX(WP), INTENT(OUT) :: TAUQ1(:)
    COMPLEX(WP), INTENT(OUT) :: TAUQ2(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'UNBDB'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_TRANS
    CHARACTER(LEN=1) :: O_SIGNS
    INTEGER :: O_INFO
    INTEGER :: LDX11
    INTEGER :: LDX12
    INTEGER :: LDX21
    INTEGER :: LDX22
    INTEGER :: M
    INTEGER :: P
    INTEGER :: Q
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    !    Not inited: 1 scalars (special=1)
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(SIGNS)) THEN
        O_SIGNS = SIGNS
    ELSE
        O_SIGNS = 'O'
    ENDIF
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    LDX11 = MAX(1,SIZE(X11,1))
    LDX12 = MAX(1,SIZE(X12,1))
    LDX21 = MAX(1,SIZE(X21,1))
    LDX22 = MAX(1,SIZE(X22,1))
    M = SIZE(TAUP1) + SIZE(TAUP2)
    P = SIZE(TAUP1)
    Q = SIZE(TAUQ1)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_UNBDB(O_TRANS,O_SIGNS,M,P,Q,X11,LDX11,X12,LDX12,X21,LDX21, &
     &  X22,LDX22,THETA,PHI,TAUP1,TAUP2,TAUQ1,TAUQ2,S_WORK,LWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ! Error while build wrapper: CUNBDB
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_UNBDB(O_TRANS,O_SIGNS,M,P,Q,X11,LDX11,X12,LDX12,X21,   &
     &    LDX21,X22,LDX22,THETA,PHI,TAUP1,TAUP2,TAUQ1,TAUQ2,WORK,LWORK, &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
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
END SUBROUTINE CUNBDB_F95
