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

PURE SUBROUTINE CGEEV_F95(A,W,VL,VR,INFO)
    ! Fortran77 call:
    ! CGEEV(JOBVL,JOBVR,N,A,LDA,W,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEEV, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(OUT) :: W(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VR(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=4), PARAMETER :: SRNAME = 'GEEV'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: JOBVL
    CHARACTER(LEN=1) :: JOBVR
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_VL(:,:)
    COMPLEX(WP), POINTER :: O_VR(:,:)
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(VL)) THEN
        JOBVL = 'V'
    ELSE
        JOBVL = 'N'
    ENDIF
    IF(PRESENT(VR)) THEN
        JOBVR = 'V'
    ELSE
        JOBVR = 'N'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    IF(PRESENT(VL)) THEN
        LDVL = MAX(1,SIZE(VL,1))
    ELSE
        LDVL = 1
    ENDIF
    IF(PRESENT(VR)) THEN
        LDVR = MAX(1,SIZE(VR,1))
    ELSE
        LDVR = 1
    ENDIF
    N = SIZE(A,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(VL)) THEN
        O_VL => VL
    ELSE
        O_VL => L_A2_COMP
    ENDIF
    IF(PRESENT(VR)) THEN
        O_VR => VR
    ELSE
        O_VR => L_A2_COMP
    ENDIF
    ALLOCATE(RWORK(2*N), STAT=L_STAT_ALLOC)
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_GEEV(JOBVL,JOBVR,N,A,LDA,W,O_VL,LDVL,O_VR,LDVR,S_WORK,     &
     &                                               LWORK,RWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GEEV(JOBVL,JOBVR,N,A,LDA,W,O_VL,LDVL,O_VR,LDVR,WORK,   &
     &                                               LWORK,RWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CGEEV_F95
