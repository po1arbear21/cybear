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

PURE SUBROUTINE ZGESDD_F95(A,S,U,VT,JOBZ,INFO)
    ! Fortran77 call:
    ! ZGESDD(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,RWORK,IWORK,INFO)
    ! JOBZ='N','A','S','O'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GESDD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(OUT) :: S(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VT(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GESDD'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_JOBZ
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDU
    INTEGER :: LDVT
    INTEGER :: LWORK
    INTEGER :: LRWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_U(:,:)
    COMPLEX(WP), POINTER :: O_VT(:,:)
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, MIN, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(JOBZ)) THEN
        O_JOBZ = JOBZ
    ELSE
        O_JOBZ = 'N'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    IF(PRESENT(U)) THEN
        LDU = MAX(1,SIZE(U,1))
    ELSE
        LDU = 1
    ENDIF
    IF(PRESENT(VT)) THEN
        LDVT = MAX(1,SIZE(VT,1))
    ELSE
        LDVT = 1
    ENDIF
    M = SIZE(A,1)
    N = SIZE(A,2)
    IF((O_JOBZ.EQ.'N'.OR.O_JOBZ.EQ.'n')) THEN
        LRWORK = MAX(1,5*MIN(M,N))
    ELSE
        LRWORK = MAX(1,5*MIN(M,N)*MIN(M,N)+7*MIN(M,N))
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(U)) THEN
        O_U => U
    ELSE
        O_U => L_A2_COMP
    ENDIF
    IF(PRESENT(VT)) THEN
        O_VT => VT
    ELSE
        O_VT => L_A2_COMP
    ENDIF
    ALLOCATE(IWORK(8*MIN(M,N)), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(RWORK(LRWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_GESDD(O_JOBZ,M,N,A,LDA,S,O_U,LDU,O_VT,LDVT,S_WORK,LWORK,   &
     &                                               RWORK,IWORK,O_INFO)
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
        CALL F77_GESDD(O_JOBZ,M,N,A,LDA,S,O_U,LDU,O_VT,LDVT,WORK,LWORK, &
     &                                               RWORK,IWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZGESDD_F95
