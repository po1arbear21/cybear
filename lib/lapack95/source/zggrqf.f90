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

PURE SUBROUTINE ZGGRQF_F95(A,B,TAUA,TAUB,INFO)
    ! Fortran77 call:
    ! ZGGRQF(M,P,N,A,LDA,TAUA,B,LDB,TAUB,WORK,LWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GGRQF, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: TAUA(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: TAUB(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GGRQF'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: P
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_TAUA(:)
    COMPLEX(WP), POINTER :: O_TAUB(:)
    COMPLEX(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    M = SIZE(A,1)
    N = SIZE(A,2)
    P = SIZE(B,1)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(TAUA)) THEN
        O_TAUA => TAUA
    ELSE
        ALLOCATE(O_TAUA(MIN(M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(TAUB)) THEN
            O_TAUB => TAUB
        ELSE
            ALLOCATE(O_TAUB(MIN(P,N)), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_GGRQF(M,P,N,A,LDA,O_TAUA,B,LDB,O_TAUB,S_WORK,LWORK,O_INFO)
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
        CALL F77_GGRQF(M,P,N,A,LDA,O_TAUA,B,LDB,O_TAUB,WORK,LWORK,      &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(TAUA)) THEN
        DEALLOCATE(O_TAUA, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(TAUB)) THEN
        DEALLOCATE(O_TAUB, STAT=L_STAT_DEALLOC)
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZGGRQF_F95
