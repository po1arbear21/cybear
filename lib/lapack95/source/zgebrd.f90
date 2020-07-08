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

PURE SUBROUTINE ZGEBRD_F95(A,D,E,TAUQ,TAUP,INFO)
    ! Fortran77 call:
    ! ZGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEBRD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: D(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: E(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: TAUQ(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: TAUP(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GEBRD'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_D(:)
    REAL(WP), POINTER :: O_E(:)
    COMPLEX(WP), POINTER :: O_TAUQ(:)
    COMPLEX(WP), POINTER :: O_TAUP(:)
    COMPLEX(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDA = MAX(1,SIZE(A,1))
    M = SIZE(A,1)
    N = SIZE(A,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(D)) THEN
        O_D => D
    ELSE
        ALLOCATE(O_D(MIN(M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(E)) THEN
            O_E => E
        ELSE
            ALLOCATE(O_E(MIN(M,N)-1), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(TAUP)) THEN
            O_TAUP => TAUP
        ELSE
            ALLOCATE(O_TAUP(MIN(M,N)), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(TAUQ)) THEN
            O_TAUQ => TAUQ
        ELSE
            ALLOCATE(O_TAUQ(MIN(M,N)), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_GEBRD(M,N,A,LDA,O_D,O_E,O_TAUQ,O_TAUP,S_WORK,LWORK,O_INFO)
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
        CALL F77_GEBRD(M,N,A,LDA,O_D,O_E,O_TAUQ,O_TAUP,WORK,LWORK,      &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(D)) THEN
        DEALLOCATE(O_D, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(E)) THEN
        DEALLOCATE(O_E, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(TAUP)) THEN
        DEALLOCATE(O_TAUP, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(TAUQ)) THEN
        DEALLOCATE(O_TAUQ, STAT=L_STAT_DEALLOC)
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZGEBRD_F95
