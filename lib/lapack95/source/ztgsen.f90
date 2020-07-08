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

PURE SUBROUTINE ZTGSEN_F95(A,B,SELECT,ALPHA,BETA,IJOB,Q,Z,PL,PR,DIF,M,  &
     &                                                             INFO)
    ! Fortran77 call:
    ! ZTGSEN(IJOB,WANTQ,WANTZ,SELECT,N,A,LDA,B,LDB,ALPHA,BETA,Q,LDQ,Z,
    !   LDZ,M,PL,PR,DIF,WORK,LWORK,IWORK,LIWORK,INFO)
    ! IJOB=0,1,2,3,4,5; default: 0
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TGSEN, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: IJOB
    REAL(WP), INTENT(OUT), OPTIONAL :: PL
    REAL(WP), INTENT(OUT), OPTIONAL :: PR
    INTEGER, INTENT(OUT), OPTIONAL :: M
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    LOGICAL, INTENT(IN) :: SELECT(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: ALPHA(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: BETA(:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: Q(:,:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: Z(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: DIF(2)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TGSEN'
    ! <<< Local scalars >>>
    INTEGER :: O_IJOB
    REAL(WP) :: O_PL
    REAL(WP) :: O_PR
    INTEGER :: O_M
    INTEGER :: O_INFO
    LOGICAL :: WANTQ
    LOGICAL :: WANTZ
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDQ
    INTEGER :: LDZ
    INTEGER :: LWORK
    INTEGER :: LIWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_ALPHA(:)
    COMPLEX(WP), POINTER :: O_BETA(:)
    COMPLEX(WP), POINTER :: O_Q(:,:)
    COMPLEX(WP), POINTER :: O_Z(:,:)
    REAL(WP), TARGET :: O_O_DIF(2)
    REAL(WP), POINTER :: O_DIF(:)
    COMPLEX(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    INTEGER :: S_IWORK(1)
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(IJOB)) THEN
        O_IJOB = IJOB
    ELSE
        O_IJOB = 0
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
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(ALPHA)) THEN
        O_ALPHA => ALPHA
    ELSE
        ALLOCATE(O_ALPHA(N), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(BETA)) THEN
            O_BETA => BETA
        ELSE
            ALLOCATE(O_BETA(N), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(PRESENT(DIF)) THEN
        O_DIF => DIF
    ELSE
        O_DIF => O_O_DIF
    ENDIF
    IF(PRESENT(Q)) THEN
        O_Q => Q
    ELSE
        O_Q => L_A2_COMP
    ENDIF
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_COMP
    ENDIF
    ! <<< Request work array(s) size >>>
    LIWORK = -1
    LWORK = -1
    CALL F77_TGSEN(O_IJOB,WANTQ,WANTZ,SELECT,N,A,LDA,B,LDB,O_ALPHA,     &
     & O_BETA,O_Q,LDQ,O_Z,LDZ,O_M,O_PL,O_PR,O_DIF,S_WORK,LWORK,S_IWORK, &
     &                                                    LIWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LIWORK = S_IWORK(1)
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(IWORK(LIWORK), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TGSEN(O_IJOB,WANTQ,WANTZ,SELECT,N,A,LDA,B,LDB,O_ALPHA, &
     &     O_BETA,O_Q,LDQ,O_Z,LDZ,O_M,O_PL,O_PR,O_DIF,WORK,LWORK,IWORK, &
     &                                                    LIWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(M)) THEN
        M = O_M
    ENDIF
    IF(PRESENT(PL)) THEN
        PL = O_PL
    ENDIF
    IF(PRESENT(PR)) THEN
        PR = O_PR
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(ALPHA)) THEN
        DEALLOCATE(O_ALPHA, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(BETA)) THEN
        DEALLOCATE(O_BETA, STAT=L_STAT_DEALLOC)
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZTGSEN_F95
