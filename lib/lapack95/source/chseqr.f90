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

PURE SUBROUTINE CHSEQR_F95(H,W,ILO,IHI,Z,JOB,COMPZ,INFO)
    ! Fortran77 call:
    ! CHSEQR(JOB,COMPZ,N,ILO,IHI,H,LDH,W,Z,LDZ,WORK,LWORK,INFO)
    ! JOB='E','S'; default: 'E'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_HSEQR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ILO
    INTEGER, INTENT(IN), OPTIONAL :: IHI
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: COMPZ
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: H(:,:)
    COMPLEX(WP), INTENT(OUT) :: W(:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: Z(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'HSEQR'
    ! <<< Local scalars >>>
    INTEGER :: O_ILO
    INTEGER :: O_IHI
    CHARACTER(LEN=1) :: O_JOB
    CHARACTER(LEN=1) :: O_COMPZ
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDH
    INTEGER :: LDZ
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_Z(:,:)
    COMPLEX(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(Z).AND.PRESENT(COMPZ)) THEN
        IF((COMPZ.EQ.'I'.OR.COMPZ.EQ.'i').OR.                           &
     &    (COMPZ.EQ.'V'.OR.COMPZ.EQ.'v')) THEN
            O_COMPZ = COMPZ
        ELSE
            O_INFO=-1001; GOTO 1001
        ENDIF
    ELSEIF(PRESENT(Z)) THEN
        O_COMPZ = 'I'
    ELSEIF(PRESENT(COMPZ)) THEN
        O_INFO=-1001; GOTO 1001
    ELSE
        O_COMPZ = 'N'
    ENDIF
    IF(PRESENT(ILO)) THEN
        O_ILO = ILO
    ELSE
        O_ILO = 1
    ENDIF
    IF(PRESENT(JOB)) THEN
        O_JOB = JOB
    ELSE
        O_JOB = 'E'
    ENDIF
    LDH = MAX(1,SIZE(H,1))
    IF(PRESENT(Z)) THEN
        LDZ = MAX(1,SIZE(Z,1))
    ELSE
        LDZ = 1
    ENDIF
    N = SIZE(H,2)
    IF(PRESENT(IHI)) THEN
        O_IHI = IHI
    ELSE
        O_IHI = N
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_COMP
    ENDIF
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_HSEQR(O_JOB,O_COMPZ,N,O_ILO,O_IHI,H,LDH,W,O_Z,LDZ,S_WORK,  &
     &                                                     LWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_HSEQR(O_JOB,O_COMPZ,N,O_ILO,O_IHI,H,LDH,W,O_Z,LDZ,WORK,&
     &                                                     LWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CHSEQR_F95
