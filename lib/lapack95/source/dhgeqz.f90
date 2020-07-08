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

PURE SUBROUTINE DHGEQZ_F95(H,T,ILO,IHI,ALPHAR,ALPHAI,BETA,Q,Z,JOB,COMPQ,&
     &                                                       COMPZ,INFO)
    ! Fortran77 call:
    ! DHGEQZ(JOB,COMPQ,COMPZ,N,ILO,IHI,H,LDH,T,LDT,ALPHAR,ALPHAI,BETA,Q,
    !   LDQ,Z,LDZ,WORK,LWORK,INFO)
    ! Default ILO=1
    ! Default IHI=N
    ! JOB='E','S'; default: 'E'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_HGEQZ, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ILO
    INTEGER, INTENT(IN), OPTIONAL :: IHI
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: COMPQ
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: COMPZ
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: H(:,:)
    REAL(WP), INTENT(INOUT) :: T(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: ALPHAR(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: ALPHAI(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: BETA(:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: Q(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: Z(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'HGEQZ'
    ! <<< Local scalars >>>
    INTEGER :: O_ILO
    INTEGER :: O_IHI
    CHARACTER(LEN=1) :: O_JOB
    CHARACTER(LEN=1) :: O_COMPQ
    CHARACTER(LEN=1) :: O_COMPZ
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDH
    INTEGER :: LDT
    INTEGER :: LDQ
    INTEGER :: LDZ
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_ALPHAR(:)
    REAL(WP), POINTER :: O_ALPHAI(:)
    REAL(WP), POINTER :: O_BETA(:)
    REAL(WP), POINTER :: O_Q(:,:)
    REAL(WP), POINTER :: O_Z(:,:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    REAL(WP) :: S_WORK(1)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(Q).AND.PRESENT(COMPQ)) THEN
        IF((COMPQ.EQ.'I'.OR.COMPQ.EQ.'i').OR.                           &
     &    (COMPQ.EQ.'V'.OR.COMPQ.EQ.'v')) THEN
            O_COMPQ = COMPQ
        ELSE
            O_INFO=-1001; GOTO 1001
        ENDIF
    ELSEIF(PRESENT(Q)) THEN
        O_COMPQ = 'V'
    ELSEIF(PRESENT(COMPQ)) THEN
        O_INFO=-1001; GOTO 1001
    ELSE
        O_COMPQ = 'N'
    ENDIF
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
    IF(PRESENT(Q)) THEN
        LDQ = MAX(1,SIZE(Q,1))
    ELSE
        LDQ = 1
    ENDIF
    LDT = MAX(1,SIZE(T,1))
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
    IF(PRESENT(ALPHAI)) THEN
        O_ALPHAI => ALPHAI
    ELSE
        ALLOCATE(O_ALPHAI(N), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(ALPHAR)) THEN
            O_ALPHAR => ALPHAR
        ELSE
            ALLOCATE(O_ALPHAR(N), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(BETA)) THEN
            O_BETA => BETA
        ELSE
            ALLOCATE(O_BETA(N), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
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
    CALL F77_HGEQZ(O_JOB,O_COMPQ,O_COMPZ,N,O_ILO,O_IHI,H,LDH,T,LDT,     &
     &     O_ALPHAR,O_ALPHAI,O_BETA,O_Q,LDQ,O_Z,LDZ,S_WORK,LWORK,O_INFO)
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
        CALL F77_HGEQZ(O_JOB,O_COMPQ,O_COMPZ,N,O_ILO,O_IHI,H,LDH,T,LDT, &
     &       O_ALPHAR,O_ALPHAI,O_BETA,O_Q,LDQ,O_Z,LDZ,WORK,LWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(ALPHAI)) THEN
        DEALLOCATE(O_ALPHAI, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(ALPHAR)) THEN
        DEALLOCATE(O_ALPHAR, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(BETA)) THEN
        DEALLOCATE(O_BETA, STAT=L_STAT_DEALLOC)
    ENDIF
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DHGEQZ_F95
