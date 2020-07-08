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

PURE SUBROUTINE ZTREVC_F95(T,HOWMNY,SELECT,VL,VR,M,INFO)
    ! Fortran77 call:
    ! ZTREVC(SIDE,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,MM,M,WORK,RWORK,
    !   INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TREVC, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: HOWMNY
    INTEGER, INTENT(OUT), OPTIONAL :: M
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: T(:,:)
    ! SELECT: INOUT intent instead of IN because PURE.
    LOGICAL, INTENT(INOUT), OPTIONAL, TARGET :: SELECT(:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: VL(:,:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: VR(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TREVC'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_HOWMNY
    INTEGER :: O_M
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: SIDE
    INTEGER :: N
    INTEGER :: LDT
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: MM
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    LOGICAL, POINTER :: O_SELECT(:)
    COMPLEX(WP), POINTER :: O_VL(:,:)
    COMPLEX(WP), POINTER :: O_VR(:,:)
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    LOGICAL, TARGET :: L_A1_LOGI(1)
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(SELECT).AND.PRESENT(HOWMNY)) THEN
        O_INFO=-1001; GOTO 1001
    ELSEIF(PRESENT(SELECT)) THEN
        O_HOWMNY = 'S'
    ELSEIF(PRESENT(HOWMNY)) THEN
        IF((HOWMNY.EQ.'A'.OR.HOWMNY.EQ.'a').OR.                         &
     &    (HOWMNY.EQ.'B'.OR.HOWMNY.EQ.'b')) THEN
            O_HOWMNY = HOWMNY
        ELSE
            O_INFO=-1001; GOTO 1001
        ENDIF
    ELSE
        O_HOWMNY = 'A'
    ENDIF
    LDT = MAX(1,SIZE(T,1))
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
    IF(PRESENT(VL)) THEN
        MM = SIZE(VL,2)
    ELSE
        MM = SIZE(VR,2)
    ENDIF
    N = SIZE(T,2)
    IF(PRESENT(VL).AND.PRESENT(VR)) THEN
        SIDE = 'B'
    ELSEIF(PRESENT(VL)) THEN
        SIDE = 'L'
    ELSEIF(PRESENT(VR)) THEN
        SIDE = 'R'
    ELSE
        O_INFO=-1001; GOTO 1001
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(SELECT)) THEN
        O_SELECT => SELECT
    ELSE
        O_SELECT => L_A1_LOGI
    ENDIF
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
    ALLOCATE(RWORK(N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(2*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TREVC(SIDE,O_HOWMNY,O_SELECT,N,T,LDT,O_VL,LDVL,O_VR,   &
     &                                    LDVR,MM,O_M,WORK,RWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(M)) THEN
        M = O_M
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZTREVC_F95
