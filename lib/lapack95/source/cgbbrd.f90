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

PURE SUBROUTINE CGBBRD_F95(AB,C,D,E,Q,PT,KL,M,INFO)
    ! Fortran77 call:
    ! CGBBRD(VECT,M,N,NCC,KL,KU,AB,LDAB,D,E,Q,LDQ,PT,LDPT,C,LDC,WORK,
    !   RWORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GBBRD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: KL
    INTEGER, INTENT(IN), OPTIONAL :: M
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: AB(:,:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: C(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: D(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: E(:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: Q(:,:)
    COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: PT(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GBBRD'
    ! <<< Local scalars >>>
    INTEGER :: O_KL
    INTEGER :: O_M
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: VECT
    INTEGER :: N
    INTEGER :: NCC
    INTEGER :: KU
    INTEGER :: LDAB
    INTEGER :: LDQ
    INTEGER :: LDPT
    INTEGER :: LDC
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_C(:,:)
    REAL(WP), POINTER :: O_D(:)
    REAL(WP), POINTER :: O_E(:)
    COMPLEX(WP), POINTER :: O_Q(:,:)
    COMPLEX(WP), POINTER :: O_PT(:,:)
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    COMPLEX(WP), TARGET :: L_A2_COMP(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDAB = MAX(1,SIZE(AB,1))
    IF(PRESENT(PT)) THEN
        LDPT = MAX(1,SIZE(PT,1))
    ELSE
        LDPT = 1
    ENDIF
    IF(PRESENT(Q)) THEN
        LDQ = MAX(1,SIZE(Q,1))
    ELSE
        LDQ = 1
    ENDIF
    N = SIZE(AB,2)
    IF(PRESENT(C)) THEN
        NCC = SIZE(C,2)
    ELSE
        NCC = 0
    ENDIF
    IF(PRESENT(Q).AND.PRESENT(PT)) THEN
        VECT = 'B'
    ELSEIF(PRESENT(Q)) THEN
        VECT = 'Q'
    ELSEIF(PRESENT(PT)) THEN
        VECT = 'P'
    ELSE
        VECT = 'N'
    ENDIF
    IF(PRESENT(KL)) THEN
        O_KL = KL
    ELSE
        O_KL = (LDAB-1)/2
    ENDIF
    IF(PRESENT(M)) THEN
        O_M = M
    ELSE
        O_M = N
    ENDIF
    KU = LDAB-O_KL-1
    LDC = MAX(1,O_M)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(C)) THEN
        O_C => C
    ELSE
        O_C => L_A2_COMP
    ENDIF
    IF(PRESENT(D)) THEN
        O_D => D
    ELSE
        ALLOCATE(O_D(MIN(O_M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(E)) THEN
            O_E => E
        ELSE
            ALLOCATE(O_E(MIN(O_M,N)-1), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(PRESENT(PT)) THEN
        O_PT => PT
    ELSE
        O_PT => L_A2_COMP
    ENDIF
    IF(PRESENT(Q)) THEN
        O_Q => Q
    ELSE
        O_Q => L_A2_COMP
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(RWORK(MAX(O_M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(MAX(O_M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GBBRD(VECT,O_M,N,NCC,O_KL,KU,AB,LDAB,O_D,O_E,O_Q,LDQ,  &
     &                              O_PT,LDPT,O_C,LDC,WORK,RWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(D)) THEN
        DEALLOCATE(O_D, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(E)) THEN
        DEALLOCATE(O_E, STAT=L_STAT_DEALLOC)
    ENDIF
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CGBBRD_F95
