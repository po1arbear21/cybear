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

PURE SUBROUTINE ZSTEQR_F95(D,E,Z,COMPZ,INFO)
    ! Fortran77 call:
    ! ZSTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_STEQR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: COMPZ
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: D(:)
    REAL(WP), INTENT(INOUT) :: E(:)
    COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: Z(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'STEQR'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_COMPZ
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDZ
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: O_Z(:,:)
    REAL(WP), POINTER :: WORK(:)
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
    IF(PRESENT(Z)) THEN
        LDZ = MAX(1,SIZE(Z,1))
    ELSE
        LDZ = 1
    ENDIF
    N = SIZE(D)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(Z)) THEN
        O_Z => Z
    ELSE
        O_Z => L_A2_COMP
    ENDIF
    ALLOCATE(WORK(MAX(1,2*N-2)), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_STEQR(O_COMPZ,N,D,E,O_Z,LDZ,WORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZSTEQR_F95
