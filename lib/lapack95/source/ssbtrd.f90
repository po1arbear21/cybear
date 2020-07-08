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

PURE SUBROUTINE SSBTRD_F95(AB,Q,VECT,UPLO,INFO)
    ! Fortran77 call:
    ! SSBTRD(VECT,UPLO,N,KD,AB,LDAB,D,E,Q,LDQ,WORK,INFO)
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_SBTRD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: VECT
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: AB(:,:)
    REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: Q(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'SBTRD'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_VECT
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: KD
    INTEGER :: LDAB
    INTEGER :: LDQ
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_Q(:,:)
    REAL(WP), POINTER :: D(:)
    REAL(WP), POINTER :: E(:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Stubs to "allocate" optional arrays >>>
    REAL(WP), TARGET :: L_A2_REAL(1,1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    IF(PRESENT(Q).AND.PRESENT(VECT)) THEN
        IF((VECT.EQ.'V'.OR.VECT.EQ.'v').OR.                             &
     &    (VECT.EQ.'U'.OR.VECT.EQ.'u')) THEN
            O_VECT = VECT
        ELSE
            O_INFO=-1001; GOTO 1001
        ENDIF
    ELSEIF(PRESENT(Q)) THEN
        O_VECT = 'V'
    ELSEIF(PRESENT(VECT)) THEN
        O_INFO=-1001; GOTO 1001
    ELSE
        O_VECT = 'N'
    ENDIF
    KD = SIZE(AB,1)-1
    LDAB = MAX(1,SIZE(AB,1))
    IF(PRESENT(Q)) THEN
        LDQ = MAX(1,SIZE(Q,1))
    ELSE
        LDQ = 1
    ENDIF
    N = SIZE(AB,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(Q)) THEN
        O_Q => Q
    ELSE
        O_Q => L_A2_REAL
    ENDIF
    ALLOCATE(D(N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(E(N-1), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_SBTRD(O_VECT,O_UPLO,N,KD,AB,LDAB,D,E,O_Q,LDQ,WORK,     &
     &                                                           O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(D, STAT=L_STAT_DEALLOC)
    DEALLOCATE(E, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SSBTRD_F95
