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

PURE SUBROUTINE SGEQPF_F95(A,JPVT,TAU,INFO)
    ! Fortran77 call:
    ! SGEQPF(M,N,A,LDA,JPVT,TAU,WORK,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEQPF, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: A(:,:)
    INTEGER, INTENT(INOUT) :: JPVT(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: TAU(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GEQPF'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_TAU(:)
    REAL(WP), POINTER :: WORK(:)
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
    IF(PRESENT(TAU)) THEN
        O_TAU => TAU
    ELSE
        ALLOCATE(O_TAU(MIN(M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(3*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GEQPF(M,N,A,LDA,JPVT,O_TAU,WORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(TAU)) THEN
        DEALLOCATE(O_TAU, STAT=L_STAT_DEALLOC)
    ENDIF
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE SGEQPF_F95
