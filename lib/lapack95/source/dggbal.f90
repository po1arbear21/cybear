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

PURE SUBROUTINE DGGBAL_F95(A,B,ILO,IHI,LSCALE,RSCALE,JOB,INFO)
    ! Fortran77 call:
    ! DGGBAL(JOB,N,A,LDA,B,LDB,ILO,IHI,LSCALE,RSCALE,WORK,INFO)
    ! Default ILO=1
    ! Default IHI=N
    ! JOB='B','S','P','N'; default: 'B'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GGBAL, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT), OPTIONAL :: ILO
    INTEGER, INTENT(OUT), OPTIONAL :: IHI
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(INOUT) :: B(:,:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: LSCALE(:)
    REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: RSCALE(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GGBAL'
    ! <<< Local scalars >>>
    INTEGER :: O_ILO
    INTEGER :: O_IHI
    CHARACTER(LEN=1) :: O_JOB
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: O_LSCALE(:)
    REAL(WP), POINTER :: O_RSCALE(:)
    REAL(WP), POINTER :: WORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(JOB)) THEN
        O_JOB = JOB
    ELSE
        O_JOB = 'B'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    N = SIZE(A,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(LSCALE)) THEN
        O_LSCALE => LSCALE
    ELSE
        ALLOCATE(O_LSCALE(N), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        IF(PRESENT(RSCALE)) THEN
            O_RSCALE => RSCALE
        ELSE
            ALLOCATE(O_RSCALE(N), STAT=L_STAT_ALLOC)
        ENDIF
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(6*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GGBAL(O_JOB,N,A,LDA,B,LDB,O_ILO,O_IHI,O_LSCALE,        &
     &                                             O_RSCALE,WORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(IHI)) THEN
        IHI = O_IHI
    ENDIF
    IF(PRESENT(ILO)) THEN
        ILO = O_ILO
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(LSCALE)) THEN
        DEALLOCATE(O_LSCALE, STAT=L_STAT_DEALLOC)
    ENDIF
    IF(.NOT. PRESENT(RSCALE)) THEN
        DEALLOCATE(O_RSCALE, STAT=L_STAT_DEALLOC)
    ENDIF
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DGGBAL_F95
