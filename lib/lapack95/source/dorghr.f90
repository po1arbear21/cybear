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

PURE SUBROUTINE DORGHR_F95(A,TAU,ILO,IHI,INFO)
    ! Fortran77 call:
    ! DORGHR(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
    ! Default ILO=1
    ! Default IHI=N
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_ORGHR, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ILO
    INTEGER, INTENT(IN), OPTIONAL :: IHI
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: A(:,:)
    REAL(WP), INTENT(IN) :: TAU(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'ORGHR'
    ! <<< Local scalars >>>
    INTEGER :: O_ILO
    INTEGER :: O_IHI
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:)
    ! <<< Arrays to request optimal sizes >>>
    REAL(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ILO)) THEN
        O_ILO = ILO
    ELSE
        O_ILO = 1
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    N = SIZE(A,2)
    IF(PRESENT(IHI)) THEN
        O_IHI = IHI
    ELSE
        O_IHI = N
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_ORGHR(N,O_ILO,O_IHI,A,LDA,TAU,S_WORK,LWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_ORGHR(N,O_ILO,O_IHI,A,LDA,TAU,WORK,LWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DORGHR_F95
