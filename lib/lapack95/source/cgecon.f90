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

PURE SUBROUTINE CGECON_F95(A,ANORM,RCOND,NORM,INFO)
    ! Fortran77 call:
    ! CGECON(NORM,N,A,LDA,ANORM,RCOND,WORK,RWORK,INFO)
    ! NORM='1','O','I'; default: '1'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GECON, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    REAL(WP), INTENT(IN) :: ANORM
    REAL(WP), INTENT(OUT) :: RCOND
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: A(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GECON'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_NORM
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(NORM)) THEN
        O_NORM = NORM
    ELSE
        O_NORM = '1'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    N = SIZE(A,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(RWORK(2*N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(2*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GECON(O_NORM,N,A,LDA,ANORM,RCOND,WORK,RWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CGECON_F95
