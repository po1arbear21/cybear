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

PURE SUBROUTINE ZHEGVD_F95(A,B,W,ITYPE,JOBZ,UPLO,INFO)
    ! Fortran77 call:
    ! ZHEGVD(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,LRWORK,
    !   IWORK,LIWORK,INFO)
    ! ITYPE=1,2,3; default: 1
    ! JOBZ='N','V'; default: 'N'
    ! UPLO='U','L'; default: 'U'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_HEGVD, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ITYPE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    REAL(WP), INTENT(OUT) :: W(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'HEGVD'
    ! <<< Local scalars >>>
    INTEGER :: O_ITYPE
    CHARACTER(LEN=1) :: O_JOBZ
    CHARACTER(LEN=1) :: O_UPLO
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LWORK
    INTEGER :: LRWORK
    INTEGER :: LIWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    INTEGER :: S_IWORK(1)
    REAL(WP) :: S_RWORK(1)
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ITYPE)) THEN
        O_ITYPE = ITYPE
    ELSE
        O_ITYPE = 1
    ENDIF
    IF(PRESENT(JOBZ)) THEN
        O_JOBZ = JOBZ
    ELSE
        O_JOBZ = 'N'
    ENDIF
    IF(PRESENT(UPLO)) THEN
        O_UPLO = UPLO
    ELSE
        O_UPLO = 'U'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    N = SIZE(A,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ! <<< Request work array(s) size >>>
    LIWORK = -1
    LRWORK = -1
    LWORK = -1
    CALL F77_HEGVD(O_ITYPE,O_JOBZ,O_UPLO,N,A,LDA,B,LDB,W,S_WORK,LWORK,  &
     &                             S_RWORK,LRWORK,S_IWORK,LIWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LIWORK = S_IWORK(1)
    LRWORK = S_RWORK(1)
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    ALLOCATE(IWORK(LIWORK), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(RWORK(LRWORK), STAT=L_STAT_ALLOC)
    ENDIF
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_HEGVD(O_ITYPE,O_JOBZ,O_UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,&
     &                                 RWORK,LRWORK,IWORK,LIWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(RWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZHEGVD_F95
