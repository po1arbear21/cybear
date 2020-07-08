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

PURE SUBROUTINE ZTPQRT_F95(A,B,T,NB,INFO)
    ! Fortran77 call:
    ! ZTPQRT( M,N,L,NB,A,LDA,B,LDB,T,LDT,WORK,INFO )
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TPQRT, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN) :: NB
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: A(:,:)
    COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    COMPLEX(WP), INTENT(OUT) :: T(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TPQRT'
    ! <<< Local scalars >>>
    INTEGER :: O_INFO
    INTEGER :: L
    INTEGER :: LDB
    INTEGER :: LDA
    INTEGER :: LDT
    INTEGER :: N
    INTEGER :: M
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: WORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    L = SIZE(B,1)
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDT = MAX(1,SIZE(T,1))
    M = SIZE(B,1)
    N = SIZE(B,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(WORK(NB*N), STAT=L_STAT_ALLOC)
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TPQRT( M,N,L,NB,A,LDA,B,LDB,T,LDT,WORK,INFO )
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZTPQRT_F95
