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

PURE SUBROUTINE CTGSYL_F95(A,B,C,D,E,F,IJOB,TRANS,SCALE,DIF,INFO)
    ! Fortran77 call:
    ! CTGSYL(TRANS,IJOB,M,N,A,LDA,B,LDB,C,LDC,D,LDD,E,LDE,F,LDF,SCALE,
    !   DIF,WORK,LWORK,IWORK,INFO)
    ! IJOB=0,1,2,3,4; default: 0
    ! TRANS='N','T'; default: 'N'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_TGSYL, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0E0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: IJOB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    REAL(WP), INTENT(OUT), OPTIONAL :: SCALE
    REAL(WP), INTENT(OUT), OPTIONAL :: DIF
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(IN) :: A(:,:)
    COMPLEX(WP), INTENT(IN) :: B(:,:)
    COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    COMPLEX(WP), INTENT(IN) :: D(:,:)
    COMPLEX(WP), INTENT(IN) :: E(:,:)
    COMPLEX(WP), INTENT(INOUT) :: F(:,:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'TGSYL'
    ! <<< Local scalars >>>
    INTEGER :: O_IJOB
    CHARACTER(LEN=1) :: O_TRANS
    REAL(WP) :: O_SCALE
    REAL(WP) :: O_DIF
    INTEGER :: O_INFO
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LDA
    INTEGER :: LDB
    INTEGER :: LDC
    INTEGER :: LDD
    INTEGER :: LDE
    INTEGER :: LDF
    INTEGER :: LWORK
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    COMPLEX(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Arrays to request optimal sizes >>>
    COMPLEX(WP) :: S_WORK(1)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(IJOB)) THEN
        O_IJOB = IJOB
    ELSE
        O_IJOB = 0
    ENDIF
    IF(PRESENT(TRANS)) THEN
        O_TRANS = TRANS
    ELSE
        O_TRANS = 'N'
    ENDIF
    LDA = MAX(1,SIZE(A,1))
    LDB = MAX(1,SIZE(B,1))
    LDC = MAX(1,SIZE(C,1))
    LDD = MAX(1,SIZE(D,1))
    LDE = MAX(1,SIZE(E,1))
    LDF = MAX(1,SIZE(F,1))
    M = SIZE(A,2)
    N = SIZE(B,2)
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(IWORK(M+N+2), STAT=L_STAT_ALLOC)
    ! <<< Request work array(s) size >>>
    LWORK = -1
    CALL F77_TGSYL(O_TRANS,O_IJOB,M,N,A,LDA,B,LDB,C,LDC,D,LDD,E,LDE,F,  &
     &                      LDF,O_SCALE,O_DIF,S_WORK,LWORK,IWORK,O_INFO)
    ! <<< Exit if error: bad parameters >>>
    IF(O_INFO /= 0) THEN
        GOTO 200
    ENDIF
    LWORK = S_WORK(1)
    ! <<< Allocate work arrays with requested sizes >>>
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(LWORK), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_TGSYL(O_TRANS,O_IJOB,M,N,A,LDA,B,LDB,C,LDC,D,LDD,E,LDE,&
     &                      F,LDF,O_SCALE,O_DIF,WORK,LWORK,IWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Set output optional scalar parameters >>>
    IF(PRESENT(DIF)) THEN
        DIF = O_DIF
    ENDIF
    IF(PRESENT(SCALE)) THEN
        SCALE = O_SCALE
    ENDIF
    ! <<< Deallocate work arrays with requested sizes >>>
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
200    CONTINUE
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE CTGSYL_F95
