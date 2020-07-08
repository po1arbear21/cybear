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

PURE SUBROUTINE DSTEBZ_F95(D,E,M,NSPLIT,W,IBLOCK,ISPLIT,ORDER,VL,VU,IL, &
     &                                                   IU,ABSTOL,INFO)
    ! Fortran77 call:
    ! DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,IBLOCK,
    !   ISPLIT,WORK,IWORK,INFO)
    ! ORDER='B','E'; default: 'B'
    ! Default VL=-HUGE(VL)
    ! Default VU=HUGE(VL)
    ! Default IL=1
    ! Default IU=N
    ! Default ABSTOL=0.0_WP
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_STEBZ, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(OUT) :: M
    INTEGER, INTENT(OUT) :: NSPLIT
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: ORDER
    REAL(WP), INTENT(IN), OPTIONAL :: VL
    REAL(WP), INTENT(IN), OPTIONAL :: VU
    INTEGER, INTENT(IN), OPTIONAL :: IL
    INTEGER, INTENT(IN), OPTIONAL :: IU
    REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(IN) :: D(:)
    REAL(WP), INTENT(IN) :: E(:)
    REAL(WP), INTENT(OUT) :: W(:)
    INTEGER, INTENT(OUT) :: IBLOCK(:)
    INTEGER, INTENT(OUT) :: ISPLIT(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'STEBZ'
    ! <<< Local scalars >>>
    CHARACTER(LEN=1) :: O_ORDER
    REAL(WP) :: O_VL
    REAL(WP) :: O_VU
    INTEGER :: O_IL
    INTEGER :: O_IU
    REAL(WP) :: O_ABSTOL
    INTEGER :: O_INFO
    CHARACTER(LEN=1) :: RANGE
    INTEGER :: N
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    REAL(WP), POINTER :: WORK(:)
    INTEGER, POINTER :: IWORK(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC HUGE, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ABSTOL)) THEN
        O_ABSTOL = ABSTOL
    ELSE
        O_ABSTOL = 0.0_WP
    ENDIF
    IF(PRESENT(IL)) THEN
        O_IL = IL
    ELSE
        O_IL = 1
    ENDIF
    IF(PRESENT(ORDER)) THEN
        O_ORDER = ORDER
    ELSE
        O_ORDER = 'B'
    ENDIF
    IF(PRESENT(VL)) THEN
        O_VL = VL
    ELSE
        O_VL = -HUGE(VL)
    ENDIF
    IF(PRESENT(VU)) THEN
        O_VU = VU
    ELSE
        O_VU = HUGE(VL)
    ENDIF
    N = SIZE(D)
    IF((PRESENT(VL).OR.PRESENT(VU)).AND.                                &
     &(PRESENT(IL).OR.PRESENT(IU))) THEN
        O_INFO=-1001; GOTO 1001
    ELSEIF((PRESENT(VL).OR.PRESENT(VU))) THEN
        RANGE = 'V'
    ELSEIF((PRESENT(IL).OR.PRESENT(IU))) THEN
        RANGE = 'I'
    ELSE
        RANGE = 'A'
    ENDIF
    IF(PRESENT(IU)) THEN
        O_IU = IU
    ELSE
        O_IU = N
    ENDIF
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    ALLOCATE(IWORK(3*N), STAT=L_STAT_ALLOC)
    IF(L_STAT_ALLOC==0) THEN
        ALLOCATE(WORK(4*N), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_STEBZ(RANGE,O_ORDER,N,O_VL,O_VU,O_IL,O_IU,O_ABSTOL,D,E,&
     &                       M,NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    DEALLOCATE(IWORK, STAT=L_STAT_DEALLOC)
    DEALLOCATE(WORK, STAT=L_STAT_DEALLOC)
1001    CONTINUE
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DSTEBZ_F95
