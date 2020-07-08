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

PURE SUBROUTINE DGEBAK_F95(V,SCALE,ILO,IHI,JOB,SIDE,INFO)
    ! Fortran77 call:
    ! DGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
    ! Default ILO=1
    ! Default IHI=N
    ! JOB='B','S','P','N'; default: 'B'
    ! SIDE='L','R'; default: 'L'
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GEBAK, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: ILO
    INTEGER, INTENT(IN), OPTIONAL :: IHI
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    REAL(WP), INTENT(INOUT) :: V(:,:)
    REAL(WP), INTENT(IN) :: SCALE(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GEBAK'
    ! <<< Local scalars >>>
    INTEGER :: O_ILO
    INTEGER :: O_IHI
    CHARACTER(LEN=1) :: O_JOB
    CHARACTER(LEN=1) :: O_SIDE
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: M
    INTEGER :: LDV
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    IF(PRESENT(ILO)) THEN
        O_ILO = ILO
    ELSE
        O_ILO = 1
    ENDIF
    IF(PRESENT(JOB)) THEN
        O_JOB = JOB
    ELSE
        O_JOB = 'B'
    ENDIF
    IF(PRESENT(SIDE)) THEN
        O_SIDE = SIDE
    ELSE
        O_SIDE = 'L'
    ENDIF
    LDV = MAX(1,SIZE(V,1))
    M = SIZE(V,2)
    N = SIZE(V,1)
    IF(PRESENT(IHI)) THEN
        O_IHI = IHI
    ELSE
        O_IHI = N
    ENDIF
    ! <<< Call lapack77 routine >>>
    CALL F77_GEBAK(O_JOB,O_SIDE,N,O_ILO,O_IHI,SCALE,M,V,LDV,O_INFO)
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE DGEBAK_F95
