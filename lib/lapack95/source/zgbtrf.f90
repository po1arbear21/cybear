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

PURE SUBROUTINE ZGBTRF_F95(AB,KL,M,IPIV,INFO)
    ! Fortran77 call:
    ! ZGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
    ! <<< Use statements >>>
    USE F77_LAPACK, ONLY: F77_GBTRF, F77_XERBLA
    ! <<< Implicit statement >>>
    IMPLICIT NONE
    ! <<< Kind parameter >>>
    INTEGER, PARAMETER :: WP = KIND(1.0D0)
    ! <<< Scalar arguments >>>
    INTEGER, INTENT(IN), OPTIONAL :: KL
    INTEGER, INTENT(IN), OPTIONAL :: M
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    ! <<< Array arguments >>>
    COMPLEX(WP), INTENT(INOUT) :: AB(:,:)
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    ! <<< Local declarations >>>
    ! <<< Parameters >>>
    CHARACTER(LEN=5), PARAMETER :: SRNAME = 'GBTRF'
    ! <<< Local scalars >>>
    INTEGER :: O_KL
    INTEGER :: O_M
    INTEGER :: O_INFO
    INTEGER :: N
    INTEGER :: KU
    INTEGER :: LDAB
    INTEGER :: L_STAT_ALLOC, L_STAT_DEALLOC
    ! <<< Local arrays >>>
    INTEGER, POINTER :: O_IPIV(:)
    ! <<< Intrinsic functions >>>
    INTRINSIC MAX, PRESENT, SIZE
    ! <<< Executable statements >>>
    ! <<< Init optional and skipped scalars >>>
    LDAB = MAX(1,SIZE(AB,1))
    N = SIZE(AB,2)
    IF(PRESENT(KL)) THEN
        O_KL = KL
    ELSE
        O_KL = (LDAB-1)/3
    ENDIF
    IF(PRESENT(M)) THEN
        O_M = M
    ELSE
        O_M = N
    ENDIF
    KU = LDAB-2*O_KL-1
    ! <<< Init allocate status >>>
    L_STAT_ALLOC = 0
    ! <<< Allocate local and work arrays >>>
    IF(PRESENT(IPIV)) THEN
        O_IPIV => IPIV
    ELSE
        ALLOCATE(O_IPIV(MIN(O_M,N)), STAT=L_STAT_ALLOC)
    ENDIF
    ! <<< Call lapack77 routine >>>
    IF(L_STAT_ALLOC==0) THEN
        CALL F77_GBTRF(O_M,N,O_KL,KU,AB,LDAB,O_IPIV,O_INFO)
    ELSE; O_INFO = -1000
    ENDIF
    ! <<< Deallocate local and work arrays >>>
    IF(.NOT. PRESENT(IPIV)) THEN
        DEALLOCATE(O_IPIV, STAT=L_STAT_DEALLOC)
    ENDIF
    ! <<< Error handler >>>
    IF(PRESENT(INFO)) THEN
        INFO = O_INFO
    ELSEIF(O_INFO <= -1000) THEN
        CALL F77_XERBLA(SRNAME,-O_INFO)
    ENDIF
END SUBROUTINE ZGBTRF_F95
