#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

#define DENSE_TYPE           CONCAT(dense,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)

#define SCHUR_DENSE          CONCAT(schur_dense,T)
#define SCHUR_HESSENBERG     CONCAT(schur_hessenberg,T)

subroutine SCHUR_DENSE(A, V, U)
  !! Computes Schur decomposition for dense matrices
  !!
  !!  A = V U V^H

  type(DENSE_TYPE),      intent(in)  :: A
  type(DENSE_TYPE),      intent(out) :: V
    !! unitary matrix: V^-1 = V^H. for real: V^-1 = V^T.
  type(HESSENBERG_TYPE), intent(out) :: U
    !! upper quasi-triangular matrix

  integer               :: n, info
  type(HESSENBERG_TYPE) :: H
  TT, allocatable       :: tau(:)
  real, allocatable     :: scale_(:)

  ASSERT(A%nrows == A%ncols)

  n = A%nrows
  call H%init(n, upper=.true.)
  H%d = A%d

  ! balance A: A2 = PAP^T
  allocate (scale_(n))
  call gebal(H%d, scale=scale_, job='P', info=info)   ! only permutation allowed
  ASSERT(info==0)

  ! reduce A to upper hessenberg form: A = V H V^H
  allocate (tau(n-1))
  call gehrd(H%d, tau=tau, info=info)
  ASSERT(info==0)

  ! form V explicitly
  call V%init(n)
  V%d = H%d
#ifdef TCMPLX
  call unghr(V%d, tau=tau, info=info)
#else
  call orghr(V%d, tau=tau, info=info)
#endif
  ASSERT(info==0)

  ! Schur decomp of Hessenberg: H = Z*U*Z^H   =>  A = (QZ)*U*(QZ)^H
  call schur(H, V, U, nonHessenberg=.true.)

  ! transform EVecs
  call gebak(V%d, scale_, job='P', side='R', info=info)
  ASSERT(info==0)
end subroutine

subroutine SCHUR_HESSENBERG(H, Z, U, nonHessenberg)
  !! Computes Schur decomposition for a Hessenberg matrices
  !!
  !!  H = Z U Z^H
  !!  U: upper quasi-triangular
  !!  Z: unitary
  !!  A = Q H Q^H

  type(HESSENBERG_TYPE), intent(in)    :: H
    !! upper Hessenberg
  type(DENSE_TYPE),      intent(inout) :: Z
    !! unitary matrix: Z^-1 = Z'. for real: Z^-1 = Z^T.
    !! might store Q on input if nonHessenberg is set
  type(HESSENBERG_TYPE), intent(out)   :: U
    !! upper quasi-triangular matrix.
    !! quasi: 2x2 blocks on diagonal representing EVals.
    !!    for TYPE=cmplx: always triangular
    !!    for TYPE=real:  triangular iff real diagonalizable
  logical, optional,     intent(in)    :: nonHessenberg
    !! Does input Z store the transformation of general A to Hessenberg H?
    !! default: false

  character :: compz
  integer   :: n, info
  logical   :: nonHessenberg_

  nonHessenberg_ = .false.
  if (present(nonHessenberg)) nonHessenberg_ = nonHessenberg

  if (nonHessenberg_) then
    ASSERT(Z%nrows == H%nrows .and. Z%ncols == H%ncols)
  end if

  ASSERT(H%upper)

  n = H%nrows
  if (nonHessenberg_) then
    compz = 'V'
  else
    compz = 'I'
    call Z%init(n)
  end if
  call U%init(n, .true.)

  U%d = H%d
#ifdef TCMPLX
  block
    complex, allocatable :: w(:)
    allocate (w(n))
    call hseqr(U%d, w,      z=Z%d, info=info, compz=compz, job='S')
  end block
#else
  block
    real, allocatable :: wr(:), wi(:)
    allocate (wr(n), wi(n))
    call hseqr(U%d, wr, wi, z=Z%d, info=info, compz=compz, job='S')
  end block
#endif
  ASSERT(info == 0)
end subroutine

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef DENSE_TYPE
#undef HESSENBERG_TYPE

#undef SCHUR_DENSE
#undef SCHUR_HESSENBERG
