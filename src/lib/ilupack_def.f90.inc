#define PASTE(X)         X
#define CONCATHELP(X, Y) PASTE(X)Y
#define CONCAT(X, Y)     CONCATHELP(X,Y)

#define INTERFACE_INIT   CONCAT(T0,gnlamginit)
#define INTERFACE_FACTOR CONCAT(T0,gnlamgfactor)
#define INTERFACE_INFO   CONCAT(T0,gnlamginfo)
#define INTERFACE_NNZ    CONCAT(T0,gnlamgnnz)
#define INTERFACE_SOL    CONCAT(T0,gnlamgsol)
#define INTERFACE_SOLVER CONCAT(T0,gnlamgsolver)
#define INTERFACE_DELETE CONCAT(T0,gnlamgdelete)

interface

  subroutine INTERFACE_INIT(n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   )
    !! init default parameters

    import :: SPARSE_IDX

    integer(SPARSE_IDX) :: n
    integer(SPARSE_IDX) :: ia(*)
    integer(SPARSE_IDX) :: ja(*)
    TT                  ::  a(*)
    integer(SPARSE_IDX) :: matching
    character(20)       :: ordering
    real                :: droptol
    real                :: droptolS
    real                :: condest
    real                :: restol
    integer(SPARSE_IDX) :: maxit
    real                :: elbow
    integer(SPARSE_IDX) :: lfil
    integer(SPARSE_IDX) :: lfilS
    integer(SPARSE_IDX) :: nrestart
    integer(SPARSE_IDX) :: mixedprecision
    integer(SPARSE_IDX) :: ind(*)
  end subroutine

  function INTERFACE_FACTOR(param,          prec,                                 &
    &                       n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   ) result(ierr)
    !! compute multilevel ILU 'prec'
    !!
    !! Note that the initial input matrix A will be rescaled by rows and
    !! by columns (powers of 2.0) and that the order in the array might have
    !! been altered
    !! if you do need the original matrix (ia,ja,a) in for different purposes,
    !! you should use a copy (ib,jb,b) instead

    import :: int64, SPARSE_IDX

    integer(int64)      :: param
    integer(int64)      :: prec
    integer(SPARSE_IDX) :: n
    integer(SPARSE_IDX) :: ia(*)
    integer(SPARSE_IDX) :: ja(*)
    TT                  ::  a(*)
    integer(SPARSE_IDX) :: matching
    character(20)       :: ordering
    real                :: droptol
    real                :: droptolS
    real                :: condest
    real                :: restol
    integer(SPARSE_IDX) :: maxit
    real                :: elbow
    integer(SPARSE_IDX) :: lfil
    integer(SPARSE_IDX) :: lfilS
    integer(SPARSE_IDX) :: nrestart
    integer(SPARSE_IDX) :: mixedprecision
    integer(SPARSE_IDX) :: ind(*)
    integer(SPARSE_IDX) :: ierr
  end function

  subroutine INTERFACE_INFO(param, prec, n, ia, ja, a)
    !! displaying the multilevel structure

    import :: int64, SPARSE_IDX

    integer(int64)      :: param
    integer(int64)      :: prec
    integer(SPARSE_IDX) :: n
    integer(SPARSE_IDX) :: ia(*)
    integer(SPARSE_IDX) :: ja(*)
    TT                  ::  a(*)
  end subroutine

  function INTERFACE_NNZ(param, prec) result(nnz)
    !! the logical number of nonzeros only

    import :: int64, SPARSE_IDX

    integer(int64)      :: param
    integer(int64)      :: prec
    integer(SPARSE_IDX) :: nnz
  end function

  subroutine INTERFACE_SOL(param, prec, rhs, sol, n)
    !! solve a single linear system with `PREC'

    import :: int64, SPARSE_IDX

    integer(int64)      :: param
    integer(int64)      :: prec
    TT                  :: rhs(*)
    TT                  :: sol(*)
    integer(SPARSE_IDX) :: n
  end subroutine

  function INTERFACE_SOLVER(param,          prec,    rhs,      sol,               &
    &                       n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   ) result(ierr)
    !! solve Ax=b iteratively until the desired accuracy is achieved

    import :: int64, SPARSE_IDX

    integer(int64)      :: param
    integer(int64)      :: prec
    TT                  :: rhs(*)
    TT                  :: sol(*)
    integer(SPARSE_IDX) :: n
    integer(SPARSE_IDX) :: ia(*)
    integer(SPARSE_IDX) :: ja(*)
    TT                  ::  a(*)
    integer(SPARSE_IDX) :: matching
    character(20)       :: ordering
    real                :: droptol
    real                :: droptolS
    real                :: condest
    real                :: restol
    integer(SPARSE_IDX) :: maxit
    real                :: elbow
    integer(SPARSE_IDX) :: lfil
    integer(SPARSE_IDX) :: lfilS
    integer(SPARSE_IDX) :: nrestart
    integer(SPARSE_IDX) :: mixedprecision
    integer(SPARSE_IDX) :: ind(*)
    integer(SPARSE_IDX) :: ierr
  end function

  subroutine INTERFACE_DELETE(param, prec)
    !! Finally release memory

    import :: int64, SPARSE_IDX

    integer(int64) :: param
    integer(int64) :: prec
  end subroutine

end interface


#undef T0
#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef CONCATHELP
#undef CONCAT

#undef INTERFACE_INIT
#undef INTERFACE_FACTOR
#undef INTERFACE_INFO
#undef INTERFACE_NNZ
#undef INTERFACE_SOL
#undef INTERFACE_SOLVER
#undef INTERFACE_DELETE
