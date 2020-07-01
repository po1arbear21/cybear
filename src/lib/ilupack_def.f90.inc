#define PASTE(X)         X
#define CONCATHELP(X, Y) PASTE(X)Y
#define CONCAT(X, Y)     CONCATHELP(X,Y)

#define INTERFACE_INIT   CONCAT(T0,gnlamginit)
#define INTERFACE_FACTOR CONCAT(T0,gnlamgfactor)
#define INTERFACE_INFO   CONCAT(T0,gnlamginfo)
#define INTERFACE_NNZ    CONCAT(T0,gnlamgnnz)
#define INTERFACE_SOLVER CONCAT(T0,gnlamgsolver)
#define INTERFACE_DELETE CONCAT(T0,gnlamgdelete)

interface
  subroutine INTERFACE_INIT(n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   )

    integer       :: n
    integer       :: ia(*)
    integer       :: ja(*)
    TT            ::  a(*)
    integer       :: matching
    character(20) :: ordering
    real          :: droptol
    real          :: droptolS
    real          :: condest
    real          :: restol
    integer       :: maxit
    real          :: elbow
    integer       :: lfil
    integer       :: lfilS
    integer       :: nrestart
    integer       :: mixedprecision
    integer       :: ind(*)
  end subroutine

  function INTERFACE_FACTOR(param,          prec,                                 &
    &                       n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   ) result(ierr)

    integer       :: param
    integer       :: prec
    integer       :: n
    integer       :: ia(*)
    integer       :: ja(*)
    TT            ::  a(*)
    integer       :: matching
    character(20) :: ordering
    real          :: droptol
    real          :: droptolS
    real          :: condest
    real          :: restol
    integer       :: maxit
    real          :: elbow
    integer       :: lfil
    integer       :: lfilS
    integer       :: nrestart
    integer       :: mixedprecision
    integer       :: ind(*)
    integer       :: ierr
  end function

  subroutine INTERFACE_INFO(param, prec, n, ia, ja, a)
    !! displaying the multilevel structure

    integer :: param
    integer :: prec
    integer :: n
    integer :: ia(*)
    integer :: ja(*)
    TT      ::  a(*)
  end subroutine

  function INTERFACE_NNZ(param, prec) result(nnz)
    !! the logical number of nonzeros only

    integer :: param
    integer :: prec
    integer :: nnz
  end function

  function INTERFACE_SOLVER(param,          prec,    rhs,      sol,               &
    &                       n,              ia,      ja,       a,       matching, &
    &                       ordering,       droptol, droptolS, condest, restol,   &
    &                       maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                       mixedprecision, ind                                   ) result(ierr)

    integer       :: param
    integer       :: prec
    TT            :: rhs(*)
    TT            :: sol(*)
    integer       :: n
    integer       :: ia(*)
    integer       :: ja(*)
    TT            ::  a(*)
    integer       :: matching
    character(20) :: ordering
    real          :: droptol
    real          :: droptolS
    real          :: condest
    real          :: restol
    integer       :: maxit
    real          :: elbow
    integer       :: lfil
    integer       :: lfilS
    integer       :: nrestart
    integer       :: mixedprecision
    integer       :: ind(*)
    integer       :: ierr
  end function

  subroutine INTERFACE_DELETE(param, prec)
    integer :: param
    integer :: prec
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
#undef INTERFACE_SOLVER
#undef INTERFACE_DELETE
