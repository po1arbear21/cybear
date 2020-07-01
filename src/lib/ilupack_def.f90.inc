#define PASTE1a(X)           X
#define PASTE1b(X)           PASTE1a(X)_
#define CONCATHELP2a(X, Y)   PASTE1a(X)Y
#define CONCATHELP2b(X, Y)   PASTE1b(X)Y
#define CONCAT2a(X, Y)       CONCATHELP2a(X,Y)
#define CONCAT2b(X, Y)       CONCATHELP2b(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE1b(CONCATHELP2b(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ILUPACK_TYPE         CONCAT2b(ilupack,T)

#define INTERFACE_INIT       CONCAT2a(T0,gnlamginit)
#define INTERFACE_FACTOR     CONCAT2a(T0,gnlamgfactor)
#define INTERFACE_INFO       CONCAT2a(T0,gnlamginfo)
#define INTERFACE_NNZ        CONCAT2a(T0,gnlamgnnz)
#define INTERFACE_SOLVER     CONCAT2a(T0,gnlamgsolver)
#define INTERFACE_DELETE     CONCAT2a(T0,gnlamgdelete)

#define TYPE_INIT            CONCAT3(ilupack,T,init)
#define TYPE_FACTOR          CONCAT3(ilupack,T,factor)
#define TYPE_INFO            CONCAT3(ilupack,T,info)
#define TYPE_NNZ             CONCAT3(ilupack,T,nnz)
#define TYPE_SOLVE           CONCAT3(ilupack,T,solve)
#define TYPE_DELETE          CONCAT3(ilupack,T,delete)
#define TYPE_PRINT_OPTS      CONCAT3(ilupack,T,print_opts)

type ILUPACK_TYPE
  !! ilupack wrapper.
  !! usage:
  !!    1. call init: sets default options
  !!    [1b. overwrite option values manually]
  !!    2. call factor
  !!    [2a. call info: writes info about multilevel structure]
  !!    [2b. get fill-in-factor]
  !!    3. call solve
  !!    [3b. call solver for multiple rhs]
  !!    4. call delete

  !
  ! ilupack options (can be set manually)
  !
  character(20) :: ordering
    !! string indicating which reordering to apply
    !! multilevel orderings
    !! 'amd'           Approximate Minimum Degree
    !! 'mmd'           Minimum Degree
    !! 'rcm'           Reverse Cuthill-McKee
    !! 'metisn'        Metis multilevel nested dissection by nodes
    !! 'metise'        Metis multilevel nested dissection by edges
    !! 'pq'            ddPQ strategy by Saad
    !! default: 'amd
  real          :: elbow
    !! elbow space factor for the fill computed during the ILU
    !! default: 10
  integer       :: lfil
    !! number of fill elements per column in L resp. rows in U
    !! default: n+1
  integer       :: lfilS
    !! maximum number of nonzeros per row in the approximate Schur complement
    !! default: n+1
  integer       :: matching
    !! maximum weight matching
    !! default: 1 == is different from zero == matching turned on
  integer       :: maxit
    !! maximum number of iterative steps
    !! default: (500??)
  integer       :: nrestart
    !! restart length for GMRES
    !! default: 30
  real          :: droptol
    !! threshold for ILU
    !! default: 1e-2
  real          :: droptolS
    !! threshold for the approximate Schur complements
    !! info: should be a smaller than droptol (~10%?)
    !! default: 1e-3
  real          :: condest
    !! norm bound for the inverse factors L^{-1}, U^{-1}
    !! default: 1e2
  real          :: restol
    !! relative error for the backward error (SPD case: relative energy norm) used during the iterative solver
    !! default: sqrt(eps)
  integer       :: mixedprecision
    !! use single precision preconditioner if == 1
    !! default: 0

  !
  ! parameters (should not be set manually)
  !
  integer :: n
    !! size of the system. matrix is of (nxn) size.
  integer :: it
    !! number of iterations at solver step


  !
  ! ilupack workspace variables
  !
  integer, allocatable :: ind(:)
    !! array of size n
    !! default: 0
  integer              :: param
    !! parameter C-pointer casted to integer
  integer              :: prec
    !! preconditioner C-pointer casted to integer

  !
  ! copy of ia,ja,a. factor rescales them...
  !
  integer, allocatable :: ia(:), ja(:)
  TT,      allocatable :: a(:)

contains
  procedure :: init       => TYPE_INIT
  procedure :: factor     => TYPE_FACTOR
  procedure :: info       => TYPE_INFO
  procedure :: nnz        => TYPE_NNZ
  procedure :: solve      => TYPE_SOLVE
  procedure :: delete     => TYPE_DELETE
  procedure :: print_opts => TYPE_PRINT_OPTS
end type

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

  function INTERFACE_SOLVER(param,        prec,    rhs,      sol,               &
    &                     n,              ia,      ja,       a,       matching, &
    &                     ordering,       droptol, droptolS, condest, restol,   &
    &                     maxit,          elbow,   lfil,     lfilS,   nrestart, &
    &                     mixedprecision, ind                                   ) result(ierr)

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

#undef PASTE1a
#undef PASTE1b
#undef CONCATHELP2a
#undef CONCATHELP2b
#undef CONCAT2a
#undef CONCAT2b
#undef CONCATHELP3
#undef CONCAT3

#undef ILUPACK_TYPE

#undef INTERFACE_INIT
#undef INTERFACE_FACTOR
#undef INTERFACE_INFO
#undef INTERFACE_NNZ
#undef INTERFACE_SOLVER
#undef INTERFACE_DELETE

#undef TYPE_INIT
#undef TYPE_FACTOR
#undef TYPE_INFO
#undef TYPE_NNZ
#undef TYPE_SOLVE
#undef TYPE_DELETE
#undef TYPE_PRINT_OPTS
