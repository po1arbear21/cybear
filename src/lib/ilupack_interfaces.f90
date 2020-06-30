module ilupack_interfaces_m
  !! creates interfaces for the external routines of ILUPACK
  !!
  !! for parameter documentation see ilupack.f90: type(ilupack_opt).

  implicit none

  private
  public :: dgnlamginit
  public :: dgnlamgfactor
  public :: dgnlamginfo
  public :: dgnlamgnnz
  public :: dgnlamgsolver
  public :: dgnlamgdelete

  interface
    subroutine dgnlamginit(n,        ia,      ja,       a,       matching, &
                           ordering, droptol, droptolS, condest, restol,   &
                           maxit,    elbow,   lfil,     lfilS,   nrestart, &
                           mixedprecision, ind)
      integer       :: n
      integer       :: ia(*)
      integer       :: ja(*)
      real          ::  a(*)
      integer       :: matching
      character(20) :: ordering
      real          :: droptol
      real          :: droptolS
      real          :: condest
      real          :: restol
      integer       :: maxit
      integer       :: elbow
      integer       :: lfil
      integer       :: lfilS
      integer       :: nrestart
      integer       :: mixedprecision
      integer       :: ind(*)
    end subroutine

    function dgnlamgfactor(param,    prec,                                 &
                           n,        ia,      ja,       a,       matching, &
                           ordering, droptol, droptolS, condest, restol,   &
                           maxit,    elbow,   lfil,     lfilS,   nrestart, &
                           mixedprecision, ind) result(ierr)

      integer       :: param
      integer       :: prec
      integer       :: n
      integer       :: ia(*)
      integer       :: ja(*)
      real          ::  a(*)
      integer       :: matching
      character(20) :: ordering
      real          :: droptol
      real          :: droptolS
      real          :: condest
      real          :: restol
      integer       :: maxit
      integer       :: elbow
      integer       :: lfil
      integer       :: lfilS
      integer       :: nrestart
      integer       :: mixedprecision
      integer       :: ind(*)
      integer       :: ierr
    end function

    subroutine dgnlamginfo(param, prec, n, ia, ja, a)
      !! displaying the multilevel structure

      integer :: param
      integer :: prec
      integer :: n
      integer :: ia(*)
      integer :: ja(*)
      real    ::  a(*)
    end subroutine

    function dgnlamgnnz(param, prec) result(nnz)
      !! the logical number of nonzeros only

      integer :: param
      integer :: prec
      integer :: nnz
    end function

    function dgnlamgsolver(param,    prec,    rhs,      sol,               &
                           n,        ia,      ja,       a,       matching, &
                           ordering, droptol, droptolS, condest, restol,   &
                           maxit,    elbow,   lfil,     lfilS,   nrestart, &
                           mixedprecision, ind) result(ierr)

      integer       :: param
      integer       :: prec
      real          :: rhs(*)
      real          :: sol(*)
      integer       :: n
      integer       :: ia(*)
      integer       :: ja(*)
      real          ::  a(*)
      integer       :: matching
      character(20) :: ordering
      real          :: droptol
      real          :: droptolS
      real          :: condest
      real          :: restol
      integer       :: maxit
      integer       :: elbow
      integer       :: lfil
      integer       :: lfilS
      integer       :: nrestart
      integer       :: mixedprecision
      integer       :: ind(*)
      integer       :: ierr
    end function

    subroutine dgnlamgdelete(param, prec)
      integer :: param
      integer :: prec
    end subroutine
  end interface

end module
