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

subroutine TYPE_INIT(this, ia, ja, a)
  !! init ilupack options.

  class(ILUPACK_TYPE), intent(out) :: this
  integer,             intent(in)  :: ia(:)
  integer,             intent(in)  :: ja(:)
  TT,                  intent(in)  ::  a(:)

  ASSERT(size(a) == size(ja))

  this%n = size(ia)-1
  allocate (this%ind(this%n), source=0)

  this%ia = ia
  this%ja = ja
  this%a  = a

  call INTERFACE_INIT(this%n,              this%ia,      this%ja,       this%a,       this%matching, &
    &                 this%ordering,       this%droptol, this%droptolS, this%condest, this%restol,   &
    &                 this%maxit,          this%elbow,   this%lfil,     this%lfilS,   this%nrestart, &
    &                 this%mixedprecision, this%ind                                                  )

end subroutine

subroutine TYPE_PRINT_OPTS(this)
  !! print ilupack options.

  class(ILUPACK_TYPE), intent(in) :: this

  print '(A)',        'ordering:      "' // this%ordering // '"'
  print '(A,E24.16)', 'elbow:          ',   this%elbow
  print '(A,I24)',    'lfil:           ',   this%lfil
  print '(A,I24)',    'lfilS:          ',   this%lfilS
  print '(A,I24)',    'matching:       ',   this%matching
  print '(A,I24)',    'maxit:          ',   this%maxit
  print '(A,I24)',    'nrestart:       ',   this%nrestart
  print '(A,E24.16)', 'droptol:        ',   this%droptol
  print '(A,E24.16)', 'droptolS:       ',   this%droptolS
  print '(A,E24.16)', 'condest:        ',   this%condest
  print '(A,E24.16)', 'restol:         ',   this%restol
  print '(A,I24)',    'mixedprecision: ',   this%mixedprecision
end subroutine

subroutine TYPE_FACTOR(this)
  class(ILUPACK_TYPE), intent(inout) :: this

  integer :: ierr

  ierr = INTERFACE_FACTOR(this%param,          this%prec,                                                &
    &                     this%n,              this%ia,      this%ja,       this%a,       this%matching, &
    &                     this%ordering,       this%droptol, this%droptolS, this%condest, this%restol,   &
    &                     this%maxit,          this%elbow,   this%lfil,     this%lfilS,   this%nrestart, &
    &                     this%mixedprecision, this%ind                                                  )

  if      ((ierr > -8) .and. (ierr < 0)) then
    call program_error("Error in ilupack: " // trim(ILUPACK_FACTOR_ERROR(ierr)))
  else if ( ierr /= 0                  ) then
    call program_error("Error in ilupack: zero pivot encountered at step number " // int2str(ierr))
  end if
end subroutine

subroutine TYPE_INFO(this)
  class(ILUPACK_TYPE), intent(in) :: this

  call INTERFACE_INFO(this%param, this%prec, this%n, this%ia, this%ja, this%a)
end subroutine

integer function TYPE_NNZ(this)
  class(ILUPACK_TYPE), intent(in) :: this

  TYPE_NNZ = INTERFACE_NNZ(this%param, this%prec)
end function

subroutine TYPE_SOLVE(this, rhs, sol)
  class(ILUPACK_TYPE), intent(inout) :: this
  TT,                  intent(in)    :: rhs(:)
    !! right hand side  of length n
  TT,                  intent(inout) :: sol(:)
    !! solution vector of length n, must be initialized on input with an initial guess

  integer :: ierr

  ASSERT(size(rhs) == this%n)
  ASSERT(size(sol) == this%n)

  this%it = this%maxit

  ierr = INTERFACE_SOLVER(this%param,          this%prec,    rhs,           sol,                         &
    &                    this%n,              this%ia,      this%ja,       this%a,       this%matching, &
    &                    this%ordering,       this%droptol, this%droptolS, this%condest, this%restol,   &
    &                    this%it,             this%elbow,   this%lfil,     this%lfilS,   this%nrestart, &
    &                    this%mixedprecision, this%ind                                                  )

  if      ((ierr >= -3) .and. (ierr < 0)) then
    call program_error("Error in ilupack: " // trim(ILUPACK_SOLVER_ERROR(ierr)))
  else if ( ierr /=  0                  ) then
    call program_error("Error in ilupack: solver exited with error code " // int2str(ierr))
  end if
end subroutine

subroutine TYPE_DELETE(this)
  class(ILUPACK_TYPE), intent(inout) :: this

  call INTERFACE_DELETE(this%param, this%prec)
end subroutine


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
