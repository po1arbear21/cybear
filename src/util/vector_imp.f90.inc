#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_TYPE          CONCAT(vector,T)
#define VECTOR_INIT          CONCAT3(vector,T,init)
#define VECTOR_DESTRUCT      CONCAT3(vector,T,destruct)
#define VECTOR_RESET         CONCAT3(vector,T,reset)
#define VECTOR_RESERVE       CONCAT3(vector,T,reserve)
#define VECTOR_RESIZE        CONCAT3(vector,T,resize)
#define VECTOR_PUSH          CONCAT3(vector,T,push)
#define VECTOR_POP           CONCAT3(vector,T,pop)
#define VECTOR_SHRINK        CONCAT3(vector,T,shrink)
#define VECTOR_GET           CONCAT3(vector,T,get)
#define VECTOR_SET           CONCAT3(vector,T,set)
#define VECTOR_FIRST         CONCAT3(vector,T,first)
#define VECTOR_LAST          CONCAT3(vector,T,last)

subroutine VECTOR_INIT(this, n, c, x)
  class(VECTOR_TYPE), intent(out) :: this
  integer,            intent(in)  :: n
    !! initial size
  integer, optional,  intent(in)  :: c
    !! initial capacity (default: c = 1)
  TT,      optional,  intent(in)  :: x(n)
     ! initial values

  ! local variables
  integer :: c_

  ! resize
  this%n = n
  if (this%n < 0) this%n = 0

  ! capacity
  c_ = this%n
  if (present(c)) then
    if (c > this%n) c_ = c
  end if
  c_ = max(c_, 1)

  ! reserve memory
  call this%reserve(c_)

  ! set initial values
  if (present(x)) then
    this%d(1:this%n) = x
  end if
end subroutine

subroutine VECTOR_DESTRUCT(this)
  class(VECTOR_TYPE), intent(inout) :: this

  if (allocated(this%d)) deallocate(this%d)
end subroutine

subroutine VECTOR_RESET(this)
  class(VECTOR_TYPE), intent(inout) :: this

  this%n = 0
end subroutine

subroutine VECTOR_RESERVE(this, c)
  class(VECTOR_TYPE), intent(inout) :: this
  integer,            intent(in)    :: c
    !! new capacity
  TT, allocatable                   :: tmp(:)

  if (c > this%c) then
    if (this%c == 0) then
      allocate (this%d(c))
    else
      allocate (tmp(c))
      tmp(1:this%c) = this%d
      call move_alloc(tmp, this%d)
    end if
    this%c = c
  end if
end subroutine

subroutine VECTOR_RESIZE(this, n)
  class(VECTOR_TYPE), intent(inout) :: this
  integer,            intent(in)    :: n
    !! new size

  ! reserve memory
  call this%reserve(n)

  ! change size
  this%n = n
  if (this%n < 0) this%n = 0
end subroutine

subroutine VECTOR_PUSH(this, x)
  !! append element to end of vector
  class(VECTOR_TYPE), intent(inout) :: this
  TT,                 intent(in)    :: x
    !! element to add

  ! increase size by 1
  this%n = this%n + 1

  ! double capacity if necessary
  if (this%n > this%c) call this%reserve(2 * this%n)

  ! save element
  this%d(this%n) = x
end subroutine

subroutine VECTOR_POP(this)
  !! remove last element from vector
  class(VECTOR_TYPE), intent(inout) :: this

  ! decrease size by 1
  if (this%n > 0) this%n = this%n - 1
end subroutine

subroutine VECTOR_SHRINK(this)
  !! remove excess capacity
  class(VECTOR_TYPE), intent(inout) :: this

  ! local variables
  TT, allocatable :: tmp(:)
  integer         :: c

  c = this%n
  if (this%n == 0) c = 1

  if (c < this%c) then
    allocate (tmp(c), source = this%d(1:c))
    deallocate (this%d)
    call move_alloc(tmp, this%d)
    this%c = c
  end if
end subroutine

function VECTOR_GET(this) result(x)
  class(VECTOR_TYPE), intent(in) :: this
  TT                             :: x(this%n)

  x = this%d(1:this%n)
end function

subroutine VECTOR_SET(this, x)
  class(VECTOR_TYPE), intent(inout) :: this
  TT                , intent(in)    :: x(this%n)

  this%d(1:this%n) = x
end subroutine

function VECTOR_FIRST(this) result(x)
  class(VECTOR_TYPE), intent(inout) :: this
  TT                                :: x

  x = this%d(1)
end function

function VECTOR_LAST(this) result(x)
  class(VECTOR_TYPE), intent(inout) :: this
  TT                                :: x

  x = this%d(this%n)
end function

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef VECTOR_TYPE
#undef VECTOR_INIT
#undef VECTOR_DESTRUCT
#undef VECTOR_RESET
#undef VECTOR_RESERVE
#undef VECTOR_RESIZE
#undef VECTOR_PUSH
#undef VECTOR_POP
#undef VECTOR_SHRINK
#undef VECTOR_GET
#undef VECTOR_SET
#undef VECTOR_FIRST
#undef VECTOR_LAST