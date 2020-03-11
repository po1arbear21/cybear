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
#define VECTOR_FRONT         CONCAT3(vector,T,front)
#define VECTOR_BACK          CONCAT3(vector,T,back)
#define VECTOR_PUSH_ELEM     CONCAT3(vector,T,push_elem)
#define VECTOR_PUSH_ELEMS    CONCAT3(vector,T,push_elems)
#define VECTOR_POP           CONCAT3(vector,T,pop)
#define VECTOR_SHRINK        CONCAT3(vector,T,shrink)
#define VECTOR_TO_ARRAY      CONCAT3(vector,T,to_array)
#define VECTOR_FROM_ARRAY    CONCAT3(vector,T,from_array)

subroutine VECTOR_INIT(this, n, c, x)
  !! initialize vector
  class(VECTOR_TYPE), intent(out) :: this
  integer,            intent(in)  :: n
    !! initial size
  integer, optional,  intent(in)  :: c
    !! initial capacity (default: c = 1)
  TT,      optional,  intent(in)  :: x(n)
    !! initial values

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
  !! destruct vector (deallocate memory)
  class(VECTOR_TYPE), intent(inout) :: this

  if (allocated(this%d)) deallocate (this%d)
end subroutine

subroutine VECTOR_RESET(this)
  !! reset vector to size zero
  class(VECTOR_TYPE), intent(inout) :: this

  this%n = 0
end subroutine

subroutine VECTOR_RESERVE(this, c)
  !! reserve memory (reallocate)
  class(VECTOR_TYPE), intent(inout) :: this
  integer,            intent(in)    :: c
    !! new capacity

  ! local variables
  TT, allocatable :: tmp(:)

  if (allocated(this%d)) then
    ! only reallocate if new capacity is actually larger than old
    if (c > size(this%d)) then
      ! reallocate and copy old data
      allocate (tmp(c))
      tmp(1:size(this%d)) = this%d
      call move_alloc(tmp, this%d)
    end if
  else
    allocate (this%d(c))
  end if
end subroutine

subroutine VECTOR_RESIZE(this, n)
  !! resize vector
  class(VECTOR_TYPE), intent(inout) :: this
  integer,            intent(in)    :: n
    !! new size

  ! reserve memory
  if (n > 0) call this%reserve(2**(ceiling(log(real(n)) / log(2.0))))

  ! change size
  this%n = n
  if (this%n < 0) this%n = 0
end subroutine

function VECTOR_FRONT(this) result(x)
  !! get first element
  class(VECTOR_TYPE), intent(in) :: this
  TT                             :: x
    !! return first element

  x = this%d(1)
end function

function VECTOR_BACK(this) result(x)
  !! get last element
  class(VECTOR_TYPE), intent(in) :: this
  TT                             :: x
    !! return last element

  x = this%d(this%n)
end function

subroutine VECTOR_PUSH_ELEM(this, x)
  !! append element to end of vector
  class(VECTOR_TYPE), intent(inout) :: this
  TT,                 intent(in)    :: x
    !! element to add

  ! double capacity if necessary
  if (this%n + 1 > size(this%d)) call this%reserve(2 * (this%n + 1))

  ! increase size by 1
  this%n = this%n + 1

  ! save element
  this%d(this%n) = x
end subroutine

subroutine VECTOR_PUSH_ELEMS(this, x)
  !! append elements to end of vector
  class(VECTOR_TYPE), intent(inout) :: this
  TT,                 intent(in)    :: x(:)
    !! elements to add

  ! increase size by length of given elements
  this%n = this%n + size(x)

  ! increase capacity to double as much as currently needed (if necessary)
  if (this%n > size(this%d)) call this%reserve(2 * this%n)

  ! save element
  this%d(this%n-size(x)+1:this%n) = x
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

  if (c < size(this%d)) then
    allocate (tmp(c), source = this%d(1:c))
    call move_alloc(tmp, this%d)
  end if
end subroutine

function VECTOR_TO_ARRAY(this) result(x)
  !! convert vector to array
  class(VECTOR_TYPE), intent(in) :: this
  TT                             :: x(this%n)
    !! return elements

  x = this%d(1:this%n)
end function

subroutine VECTOR_FROM_ARRAY(this, x)
  !! convert array to vector (old contents are overridden)
  class(VECTOR_TYPE), intent(inout) :: this
  TT                , intent(in)    :: x(:)
    !! new elements

  call this%resize(size(x))
  this%d(1:this%n) = x
end subroutine

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
#undef VECTOR_FRONT
#undef VECTOR_BACK
#undef VECTOR_PUSH_ELEM
#undef VECTOR_PUSH_ELEMS
#undef VECTOR_POP
#undef VECTOR_SHRINK
#undef VECTOR_TO_ARRAY
#undef VECTOR_FROM_ARRAY
