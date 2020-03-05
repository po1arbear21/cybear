#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DEQUE_TYPE           CONCAT(deque,T)
#define DEQUE_INIT           CONCAT3(deque,T,init)
#define DEQUE_DESTRUCT       CONCAT3(deque,T,destruct)
#define DEQUE_RESET          CONCAT3(deque,T,reset)
#define DEQUE_RESERVE        CONCAT3(deque,T,reserve)
#define DEQUE_RESIZE         CONCAT3(deque,T,resize)
#define DEQUE_FRONT          CONCAT3(deque,T,first)
#define DEQUE_BACK           CONCAT3(deque,T,last)
#define DEQUE_PUSH_FRONT     CONCAT3(deque,T,push_front)
#define DEQUE_PUSH_BACK      CONCAT3(deque,T,push_back)
#define DEQUE_POP_FRONT      CONCAT3(deque,T,pop_front)
#define DEQUE_POP_BACK       CONCAT3(deque,T,pop_back)
#define DEQUE_SHRINK         CONCAT3(deque,T,shrink)
#define DEQUE_TO_ARRAY       CONCAT3(deque,T,to_array)
#define DEQUE_FROM_ARRAY     CONCAT3(deque,T,from_array)

subroutine DEQUE_INIT(this, n, c, x)
  !! initialize deque
  class(DEQUE_TYPE), intent(out) :: this
  integer,           intent(in)  :: n
    !! initial size
  integer, optional, intent(in)  :: c
    !! initial capacity (default: c = 1)
  TT,      optional, intent(in)  :: x(n)
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

  ! set front an back pointers
  this%ifront = 1
  this%iback  = this%n

  ! set initial values
  if (present(x)) then
    this%d(1:this%n) = x
  end if
end subroutine

subroutine DEQUE_DESTRUCT(this)
  !! destruct deque (deallocate memory)
  class(DEQUE_TYPE), intent(inout) :: this

  if (allocated(this%d)) deallocate (this%d)
end subroutine

subroutine DEQUE_RESET(this)
  !! reset deque size to zero
  class(DEQUE_TYPE), intent(inout) :: this

  this%n      = 0
  this%ifront = 1
  this%iback  = 0
end subroutine

subroutine DEQUE_RESERVE(this, c)
  !! reserve memory (reallocate)
  class(DEQUE_TYPE), intent(inout) :: this
  integer,           intent(in)    :: c
    !! new capacity

  ! local variables
  TT, allocatable :: tmp(:)

  if (allocated(this%d)) then
    ! only reallocate if new capacity is actually larger than old
    if (c > size(this%d)) then
      ! reallocate and copy old data
      allocate (tmp(c))
      if (this%n > 0) then
        tmp(1:this%n) = this%to_array()
      end if

      ! move allocation and update front an back pointers
      call move_alloc(tmp, this%d)
      this%ifront = 1
      this%iback  = this%n
    end if
  else
    allocate (this%d(c))
  end if
end subroutine

subroutine DEQUE_RESIZE(this, n)
  !! resize deque
  class(DEQUE_TYPE), intent(inout) :: this
  integer,           intent(in)    :: n
    !! new size

  ! reserve memory
  call this%reserve(n)

  ! change size
  this%n = n
  if (this%n < 0) this%n = 0
  this%iback = this%ifront + this%n - 1
  if (this%iback > size(this%d)) this%iback = this%iback - size(this%d)
end subroutine

function DEQUE_FRONT(this) result(x)
  !! get first element
  class(DEQUE_TYPE), intent(in) :: this
  TT                            :: x
    !! return first element

  x = this%d(this%ifront)
end function

function DEQUE_BACK(this) result(x)
  !! get last element
  class(DEQUE_TYPE), intent(in) :: this
  TT                            :: x
    !! return last element

  x = this%d(this%iback)
end function

subroutine DEQUE_PUSH_FRONT(this, x)
  !! prepend element to front of deque
  class(DEQUE_TYPE), intent(inout) :: this
  TT,                intent(in)    :: x
    !! element to add

  ! double capacity if necessary
  if (this%n + 1 > size(this%d)) call this%reserve(2 * (this%n + 1))

  if (this%n == 0) then
    this%n = 1
    this%ifront = 1
    this%iback  = 1
  else
    this%n = this%n + 1
    this%ifront = this%ifront - 1
    if (this%ifront < 1) this%ifront = this%ifront + size(this%d)
  end if

  ! save element
  this%d(this%ifront) = x
end subroutine

subroutine DEQUE_PUSH_BACK(this, x)
  !! append element to back of deque
  class(DEQUE_TYPE), intent(inout) :: this
  TT,                intent(in)    :: x
    !! element to add

  ! double capacity if necessary
  if (this%n + 1 > size(this%d)) call this%reserve(2 * (this%n + 1))

  if (this%n == 0) then
    this%n = 1
    this%ifront = 1
    this%iback  = 1
  else
    this%n = this%n + 1
    this%iback = this%iback + 1
    if (this%iback > size(this%d)) this%iback = this%iback - size(this%d)
  end if

  ! save element
  this%d(this%iback) = x
end subroutine

subroutine DEQUE_POP_FRONT(this)
  !! remove first element from deque
  class(DEQUE_TYPE), intent(inout) :: this

  if (this%n > 0) then
    this%n = this%n - 1
    this%ifront = this%ifront + 1
    if (this%ifront > size(this%d)) this%ifront = this%ifront - size(this%d)
  end if
end subroutine

subroutine DEQUE_POP_BACK(this)
  !! remove last element from deque
  class(DEQUE_TYPE), intent(inout) :: this

  if (this%n > 0) then
    this%n = this%n - 1
    this%iback = this%iback - 1
    if (this%iback < 1) this%iback = this%iback + size(this%d)
  end if
end subroutine

subroutine DEQUE_SHRINK(this)
  !! remove excess capacity
  class(DEQUE_TYPE), intent(inout) :: this

  ! local variables
  TT, allocatable :: tmp(:)
  integer         :: c

  c = this%n
  if (this%n == 0) c = 1

  if (c < size(this%d)) then
    allocate (tmp(c))
    if (this%n > 0) then
      tmp(1:this%n) = this%to_array()

      call move_alloc(tmp, this%d)
      this%ifront = 1
      this%iback  = this%n
    end if
  end if
end subroutine

function DEQUE_TO_ARRAY(this) result(x)
  !! convert deque to array
  class(DEQUE_TYPE), intent(in) :: this
  TT                            :: x(this%n)
    !! return elements in deque

  if (this%n > 0) then
    if (this%ifront <= this%iback) then
      ! one continuous block of data
      x(1:this%n) = this%d(this%ifront:this%iback)
    else
      ! two data blocks
      x(1:(size(this%d)-this%ifront+1)) = this%d(this%ifront:size(this%d))
      x((size(this%d)-this%ifront+2):this%n) = this%d(1:this%iback)
    end if
  end if
end function

subroutine DEQUE_FROM_ARRAY(this, x)
  !! convert array to deque (old contents are overridden)
  class(DEQUE_TYPE), intent(inout) :: this
  TT               , intent(in)    :: x(:)
    !! new elements

  if (allocated(this%d)) then
    if (size(this%d) < size(x)) then
      deallocate (this%d)
      allocate (this%d(1:size(x)), source = x)
    else
      this%d(1:size(x)) = x
    end if
  else
    allocate (this%d(1:size(x)), source = x)
  end if

  this%n      = size(x)
  this%ifront = 1
  this%iback  = this%n
end subroutine

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef DEQUE_TYPE
#undef DEQUE_INIT
#undef DEQUE_DESTRUCT
#undef DEQUE_RESET
#undef DEQUE_RESERVE
#undef DEQUE_RESIZE
#undef DEQUE_FRONT
#undef DEQUE_BACK
#undef DEQUE_PUSH_FRONT
#undef DEQUE_PUSH_BACK
#undef DEQUE_POP_FRONT
#undef DEQUE_POP_BACK
#undef DEQUE_SHRINK
#undef DEQUE_TO_ARRAY
#undef DEQUE_FROM_ARRAY