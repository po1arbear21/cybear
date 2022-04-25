subroutine deque_{}T{}_init(this, n, c, x)
  !! initialize deque
  class(deque_{}T),     intent(out) :: this
  integer,              intent(in)  :: n
    !! initial size
  integer,    optional, intent(in)  :: c
    !! initial capacity (default: c = 1)
  m4_type(T), optional, intent(in)  :: x(n)
    !! initial values

  ! local variables
  integer :: c_

  ! resize
  this%n = max(n, 0)

  ! capacity
  c_ = max(this%n, 1)
  if (present(c)) c_ = max(c_, c)

  ! reserve memory
  call this%reserve(c_)

  ! set front and back pointers
  this%ifront = 1
  this%iback  = this%n

  ! set initial values
  if (present(x)) then
    this%d(1:this%n) = x
  end if
end subroutine

subroutine deque_{}T{}_destruct(this)
  !! destruct deque (deallocate memory)
  class(deque_{}T), intent(inout) :: this

  if (allocated(this%d)) deallocate (this%d)
end subroutine

subroutine deque_{}T{}_reset(this)
  !! reset deque size to zero
  class(deque_{}T), intent(inout) :: this

  this%n      = 0
  this%ifront = 1
  this%iback  = 0
end subroutine

subroutine deque_{}T{}_reserve(this, c)
  !! reserve memory (reallocate)
  class(deque_{}T), intent(inout) :: this
  integer,          intent(in)    :: c
    !! new capacity

  m4_type(T), allocatable :: tmp(:)

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

subroutine deque_{}T{}_resize(this, n)
  !! resize deque
  class(deque_{}T), intent(inout) :: this
  integer,          intent(in)    :: n
    !! new size

  ! reserve memory
  call this%reserve(n)

  ! change size
  this%n = n
  if (this%n < 0) this%n = 0
  this%iback = this%ifront + this%n - 1
  if (this%iback > size(this%d)) this%iback = this%iback - size(this%d)
end subroutine

function deque_{}T{}_front(this) result(x)
  !! get first element
  class(deque_{}T), intent(in) :: this
  m4_type(T)                   :: x
    !! return first element

  x = this%d(this%ifront)
end function

function deque_{}T{}_back(this) result(x)
  !! get last element
  class(deque_{}T), intent(in) :: this
  m4_type(T)                   :: x
    !! return last element

  x = this%d(this%iback)
end function

subroutine deque_{}T{}_push_front(this, x)
  !! prepend element to front of deque
  class(deque_{}T), intent(inout) :: this
  m4_type(T),       intent(in)    :: x
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

subroutine deque_{}T{}_push_back(this, x)
  !! append element to back of deque
  class(deque_{}T), intent(inout) :: this
  m4_type(T),       intent(in)    :: x
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

subroutine deque_{}T{}_pop_front(this)
  !! remove first element from deque
  class(deque_{}T), intent(inout) :: this

  if (this%n > 0) then
    this%n = this%n - 1
    this%ifront = this%ifront + 1
    if (this%ifront > size(this%d)) this%ifront = this%ifront - size(this%d)
  end if
end subroutine

subroutine deque_{}T{}_pop_back(this)
  !! remove last element from deque
  class(deque_{}T), intent(inout) :: this

  if (this%n > 0) then
    this%n = this%n - 1
    this%iback = this%iback - 1
    if (this%iback < 1) this%iback = this%iback + size(this%d)
  end if
end subroutine

subroutine deque_{}T{}_shrink(this)
  !! remove excess capacity
  class(deque_{}T), intent(inout) :: this

  integer                 :: c
  m4_type(T), allocatable :: tmp(:)

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

function deque_{}T{}_to_array(this) result(x)
  !! convert deque to array
  class(deque_{}T), intent(in) :: this
  m4_type(T)                   :: x(this%n)
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

subroutine deque_{}T{}_from_array(this, x)
  !! convert array to deque (old contents are overridden)
  class(deque_{}T), intent(inout) :: this
  m4_type(T),       intent(in)    :: x(:)
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

m4_undefine({T})
