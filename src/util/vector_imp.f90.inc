subroutine vector_{}T{}_init(this, n, c, x)
  !! initialize vector
  class(vector_{}T),    intent(out) :: this
  integer,              intent(in)  :: n
    !! initial size
  integer,    optional, intent(in)  :: c
    !! initial capacity (default: c = 1)
  m4_type(T), optional, intent(in)  :: x(n)
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

subroutine vector_{}T{}_reset(this)
  !! reset vector to size zero
  class(vector_{}T), intent(inout) :: this

  this%n = 0
end subroutine

subroutine vector_{}T{}_reserve(this, c)
  !! reserve memory (reallocate)
  class(vector_{}T), intent(inout) :: this
  integer,           intent(in)    :: c
    !! new capacity

  m4_type(T), allocatable :: tmp(:)

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

subroutine vector_{}T{}_resize(this, n)
  !! resize vector
  class(vector_{}T), intent(inout) :: this
  integer,           intent(in)    :: n
    !! new size

  ! reserve memory
  if (n > 0) call this%reserve(2**(ceiling(log(real(n)) / log(2.0))))

  ! change size
  this%n = n
  if (this%n < 0) this%n = 0
end subroutine

function vector_{}T{}_front(this) result(x)
  !! get first element
  class(vector_{}T), intent(in) :: this
  m4_type(T)                    :: x
    !! return first element

  x = this%d(1)
end function

function vector_{}T{}_back(this) result(x)
  !! get last element
  class(vector_{}T), intent(in) :: this
  m4_type(T)                    :: x
    !! return last element

  x = this%d(this%n)
end function

subroutine vector_{}T{}_pop(this)
  !! remove last element from vector
  class(vector_{}T), intent(inout) :: this

  ! decrease size by 1
  if (this%n > 0) this%n = this%n - 1
end subroutine

subroutine vector_{}T{}_shrink(this)
  !! remove excess capacity
  class(vector_{}T), intent(inout) :: this

  integer                 :: c
  m4_type(T), allocatable :: tmp(:)

  c = this%n
  if (this%n == 0) c = 1

  if (c < size(this%d)) then
    allocate (tmp(c), source = this%d(1:c))
    call move_alloc(tmp, this%d)
  end if
end subroutine

function vector_{}T{}_to_array(this) result(x)
  !! convert vector to array
  class(vector_{}T), intent(in) :: this
  m4_type(T)                    :: x(this%n)
    !! return elements

  x = this%d(1:this%n)
end function

subroutine vector_{}T{}_from_array(this, x)
  !! convert array to vector (old contents are overridden)
  class(vector_{}T), intent(inout) :: this
  m4_type(T),        intent(in)    :: x(:)
    !! new elements

  call this%resize(size(x))
  this%d(1:this%n) = x
end subroutine

subroutine vector_{}T{}_push(this, x)
  !! append element to end of vector
  class(vector_{}T), intent(inout) :: this
  m4_type(T),        intent(in)    :: x
    !! element to add

  ! double capacity if necessary
  if (this%n + 1 > size(this%d)) call this%reserve(2 * (this%n + 1))

  ! increase size by 1
  this%n = this%n + 1

  ! save element
  this%d(this%n) = x
end subroutine

subroutine vector_{}T{}_push_n(this, x)
  !! append elements to end of vector
  class(vector_{}T), intent(inout) :: this
  m4_type(T),        intent(in)    :: x(:)
    !! elements to add

  ! increase size by length of given elements
  this%n = this%n + size(x)

  ! increase capacity to double as much as currently needed (if necessary)
  if (this%n > size(this%d)) call this%reserve(2 * this%n)

  ! save element
  this%d(this%n-size(x)+1:this%n) = x
end subroutine

m4_undefine({T})
