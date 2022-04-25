subroutine hashmap_{}T{}_init(this, c)
  !! initialize hashmap
  class(hashmap_{}T), intent(out) :: this
  integer, optional,  intent(in)  :: c
    !! initial capacity (default: c = 16)

  integer :: c_

  ! optional args
  c_ = 16
  if (present(c)) c_ = c
  c_ = max(c_, 1)

  ! init vectors
  call this%keys%init(   0, c = c_)
  call this%keys_i0%init(0, c = c_+1)
  call this%keys_i0%push(1)
  call this%values%init( 0, c = c_)

  ! reserve
  call this%reserve(c_)
end subroutine

subroutine hashmap_{}T{}_reset(this)
  !! reset hashmap
  class(hashmap_{}T), intent(inout) :: this

  call this%keys%reset()
  call this%values%reset()
  if (allocated(this%table)) this%table = -1
end subroutine

subroutine hashmap_{}T{}_reserve(this, c)
  !! reserve memory
  class(hashmap_{}T), intent(inout) :: this
  integer,            intent(in)    :: c
    !! new capacity

  integer                 :: i
  integer,    allocatable :: keys(:), keys_i0(:)
  m4_type(T), allocatable :: values(:)

  if (allocated(this%table)) then
    ! only reallocate if new capacity is actually larger than old
    if (4*c > size(this%table)) then
      ! copy old keys and values
      allocate (keys(   this%keys%n   ), source = this%keys%d(   1:this%keys%n   ))
      allocate (keys_i0(this%keys_i0%n), source = this%keys_i0%d(1:this%keys_i0%n))
      allocate (values( this%values%n ), source = this%values%d( 1:this%values%n ))

      ! reinit vectors
      call this%keys%init(   n = 0, c = this%keys%n * c / this%values%n)
      call this%keys_i0%init(n = 0, c = c+1)
      call this%keys_i0%push(1)
      call this%values%init( n = 0, c = c)

      ! reallocate table
      deallocate (this%table)
      allocate (this%table(4*c), source = -1)

      ! reinsert values
      do i = 1, size(values)
        call this%set(keys(keys_i0(i):keys_i0(i+1)-1), values(i))
      end do
    end if
  else
    ! allocate table
    allocate (this%table(4*c), source = -1)
  end if
end subroutine

function hashmap_{}T{}_lookup(this, key) result(i)
  !! lookup key and return index
  class(hashmap_{}T), intent(in) :: this
  integer,            intent(in) :: key(:)
    !! key to lookup
  integer                        :: i
    !! return index

  integer :: j, i0, i1, h, failsafe

  m4_assert(size(key) > 0)

  failsafe = 0

  h = hash(key)
  do while (failsafe <= size(this%table))
    i = modulo(h, size(this%table)) + 1
    j = this%table(i)

    ! found empty spot?
    if (j == -1) exit

    ! compare keys (size and content)
    i0 = this%keys_i0%d(j)
    i1 = this%keys_i0%d(j+1)-1
    if (size(key) == (i1 - i0 + 1)) then
      if (all(key == this%keys%d(i0:i1))) exit
    end if

    ! collision
    h        = hash(h + 1)
    failsafe = failsafe + 1
  end do

  ! too many collisions
  if (failsafe > size(this%table)) i = -1
end function

subroutine hashmap_{}T{}_get_key_1(this, key, value, status)
  !! lookup key and return value
  class(hashmap_{}T), intent(in)  :: this
  integer,            intent(in)  :: key
  m4_type(T),         intent(out) :: value
  logical, optional,  intent(out) :: status

  call hashmap_{}T{}_get_key_n(this, [key], value, status=status)
end subroutine

subroutine hashmap_{}T{}_get_key_n(this, key, value, status)
  !! lookup key and return value
  class(hashmap_{}T), intent(in)  :: this
  integer,            intent(in)  :: key(:)
  m4_type(T),         intent(out) :: value
  logical, optional,  intent(out) :: status

  integer :: i, j

  if (present(status)) status = .false.

  i = this%lookup(key)
  if (i <= 0) return

  j = this%table(i)
  if (j <= 0) return

  value = this%values%d(j)
  if (present(status)) status = .true.
end subroutine

subroutine hashmap_{}T{}_set_key_1(this, key, value, status)
  !! lookup key and insert/update value
  class(hashmap_{}T), intent(inout) :: this
  integer,            intent(in)    :: key
  m4_type(T),         intent(in)    :: value
  logical, optional,  intent(out)   :: status

  call hashmap_{}T{}_set_key_n(this, [key], value, status=status)
end subroutine

recursive subroutine hashmap_{}T{}_set_key_n(this, key, value, status)
  !! lookup key and insert/update value
  class(hashmap_{}T), intent(inout) :: this
  integer,            intent(in)    :: key(:)
  m4_type(T),         intent(in)    :: value
  logical, optional,  intent(out)   :: status

  integer :: i, j

  if (present(status)) status = .false.

  i = this%lookup(key)
  if (i <= 0) return

  j = this%table(i)
  if (j <= 0) then
    ! new key
    call this%keys%push(key)
    call this%keys_i0%push(this%keys%n + 1)
    call this%values%push(value)
    this%table(i) = this%values%n

    ! reserve more memory if table is about 3/4 full
    if (4 * this%values%n > (3 * size(this%table))) call this%reserve(2 * this%values%n)
  else
    ! update value
    this%values%d(j) = value
  end if

  if (present(status)) status = .true.
end subroutine

m4_undefine({T})
