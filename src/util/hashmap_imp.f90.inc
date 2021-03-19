#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_KEY       CONCAT(vector,TKEY)
#define VECTOR_VALUE     CONCAT(vector,TVALUE)
#define HASHMAP_TYPE     CONCAT3(hashmap,TKEY,TVALUE)
#define HASHMAP_INIT     CONCAT(HASHMAP_TYPE,init)
#define HASHMAP_DESTRUCT CONCAT(HASHMAP_TYPE,destruct)
#define HASHMAP_RESET    CONCAT(HASHMAP_TYPE,reset)
#define HASHMAP_RESERVE  CONCAT(HASHMAP_TYPE,reserve)
#define HASHMAP_LOOKUP   CONCAT(HASHMAP_TYPE,lookup)
#define HASHMAP_GET      CONCAT(HASHMAP_TYPE,get)
#define HASHMAP_SET      CONCAT(HASHMAP_TYPE,set)

subroutine HASHMAP_INIT(this, c)
  !! initialize hashmap
  class(HASHMAP_TYPE), intent(out) :: this
  integer, optional,   intent(in)  :: c
    !! initial capacity (default: c = 16)

  integer :: c_

  ! optional args
  c_ = 16
  if (present(c)) c_ = c
  c_ = max(c_, 1)

  ! init vectors
  call this%keys%init(0, c = c_)
  call this%values%init(0, c = c_)

  ! reserve
  call this%reserve(c_)
end subroutine

subroutine HASHMAP_DESTRUCT(this)
  !! destruct hashmap
  class(HASHMAP_TYPE), intent(inout) :: this

  call this%keys%destruct()
  call this%values%destruct()
  if (allocated(this%table)) deallocate (this%table)
end subroutine

subroutine HASHMAP_RESET(this)
  !! reset hashmap
  class(HASHMAP_TYPE), intent(inout) :: this

  call this%keys%reset()
  call this%values%reset()
  if (allocated(this%table)) this%table = -1
end subroutine

subroutine HASHMAP_RESERVE(this, c)
  !! reserve memory
  class(HASHMAP_TYPE), intent(inout) :: this
  integer,             intent(in)    :: c
    !! new capacity

  integer              :: i
  TTKEY,   allocatable :: keys(:)
  TTVALUE, allocatable :: values(:)

  if (allocated(this%table)) then
    ! only reallocate if new capacity is actually larger than old
    if (4*c > size(this%table)) then
      ! copy old keys and values
      allocate (  keys(this%keys%n  ), source = this%keys%d(  1:this%keys%n  ))
      allocate (values(this%values%n), source = this%values%d(1:this%values%n))

      ! re-init vectors
      call this%keys%init(  n = 0, c = c)
      call this%values%init(n = 0, c = c)

      ! reallocate table
      deallocate (this%table)
      allocate (this%table(4*c), source = -1)

      ! insert values
      do i = 1, size(keys)
        call this%set(keys(i), values(i))
      end do
    end if
  else
    ! allocate table
    allocate (this%table(4*c), source = -1)
  end if
end subroutine

function HASHMAP_LOOKUP(this, key) result(i)
  !! lookup key and return index
  class(HASHMAP_TYPE), intent(in) :: this
  TTKEY,               intent(in) :: key
    !! key to lookup
  integer                         :: i
    !! return index

  integer :: h, failsafe

  failsafe = 0

  h = hash(key)
  do while (failsafe <= size(this%table))
    i = modulo(h, size(this%table)) + 1
    if (this%table(i) == -1) then
      exit
    elseif (key_compare(key, this%keys%d(this%table(i)))) then
      exit
    end if

    ! collision
    h = hash(h + 1)
    failsafe = failsafe + 1
  end do

  ! too many collisions
  if (failsafe > size(this%table)) i = -1
end function

subroutine HASHMAP_GET(this, key, value, status)
  !! lookup key and return value
  class(HASHMAP_TYPE), intent(in)  :: this
  TTKEY,               intent(in)  :: key
  TTVALUE,             intent(out) :: value
  logical, optional,   intent(out) :: status

  integer :: i, j

  if (present(status)) status = .false.

  i = this%lookup(key)
  if (i <= 0) return

  j = this%table(i)
  if (j <= 0) return

  value = this%values%d(j)
  if (present(status)) status = .true.
end subroutine

subroutine HASHMAP_SET(this, key, value, status)
  !! lookup key and insert/update value
  class(HASHMAP_TYPE), intent(inout) :: this
  TTKEY,               intent(in)    :: key
  TTVALUE,             intent(in)    :: value
  logical, optional,   intent(out)   :: status

  integer :: i, j

  if (present(status)) status = .false.

  i = this%lookup(key)
  if (i <= 0) return

  j = this%table(i)
  if (j <= 0) then
    ! new key
    call this%keys%push(key)
    call this%values%push(value)
    this%table(i) = this%keys%n

    ! reserve more memory if table is about 3/4 full
    if (4 * this%keys%n > (3 * size(this%table))) call this%reserve(2 * this%keys%n)
  else
    ! update value
    this%values%d(j) = value
  end if

  if (present(status)) status = .true.
end subroutine

#undef TKEY
#undef TTKEY
#undef TVALUE
#undef TTVALUE

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef VECTOR_KEY
#undef VECTOR_VALUE
#undef HASHMAP_TYPE
#undef HASHMAP_INIT
#undef HASHMAP_DESTRUCT
#undef HASHMAP_RESET
#undef HASHMAP_RESERVE
#undef HASHMAP_LOOKUP
#undef HASHMAP_GET
#undef HASHMAP_SET
