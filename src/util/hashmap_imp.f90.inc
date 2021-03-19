#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_VALUE      CONCAT(vector,T)
#define HASHMAP_TYPE      CONCAT(hashmap,T)
#define HASHMAP_INIT      CONCAT3(hashmap,T,init)
#define HASHMAP_DESTRUCT  CONCAT3(hashmap,T,destruct)
#define HASHMAP_RESET     CONCAT3(hashmap,T,reset)
#define HASHMAP_RESERVE   CONCAT3(hashmap,T,reserve)
#define HASHMAP_LOOKUP    CONCAT3(hashmap,T,lookup)
#define HASHMAP_GET_KEY_1 CONCAT3(hashmap,T,get_key_1)
#define HASHMAP_GET_KEY_N CONCAT3(hashmap,T,get_key_n)
#define HASHMAP_SET_KEY_1 CONCAT3(hashmap,T,set_key_1)
#define HASHMAP_SET_KEY_N CONCAT3(hashmap,T,set_key_n)

subroutine HASHMAP_INIT(this, keysize, c)
  !! initialize hashmap
  class(HASHMAP_TYPE), intent(out) :: this
  integer, optional,   intent(in)  :: keysize
    !! number of integers each key consists of (default: keysize = 1)
  integer, optional,   intent(in)  :: c
    !! initial capacity (default: c = 16)

  integer :: c_

  ! optional args
  this%keysize = 1
  if (present(keysize)) this%keysize = keysize
  c_ = 16
  if (present(c)) c_ = c
  c_ = max(c_, 1)

  ! init vectors
  call this%keys%init(0, c = this%keysize * c_)
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

  integer              :: i, i0, i1
  integer, allocatable :: keys(:)
  TT,      allocatable :: values(:)

  if (allocated(this%table)) then
    ! only reallocate if new capacity is actually larger than old
    if (4*c > size(this%table)) then
      ! copy old keys and values
      allocate (  keys(this%keys%n  ), source = this%keys%d(  1:this%keys%n  ))
      allocate (values(this%values%n), source = this%values%d(1:this%values%n))

      ! reinit vectors
      call this%keys%init(  n = 0, c = this%keysize * c)
      call this%values%init(n = 0, c = c)

      ! reallocate table
      deallocate (this%table)
      allocate (this%table(4*c), source = -1)

      ! reinsert values
      i1 = 0
      do i = 1, size(values)
        i0 = i1 + 1
        i1 = i1 + this%keysize
        call this%set(keys(i0:i1), values(i))
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
  integer,             intent(in) :: key(:)
    !! key to lookup
  integer                         :: i
    !! return index

  integer :: j, j0, j1, h, failsafe

  ASSERT(size(key) == this%keysize)

  failsafe = 0

  h = hash(key)
  do while (failsafe <= size(this%table))
    i = modulo(h, size(this%table)) + 1
    j = this%table(i)
    if (j == -1) exit
    j1 = this%table(i) * this%keysize
    j0 = j1 - this%keysize + 1
    if (all(key == this%keys%d(j0:j1))) exit

    ! collision
    h = hash(h + 1)
    failsafe = failsafe + 1
  end do

  ! too many collisions
  if (failsafe > size(this%table)) i = -1
end function

subroutine HASHMAP_GET_KEY_1(this, key, value, status)
  !! lookup key and return value
  class(HASHMAP_TYPE), intent(in)  :: this
  integer,             intent(in)  :: key
  TT,                  intent(out) :: value
  logical, optional,   intent(out) :: status

  ASSERT(this%keysize == 1)

  call HASHMAP_GET_KEY_N(this, [key], value, status=status)
end subroutine

subroutine HASHMAP_GET_KEY_N(this, key, value, status)
  !! lookup key and return value
  class(HASHMAP_TYPE), intent(in)  :: this
  integer,             intent(in)  :: key(:)
  TT,                  intent(out) :: value
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

subroutine HASHMAP_SET_KEY_1(this, key, value, status)
  !! lookup key and insert/update value
  class(HASHMAP_TYPE), intent(inout) :: this
  integer,             intent(in)    :: key
  TT,                  intent(in)    :: value
  logical, optional,   intent(out)   :: status

  ASSERT(this%keysize == 1)

  call HASHMAP_SET_KEY_N(this, [key], value, status=status)
end subroutine

subroutine HASHMAP_SET_KEY_N(this, key, value, status)
  !! lookup key and insert/update value
  class(HASHMAP_TYPE), intent(inout) :: this
  integer,             intent(in)    :: key(:)
  TT,                  intent(in)    :: value
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
    this%table(i) = this%values%n

    ! reserve more memory if table is about 3/4 full
    if (4 * this%values%n > (3 * size(this%table))) call this%reserve(2 * this%values%n)
  else
    ! update value
    this%values%d(j) = value
  end if

  if (present(status)) status = .true.
end subroutine

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef VECTOR_VALUE
#undef HASHMAP_TYPE
#undef HASHMAP_INIT
#undef HASHMAP_DESTRUCT
#undef HASHMAP_RESET
#undef HASHMAP_RESERVE
#undef HASHMAP_LOOKUP
#undef HASHMAP_GET_KEY_1
#undef HASHMAP_GET_KEY_N
#undef HASHMAP_SET_KEY_1
#undef HASHMAP_SET_KEY_N
