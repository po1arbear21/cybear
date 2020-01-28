#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_TYPE          CONCAT(vector,T)
#define HASHMAP_TYPE         CONCAT(hashmap,T)
#define HASHMAP_INIT         CONCAT3(hashmap,T,init)
#define HASHMAP_RESIZE       CONCAT3(hashmap,T,resize)
#define HASHMAP_GET_INDEX    CONCAT3(hashmap,T,get_index)
#define HASHMAP_GET          CONCAT3(hashmap,T,get)
#define HASHMAP_INSERT       CONCAT3(hashmap,T,insert)

subroutine HASHMAP_INIT(this, n)
  !! initialize hashmap
  class(HASHMAP_TYPE), intent(out) :: this
  integer,             intent(in)  :: n
    !! size of hashmap

  ! init vectors
  call this%keys%init(0, c = n / 4)
  call this%values%init(0, c = n / 4)

  ! allocate table
  allocate (this%table(1:n), source = -1)
end subroutine

subroutine HASHMAP_RESIZE(this, n)
  !! resize hashmap
  class(HASHMAP_TYPE), intent(inout) :: this
  integer,             intent(in)    :: n
    !! new size

  ! local variables
  integer              :: i
  integer, allocatable :: keys(:)
  TT,      allocatable :: values(:)

  ! copy old values
  allocate(keys(this%keys%n), source = this%keys%d(1:this%keys%n))
  allocate(values(this%values%n), source = this%values%d(1:this%values%n))

  ! re-init vectors
  call this%keys%init(2 * size(keys))
  call this%values%init(2 * size(values))

  ! reallocate new table
  deallocate (this%table)
  allocate (this%table(1:n), source = -1)

  ! insert values
  do i = 1, size(keys)
    call this%insert(keys(i), values(i))
  end do
end subroutine

function HASHMAP_GET_INDEX(this, key) result(i)
  !! lookup key and return index
  class(HASHMAP_TYPE), intent(in) :: this
  integer,             intent(in) :: key
  integer                         :: i

  ! local variables
  integer :: j

  i = 0
  do j = 0, size(this%table)
    i = modulo(hash_int(key + j), size(this%table)) + 1
    if (this%table(i) == -1) then
      exit
    else if (this%keys%d(this%table(i)) == key) then
      exit
    end if
  end do
end function

function HASHMAP_GET(this, key) result(v)
  !! lookup key and return pointer to value
  class(HASHMAP_TYPE), target, intent(in) :: this
  integer,                     intent(in) :: key
  TT, pointer                             :: v

  ! local variables
  integer :: i

  nullify(v)

  i = this%get_index(key)
  if (this%table(i) /= -1) then
    v => this%values%d(this%table(i))
  end if
end function

subroutine HASHMAP_INSERT(this, key, value, pos)
  !! insert key value pair
  class(HASHMAP_TYPE), intent(inout) :: this
  integer,             intent(in)    :: key
  TT,                  intent(in)    :: value
  integer, optional,   intent(out)   :: pos
    !! optional output position in keys/values vector

  ! local variables
  integer :: i

  i = this%get_index(key)
  if (this%table(i) == -1) then
    ! new key
    call this%keys%push(key)
    call this%values%push(value)
    this%table(i) = this%keys%n
  else if (this%keys%d(this%table(i)) == key) then
    print *, "key = ", key
    call program_error("key already exists!")
  else
    call program_error("table full!")
  end if
  if (present(pos)) pos = this%table(i)
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
#undef HASHMAP_TYPE
#undef HASHMAP_INIT
#undef HASHMAP_RESIZE
#undef HASHMAP_GET_INDEX
#undef HASHMAP_GET
#undef HASHMAP_INSERT
