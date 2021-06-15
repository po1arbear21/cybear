#define PASTE(X)                X
#define PASTE2(X)               PASTE(X)_
#define CONCATHELP(X, Y)        PASTE2(X)Y
#define CONCAT(X, Y)            CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z)    PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)        CONCATHELP3(X,Y,Z)
#define CONCATHELP4(X, Y, Z, A) PASTE2(CONCATHELP3(X,Y,Z))A
#define CONCAT4(X, Y, Z, A)     CONCATHELP4(X,Y,Z,A)

#define RB_RED                  1
#define RB_BLACK                2
#define MAPNODE_TYPE            CONCAT3(mapnode,T,U)
#define MAPNODE_INIT            CONCAT4(mapnode,T,U,init)
#define MAPNODE_DESTRUCT        CONCAT4(mapnode,T,U,destruct)
#define MAPNODE_COPY            CONCAT4(mapnode,T,U,copy)
#define MAPNODE_GRANDPARENT     CONCAT4(mapnode,T,U,grandparent)
#define MAPNODE_SIBLING         CONCAT4(mapnode,T,U,sibling)
#define MAPNODE_UNCLE           CONCAT4(mapnode,T,U,uncle)
#define MAPNODE_NUM_CHILDREN    CONCAT4(mapnode,T,U,num_children)
#define MAPNODE_NEXT            CONCAT4(mapnode,T,U,next)
#define MAPNODE_PREVIOUS        CONCAT4(mapnode,T,U,previous)
#define MAPNODE_ROTATE_LEFT     CONCAT4(mapnode,T,U,rotate_left)
#define MAPNODE_ROTATE_RIGHT    CONCAT4(mapnode,T,U,rotate_right)
#define MAP_TYPE                CONCAT3(map,T,U)
#define MAP_INIT                CONCAT4(map,T,U,init)
#define MAP_DESTRUCT            CONCAT4(map,T,U,destruct)
#define MAP_COPY                CONCAT4(map,T,U,copy)
#define MAP_FIND                CONCAT4(map,T,U,find)
#define MAP_FIND_POS            CONCAT4(map,T,U,find_pos)
#define MAP_GET                 CONCAT4(map,T,U,get)
#define MAP_SET                 CONCAT4(map,T,U,set)
#define MAP_LEFTMOST            CONCAT4(map,T,U,leftmost)
#define MAP_RIGHTMOST           CONCAT4(map,T,U,rightmost)
#define MAP_INSERT              CONCAT4(map,T,U,insert)
#define MAP_TO_ARRAY            CONCAT4(map,T,U,to_array)

subroutine MAPNODE_INIT(this, key, value, color, left, right, parent)
  !! initialize map-node
  class(MAPNODE_TYPE),         intent(out) :: this
  TT,                          intent(in)  :: key
    !! element key
  UU,                          intent(in)  :: value
    !! element value
  integer,                     intent(in)  :: color
    !! node color
  type(MAPNODE_TYPE), pointer, intent(in)  :: left
    !! pointer to left child
  type(MAPNODE_TYPE), pointer, intent(in)  :: right
    !! pointer to right child
  type(MAPNODE_TYPE), pointer, intent(in)  :: parent
    !! pointer to parent

  this%key    =  key
  this%value  =  value
  this%color  =  color
  this%left   => left
  this%right  => right
  this%parent => parent
end subroutine

recursive subroutine MAPNODE_DESTRUCT(this)
  !! destruct map-node (and children)
  class(MAPNODE_TYPE), intent(inout) :: this

  ! destruct left child
  if (associated(this%left)) then
    call this%left%destruct()
    deallocate(this%left)
  end if

  ! destruct right child
  if (associated(this%right)) then
    call this%right%destruct()
    deallocate(this%right)
  end if
end subroutine

recursive subroutine MAPNODE_COPY(this, n)
  !! copy map-node from source
  class(MAPNODE_TYPE), target, intent(inout) :: this
  type(MAPNODE_TYPE),          intent(in)    :: n
    !! source mapnode

  ! destruct existing children
  call this%destruct()

  ! copy key, value and color
  this%key   = n%key
  this%value = n%value
  this%color = n%color

  ! copy left child
  if (associated(n%left)) then
    allocate (this%left)
    this%left%parent => this
    call this%left%copy(n%left)
  end if

  ! copy right child
  if (associated(n%right)) then
    allocate (this%right)
    this%right%parent => this
    call this%right%copy(n%right)
  end if
end subroutine

function MAPNODE_GRANDPARENT(this) result(p)
  !! get map-node grandparent (parent of parent)
  class(MAPNODE_TYPE), intent(in) :: this
  type(MAPNODE_TYPE), pointer     :: p
    !! return pointer to grandparent

  if (.not. associated(this%parent)) then
    nullify(p)
  else
    p => this%parent%parent
  end if
end function

function MAPNODE_SIBLING(this) result(p)
  !! get map-node sibling (the other child of this%parent)
  class(MAPNODE_TYPE), target, intent(in) :: this
  type(MAPNODE_TYPE), pointer             :: p
    !! return pointer to sibling

  if (.not. associated(this%parent)) then
    nullify(p)
  else
    if (associated(this%parent%left, this)) then
      p => this%parent%right
    else
      p => this%parent%left
    end if
  end if
end function

function MAPNODE_UNCLE(this) result(p)
  !! get map-node uncle (the other child of this%parent%parent)
  class(MAPNODE_TYPE), intent(in) :: this
  type(MAPNODE_TYPE), pointer     :: p
    !! return pointer to uncle

  if (.not. associated(this%parent)) then
    nullify(p)
  else
    p => this%parent%sibling()
  end if
end function

function MAPNODE_NUM_CHILDREN(this) result(n)
  !! get number of children
  class(MAPNODE_TYPE), intent(in) :: this
  integer                         :: n
    !! return number of children

  n = 0
  if (associated(this%left)) then
    n = n + 1
  end if
  if (associated(this%right)) then
    n = n + 1
  end if
end function

function MAPNODE_NEXT(this) result(next)
  !! get next map-node
  class(MAPNODE_TYPE), target, intent(in) :: this
  type(MAPNODE_TYPE), pointer             :: next
    !! return next map-node

  ! local variables
  type(MAPNODE_TYPE), pointer :: p

  if (associated(this%right)) then
    next => this%right
    do while (associated(next%left))
      next => next%left
    end do
  else
    p => this
    next => p%parent
    if (.not. associated(next)) return
    do while (associated(next%right, p))
      p => next
      next => p%parent
      if (.not. associated(next)) return
    end do
  end if
end function

function MAPNODE_PREVIOUS(this) result(prev)
  !! get previous map-node
  class(MAPNODE_TYPE), target, intent(in) :: this
  type(MAPNODE_TYPE), pointer           :: prev
    !! return previous map-node

  ! local variables
  type(MAPNODE_TYPE), pointer :: p

  if (associated(this%left)) then
    prev => this%left
    do while (associated(prev%left))
      prev => prev%left
    end do
  else
    p => this
    prev => p%parent
    if (.not. associated(prev)) return
    do while (associated(prev%left, p))
      p => prev
      prev => p%parent
      if (.not. associated(prev)) return
    end do
  end if
end function

subroutine MAPNODE_ROTATE_LEFT(this)
  !! rotate tree left
  class(MAPNODE_TYPE), target, intent(inout) :: this

  ! local variables
  type(MAPNODE_TYPE), pointer :: n, g

  g => this%parent
  n => this%right ! must be associated!

  if (associated(g)) then
    if (associated(g%left, this)) then
      g%left => n
    else
      g%right => n
    end if
  end if
  n%parent    => g
  this%parent => n
  this%right  => n%left
  if (associated(n%left)) then
    n%left%parent => this
  end if
  n%left => this
end subroutine

subroutine MAPNODE_ROTATE_RIGHT(this)
  !! rotate tree left
  class(MAPNODE_TYPE), target, intent(inout) :: this

  ! local variables
  type(MAPNODE_TYPE), pointer :: p, g

  g => this%parent
  p => this%left ! must be associated

  if (associated(g)) then
    if (associated(g%left, this)) then
      g%left => p
    else
      g%right => p
    end if
  end if
  p%parent    => g
  this%parent => p
  this%left   => p%right
  if (associated(p%right)) then
    p%right%parent => this
  end if
  p%right => this
end subroutine

subroutine MAP_INIT(this)
  !! initialize map
  class(MAP_TYPE), intent(out) :: this

  ! no elements
  this%n = 0
end subroutine

subroutine MAP_DESTRUCT(this)
  !! destruct map
  class(MAP_TYPE), intent(inout) :: this

  ! destruct root node
  if (associated(this%root)) then
    call this%root%destruct()
    deallocate(this%root)
  end if

  ! clear number of elements
  this%n = 0
end subroutine

subroutine MAP_COPY(this, m)
  !! copy map
  class(MAP_TYPE), intent(inout) :: this
    !! target
  type(MAP_TYPE),  intent(in)    :: m
    !! source

  ! delete old data
  call this%destruct()

  ! copy root
  if (associated(m%root)) then
    allocate (this%root)
    call this%root%copy(m%root)
  end if

  ! set number of elements
  this%n = m%n
end subroutine

function MAP_FIND(this, key) result(node)
  !! find using the key in map
  class(MAP_TYPE), intent(in) :: this
  TT,              intent(in) :: key
    !! key to lookup
  type(MAPNODE_TYPE), pointer :: node
    !! return found node
    !! (null if not found)

  node => this%root
  do while (associated(node))
    if (node%key == key) then
      ! key found
      return
    elseif (key < node%key) then
      ! goto left subtree
      node => node%left
    else
      ! goto right subtree
      node => node%right
    end if
  end do
end function

function MAP_FIND_POS(this, key) result(node)
  ! find nearest key in map
  class(MAP_TYPE),    intent(in) :: this
  TT,                 intent(in) :: key
    !! key to lookup
  type(MAPNODE_TYPE), pointer    :: node
    !! return node which is nearest to key

  ! local variables
  type(MAPNODE_TYPE), pointer :: tmp

  tmp => this%root
  nullify(node)
  do while (associated(tmp))
    node => tmp
    if (tmp%key == key) then
      ! key found
      return
    elseif (key < tmp%key) then
      ! goto left subtree
      tmp => tmp%left
    else
      ! goto right subtree
      tmp => tmp%right
    end if
  end do
end function

function MAP_GET(this, key) result(val)
  !! lookup key in map and return value
  class(MAP_TYPE), intent(in) :: this
  TT,              intent(in) :: key
    !! key to lookup
  UU                          :: val
    !! return value

  ! local variables
  type(MAPNODE_TYPE), pointer :: n

  n => this%find(key)
  val = n%value
end function

subroutine MAP_SET(this, key, value)
  !! insert new key-value-pair into map or update existing one
  class(MAP_TYPE), intent(inout) :: this
  TT,              intent(in)    :: key
  UU,              intent(in)    :: value

  ! local variables
  type(MAPNODE_TYPE), pointer :: node
  logical                     :: status

  call this%insert(key, value, node = node, status = status)

  if (.not. status) node%value = value
end subroutine

function MAP_LEFTMOST(this) result(node)
  !! get leftmost node in map
  class(MAP_TYPE), intent(in) :: this
  type(MAPNODE_TYPE), pointer :: node
    !! return leftmost node

  node => this%root
  if (.not. associated(node)) return
  do while (associated(node%left))
    node => node%left
  end do
end function

function MAP_RIGHTMOST(this) result(node)
  !! get rightmost node in map
  class(MAP_TYPE), intent(in) :: this
  type(MAPNODE_TYPE), pointer :: node
    !! return rightmost node

  node => this%root
  if (.not. associated(node)) return
  do while (associated(node%right))
    node => node%right
  end do
end function

subroutine MAP_INSERT(this, key, value, node, status)
  !! insert key-value-pair into map
    class(MAP_TYPE),                       intent(inout) :: this
    TT,                                    intent(in)    :: key
      !! key to insert
    UU,                                    intent(in)    :: value
      !! value to insert
    type(MAPNODE_TYPE), optional, pointer, intent(out)   :: node
      !! optional output pointer to newly created node (or to already existing one with the same key)
    logical,            optional,          intent(out)   :: status
      !! optional output status of insertion: success => true, fail => false

    ! local variables
    type(MAPNODE_TYPE), pointer :: n, p, g, u

    ! check for empty map
    if (.not. associated(this%root)) then
      allocate (this%root)
      call this%root%init(key, value, RB_BLACK, null(), null(), null())
      this%n = 1
      if (present(node)) node => this%root
      if (present(status)) status = .true.
      return
    end if

    ! find insertion position
    p => this%find_pos(key)

    ! make sure key does not exist already
    if (key == p%key) then
      ! insertion failed: duplicate key
      if (present(node)) node => p
      if (present(status)) status = .false.
      return
    end if

    ! allocate new node and insert it into the map
    allocate (n)
    if (present(node)) node => n
    call n%init(key, value, RB_RED, null(), null(), p)
    if (key < p%key) then
      p%left => n
    else
      p%right => n
    end if

    ! rebalance map
    do while (.true.)
      ! condition 1: n has no parent
      if (.not. associated(p)) then
        this%root => n ! new root
        n%color = RB_BLACK
        exit
      end if

      ! condition 2: parent is black
      if (p%color == RB_BLACK) then
        exit
      end if

      g => p%parent  ! grandparent
      u => n%uncle() ! uncle

      ! condition 3: parent and uncle are red
      if (associated(u)) then
        if (u%color == RB_RED) then
          p%color = RB_BLACK
          u%color = RB_BLACK
          g%color = RB_RED
          n => g
          p => n%parent
          cycle
        end if
      end if

      ! condition 4a: uncle is black and n is the right child of p which is the left child of g
      if (associated(n, p%right) .and. associated(p, g%left)) then
        call p%rotate_left()
        n => p
        p => g%left
      ! condition 4b: uncle is black and n is the left child of p which is the right child of g
      else if (associated(n, p%left) .and. associated(p, g%right)) then
        call p%rotate_right()
        n => p
        p => g%right
      end if

      ! condition 5a: uncle is black and n is the left child of p which is the left child of g
      if (associated(n, p%left) .and. associated(p, g%left)) then
        call g%rotate_right()
        if (associated(g, this%root)) this%root => p
      ! condition 5b: uncle is black and n is the right child of p which is the right child of g
      else
        call g%rotate_left()
        if (associated(g, this%root)) this%root => p
      end if

      p%color = RB_BLACK
      g%color = RB_RED
      exit
    end do

    ! update number of nodes
    this%n = this%n + 1

    if (present(status)) then
      status = .true.
    end if
end subroutine

subroutine MAP_TO_ARRAY(this, keys, values)
  !! convert map to array
  class(MAP_TYPE), intent(in)  :: this
  TT, optional,    intent(out) :: keys(:)
    !! optional output array with keys
  UU, optional,    intent(out) :: values(:)
    !! optional output array with values

  ! local variables
  type(MAPNODE_TYPE), pointer  :: n
  integer                      :: i

  if (present(keys)) then
    ASSERT(size(keys) == this%n)
  end if

  if (present(values)) then
    ASSERT(size(values) == this%n)
  end if

  ! goto leftmost node
  n => this%leftmost()

  ! iterate through map
  i = 1
  do while (associated(n))
    if (present(keys  )) keys(i)   = n%key
    if (present(values)) values(i) = n%value
    n => n%next()
    i = i + 1
  end do
end subroutine

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3
#undef CONCATHELP4
#undef CONCAT4

#undef RB_RED
#undef RB_BLACK
#undef MAPNODE_TYPE
#undef MAPNODE_INIT
#undef MAPNODE_DESTRUCT
#undef MAPNODE_COPY
#undef MAPNODE_GRANDPARENT
#undef MAPNODE_SIBLING
#undef MAPNODE_UNCLE
#undef MAPNODE_NUM_CHILDREN
#undef MAPNODE_NEXT
#undef MAPNODE_PREVIOUS
#undef MAPNODE_ROTATE_LEFT
#undef MAPNODE_ROTATE_RIGHT
#undef MAP_TYPE
#undef MAP_INIT
#undef MAP_DESTRUCT
#undef MAP_COPY
#undef MAP_FIND
#undef MAP_FIND_POS
#undef MAP_GET
#undef MAP_SET
#undef MAP_LEFTMOST
#undef MAP_RIGHTMOST
#undef MAP_INSERT
#undef MAP_TO_ARRAY

#undef T
#undef TT
#undef U
#undef UU
