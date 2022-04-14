subroutine mapnode_{}T{}_{}U{}_init(this, key, value, color, left, right, parent)
  !! initialize map-node
  class(mapnode_{}T{}_{}U),         intent(out) :: this
  m4_type(T),                       intent(in)  :: key
    !! element key
  m4_type(U),                       intent(in)  :: value
    !! element value
  integer,                          intent(in)  :: color
    !! node color
  type(mapnode_{}T{}_{}U), pointer, intent(in)  :: left
    !! pointer to left child
  type(mapnode_{}T{}_{}U), pointer, intent(in)  :: right
    !! pointer to right child
  type(mapnode_{}T{}_{}U), pointer, intent(in)  :: parent
    !! pointer to parent

  this%key    =  key
  this%value  =  value
  this%color  =  color
  this%left   => left
  this%right  => right
  this%parent => parent
end subroutine

recursive subroutine mapnode_{}T{}_{}U{}_destruct(this)
  !! destruct map-node (and children)
  class(mapnode_{}T{}_{}U), intent(inout) :: this

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

recursive subroutine mapnode_{}T{}_{}U{}_copy(this, n)
  !! copy map-node from source
  class(mapnode_{}T{}_{}U), target, intent(inout) :: this
  type(mapnode_{}T{}_{}U),          intent(in)    :: n
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

function mapnode_{}T{}_{}U{}_grandparent(this) result(p)
  !! get map-node grandparent (parent of parent)
  class(mapnode_{}T{}_{}U), intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer     :: p
    !! return pointer to grandparent

  if (.not. associated(this%parent)) then
    nullify(p)
  else
    p => this%parent%parent
  end if
end function

function mapnode_{}T{}_{}U{}_sibling(this) result(p)
  !! get map-node sibling (the other child of this%parent)
  class(mapnode_{}T{}_{}U), target, intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer             :: p
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

function mapnode_{}T{}_{}U{}_uncle(this) result(p)
  !! get map-node uncle (the other child of this%parent%parent)
  class(mapnode_{}T{}_{}U), intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer     :: p
    !! return pointer to uncle

  if (.not. associated(this%parent)) then
    nullify(p)
  else
    p => this%parent%sibling()
  end if
end function

function mapnode_{}T{}_{}U{}_num_children(this) result(n)
  !! get number of children
  class(mapnode_{}T{}_{}U), intent(in) :: this
  integer                              :: n
    !! return number of children

  n = 0
  if (associated(this%left)) then
    n = n + 1
  end if
  if (associated(this%right)) then
    n = n + 1
  end if
end function

function mapnode_{}T{}_{}U{}_next(this) result(next)
  !! get next map-node
  class(mapnode_{}T{}_{}U), target, intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer             :: next
    !! return next map-node

  type(mapnode_{}T{}_{}U), pointer :: p

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

function mapnode_{}T{}_{}U{}_previous(this) result(prev)
  !! get previous map-node
  class(mapnode_{}T{}_{}U), target, intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer             :: prev
    !! return previous map-node

  type(mapnode_{}T{}_{}U), pointer :: p

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

subroutine mapnode_{}T{}_{}U{}_rotate_left(this)
  !! rotate tree left
  class(mapnode_{}T{}_{}U), target, intent(inout) :: this

  type(mapnode_{}T{}_{}U), pointer :: n, g

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

subroutine mapnode_{}T{}_{}U{}_rotate_right(this)
  !! rotate tree left
  class(mapnode_{}T{}_{}U), target, intent(inout) :: this

  type(mapnode_{}T{}_{}U), pointer :: p, g

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

subroutine map_{}T{}_{}U{}_init(this)
  !! initialize map
  class(map_{}T{}_{}U), intent(out) :: this

  ! no elements
  this%n = 0
end subroutine

subroutine map_{}T{}_{}U{}_destruct(this)
  !! destruct map
  class(map_{}T{}_{}U), intent(inout) :: this

  ! destruct root node
  if (associated(this%root)) then
    call this%root%destruct()
    deallocate(this%root)
  end if

  ! clear number of elements
  this%n = 0
end subroutine

subroutine map_{}T{}_{}U{}_copy(this, m)
  !! copy map
  class(map_{}T{}_{}U), intent(inout) :: this
    !! target
  type(map_{}T{}_{}U),  intent(in)    :: m
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

function map_{}T{}_{}U{}_find(this, key) result(node)
  !! find using the key in map
  class(map_{}T{}_{}U), intent(in) :: this
  m4_type(T),           intent(in) :: key
    !! key to lookup
  type(mapnode_{}T{}_{}U), pointer :: node
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

function map_{}T{}_{}U{}_find_pos(this, key) result(node)
  ! find nearest key in map
  class(map_{}T{}_{}U),    intent(in) :: this
  m4_type(T),              intent(in) :: key
    !! key to lookup
  type(mapnode_{}T{}_{}U), pointer    :: node
    !! return node which is nearest to key

  type(mapnode_{}T{}_{}U), pointer :: tmp

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

function map_{}T{}_{}U{}_get(this, key) result(val)
  !! lookup key in map and return value
  class(map_{}T{}_{}U), intent(in) :: this
  m4_type(T),           intent(in) :: key
    !! key to lookup
  m4_type(U)                       :: val
    !! return value

  type(mapnode_{}T{}_{}U), pointer :: n

  n => this%find(key)
  val = n%value
end function

subroutine map_{}T{}_{}U{}_set(this, key, value)
  !! insert new key-value-pair into map or update existing one
  class(map_{}T{}_{}U), intent(inout) :: this
  m4_type(T),           intent(in)    :: key
  m4_type(U),           intent(in)    :: value

  type(mapnode_{}T{}_{}U), pointer :: node
  logical                     :: status

  call this%insert(key, value, node = node, status = status)

  if (.not. status) node%value = value
end subroutine

function map_{}T{}_{}U{}_leftmost(this) result(node)
  !! get leftmost node in map
  class(map_{}T{}_{}U), intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer :: node
    !! return leftmost node

  node => this%root
  if (.not. associated(node)) return
  do while (associated(node%left))
    node => node%left
  end do
end function

function map_{}T{}_{}U{}_rightmost(this) result(node)
  !! get rightmost node in map
  class(map_{}T{}_{}U), intent(in) :: this
  type(mapnode_{}T{}_{}U), pointer :: node
    !! return rightmost node

  node => this%root
  if (.not. associated(node)) return
  do while (associated(node%right))
    node => node%right
  end do
end function

subroutine map_{}T{}_{}U{}_insert(this, key, value, node, status)
  !! insert key-value-pair into map
  class(map_{}T{}_{}U),                       intent(inout) :: this
  m4_type(T),                                 intent(in)    :: key
    !! key to insert
  m4_type(U),                                 intent(in)    :: value
    !! value to insert
  type(mapnode_{}T{}_{}U), optional, pointer, intent(out)   :: node
    !! optional output pointer to newly created node (or to already existing one with the same key)
  logical,                 optional,          intent(out)   :: status
    !! optional output status of insertion: success => true, fail => false

  integer, parameter :: RB_RED = 1, RB_BLACK = 2

  type(mapnode_{}T{}_{}U), pointer :: n, p, g, u

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

subroutine map_{}T{}_{}U{}_to_array(this, keys, values)
  !! convert map to array
  class(map_{}T{}_{}U), intent(in)  :: this
  m4_type(T), optional, intent(out) :: keys(:)
    !! optional output array with keys
  m4_type(U), optional, intent(out) :: values(:)
    !! optional output array with values

  integer :: i
  type(mapnode_{}T{}_{}U), pointer :: n

  if (present(keys)) then
    m4_assert(size(keys) == this%n)
  end if

  if (present(values)) then
    m4_assert(size(values) == this%n)
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

m4_undefine({T},{U})
