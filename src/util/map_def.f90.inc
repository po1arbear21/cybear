type mapnode_{}T{}_{}U
  m4_type(T) :: key
    !! element key
  m4_type(U) :: value
    !! element value
  integer    :: color
    !! node color

  type(mapnode_{}T{}_{}U), pointer :: left   => null()
    !! pointer to left child
  type(mapnode_{}T{}_{}U), pointer :: right  => null()
    !! pointer to right child
  type(mapnode_{}T{}_{}U), pointer :: parent => null()
    !! pointer to parent
contains
  procedure :: init         => mapnode_{}T{}_{}U{}_init
  procedure :: copy         => mapnode_{}T{}_{}U{}_copy
  procedure :: grandparent  => mapnode_{}T{}_{}U{}_grandparent
  procedure :: sibling      => mapnode_{}T{}_{}U{}_sibling
  procedure :: uncle        => mapnode_{}T{}_{}U{}_uncle
  procedure :: num_children => mapnode_{}T{}_{}U{}_num_children
  procedure :: next         => mapnode_{}T{}_{}U{}_next
  procedure :: previous     => mapnode_{}T{}_{}U{}_previous
  procedure :: rotate_left  => mapnode_{}T{}_{}U{}_rotate_left
  procedure :: rotate_right => mapnode_{}T{}_{}U{}_rotate_right

  final :: mapnode_{}T{}_{}U{}_destruct
end type

type map_{}T{}_{}U
  type(mapnode_{}T{}_{}U), pointer :: root => null()
    !! pointer to root node
  integer :: n
    !! number of elements
contains
  procedure :: init     => map_{}T{}_{}U{}_init
  procedure :: copy     => map_{}T{}_{}U{}_copy
  procedure :: find     => map_{}T{}_{}U{}_find
  procedure :: get      => map_{}T{}_{}U{}_get
  procedure :: set      => map_{}T{}_{}U{}_set
  procedure :: insert   => map_{}T{}_{}U{}_insert
  procedure :: to_array => map_{}T{}_{}U{}_to_array

  procedure, private :: find_pos  => map_{}T{}_{}U{}_find_pos
  procedure, private :: leftmost  => map_{}T{}_{}U{}_leftmost
  procedure, private :: rightmost => map_{}T{}_{}U{}_rightmost

  final :: map_{}T{}_{}U{}_destruct
end type

m4_undefine({T},{U})
