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

type MAPNODE_TYPE
  TT                          :: key
    !! element key
  UU                          :: value
    !! element value
  integer                     :: color
    !! node color
  type(MAPNODE_TYPE), pointer :: left   => null()
    !! pointer to left child
  type(MAPNODE_TYPE), pointer :: right  => null()
    !! pointer to right child
  type(MAPNODE_TYPE), pointer :: parent => null()
    !! pointer to parent
contains
  procedure :: init         => MAPNODE_INIT
  procedure :: destruct     => MAPNODE_DESTRUCT
  procedure :: copy         => MAPNODE_COPY
  procedure :: grandparent  => MAPNODE_GRANDPARENT
  procedure :: sibling      => MAPNODE_SIBLING
  procedure :: uncle        => MAPNODE_UNCLE
  procedure :: num_children => MAPNODE_NUM_CHILDREN
  procedure :: next         => MAPNODE_NEXT
  procedure :: previous     => MAPNODE_PREVIOUS
  procedure :: rotate_left  => MAPNODE_ROTATE_LEFT
  procedure :: rotate_right => MAPNODE_ROTATE_RIGHT
end type

type MAP_TYPE
  type(MAPNODE_TYPE), pointer :: root => null()
    !! pointer to root node
  integer                     :: n
    !! number of elements
contains
  procedure :: init      => MAP_INIT
  procedure :: destruct  => MAP_DESTRUCT
  procedure :: copy      => MAP_COPY
  procedure :: find      => MAP_FIND
  procedure :: find_pos  => MAP_FIND_POS
  procedure :: get       => MAP_GET
  procedure :: set       => MAP_SET
  procedure :: leftmost  => MAP_LEFTMOST
  procedure :: rightmost => MAP_RIGHTMOST
  procedure :: insert    => MAP_INSERT
  procedure :: to_array  => MAP_TO_ARRAY
end type

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
