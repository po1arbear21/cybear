#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DEQUE_TYPE           CONCAT(deque,T)
#define DEQUE_INIT           CONCAT3(deque,T,init)
#define DEQUE_DESTRUCT       CONCAT3(deque,T,destruct)
#define DEQUE_RESET          CONCAT3(deque,T,reset)
#define DEQUE_RESERVE        CONCAT3(deque,T,reserve)
#define DEQUE_RESIZE         CONCAT3(deque,T,resize)
#define DEQUE_FRONT          CONCAT3(deque,T,first)
#define DEQUE_BACK           CONCAT3(deque,T,last)
#define DEQUE_PUSH_FRONT     CONCAT3(deque,T,push_front)
#define DEQUE_PUSH_BACK      CONCAT3(deque,T,push_back)
#define DEQUE_POP_FRONT      CONCAT3(deque,T,pop_front)
#define DEQUE_POP_BACK       CONCAT3(deque,T,pop_back)
#define DEQUE_SHRINK         CONCAT3(deque,T,shrink)
#define DEQUE_TO_ARRAY       CONCAT3(deque,T,to_array)
#define DEQUE_FROM_ARRAY     CONCAT3(deque,T,from_array)

type DEQUE_TYPE
  !! double-ended queue
  TT, allocatable :: d(:)
    !! data array
  integer :: n
    !! number of elements
  integer :: ifront
    !! front index
  integer :: iback
    !! back index
contains
  procedure :: init       => DEQUE_INIT
  procedure :: destruct   => DEQUE_DESTRUCT
  procedure :: reset      => DEQUE_RESET
  procedure :: reserve    => DEQUE_RESERVE
  procedure :: resize     => DEQUE_RESIZE
  procedure :: front      => DEQUE_FRONT
  procedure :: back       => DEQUE_BACK
  procedure :: push_front => DEQUE_PUSH_FRONT
  procedure :: push_back  => DEQUE_PUSH_BACK
  procedure :: pop_front  => DEQUE_POP_FRONT
  procedure :: pop_back   => DEQUE_POP_BACK
  procedure :: shrink     => DEQUE_SHRINK
  procedure :: to_array   => DEQUE_TO_ARRAY
  procedure :: from_array => DEQUE_FROM_ARRAY
end type

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef DEQUE_TYPE
#undef DEQUE_INIT
#undef DEQUE_DESTRUCT
#undef DEQUE_RESET
#undef DEQUE_RESERVE
#undef DEQUE_RESIZE
#undef DEQUE_FRONT
#undef DEQUE_BACK
#undef DEQUE_PUSH_FRONT
#undef DEQUE_PUSH_BACK
#undef DEQUE_POP_FRONT
#undef DEQUE_POP_BACK
#undef DEQUE_SHRINK
#undef DEQUE_TO_ARRAY
#undef DEQUE_FROM_ARRAY