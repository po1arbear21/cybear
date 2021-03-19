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

type HASHMAP_TYPE
  !! hashmap using integer (array) keys
  type(vector_int)     :: keys
    !! key data
  integer              :: keysize
    !! how many entries in keys%d yield a key
  type(VECTOR_VALUE)   :: values
    !! value data
  integer, allocatable :: table(:)
    !! lookup table (used internally)
contains
  procedure :: init     => HASHMAP_INIT
  procedure :: destruct => HASHMAP_DESTRUCT
  procedure :: reset    => HASHMAP_RESET
  procedure :: reserve  => HASHMAP_RESERVE
  generic   :: get      => HASHMAP_GET_KEY_1, HASHMAP_GET_KEY_N
  generic   :: set      => HASHMAP_SET_KEY_1, HASHMAP_SET_KEY_N

  procedure, private :: lookup => HASHMAP_LOOKUP
  procedure, private :: HASHMAP_GET_KEY_1
  procedure, private :: HASHMAP_GET_KEY_N
  procedure, private :: HASHMAP_SET_KEY_1
  procedure, private :: HASHMAP_SET_KEY_N
end type

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
