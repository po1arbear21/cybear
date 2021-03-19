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

type HASHMAP_TYPE
  type(VECTOR_KEY)     :: keys
  type(VECTOR_VALUE)   :: values
  integer, allocatable :: table(:)
contains
  procedure :: init     => HASHMAP_INIT
  procedure :: destruct => HASHMAP_DESTRUCT
  procedure :: reset    => HASHMAP_RESET
  procedure :: reserve  => HASHMAP_RESERVE
  procedure :: lookup   => HASHMAP_LOOKUP
  procedure :: get      => HASHMAP_GET
  procedure :: set      => HASHMAP_SET
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
