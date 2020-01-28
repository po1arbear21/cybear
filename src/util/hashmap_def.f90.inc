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

type HASHMAP_TYPE
  type(vector_int)     :: keys
  type(VECTOR_TYPE)    :: values
  integer, allocatable :: table(:)
contains
  procedure :: init      => HASHMAP_INIT
  procedure :: resize    => HASHMAP_RESIZE
  procedure :: get_index => HASHMAP_GET_INDEX
  procedure :: get       => HASHMAP_GET
  procedure :: insert    => HASHMAP_INSERT
end type

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
