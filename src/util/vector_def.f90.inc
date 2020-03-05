#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define VECTOR_TYPE          CONCAT(vector,T)
#define VECTOR_INIT          CONCAT3(vector,T,init)
#define VECTOR_DESTRUCT      CONCAT3(vector,T,destruct)
#define VECTOR_RESET         CONCAT3(vector,T,reset)
#define VECTOR_RESERVE       CONCAT3(vector,T,reserve)
#define VECTOR_RESIZE        CONCAT3(vector,T,resize)
#define VECTOR_FRONT         CONCAT3(vector,T,front)
#define VECTOR_BACK          CONCAT3(vector,T,back)
#define VECTOR_PUSH_ELEM     CONCAT3(vector,T,push_elem)
#define VECTOR_PUSH_ELEMS    CONCAT3(vector,T,push_elems)
#define VECTOR_POP           CONCAT3(vector,T,pop)
#define VECTOR_SHRINK        CONCAT3(vector,T,shrink)
#define VECTOR_TO_ARRAY      CONCAT3(vector,T,to_array)
#define VECTOR_FROM_ARRAY    CONCAT3(vector,T,from_array)

type VECTOR_TYPE
  TT, allocatable :: d(:)
    !! data array
  integer         :: n = 0
    !! size of vector
contains
  procedure :: init       => VECTOR_INIT
  procedure :: destruct   => VECTOR_DESTRUCT
  procedure :: reset      => VECTOR_RESET
  procedure :: reserve    => VECTOR_RESERVE
  procedure :: resize     => VECTOR_RESIZE
  procedure :: front      => VECTOR_FRONT
  procedure :: back       => VECTOR_BACK
  generic   :: push       => VECTOR_PUSH_ELEM, VECTOR_PUSH_ELEMS
  procedure ::               VECTOR_PUSH_ELEM       ! append element to end of vector
  procedure ::               VECTOR_PUSH_ELEMS      ! append elements/array to end of vector
  procedure :: pop        => VECTOR_POP             ! removes last element
  procedure :: shrink     => VECTOR_SHRINK
  procedure :: to_array   => VECTOR_TO_ARRAY        ! returns all data in an array
  procedure :: from_array => VECTOR_FROM_ARRAY      ! sets all data from an array
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
#undef VECTOR_INIT
#undef VECTOR_DESTRUCT
#undef VECTOR_RESET
#undef VECTOR_RESERVE
#undef VECTOR_RESIZE
#undef VECTOR_FRONT
#undef VECTOR_BACK
#undef VECTOR_PUSH_ELEM
#undef VECTOR_PUSH_ELEMS
#undef VECTOR_POP
#undef VECTOR_SHRINK
#undef VECTOR_TO_ARRAY
#undef VECTOR_FROM_ARRAY