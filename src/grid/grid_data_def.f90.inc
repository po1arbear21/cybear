#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define GRID_DATA_TYPE       CONCAT(grid_data,T)
#define GRID_DATA_INIT       CONCAT3(grid_data,T,init)
#define GRID_DATA_DESTRUCT   CONCAT3(grid_data,T,destruct)
#define GRID_DATA_GET        CONCAT3(grid_data,T,get)
#define GRID_DATA_SET        CONCAT3(grid_data,T,set)

type GRID_DATA_TYPE
  !! grid based array with bounds check
  integer, allocatable :: idx_bnd(:)
    !! index bounds
  TT,      allocatable :: data(:)
    !! data array
contains
  procedure :: init     => GRID_DATA_INIT
  procedure :: destruct => GRID_DATA_DESTRUCT
  procedure :: get      => GRID_DATA_GET
  procedure :: set      => GRID_DATA_SET
end type

#undef T
#undef TT
#undef TLOG

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef GRID_DATA_TYPE
#undef GRID_DATA_INIT
#undef GRID_DATA_DESTRUCT
#undef GRID_DATA_GET
#undef GRID_DATA_SET
