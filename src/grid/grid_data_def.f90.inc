#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define GRID_DATA_TYPE         CONCAT(grid_data,T)
#define GRID_DATA1_TYPE        CONCAT(grid_data1,T)
#define GRID_DATA2_TYPE        CONCAT(grid_data2,T)
#define GRID_DATA3_TYPE        CONCAT(grid_data3,T)
#define GRID_DATA4_TYPE        CONCAT(grid_data4,T)
#define GRID_DATA5_TYPE        CONCAT(grid_data5,T)
#define GRID_DATA6_TYPE        CONCAT(grid_data6,T)
#define GRID_DATA7_TYPE        CONCAT(grid_data7,T)
#define GRID_DATA8_TYPE        CONCAT(grid_data8,T)

#define GRID_DATA_INIT         CONCAT3(grid_data,T,init)
#define GRID_DATA_DESTRUCT     CONCAT3(grid_data,T,destruct)
#define GRID_DATA_GET_PTR1     CONCAT3(grid_data,T,get_ptr1)
#define GRID_DATA_GET_PTR2     CONCAT3(grid_data,T,get_ptr2)
#define GRID_DATA_GET_PTR3     CONCAT3(grid_data,T,get_ptr3)
#define GRID_DATA_GET_PTR4     CONCAT3(grid_data,T,get_ptr4)
#define GRID_DATA_GET_PTR5     CONCAT3(grid_data,T,get_ptr5)
#define GRID_DATA_GET_PTR6     CONCAT3(grid_data,T,get_ptr6)
#define GRID_DATA_GET_PTR7     CONCAT3(grid_data,T,get_ptr7)
#define GRID_DATA_GET_PTR8     CONCAT3(grid_data,T,get_ptr8)
#define GRID_DATA_GET_POINT    CONCAT3(grid_data,T,get_point)
#define GRID_DATA_GET_ALL      CONCAT3(grid_data,T,get_all)
#define GRID_DATA_SET_POINT    CONCAT3(grid_data,T,set_point)
#define GRID_DATA_SET_ALL      CONCAT3(grid_data,T,set_all)
#define GRID_DATA_UPDATE_POINT CONCAT3(grid_data,T,update_point)
#define GRID_DATA_UPDATE_ALL   CONCAT3(grid_data,T,update_all)

#define ALLOCATE_GRID_DATA     CONCAT(allocate_grid_data,T)

interface allocate_grid_data
  module procedure :: ALLOCATE_GRID_DATA
end interface

type, abstract :: GRID_DATA_TYPE
  !! grid based array with bounds check
  integer, allocatable :: idx_bnd(:)
    !! index bounds
  integer              :: n
    !! total number of values
contains
  procedure :: init     => GRID_DATA_INIT
  procedure :: destruct => GRID_DATA_DESTRUCT

  procedure :: get_ptr1 => GRID_DATA_GET_PTR1
  procedure :: get_ptr2 => GRID_DATA_GET_PTR2
  procedure :: get_ptr3 => GRID_DATA_GET_PTR3
  procedure :: get_ptr4 => GRID_DATA_GET_PTR4
  procedure :: get_ptr5 => GRID_DATA_GET_PTR5
  procedure :: get_ptr6 => GRID_DATA_GET_PTR6
  procedure :: get_ptr7 => GRID_DATA_GET_PTR7
  procedure :: get_ptr8 => GRID_DATA_GET_PTR8

  procedure :: GRID_DATA_GET_POINT
  procedure :: GRID_DATA_GET_ALL
  generic   :: get => GRID_DATA_GET_POINT, GRID_DATA_GET_ALL

  procedure :: GRID_DATA_SET_POINT
  procedure :: GRID_DATA_SET_ALL
  generic   :: set => GRID_DATA_SET_POINT, GRID_DATA_SET_ALL

#ifndef TLOG
  procedure :: GRID_DATA_UPDATE_POINT
  procedure :: GRID_DATA_UPDATE_ALL
  generic   :: update => GRID_DATA_UPDATE_POINT, GRID_DATA_UPDATE_ALL
#endif
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA1_TYPE
  TT, allocatable :: data(:)
    !! 1D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA2_TYPE
  TT, allocatable :: data(:,:)
    !! 2D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA3_TYPE
  TT, allocatable :: data(:,:,:)
    !! 3D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA4_TYPE
  TT, allocatable :: data(:,:,:,:)
    !! 4D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA5_TYPE
  TT, allocatable :: data(:,:,:,:,:)
    !! 5D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA6_TYPE
  TT, allocatable :: data(:,:,:,:,:,:)
    !! 6D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA7_TYPE
  TT, allocatable :: data(:,:,:,:,:,:,:)
    !! 7D data
end type

type, extends(GRID_DATA_TYPE) :: GRID_DATA8_TYPE
  TT, allocatable :: data(:,:,:,:,:,:,:,:)
    !! 8D data
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
#undef GRID_DATA1_TYPE
#undef GRID_DATA2_TYPE
#undef GRID_DATA3_TYPE
#undef GRID_DATA4_TYPE
#undef GRID_DATA5_TYPE
#undef GRID_DATA6_TYPE
#undef GRID_DATA7_TYPE
#undef GRID_DATA8_TYPE

#undef GRID_DATA_INIT
#undef GRID_DATA_DESTRUCT
#undef GRID_DATA_GET_PTR1
#undef GRID_DATA_GET_PTR2
#undef GRID_DATA_GET_PTR3
#undef GRID_DATA_GET_PTR4
#undef GRID_DATA_GET_PTR5
#undef GRID_DATA_GET_PTR6
#undef GRID_DATA_GET_PTR7
#undef GRID_DATA_GET_PTR8
#undef GRID_DATA_GET_POINT
#undef GRID_DATA_GET_ALL
#undef GRID_DATA_SET_POINT
#undef GRID_DATA_SET_ALL
#undef GRID_DATA_UPDATE_POINT
#undef GRID_DATA_UPDATE_ALL

#undef ALLOCATE_GRID_DATA
