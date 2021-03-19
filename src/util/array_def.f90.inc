#define PASTE(X)        X
#define PASTE2(X)       PASTE(X)_
#define CONCATHELP(X,Y) PASTE2(X)Y
#define CONCAT(X,Y)     CONCATHELP(X,Y)

#define ARRAY_TYPE      CONCAT(array,T)
#define ARRAY2_TYPE     CONCAT(array2,T)
#define ARRAY3_TYPE     CONCAT(array3,T)
#define ARRAY4_TYPE     CONCAT(array4,T)
#define ARRAY5_TYPE     CONCAT(array5,T)
#define ARRAY6_TYPE     CONCAT(array6,T)
#define ARRAY7_TYPE     CONCAT(array7,T)

#define PTR_TYPE        CONCAT(ptr,T)
#define PTR1_TYPE       CONCAT(ptr1,T)
#define PTR2_TYPE       CONCAT(ptr2,T)
#define PTR3_TYPE       CONCAT(ptr3,T)
#define PTR4_TYPE       CONCAT(ptr4,T)
#define PTR5_TYPE       CONCAT(ptr5,T)
#define PTR6_TYPE       CONCAT(ptr6,T)
#define PTR7_TYPE       CONCAT(ptr7,T)

public ARRAY_TYPE, ARRAY2_TYPE, ARRAY3_TYPE, ARRAY4_TYPE, ARRAY5_TYPE, ARRAY6_TYPE, ARRAY7_TYPE
public PTR_TYPE,   PTR2_TYPE,   PTR3_TYPE,   PTR4_TYPE,   PTR5_TYPE,   PTR6_TYPE,   PTR7_TYPE

type ARRAY_TYPE
  TT, allocatable :: d(:)
end type

type ARRAY2_TYPE
  TT, allocatable :: d(:,:)
end type

type ARRAY3_TYPE
  TT, allocatable :: d(:,:,:)
end type

type ARRAY4_TYPE
  TT, allocatable :: d(:,:,:,:)
end type

type ARRAY5_TYPE
  TT, allocatable :: d(:,:,:,:,:)
end type

type ARRAY6_TYPE
  TT, allocatable :: d(:,:,:,:,:,:)
end type

type ARRAY7_TYPE
  TT, allocatable :: d(:,:,:,:,:,:,:)
end type

type PTR_TYPE
  TT, pointer :: p => null()
end type

type PTR1_TYPE
  TT, pointer :: p(:) => null()
end type

type PTR2_TYPE
  TT, pointer :: p(:,:) => null()
end type

type PTR3_TYPE
  TT, pointer :: p(:,:,:) => null()
end type

type PTR4_TYPE
  TT, pointer :: p(:,:,:,:) => null()
end type

type PTR5_TYPE
  TT, pointer :: p(:,:,:,:,:) => null()
end type

type PTR6_TYPE
  TT, pointer :: p(:,:,:,:,:,:) => null()
end type

type PTR7_TYPE
  TT, pointer :: p(:,:,:,:,:,:,:) => null()
end type

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef ARRAY_TYPE
#undef ARRAY2_TYPE
#undef ARRAY3_TYPE
#undef ARRAY4_TYPE
#undef ARRAY5_TYPE
#undef ARRAY6_TYPE
#undef ARRAY7_TYPE

#undef PTR_TYPE
#undef PTR1_TYPE
#undef PTR2_TYPE
#undef PTR3_TYPE
#undef PTR4_TYPE
#undef PTR5_TYPE
#undef PTR6_TYPE
#undef PTR7_TYPE
