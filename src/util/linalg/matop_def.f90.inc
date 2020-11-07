#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_PTR_TYPE      CONCAT(matrix_ptr,T)

#define MATOP_TYPE           CONCAT( matop,T)
#define MATOP_INIT           CONCAT3(matop,T,init)
#define MATOP_EXEC1          CONCAT3(matop,T,exec1)
#define MATOP_EXEC2          CONCAT3(matop,T,exec2)

#define SINGLE_MATOP_TYPE    CONCAT( single_matop,T)
#define SINGLE_MATOP_INIT    CONCAT3(single_matop,T,init)
#define SINGLE_MATOP_EXEC1   CONCAT3(single_matop,T,exec1)
#define SINGLE_MATOP_EXEC2   CONCAT3(single_matop,T,exec2)

#define CHAIN_MATOP_TYPE     CONCAT( chain_matop,T)
#define CHAIN_MATOP_INIT     CONCAT3(chain_matop,T,init)
#define CHAIN_MATOP_SET_MAT  CONCAT3(chain_matop,T,set_mat)
#define CHAIN_MATOP_EXEC1    CONCAT3(chain_matop,T,exec1)
#define CHAIN_MATOP_EXEC2    CONCAT3(chain_matop,T,exec2)

type, abstract :: MATOP_TYPE
  !! Base type for matrix operations on vectors

  integer :: nrows
    !! Number of rows (size of output vector)
  integer :: ncols
    !! Number of columns (size of input vector)
contains
  procedure :: MATOP_INIT
  generic   :: init  => MATOP_INIT
  procedure :: exec1 => MATOP_EXEC1
  procedure :: exec2 => MATOP_EXEC2
  generic   :: exec  => exec1, exec2
end type

type, extends(MATOP_TYPE) :: SINGLE_MATOP_TYPE
  !! matrix operation using a single matrix

  class(MATRIX_TYPE), pointer :: mat => null()
    !! pointer to matrix
  logical                     :: inv
    !! inverse flag
contains
  procedure :: SINGLE_MATOP_INIT
  generic   :: init  => SINGLE_MATOP_INIT
  procedure :: exec1 => SINGLE_MATOP_EXEC1
  procedure :: exec2 => SINGLE_MATOP_EXEC2
end type

type, extends(MATOP_TYPE) :: CHAIN_MATOP_TYPE
  !! multiple matrices chainged together

  type(MATRIX_PTR_TYPE), allocatable :: mat(:)
    !! pointers to matrices
  logical,               allocatable :: inv(:)
    !! inverse flags
contains
  procedure :: CHAIN_MATOP_INIT
  generic   :: init    => CHAIN_MATOP_INIT
  procedure :: set_mat => CHAIN_MATOP_SET_MAT
  procedure :: exec1   => CHAIN_MATOP_EXEC1
  procedure :: exec2   => CHAIN_MATOP_EXEC2
end type

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef MATRIX_TYPE
#undef MATRIX_PTR_TYPE

#undef MATOP_TYPE
#undef MATOP_INIT
#undef MATOP_EXEC1
#undef MATOP_EXEC2

#undef SINGLE_MATOP_TYPE
#undef SINGLE_MATOP_INIT
#undef SINGLE_MATOP_EXEC1
#undef SINGLE_MATOP_EXEC2

#undef CHAIN_MATOP_TYPE
#undef CHAIN_MATOP_INIT
#undef CHAIN_MATOP_SET_MAT
#undef CHAIN_MATOP_EXEC1
#undef CHAIN_MATOP_EXEC2
