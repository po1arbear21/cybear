#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_PTR_TYPE        CONCAT(band_ptr,T)
#define BLOCK_TYPE           CONCAT(block,T)
#define BLOCK_PTR_TYPE       CONCAT(block_ptr,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define DENSE_PTR_TYPE       CONCAT(dense_ptr,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define HESSENBERG_PTR_TYPE  CONCAT(hessenberg_ptr,T)
#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_PTR_TYPE      CONCAT(matrix_ptr,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPARSE_PTR_TYPE      CONCAT(sparse_ptr,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define TRIANG_TYPE          CONCAT(triang,T)
#define TRIANG_PTR_TYPE      CONCAT(triang_ptr,T)

#define MATRIX_INIT          CONCAT3(matrix,T,init)
#define MATRIX_DESTRUCT      CONCAT3(matrix,T,destruct)
#define MATRIX_RESET         CONCAT3(matrix,T,reset)
#define MATRIX_SCALE         CONCAT3(matrix,T,scale)
#define MATRIX_MUL_VEC       CONCAT3(matrix,T,mul_vec)
#define MATRIX_MUL_MAT       CONCAT3(matrix,T,mul_mat)
#define MATRIX_FACTORIZE     CONCAT3(matrix,T,factorize)
#define MATRIX_SOLVE_VEC     CONCAT3(matrix,T,solve_vec)
#define MATRIX_SOLVE_MAT     CONCAT3(matrix,T,solve_mat)

type, abstract :: MATRIX_TYPE
  !! Base type for matrices.
  !! Can be extended if more specialized matrices are needed besides the built-in types.
  integer :: nrows
    !! Number of rows
  integer :: ncols
    !! Number of columns
  logical :: factorized
    !! Is matrix factorized and ready to solve?
contains
  procedure :: MATRIX_INIT
  generic   :: init => MATRIX_INIT

  procedure(MATRIX_DESTRUCT),  deferred :: destruct
  procedure(MATRIX_RESET),     deferred :: reset
  procedure(MATRIX_SCALE),     deferred :: scale
  procedure(MATRIX_MUL_VEC),   deferred :: mul_vec
  procedure(MATRIX_MUL_MAT),   deferred :: mul_mat
  procedure(MATRIX_FACTORIZE), deferred :: factorize
  procedure(MATRIX_SOLVE_VEC), deferred :: solve_vec
  procedure(MATRIX_SOLVE_MAT), deferred :: solve_mat
end type

type MATRIX_PTR_TYPE
  class(MATRIX_TYPE), pointer :: p => null()
end type

interface
  module subroutine MATRIX_DESTRUCT(this)
    class(MATRIX_TYPE), intent(inout) :: this
  end subroutine

  module subroutine MATRIX_RESET(this, only_factorization)
    class(MATRIX_TYPE), intent(inout) :: this
    logical, optional,  intent(in)    :: only_factorization
  end subroutine

  module subroutine MATRIX_SCALE(this, fact)
    class(MATRIX_TYPE), intent(inout) :: this
    TT,                 intent(in)    :: fact
  end subroutine

  module subroutine MATRIX_MUL_VEC(this, x, y, fact_y, trans)
    class(MATRIX_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine MATRIX_MUL_MAT(this, x, y, fact_y, trans)
    class(MATRIX_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine MATRIX_FACTORIZE(this)
    class(MATRIX_TYPE), intent(inout) :: this
  end subroutine

  module subroutine MATRIX_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(MATRIX_TYPE),  intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine MATRIX_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(MATRIX_TYPE),  intent(in)  :: this
    TT,                  intent(in)  :: rhs(:,:)
    TT,                  intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine
end interface

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef ARRAY_TYPE
#undef VECTOR_TYPE

#undef BAND_TYPE
#undef BAND_PTR_TYPE
#undef BLOCK_TYPE
#undef BLOCK_PTR_TYPE
#undef DENSE_TYPE
#undef DENSE_PTR_TYPE
#undef HESSENBERG_TYPE
#undef HESSENBERG_PTR_TYPE
#undef MATRIX_TYPE
#undef MATRIX_PTR_TYPE
#undef SPARSE_TYPE
#undef SPARSE_PTR_TYPE
#undef SPBUILD_TYPE
#undef TRIANG_TYPE
#undef TRIANG_PTR_TYPE

#undef MATRIX_INIT
#undef MATRIX_DESTRUCT
#undef MATRIX_RESET
#undef MATRIX_SCALE
#undef MATRIX_MUL_VEC
#undef MATRIX_MUL_MAT
#undef MATRIX_FACTORIZE
#undef MATRIX_SOLVE_VEC
#undef MATRIX_SOLVE_MAT
