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

#define TRIANG_INIT          CONCAT3(triang,T,init)
#define TRIANG_DESTRUCT      CONCAT3(triang,T,destruct)
#define TRIANG_RESET         CONCAT3(triang,T,reset)
#define TRIANG_SCALE         CONCAT3(triang,T,scale)
#define TRIANG_MUL_VEC       CONCAT3(triang,T,mul_vec)
#define TRIANG_MUL_MAT_ADD   CONCAT3(triang,T,mul_mat_add)
#define TRIANG_MUL_MAT_SIDES CONCAT3(triang,T,mul_mat_sides)
#define TRIANG_FACTORIZE     CONCAT3(triang,T,factorize)
#define TRIANG_SOLVE_VEC     CONCAT3(triang,T,solve_vec)
#define TRIANG_SOLVE_MAT     CONCAT3(triang,T,solve_mat)

type, extends(MATRIX_TYPE) :: TRIANG_TYPE
  logical         :: upper
    !! upper or lower triangular matrix

  TT, allocatable :: d(:,:)
    !! unpacked data
contains
  procedure :: TRIANG_INIT
  generic   :: init      => TRIANG_INIT
  procedure :: destruct  => TRIANG_DESTRUCT
  procedure :: reset     => TRIANG_RESET
  procedure :: scale     => TRIANG_SCALE
  procedure :: mul_vec   => TRIANG_MUL_VEC
  procedure :: mul_mat   => TRIANG_MUL_MAT_ADD          ! mat2 <- thisMat * mat1 + fact * mat2
  procedure :: mul_mat2  => TRIANG_MUL_MAT_SIDES        ! mat2 <- fact * thisMat * mat1
  procedure :: factorize => TRIANG_FACTORIZE            ! unneeded => throws error
  procedure :: solve_vec => TRIANG_SOLVE_VEC
  procedure :: solve_mat => TRIANG_SOLVE_MAT
end type

type TRIANG_PTR_TYPE
  type(TRIANG_TYPE), pointer :: p => null()
end type

interface
  module subroutine TRIANG_INIT(this, n, upper)
    class(TRIANG_TYPE), intent(out) :: this
    integer,            intent(in)  :: n
    logical,            intent(in)  :: upper
  end subroutine

  module subroutine TRIANG_DESTRUCT(this)
    class(TRIANG_TYPE), intent(inout) :: this
  end subroutine

  module subroutine TRIANG_RESET(this, only_factorization)
    class(TRIANG_TYPE), intent(inout) :: this
    logical, optional,  intent(in)    :: only_factorization
  end subroutine

  module subroutine TRIANG_SCALE(this, fact)
    class(TRIANG_TYPE), intent(inout) :: this
    TT,                 intent(in)    :: fact
  end subroutine

  module subroutine TRIANG_MUL_VEC(this, x, y, fact_y, trans)
    class(TRIANG_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine TRIANG_MUL_MAT_ADD(this, x, y, fact_y, trans)
    class(TRIANG_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine TRIANG_MUL_MAT_SIDES(this, y, fact, trans, right)
    class(TRIANG_TYPE),  intent(in)    :: this
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact
    character, optional, intent(in)    :: trans
    logical,   optional, intent(in)    :: right
  end subroutine

  module subroutine TRIANG_FACTORIZE(this)
    class(TRIANG_TYPE), intent(inout) :: this
  end subroutine

  module subroutine TRIANG_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(TRIANG_TYPE),  intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine TRIANG_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(TRIANG_TYPE),  intent(in)  :: this
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

#undef TRIANG_INIT
#undef TRIANG_DESTRUCT
#undef TRIANG_RESET
#undef TRIANG_SCALE
#undef TRIANG_MUL_VEC
#undef TRIANG_MUL_MAT_ADD
#undef TRIANG_MUL_MAT_SIDES
#undef TRIANG_FACTORIZE
#undef TRIANG_SOLVE_VEC
#undef TRIANG_SOLVE_MAT
