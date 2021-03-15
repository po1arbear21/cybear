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

#define HESSENBERG_INIT      CONCAT3(hessenberg,T,init)
#define HESSENBERG_DESTRUCT  CONCAT3(hessenberg,T,destruct)
#define HESSENBERG_RESET     CONCAT3(hessenberg,T,reset)
#define HESSENBERG_SCALE     CONCAT3(hessenberg,T,scale)
#define HESSENBERG_MUL_VEC   CONCAT3(hessenberg,T,mul_vec)
#define HESSENBERG_MUL_MAT   CONCAT3(hessenberg,T,mul_mat)
#define HESSENBERG_FACTORIZE CONCAT3(hessenberg,T,factorize)
#define HESSENBERG_SOLVE_VEC CONCAT3(hessenberg,T,solve_vec)
#define HESSENBERG_SOLVE_MAT CONCAT3(hessenberg,T,solve_mat)
#define HESSENBERG_IS_TRIANG CONCAT3(hessenberg,T,is_triang)

type, extends(MATRIX_TYPE) :: HESSENBERG_TYPE
  ! upper or lower hessenberg matrix
  logical :: upper

  ! unpacked data
  TT, allocatable :: d(:,:)

  ! factorization
  TT,      allocatable :: f(:,:)
  integer, allocatable :: ipiv(:)
contains
  procedure :: HESSENBERG_INIT
  generic   :: init      => HESSENBERG_INIT
  procedure :: destruct  => HESSENBERG_DESTRUCT
  procedure :: reset     => HESSENBERG_RESET
  procedure :: scale     => HESSENBERG_SCALE
  procedure :: mul_vec   => HESSENBERG_MUL_VEC
  procedure :: mul_mat   => HESSENBERG_MUL_MAT
  procedure :: factorize => HESSENBERG_FACTORIZE
  procedure :: solve_vec => HESSENBERG_SOLVE_VEC
  procedure :: solve_mat => HESSENBERG_SOLVE_MAT
  procedure :: is_triang => HESSENBERG_IS_TRIANG
end type

type HESSENBERG_PTR_TYPE
  type(HESSENBERG_TYPE), pointer :: p => null()
end type

interface
  module subroutine HESSENBERG_INIT(this, nrows, upper)
    class(HESSENBERG_TYPE), intent(out) :: this
    integer,                intent(in)  :: nrows
    logical,                intent(in)  :: upper
  end subroutine

  module subroutine HESSENBERG_DESTRUCT(this)
    class(HESSENBERG_TYPE), intent(inout) :: this
  end subroutine

  module subroutine HESSENBERG_RESET(this, only_factorization)
    class(HESSENBERG_TYPE), intent(inout) :: this
    logical, optional,      intent(in)    :: only_factorization
  end subroutine

  module subroutine HESSENBERG_SCALE(this, fact)
    class(HESSENBERG_TYPE), intent(inout) :: this
    TT,                     intent(in)    :: fact
  end subroutine

  module subroutine HESSENBERG_MUL_VEC(this, x, y, fact_y, trans)
    class(HESSENBERG_TYPE), intent(in)    :: this
    TT,                     intent(in)    :: x(:)
    TT,                     intent(inout) :: y(:)
    TT,        optional,    intent(in)    :: fact_y
    character, optional,    intent(in)    :: trans
  end subroutine

  module subroutine HESSENBERG_MUL_MAT(this, x, y, fact_y, trans)
    class(HESSENBERG_TYPE), intent(in)    :: this
    TT,                     intent(in)    :: x(:,:)
    TT,                     intent(inout) :: y(:,:)
    TT,        optional,    intent(in)    :: fact_y
    character, optional,    intent(in)    :: trans
  end subroutine

  module subroutine HESSENBERG_FACTORIZE(this)
    class(HESSENBERG_TYPE), intent(inout) :: this
  end subroutine

  module subroutine HESSENBERG_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(HESSENBERG_TYPE), intent(in)  :: this
    TT,                     intent(in)  :: rhs(:)
    TT,                     intent(out) :: x(:)
    character, optional,    intent(in)  :: trans
    logical,   optional,    intent(in)  :: iter_refine
  end subroutine

  module subroutine HESSENBERG_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(HESSENBERG_TYPE), intent(in)  :: this
    TT,                     intent(in)  :: rhs(:,:)
    TT,                     intent(out) :: x(:,:)
    character, optional,    intent(in)  :: trans
    logical,   optional,    intent(in)  :: iter_refine
  end subroutine

  logical module function HESSENBERG_IS_TRIANG(this, abs_tol)
    class(HESSENBERG_TYPE), intent(in) :: this
    real, optional,         intent(in) :: abs_tol
  end function
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

#undef HESSENBERG_INIT
#undef HESSENBERG_DESTRUCT
#undef HESSENBERG_RESET
#undef HESSENBERG_SCALE
#undef HESSENBERG_MUL_VEC
#undef HESSENBERG_MUL_MAT
#undef HESSENBERG_FACTORIZE
#undef HESSENBERG_SOLVE_VEC
#undef HESSENBERG_SOLVE_MAT
#undef HESSENBERG_IS_TRIANG
