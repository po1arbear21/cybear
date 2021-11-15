#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define DENSE_TYPE           CONCAT(dense,T)
#define DENSE_PTR_TYPE       CONCAT(dense_ptr,T)
#define MATRIX_TYPE          CONCAT(matrix,T)

#define DENSE_INIT_BY_DATA   CONCAT3(dense,T,init_by_data)
#define DENSE_INIT_BY_SIZE   CONCAT3(dense,T,init_by_size)
#define DENSE_DESTRUCT       CONCAT3(dense,T,destruct)
#define DENSE_RESET          CONCAT3(dense,T,reset)
#define DENSE_SCALE          CONCAT3(dense,T,scale)
#define DENSE_MUL_VEC        CONCAT3(dense,T,mul_vec)
#define DENSE_MUL_MAT        CONCAT3(dense,T,mul_mat)
#define DENSE_FACTORIZE      CONCAT3(dense,T,factorize)
#define DENSE_SOLVE_VEC      CONCAT3(dense,T,solve_vec)
#define DENSE_SOLVE_MAT      CONCAT3(dense,T,solve_mat)
#define DENSE_TRANSPOSE      CONCAT3(dense,T,transpose)
#define DENSE_TRANSPOSE2     CONCAT3(dense,T,transpose2)
#define DENSE_MUL_DENSE      CONCAT3(dense,T,mul_dense)
#define DENSE_EIG            CONCAT3(dense,T,eig)
#define DENSE_OUTPUT         CONCAT3(dense,T,output)
#define DENSE_EYE            CONCAT3(dense,eye,T)

type, extends(MATRIX_TYPE) :: DENSE_TYPE
  !! Dense matrix type.

  TT,      allocatable :: d(:,:)
    !! Raw matrix data
  TT,      allocatable :: lu(:,:)
    !! LU decomposition data
  integer, allocatable :: ipiv(:)
    !! Permutation used in conjunction with lu
contains
  procedure :: DENSE_INIT_BY_SIZE, DENSE_INIT_BY_DATA
  generic   :: init       => DENSE_INIT_BY_SIZE, DENSE_INIT_BY_DATA
  procedure :: destruct   => DENSE_DESTRUCT
  procedure :: reset      => DENSE_RESET
  procedure :: scale      => DENSE_SCALE
  procedure :: mul_vec    => DENSE_MUL_VEC
  procedure :: mul_mat    => DENSE_MUL_MAT
  procedure :: factorize  => DENSE_FACTORIZE
  procedure :: solve_vec  => DENSE_SOLVE_VEC
  procedure :: solve_mat  => DENSE_SOLVE_MAT
  procedure :: DENSE_TRANSPOSE                        ! transposes this matrix inplace
  procedure :: DENSE_TRANSPOSE2                       ! transpose this matrix into another
  generic   :: transpose  => DENSE_TRANSPOSE, DENSE_TRANSPOSE2
  procedure :: mul_dense  => DENSE_MUL_DENSE
  procedure :: eig        => DENSE_EIG
  procedure :: output     => DENSE_OUTPUT
end type

type DENSE_PTR_TYPE
  type(DENSE_TYPE), pointer :: p => null()
end type

interface
  module subroutine DENSE_INIT_BY_DATA(this, d)
    class(DENSE_TYPE), intent(out) :: this
    TT,                intent(in)  :: d(:,:)
  end subroutine

  module subroutine DENSE_INIT_BY_SIZE(this, nrows, ncols)
    class(DENSE_TYPE), intent(out) :: this
    integer,           intent(in)  :: nrows
    integer, optional, intent(in)  :: ncols
  end subroutine

  module subroutine DENSE_DESTRUCT(this)
    class(DENSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine DENSE_RESET(this, only_factorization)
    class(DENSE_TYPE), intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine DENSE_SCALE(this, fact)
    class(DENSE_TYPE), intent(inout) :: this
    TT,                intent(in)    :: fact
  end subroutine

  module subroutine DENSE_MUL_VEC(this, x, y, fact_y, trans)
    class(DENSE_TYPE),   intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine DENSE_MUL_MAT(this, x, y, fact_y, trans)
    class(DENSE_TYPE),   intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine DENSE_FACTORIZE(this)
    class(DENSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine DENSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(DENSE_TYPE),   intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine DENSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(DENSE_TYPE),   intent(in)  :: this
    TT,                  intent(in)  :: rhs(:,:)
    TT,                  intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine DENSE_TRANSPOSE(this)
    class(DENSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine DENSE_TRANSPOSE2(this, d2)
    class(DENSE_TYPE), intent(in) :: this
    type(DENSE_TYPE), intent(out) :: d2
  end subroutine

  module subroutine DENSE_MUL_DENSE(this, d2, d3)
    class(DENSE_TYPE), intent(in)  :: this
    type(DENSE_TYPE),  intent(in)  :: d2
    type(DENSE_TYPE),  intent(out) :: d3
  end subroutine

  module subroutine DENSE_EIG(this, e, R, L, sort)
    class(DENSE_TYPE), intent(in)  :: this
    complex,           intent(out) :: e(:)
    complex, optional, intent(out) :: R(:,:)
    complex, optional, intent(out) :: L(:,:)
    logical, optional, intent(in)  :: sort
  end subroutine

  module subroutine DENSE_OUTPUT(this, file, funit, fmt)
    class(DENSE_TYPE),      intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module function DENSE_EYE(nrows, ncols) result(eye)
    integer,            intent(in)  :: nrows
    integer, optional,  intent(in)  :: ncols
    type(DENSE_TYPE)                :: eye
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

#undef DENSE_TYPE
#undef DENSE_PTR_TYPE
#undef MATRIX_TYPE

#undef DENSE_INIT_BY_DATA
#undef DENSE_INIT_BY_SIZE
#undef DENSE_DESTRUCT
#undef DENSE_RESET
#undef DENSE_SCALE
#undef DENSE_MUL_VEC
#undef DENSE_MUL_MAT
#undef DENSE_FACTORIZE
#undef DENSE_SOLVE_VEC
#undef DENSE_SOLVE_MAT
#undef DENSE_TRANSPOSE
#undef DENSE_TRANSPOSE2
#undef DENSE_MUL_DENSE
#undef DENSE_EIG
#undef DENSE_OUTPUT
#undef DENSE_EYE
