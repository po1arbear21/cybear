#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_PTR_TYPE        CONCAT(band_ptr,T)
#define MATRIX_TYPE          CONCAT(matrix,T)

#define BAND_INIT            CONCAT3(band,T,init)
#define BAND_DESTRUCT        CONCAT3(band,T,destruct)
#define BAND_RESET           CONCAT3(band,T,reset)
#define BAND_RESET_ROW       CONCAT3(band,T,reset_row)
#define BAND_SCALE           CONCAT3(band,T,scale)
#define BAND_MUL_VEC         CONCAT3(band,T,mul_vec)
#define BAND_MUL_MAT         CONCAT3(band,T,mul_mat)
#define BAND_FACTORIZE       CONCAT3(band,T,factorize)
#define BAND_SOLVE_VEC       CONCAT3(band,T,solve_vec)
#define BAND_SOLVE_MAT       CONCAT3(band,T,solve_mat)
#define BAND_ADD_ELEM        CONCAT3(band,T,add_elem)
#define BAND_GET_ELEM        CONCAT3(band,T,get_elem)
#define BAND_SET_ELEM        CONCAT3(band,T,set_elem)
#define BAND_SET_DIAG_ARR    CONCAT3(band,T,set_diag_arr)
#define BAND_SET_DIAG_VAL    CONCAT3(band,T,set_diag_val)
#define BAND_OUTPUT          CONCAT3(band,T,output)
#define BAND_EYE             CONCAT3(band,eye,T)

type, extends(MATRIX_TYPE) :: BAND_TYPE
  !! Band matrix.
  !! Storage example
  !!    / a_11  a_12  a_13       \
  !!  A=| a_21  a_22  a_23  a_24 |
  !!    |       a_32  a_33  a_34 |
  !!    \             a_43  a_44 /
  !!  here n=4, nlower=1, nupper=2
  !!
  !! Band storage used here (shifted rows)
  !!
  !!          /             a_13  a_24 \    ! <- row -2
  !!          |       a_12  a_23  a_34 |    ! <- row -1
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 0
  !!          \ a_21  a_32  a_43       /    ! <- row 1
  !!  fixme numbering should be reversed b.c. matlab uses row>0 for upper diagonals. has effect on set_diag routines.
  !!
  !! LAPACK's band storage
  !!
  !!          /             a_13  a_24 \    ! <- row 1
  !!          |       a_12  a_23  a_34 \    ! <- row 2
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 3
  !!          \ a_21  a_32  a_43       /    ! <- row 4


  integer              :: nupper, nlower
    !! upper/lower number of diagonals (n>=0).
    !! nupper=nlower=0 <=> diagonal matrix

  TT, allocatable      :: d(:,:)
    !! matrix entries in band storage format.
    !! index1: diagonal number (-nupper:+nlower)
    !! index2: column (1:ncols)

  ! factorization data
  TT,      allocatable :: f(:,:)
  integer, allocatable :: ipiv(:)
contains
  procedure :: BAND_INIT
  generic   :: init      => BAND_INIT
  procedure :: destruct  => BAND_DESTRUCT
  procedure :: reset     => BAND_RESET
  procedure :: reset_row => BAND_RESET_ROW
  procedure :: scale     => BAND_SCALE
  procedure :: mul_vec   => BAND_MUL_VEC
  procedure :: mul_mat   => BAND_MUL_MAT
  procedure :: factorize => BAND_FACTORIZE
  procedure :: solve_vec => BAND_SOLVE_VEC
  procedure :: solve_mat => BAND_SOLVE_MAT
  procedure :: add_elem  => BAND_ADD_ELEM
  procedure :: get_elem  => BAND_GET_ELEM
  procedure :: set_elem  => BAND_SET_ELEM
  procedure :: BAND_SET_DIAG_VAL, BAND_SET_DIAG_ARR
  generic   :: set_diag  => BAND_SET_DIAG_VAL, BAND_SET_DIAG_ARR
  procedure :: output    => BAND_OUTPUT
end type

type BAND_PTR_TYPE
  type(BAND_TYPE), pointer :: p => null()
end type

interface
  module subroutine BAND_INIT(this, nrows, nlower, nupper, d0)
    class(BAND_TYPE),  intent(out) :: this
    integer,           intent(in)  :: nrows
    integer,           intent(in)  :: nlower
    integer, optional, intent(in)  :: nupper
    TT,      optional, intent(in)  :: d0(:,:)
  end subroutine

  module subroutine BAND_DESTRUCT(this)
    class(BAND_TYPE), intent(inout) :: this
  end subroutine

  module subroutine BAND_RESET(this, only_factorization)
    class(BAND_TYPE),  intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine BAND_RESET_ROW(this, i)
    class(BAND_TYPE),  intent(inout) :: this
    integer,           intent(in)    :: i
  end subroutine

  module subroutine BAND_SCALE(this, fact)
    class(BAND_TYPE), intent(inout) :: this
    TT,               intent(in)    :: fact
  end subroutine

  module subroutine BAND_MUL_VEC(this, x, y, fact_y, trans)
    class(BAND_TYPE),    intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine BAND_MUL_MAT(this, x, y, fact_y, trans)
    class(BAND_TYPE),    intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine BAND_FACTORIZE(this)
    class(BAND_TYPE), intent(inout) :: this
  end subroutine

  module subroutine BAND_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(BAND_TYPE),    intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine BAND_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(BAND_TYPE),    intent(in)  :: this
    TT,                  intent(in)  :: rhs(:,:)
    TT,                  intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine BAND_OUTPUT(this, file, funit, fmt)
    class(BAND_TYPE),       intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module subroutine BAND_ADD_ELEM(this, i, j, x)
    class(BAND_TYPE), intent(inout) :: this
    integer,          intent(in)    :: i
    integer,          intent(in)    :: j
    TT,               intent(in)    :: x
  end subroutine

  module function BAND_GET_ELEM(this, i, j) result(x)
    class(BAND_TYPE), intent(in) :: this
    integer,          intent(in) :: i
    integer,          intent(in) :: j
    TT                           :: x
  end function

  module subroutine BAND_SET_ELEM(this, i, j, x)
    class(BAND_TYPE), intent(inout) :: this
    integer,          intent(in)    :: i
    integer,          intent(in)    :: j
    TT,               intent(in)    :: x
  end subroutine

  module subroutine BAND_SET_DIAG_VAL(this, val, k)
    class(BAND_TYPE), intent(inout)        :: this
    TT,               intent(in)           :: val
    integer,          intent(in), optional :: k
  end subroutine

  module subroutine BAND_SET_DIAG_ARR(this, vals, k)
    class(BAND_TYPE), intent(inout)        :: this
    TT,               intent(in)           :: vals(:)
    integer,          intent(in), optional :: k
  end subroutine

  module function BAND_EYE(nrows) result(eye)
    integer, intent(in) :: nrows
    type(BAND_TYPE)     :: eye
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

#undef BAND_TYPE
#undef BAND_PTR_TYPE
#undef MATRIX_TYPE

#undef BAND_INIT
#undef BAND_DESTRUCT
#undef BAND_RESET
#undef BAND_RESET_ROW
#undef BAND_SCALE
#undef BAND_MUL_VEC
#undef BAND_MUL_MAT
#undef BAND_FACTORIZE
#undef BAND_SOLVE_VEC
#undef BAND_SOLVE_MAT
#undef BAND_ADD_ELEM
#undef BAND_GET_ELEM
#undef BAND_SET_ELEM
#undef BAND_SET_DIAG_ARR
#undef BAND_SET_DIAG_VAL
#undef BAND_OUTPUT
#undef BAND_EYE
