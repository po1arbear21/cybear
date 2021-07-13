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

#define SPARSE_INIT          CONCAT3(sparse,T,init)
#define SPARSE_INIT_SOLVER   CONCAT3(sparse,T,init_solver)
#define SPARSE_DESTRUCT      CONCAT3(sparse,T,destruct)
#define SPARSE_RESET         CONCAT3(sparse,T,reset)
#define SPARSE_SCALE         CONCAT3(sparse,T,scale)
#define SPARSE_MUL_VEC       CONCAT3(sparse,T,mul_vec)
#define SPARSE_MUL_VEC_SLICE CONCAT3(sparse,T,mul_vec_slice)
#define SPARSE_MUL_MAT       CONCAT3(sparse,T,mul_mat)
#define SPARSE_FACTORIZE     CONCAT3(sparse,T,factorize)
#define SPARSE_SOLVE_VEC     CONCAT3(sparse,T,solve_vec)
#define SPARSE_SOLVE_MAT     CONCAT3(sparse,T,solve_mat)
#define SPARSE_TRANSPOSE     CONCAT3(sparse,T,transpose)
#define SPARSE_TRANSPOSE2    CONCAT3(sparse,T,transpose2)
#define SPARSE_MUL_SPARSE    CONCAT3(sparse,T,mul_sparse)
#define SPARSE_INPUT         CONCAT3(sparse,T,input)
#define SPARSE_OUTPUT        CONCAT3(sparse,T,output)
#define SPARSE_NNZ           CONCAT3(sparse,T,nnz)
#define SPARSE_IS_EMPTY      CONCAT3(sparse,T,is_empty)
#define SPARSE_EYE           CONCAT3(sparse,eye,T)
#define SPARSE_ZERO          CONCAT3(sparse,zero,T)

#define SPBUILD_INIT         CONCAT3(spbuild,T,init)
#define SPBUILD_DESTRUCT     CONCAT3(spbuild,T,destruct)
#define SPBUILD_RESET        CONCAT3(spbuild,T,reset)
#define SPBUILD_RESET_ROW    CONCAT3(spbuild,T,reset_row)
#define SPBUILD_SET          CONCAT3(spbuild,T,set)
#define SPBUILD_SET_ROW      CONCAT3(spbuild,T,set_row)
#define SPBUILD_ADD          CONCAT3(spbuild,T,add)
#define SPBUILD_SORT         CONCAT3(spbuild,T,sort)
#define SPBUILD_LOAD         CONCAT3(spbuild,T,load)
#define SPBUILD_SAVE         CONCAT3(spbuild,T,save)

type, extends(MATRIX_TYPE) :: SPARSE_TYPE
  !! Sparse matrix type in CSR format.

  integer(SPARSE_IDX), allocatable :: ia(:)
  integer,             allocatable :: ja(:)
  TT,                  allocatable :: a(:)

  integer :: solver = 0
    !! SOLVER_PARDISO, SOLVER_MUMPS or SOLVER_ILUPACK
  integer :: solver_handle = 0
    !! pardiso, mumps or ilupack solver handle
contains
  procedure :: SPARSE_INIT
  generic   :: init          => SPARSE_INIT
  procedure :: init_solver   => SPARSE_INIT_SOLVER
  procedure :: destruct      => SPARSE_DESTRUCT
  procedure :: reset         => SPARSE_RESET
  procedure :: scale         => SPARSE_SCALE
  procedure :: mul_vec       => SPARSE_MUL_VEC
  procedure :: mul_vec_slice => SPARSE_MUL_VEC_SLICE
  procedure :: mul_mat       => SPARSE_MUL_MAT
  procedure :: factorize     => SPARSE_FACTORIZE
  procedure :: solve_vec     => SPARSE_SOLVE_VEC
  procedure :: solve_mat     => SPARSE_SOLVE_MAT
  procedure :: SPARSE_TRANSPOSE
  procedure :: SPARSE_TRANSPOSE2
  generic   :: transpose     => SPARSE_TRANSPOSE, SPARSE_TRANSPOSE2
  procedure :: mul_sparse    => SPARSE_MUL_SPARSE
  procedure :: input         => SPARSE_INPUT
  procedure :: output        => SPARSE_OUTPUT
  procedure :: nnz           => SPARSE_NNZ
  procedure :: is_empty      => SPARSE_IS_EMPTY
end type

type SPARSE_PTR_TYPE
  type(SPARSE_TYPE), pointer :: p => null()
end type

type SPBUILD_TYPE
  !! Sparse matrix builder which is more flexible than simple CSR.
  type(SPARSE_TYPE), pointer     :: sp
    !! Non-owning pointer to sparse matrix this builder is tied to
  type(vector_int),  allocatable :: cols(:)
    !! Column indices of non-zero elements in each row
  type(VECTOR_TYPE), allocatable :: vals(:)
    !! Values of non-zero elements in each row
  type(vector_log),  allocatable :: filled(:)
    !! Indicates whether a value is actually used (keep structure after reset)
  logical, allocatable           :: sorted(:)
    !! Indicates whether a row is sorted or not
  logical, allocatable           :: keep_struct(:)
    !! Indicates whether a row should be deleted completely on reset, or structure should be kept
  logical                        :: assume_unique
    !! user ensures no duplicate entries
contains
  procedure :: init      => SPBUILD_INIT
  procedure :: destruct  => SPBUILD_DESTRUCT
  procedure :: reset     => SPBUILD_RESET
  procedure :: reset_row => SPBUILD_RESET_ROW
  procedure :: set       => SPBUILD_SET
  procedure :: set_row   => SPBUILD_SET_ROW
  procedure :: add       => SPBUILD_ADD
  procedure :: load      => SPBUILD_LOAD
  procedure :: save      => SPBUILD_SAVE
  procedure, private :: sort => SPBUILD_SORT
end type

interface
  module subroutine SPARSE_INIT(this, nrows, ncols)
    class(SPARSE_TYPE), intent(out) :: this
    integer,            intent(in)  :: nrows
    integer, optional,  intent(in)  :: ncols
  end subroutine

  module subroutine SPARSE_DESTRUCT(this)
    class(SPARSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPARSE_RESET(this, only_factorization)
    class(SPARSE_TYPE), intent(inout) :: this
    logical, optional,  intent(in)    :: only_factorization
  end subroutine

  module subroutine SPARSE_SCALE(this, fact)
    class(SPARSE_TYPE), intent(inout) :: this
    TT,                 intent(in)    :: fact
  end subroutine

  module subroutine SPARSE_MUL_VEC_SLICE(this, x, y, fact_y, i0, i1)
    class(SPARSE_TYPE), intent(in)    :: this
    TT,                 intent(in)    :: x(:)
    TT,                 intent(inout) :: y(:)
    TT,      optional,  intent(in)    :: fact_y
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: i1
  end subroutine

  module subroutine SPARSE_MUL_VEC(this, x, y, fact_y, trans)
    class(SPARSE_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine SPARSE_MUL_MAT(this, x, y, fact_y, trans)
    class(SPARSE_TYPE),  intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine SPARSE_INIT_SOLVER(this)
    class(SPARSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPARSE_FACTORIZE(this)
    class(SPARSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPARSE_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(SPARSE_TYPE),  intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine SPARSE_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(SPARSE_TYPE),  intent(in)  :: this
    TT,                  intent(in)  :: rhs(:,:)
    TT,                  intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine SPARSE_TRANSPOSE(this)
    class(SPARSE_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPARSE_TRANSPOSE2(this, s2)
    class(SPARSE_TYPE), intent(in)  :: this
    type(SPARSE_TYPE),  intent(out) :: s2
  end subroutine

  module subroutine SPARSE_MUL_SPARSE(this, s2, s3, sort)
    class(SPARSE_TYPE), intent(in)  :: this
    type(SPARSE_TYPE),  intent(in)  :: s2
    type(SPARSE_TYPE),  intent(out) :: s3
    logical, optional,  intent(in)  :: sort
  end subroutine

  module subroutine SPARSE_INPUT(this, file, funit)
    class(SPARSE_TYPE),     intent(out) :: this
    character(*), optional, intent(in)  :: file
    integer,      optional, intent(in)  :: funit
  end subroutine

  module subroutine SPARSE_OUTPUT(this, file, funit, fmt)
    class(SPARSE_TYPE),     intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module function SPARSE_NNZ(this, only_nonzeros) result(nnz)
    class(SPARSE_TYPE), intent(in) :: this
    logical, optional,  intent(in) :: only_nonzeros
    integer(SPARSE_IDX)            :: nnz
  end function

  module function SPARSE_IS_EMPTY(this) result(empty)
    class(SPARSE_TYPE), intent(in) :: this
    logical                        :: empty
  end function

  module function SPARSE_EYE(nrows, ncols) result(eye)
    integer,            intent(in)  :: nrows
    integer, optional,  intent(in)  :: ncols
    type(SPARSE_TYPE)               :: eye
  end function

  module function SPARSE_ZERO(nrows, ncols) result(zero)
    integer,           intent(in) :: nrows
    integer, optional, intent(in) :: ncols
    type(SPARSE_TYPE)             :: zero
  end function

  module subroutine SPBUILD_INIT(this, s, n_per_row, assume_unique)
    class(SPBUILD_TYPE),       intent(out) :: this
    type(SPARSE_TYPE), target, intent(in)  :: s
    integer, optional,         intent(in)  :: n_per_row
    logical, optional,         intent(in)  :: assume_unique
  end subroutine

  module subroutine SPBUILD_DESTRUCT(this)
    class(SPBUILD_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPBUILD_RESET(this)
    class(SPBUILD_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPBUILD_RESET_ROW(this, row)
    class(SPBUILD_TYPE), intent(inout) :: this
    integer,             intent(in)    :: row
  end subroutine

  module subroutine SPBUILD_SET(this, row, col, val, fact, search)
    class(SPBUILD_TYPE), intent(inout) :: this
    integer,             intent(in)    :: row
    integer,             intent(in)    :: col
    TT,                  intent(in)    :: val
    TT,      optional,   intent(in)    :: fact
    logical, optional,   intent(in)    :: search
  end subroutine

  module subroutine SPBUILD_SET_ROW(this, row, vals, j0, j1)
    class(SPBUILD_TYPE), intent(inout) :: this
    integer,             intent(in)    :: row
    TT,                  intent(in)    :: vals(:)
    integer, optional,   intent(in)    :: j0
    integer, optional,   intent(in)    :: j1
  end subroutine

  module subroutine SPBUILD_ADD(this, row, col, val)
    class(SPBUILD_TYPE), intent(inout) :: this
    integer,             intent(in)    :: row
    integer,             intent(in)    :: col
    TT,                  intent(in)    :: val
  end subroutine

  module subroutine SPBUILD_SORT(this)
    class(SPBUILD_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPBUILD_LOAD(this)
    class(SPBUILD_TYPE), intent(inout) :: this
  end subroutine

  module subroutine SPBUILD_SAVE(this)
    class(SPBUILD_TYPE), intent(inout) :: this
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

#undef SPARSE_INIT
#undef SPARSE_INIT_SOLVER
#undef SPARSE_DESTRUCT
#undef SPARSE_RESET
#undef SPARSE_SCALE
#undef SPARSE_MUL_VEC
#undef SPARSE_MUL_VEC_SLICE
#undef SPARSE_MUL_MAT
#undef SPARSE_FACTORIZE
#undef SPARSE_SOLVE_VEC
#undef SPARSE_SOLVE_MAT
#undef SPARSE_TRANSPOSE
#undef SPARSE_TRANSPOSE2
#undef SPARSE_MUL_SPARSE
#undef SPARSE_INPUT
#undef SPARSE_OUTPUT
#undef SPARSE_NNZ
#undef SPARSE_IS_EMPTY
#undef SPARSE_EYE
#undef SPARSE_ZERO

#undef SPBUILD_INIT
#undef SPBUILD_DESTRUCT
#undef SPBUILD_RESET
#undef SPBUILD_RESET_ROW
#undef SPBUILD_SET
#undef SPBUILD_SET_ROW
#undef SPBUILD_ADD
#undef SPBUILD_SORT
#undef SPBUILD_LOAD
#undef SPBUILD_SAVE
