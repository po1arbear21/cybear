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

#define BLOCK_INIT           CONCAT3(block,T,init)
#define BLOCK_DESTRUCT       CONCAT3(block,T,destruct)
#define BLOCK_RESET          CONCAT3(block,T,reset)
#define BLOCK_DELETE_BLOCK   CONCAT3(block,T,delete_block)
#define BLOCK_SCALE          CONCAT3(block,T,scale)
#define BLOCK_MUL_VEC        CONCAT3(block,T,mul_vec)
#define BLOCK_MUL_MAT        CONCAT3(block,T,mul_mat)
#define BLOCK_FACTORIZE      CONCAT3(block,T,factorize)
#define BLOCK_SOLVE_VEC      CONCAT3(block,T,solve_vec)
#define BLOCK_SOLVE_MAT      CONCAT3(block,T,solve_mat)
#define BLOCK_GET_GLOBAL_IDX CONCAT3(block,T,get_global_idx)
#define BLOCK_GET_LOCAL_IDX  CONCAT3(block,T,get_local_idx)
#define BLOCK_SET_BLOCK_TYPE CONCAT3(block,T,set_block_type)
#define BLOCK_GET_DENSE      CONCAT3(block,T,get_dense)
#define BLOCK_GET_SPARSE     CONCAT3(block,T,get_sparse)
#define BLOCK_GET_BAND       CONCAT3(block,T,get_band)
#define BLOCK_GET_BLOCK      CONCAT3(block,T,get_block)
#define BLOCK_SET            CONCAT3(block,T,set)
#define BLOCK_SET_PTR        CONCAT3(block,T,set_ptr)
#define BLOCK_COPY           CONCAT3(block,T,copy)
#define BLOCK_BLOCK_ALLOWED  CONCAT3(block,T,block_allowed)
#define BLOCK_CHECK_SIZES    CONCAT3(block,T,check_sizes)

type, extends(MATRIX_TYPE) :: BLOCK_TYPE
  !! Block matrix
  !!
  !! Example
  !!        / A B \
  !!    M = | C D | for A=(2x3), B=(2x4), C=(3x3), D=(3x4), E=(1x3), F=(1x4)-matrix
  !!        \ E F /
  !!
  !!    => nbrows=3, nbcols=2, i0=[1,3,6], i1=[2,5,6], j0=[1,4], j1=[3,7]

  ! dimensions
  integer              :: nbrows
    !! number of block rows
  integer              :: nbcols
    !! number of block cols
  integer, allocatable :: i0(:)
    !! block row start indices
  integer, allocatable :: i1(:)
    !! block row stop  indices
  integer, allocatable :: j0(:)
    !! block col start indices
  integer, allocatable :: j1(:)
    !! block col stop  indices

  ! blocks
  type(MATRIX_PTR_TYPE), allocatable :: b(:,:)
  logical,               allocatable :: alc(:,:)
    !! internal flags to keep track of allocations of b(i,j)%p (delete b(i,j)%p in destruct if alc(i,j) == true)

  ! flags
  logical :: diag
  logical :: tridiag
contains
  procedure :: BLOCK_INIT
  generic   :: init           => BLOCK_INIT
  procedure :: destruct       => BLOCK_DESTRUCT
  procedure :: reset          => BLOCK_RESET
  procedure :: delete_block   => BLOCK_DELETE_BLOCK
  procedure :: scale          => BLOCK_SCALE
  procedure :: mul_vec        => BLOCK_MUL_VEC
  procedure :: mul_mat        => BLOCK_MUL_MAT
  procedure :: factorize      => BLOCK_FACTORIZE
  procedure :: solve_vec      => BLOCK_SOLVE_VEC
  procedure :: solve_mat      => BLOCK_SOLVE_MAT
  procedure :: get_global_idx => BLOCK_GET_GLOBAL_IDX
  procedure :: get_local_idx  => BLOCK_GET_LOCAL_IDX
  procedure :: BLOCK_GET_DENSE
  procedure :: BLOCK_GET_SPARSE
  procedure :: BLOCK_GET_BAND
  procedure :: BLOCK_GET_BLOCK
  generic   :: get            => BLOCK_GET_DENSE, BLOCK_GET_SPARSE, BLOCK_GET_BAND, BLOCK_GET_BLOCK
  procedure :: set            => BLOCK_SET
  procedure :: set_ptr        => BLOCK_SET_PTR
  procedure :: copy           => BLOCK_COPY

  procedure, private :: block_allowed => BLOCK_BLOCK_ALLOWED
  procedure, private :: check_sizes   => BLOCK_CHECK_SIZES
end type

type BLOCK_PTR_TYPE
  type(BLOCK_TYPE), pointer :: p => null()
end type

interface
  module subroutine BLOCK_INIT(this, row_dim, col_dim, diag, tridiag)
    class(BLOCK_TYPE), intent(out) :: this
    integer,           intent(in)  :: row_dim(:)
    integer, optional, intent(in)  :: col_dim(:)
    logical, optional, intent(in)  :: diag
    logical, optional, intent(in)  :: tridiag
  end subroutine

  module subroutine BLOCK_DESTRUCT(this)
    class(BLOCK_TYPE), intent(inout) :: this
  end subroutine

  module subroutine BLOCK_DELETE_BLOCK(this, ib, jb)
    class(BLOCK_TYPE), intent(inout) :: this
    integer,           intent(in)    :: ib, jb
  end subroutine

  module subroutine BLOCK_RESET(this, only_factorization)
    class(BLOCK_TYPE), intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine BLOCK_SCALE(this, fact)
    class(BLOCK_TYPE), intent(inout) :: this
    TT,                intent(in)    :: fact
  end subroutine

  module subroutine BLOCK_MUL_VEC(this, x, y, fact_y, trans)
    class(BLOCK_TYPE),   intent(in)    :: this
    TT,                  intent(in)    :: x(:)
    TT,                  intent(inout) :: y(:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine BLOCK_MUL_MAT(this, x, y, fact_y, trans)
    class(BLOCK_TYPE),   intent(in)    :: this
    TT,                  intent(in)    :: x(:,:)
    TT,                  intent(inout) :: y(:,:)
    TT,        optional, intent(in)    :: fact_y
    character, optional, intent(in)    :: trans
  end subroutine

  module subroutine BLOCK_FACTORIZE(this)
    class(BLOCK_TYPE), intent(inout) :: this
  end subroutine

  module subroutine BLOCK_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    class(BLOCK_TYPE),   intent(in)  :: this
    TT,                  intent(in)  :: rhs(:)
    TT,                  intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine BLOCK_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    class(BLOCK_TYPE),   intent(in)  :: this
    TT,                  intent(in)  :: rhs(:,:)
    TT,                  intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine BLOCK_GET_GLOBAL_IDX(this, ib, jb, il, jl, ig, jg)
    class(BLOCK_TYPE), intent(in)  :: this
    integer,           intent(in)  :: ib
    integer,           intent(in)  :: jb
    integer,           intent(in)  :: il
    integer,           intent(in)  :: jl
    integer,           intent(out) :: ig
    integer,           intent(out) :: jg
  end subroutine

  module subroutine BLOCK_GET_LOCAL_IDX(this, ig, jg, ib, jb, il, jl)
    class(BLOCK_TYPE), intent(in)  :: this
    integer,           intent(in)  :: ig
    integer,           intent(in)  :: jg
    integer,           intent(out) :: ib
    integer,           intent(out) :: jb
    integer,           intent(out) :: il
    integer,           intent(out) :: jl
  end subroutine

  module subroutine BLOCK_GET_DENSE(this, ib, jb, p)
    class(BLOCK_TYPE),         intent(inout) :: this
    integer,                   intent(in)    :: ib
    integer,                   intent(in)    :: jb
    type(DENSE_TYPE), pointer, intent(out)   :: p
  end subroutine

  module subroutine BLOCK_GET_SPARSE(this, ib, jb, p)
    class(BLOCK_TYPE),          intent(inout) :: this
    integer,                    intent(in)    :: ib
    integer,                    intent(in)    :: jb
    type(SPARSE_TYPE), pointer, intent(out)   :: p
  end subroutine

  module subroutine BLOCK_GET_BAND(this, ib, jb, p)
    class(BLOCK_TYPE),         intent(inout) :: this
    integer,                   intent(in)    :: ib
    integer,                   intent(in)    :: jb
    type(BAND_TYPE), pointer,  intent(out)   :: p
  end subroutine

  module subroutine BLOCK_GET_BLOCK(this, ib, jb, p)
    class(BLOCK_TYPE),         intent(inout) :: this
    integer,                   intent(in)    :: ib
    integer,                   intent(in)    :: jb
    type(BLOCK_TYPE), pointer, intent(out)   :: p
  end subroutine

  module subroutine BLOCK_SET(this, ib, jb, b, fact)
    class(BLOCK_TYPE),  intent(inout) :: this
    integer,            intent(in)    :: ib
    integer,            intent(in)    :: jb
    class(MATRIX_TYPE), intent(in)    :: b
    TT, optional,       intent(in)    :: fact
  end subroutine

  module subroutine BLOCK_SET_PTR(this, ib, jb, b)
    class(BLOCK_TYPE),          intent(inout) :: this
    integer,                    intent(in)    :: ib
    integer,                    intent(in)    :: jb
    class(MATRIX_TYPE), target, intent(in)    :: b
  end subroutine

  module subroutine BLOCK_COPY(this, copy)
    class(BLOCK_TYPE), intent(in)  :: this
    type(BLOCK_TYPE),  intent(out) :: copy
  end subroutine

  module function BLOCK_BLOCK_ALLOWED(this, ib, jb) result(tf)
    class(BLOCK_TYPE), intent(in) :: this
    integer,           intent(in) :: ib
    integer,           intent(in) :: jb
    logical                       :: tf
  end function

  module function BLOCK_CHECK_SIZES(this) result(tf)
    class(BLOCK_TYPE), intent(in) :: this
    logical                       :: tf
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

#undef BLOCK_INIT
#undef BLOCK_DESTRUCT
#undef BLOCK_RESET
#undef BLOCK_DELETE_BLOCK
#undef BLOCK_SCALE
#undef BLOCK_MUL_VEC
#undef BLOCK_MUL_MAT
#undef BLOCK_FACTORIZE
#undef BLOCK_SOLVE_VEC
#undef BLOCK_SOLVE_MAT
#undef BLOCK_GET_GLOBAL_IDX
#undef BLOCK_GET_LOCAL_IDX
#undef BLOCK_SET_BLOCK_TYPE
#undef BLOCK_GET_DENSE
#undef BLOCK_GET_SPARSE
#undef BLOCK_GET_BAND
#undef BLOCK_GET_BLOCK
#undef BLOCK_SET
#undef BLOCK_SET_PTR
#undef BLOCK_COPY
#undef BLOCK_BLOCK_ALLOWED
#undef BLOCK_CHECK_SIZES
