#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)

#define MATRIX_TYPE          CONCAT(matrix,T)
#define MATRIX_ALC_TYPE      CONCAT(matrix_alc,T)
#define MATRIX_PTR_TYPE      CONCAT(matrix_ptr,T)
#define MATRIX_INIT          CONCAT3(matrix,T,init)
#define MATRIX_DESTRUCT      CONCAT3(matrix,T,destruct)
#define MATRIX_RESET         CONCAT3(matrix,T,reset)
#define MATRIX_SCALE         CONCAT3(matrix,T,scale)
#define MATRIX_MUL_VEC       CONCAT3(matrix,T,mul_vec)
#define MATRIX_MUL_MAT       CONCAT3(matrix,T,mul_mat)
#define MATRIX_FACTORIZE     CONCAT3(matrix,T,factorize)
#define MATRIX_SOLVE_VEC     CONCAT3(matrix,T,solve_vec)
#define MATRIX_SOLVE_MAT     CONCAT3(matrix,T,solve_mat)
#define MATRIX_TO_DENSE      CONCAT3(matrix,T,to_dense)
#define MATRIX_TO_SPARSE     CONCAT3(matrix,T,to_sparse)

#define DENSE_TYPE           CONCAT(dense,T)
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
#define DENSE_TO_DENSE       CONCAT3(dense,T,to_dense)
#define DENSE_TO_SPARSE      CONCAT3(dense,T,to_sparse)
#ifdef TCMPLX
#define DENSE_TO_REAL        CONCAT3(dense,T,to_real)
#else
#define DENSE_TO_CMPLX       CONCAT3(dense,T,to_cmplx)
#endif
#define DENSE_TRANSPOSE      CONCAT3(dense,T,transpose)
#define DENSE_TRANSPOSE2     CONCAT3(dense,T,transpose2)
#define DENSE_ADD_DENSE      CONCAT3(dense,T,add_dense)
#define DENSE_ADD_DENSE3     CONCAT3(dense,T,add_dense3)
#define DENSE_ADD_SPARSE     CONCAT3(dense,T,add_sparse)
#define DENSE_ADD_SPARSE3    CONCAT3(dense,T,add_sparse3)
#define DENSE_ADD_BAND       CONCAT3(dense,T,add_band)
#define DENSE_ADD_BAND3      CONCAT3(dense,T,add_band3)
#define DENSE_MUL_DENSE      CONCAT3(dense,T,mul_dense)
#define DENSE_EIG            CONCAT3(dense,T,eig)
#define DENSE_OUTPUT         CONCAT3(dense,T,output)

#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPARSE_INIT          CONCAT3(sparse,T,init)
#define SPARSE_DESTRUCT      CONCAT3(sparse,T,destruct)
#define SPARSE_RESET         CONCAT3(sparse,T,reset)
#define SPARSE_SCALE         CONCAT3(sparse,T,scale)
#define SPARSE_MUL_VEC       CONCAT3(sparse,T,mul_vec)
#define SPARSE_MUL_MAT       CONCAT3(sparse,T,mul_mat)
#define SPARSE_FACTORIZE     CONCAT3(sparse,T,factorize)
#define SPARSE_SOLVE_VEC     CONCAT3(sparse,T,solve_vec)
#define SPARSE_SOLVE_MAT     CONCAT3(sparse,T,solve_mat)
#define SPARSE_TO_DENSE      CONCAT3(sparse,T,to_dense)
#define SPARSE_TO_SPARSE     CONCAT3(sparse,T,to_sparse)
#ifdef TCMPLX
#define SPARSE_TO_REAL       CONCAT3(sparse,T,to_real)
#else
#define SPARSE_TO_CMPLX      CONCAT3(sparse,T,to_cmplx)
#endif
#define SPARSE_TRANSPOSE     CONCAT3(sparse,T,transpose)
#define SPARSE_TRANSPOSE2    CONCAT3(sparse,T,transpose2)
#define SPARSE_ADD_SPARSE    CONCAT3(sparse,T,add_sparse)
#define SPARSE_ADD_SPARSE3   CONCAT3(sparse,T,add_sparse3)
#define SPARSE_ADD_BAND      CONCAT3(sparse,T,add_band)
#define SPARSE_ADD_BAND3     CONCAT3(sparse,T,add_band3)
#define SPARSE_MUL_SPARSE    CONCAT3(sparse,T,mul_sparse)
#define SPARSE_OUTPUT        CONCAT3(sparse,T,output)

#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define SPBUILD_INIT         CONCAT3(spbuild,T,init)
#define SPBUILD_DESTRUCT     CONCAT3(spbuild,T,destruct)
#define SPBUILD_RESET        CONCAT3(spbuild,T,reset)
#define SPBUILD_SET          CONCAT3(spbuild,T,set)
#define SPBUILD_SET_ROW      CONCAT3(spbuild,T,set_row)
#define SPBUILD_ADD          CONCAT3(spbuild,T,add)
#define SPBUILD_SORT         CONCAT3(spbuild,T,sort)
#define SPBUILD_LOAD         CONCAT3(spbuild,T,load)
#define SPBUILD_SAVE         CONCAT3(spbuild,T,save)

#define BAND_TYPE            CONCAT(band,T)
#define BAND_INIT            CONCAT3(band,T,init)
#define BAND_DESTRUCT        CONCAT3(band,T,destruct)
#define BAND_RESET           CONCAT3(band,T,reset)
#define BAND_SCALE           CONCAT3(band,T,scale)
#define BAND_MUL_VEC         CONCAT3(band,T,mul_vec)
#define BAND_MUL_MAT         CONCAT3(band,T,mul_mat)
#define BAND_FACTORIZE       CONCAT3(band,T,factorize)
#define BAND_SOLVE_VEC       CONCAT3(band,T,solve_vec)
#define BAND_SOLVE_MAT       CONCAT3(band,T,solve_mat)
#define BAND_TO_DENSE        CONCAT3(band,T,to_dense)
#define BAND_TO_SPARSE       CONCAT3(band,T,to_sparse)
#ifdef TCMPLX
#define BAND_TO_REAL         CONCAT3(band,T,to_real)
#else
#define BAND_TO_CMPLX        CONCAT3(band,T,to_cmplx)
#endif
#define BAND_ADD_BAND        CONCAT3(band,T,add_band)
#define BAND_ADD_BAND3       CONCAT3(band,T,add_band3)
#define BAND_OUTPUT          CONCAT3(band,T,output)

#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define HESSENBERG_INIT      CONCAT3(hessenberg,T,init)
#define HESSENBERG_DESTRUCT  CONCAT3(hessenberg,T,destruct)
#define HESSENBERG_RESET     CONCAT3(hessenberg,T,reset)
#define HESSENBERG_SCALE     CONCAT3(hessenberg,T,scale)
#define HESSENBERG_MUL_VEC   CONCAT3(hessenberg,T,mul_vec)
#define HESSENBERG_MUL_MAT   CONCAT3(hessenberg,T,mul_mat)
#define HESSENBERG_FACTORIZE CONCAT3(hessenberg,T,factorize)
#define HESSENBERG_SOLVE_VEC CONCAT3(hessenberg,T,solve_vec)
#define HESSENBERG_SOLVE_MAT CONCAT3(hessenberg,T,solve_mat)
#define HESSENBERG_TO_DENSE  CONCAT3(hessenberg,T,to_dense)
#define HESSENBERG_TO_SPARSE CONCAT3(hessenberg,T,to_sparse)
#define HESSENBERG_TO_TRIANG CONCAT3(hessenberg,T,to_triang)
#define HESSENBERG_IS_TRIANG CONCAT3(hessenberg,T,is_triang)
#ifdef TCMPLX
#define HESSENBERG_TO_REAL   CONCAT3(hessenberg,T,to_real)
#else
#define HESSENBERG_TO_CMPLX  CONCAT3(hessenberg,T,to_cmplx)
#endif

#define TRIANG_TYPE          CONCAT(triang,T)
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
#define TRIANG_TO_DENSE      CONCAT3(triang,T,to_dense)
#define TRIANG_TO_SPARSE     CONCAT3(triang,T,to_sparse)
#ifdef TCMPLX
#define TRIANG_TO_REAL       CONCAT3(triang,T,to_real)
#else
#define TRIANG_TO_CMPLX      CONCAT3(triang,T,to_cmplx)
#endif

#define BLOCK_TYPE           CONCAT(block,T)
#define BLOCK_INIT           CONCAT3(block,T,init)
#define BLOCK_DESTRUCT       CONCAT3(block,T,destruct)
#define BLOCK_RESET          CONCAT3(block,T,reset)
#define BLOCK_SCALE          CONCAT3(block,T,scale)
#define BLOCK_MUL_VEC        CONCAT3(block,T,mul_vec)
#define BLOCK_MUL_MAT        CONCAT3(block,T,mul_mat)
#define BLOCK_FACTORIZE      CONCAT3(block,T,factorize)
#define BLOCK_SOLVE_VEC      CONCAT3(block,T,solve_vec)
#define BLOCK_SOLVE_MAT      CONCAT3(block,T,solve_mat)
#define BLOCK_TO_DENSE       CONCAT3(block,T,to_dense)
#define BLOCK_TO_SPARSE      CONCAT3(block,T,to_sparse)
#define BLOCK_GET_GLOBAL_IDX CONCAT3(block,T,get_global_idx)
#define BLOCK_GET_LOCAL_IDX  CONCAT3(block,T,get_local_idx)
#define BLOCK_SET_BLOCK_TYPE CONCAT3(block,T,set_block_type)
#define BLOCK_GET_DENSE      CONCAT3(block,T,get_dense)
#define BLOCK_GET_SPARSE     CONCAT3(block,T,get_sparse)
#define BLOCK_GET_BAND       CONCAT3(block,T,get_band)
#define BLOCK_GET_BLOCK      CONCAT3(block,T,get_block)

type, abstract :: MATRIX_TYPE
  !! Base type for matrices.
  !! Can be extended if more specialized matrices are needed besides the built-in types.

  integer :: nrows
    !! Number of rows
  integer :: ncols
    !! Number of columns
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
  procedure(MATRIX_TO_DENSE),  deferred :: to_dense
  procedure(MATRIX_TO_SPARSE), deferred :: to_sparse
end type

type MATRIX_ALC_TYPE
  class(MATRIX_TYPE), allocatable :: p
end type

type MATRIX_PTR_TYPE
  class(MATRIX_TYPE), pointer :: p => null()
end type

type, extends(MATRIX_TYPE) :: DENSE_TYPE
  !! Dense matrix type.

  TT, allocatable :: d(:,:)
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
  procedure :: to_dense   => DENSE_TO_DENSE           ! insert this matrix into another dense matrix
  procedure :: to_sparse  => DENSE_TO_SPARSE          ! insert this matrix into a sparse matrix
#ifdef TCMPLX
  procedure :: to_real    => DENSE_TO_REAL
#else
  procedure :: to_cmplx   => DENSE_TO_CMPLX
#endif
  procedure :: DENSE_TRANSPOSE                        ! transposes this matrix inplace
  procedure :: DENSE_TRANSPOSE2                       ! transpose this matrix into another
  generic   :: transpose  => DENSE_TRANSPOSE, DENSE_TRANSPOSE2
  procedure :: DENSE_ADD_DENSE                        ! adds matrix to this one
  procedure :: DENSE_ADD_DENSE3                       ! adds this and a second into a third matrix
  generic   :: add_dense  => DENSE_ADD_DENSE, DENSE_ADD_DENSE3
  procedure :: DENSE_ADD_SPARSE
  procedure :: DENSE_ADD_SPARSE3
  generic   :: add_sparse => DENSE_ADD_SPARSE, DENSE_ADD_SPARSE3
  procedure :: DENSE_ADD_BAND
  procedure :: DENSE_ADD_BAND3
  generic   :: add_band   => DENSE_ADD_BAND, DENSE_ADD_BAND3
  procedure :: mul_dense  => DENSE_MUL_DENSE
  procedure :: eig        => DENSE_EIG
  procedure :: output     => DENSE_OUTPUT
end type

type, extends(MATRIX_TYPE) :: SPARSE_TYPE
  !! Sparse matrix type in CSR format.

  integer, allocatable :: ia(:)
  integer, allocatable :: ja(:)
  TT,      allocatable :: a(:)

  integer              :: pard = 0
    !! Handle to pardiso for solving
contains
  procedure :: SPARSE_INIT
  generic   :: init         => SPARSE_INIT
  procedure :: destruct     => SPARSE_DESTRUCT
  procedure :: reset        => SPARSE_RESET
  procedure :: scale        => SPARSE_SCALE
  procedure :: mul_vec      => SPARSE_MUL_VEC
  procedure :: mul_mat      => SPARSE_MUL_MAT
  procedure :: factorize    => SPARSE_FACTORIZE
  procedure :: solve_vec    => SPARSE_SOLVE_VEC
  procedure :: solve_mat    => SPARSE_SOLVE_MAT
  procedure :: to_dense     => SPARSE_TO_DENSE
  procedure :: to_sparse    => SPARSE_TO_SPARSE
#ifdef TCMPLX
  procedure :: to_real      => SPARSE_TO_REAL
#else
  procedure :: to_cmplx     => SPARSE_TO_CMPLX
#endif
  procedure :: SPARSE_TRANSPOSE
  procedure :: SPARSE_TRANSPOSE2
  generic   :: transpose    => SPARSE_TRANSPOSE, SPARSE_TRANSPOSE2
  procedure :: SPARSE_ADD_SPARSE
  procedure :: SPARSE_ADD_SPARSE3
  generic   :: add_sparse   => SPARSE_ADD_SPARSE, SPARSE_ADD_SPARSE3
  procedure :: SPARSE_ADD_BAND
  procedure :: SPARSE_ADD_BAND3
  generic   :: add_band     => SPARSE_ADD_BAND, SPARSE_ADD_BAND3
  procedure :: mul_sparse   => SPARSE_MUL_SPARSE
  procedure :: output       => SPARSE_OUTPUT
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
contains
  procedure :: init     => SPBUILD_INIT
  procedure :: destruct => SPBUILD_DESTRUCT
  procedure :: reset    => SPBUILD_RESET
  procedure :: set      => SPBUILD_SET
  procedure :: set_row  => SPBUILD_SET_ROW
  procedure :: add      => SPBUILD_ADD
  procedure :: sort     => SPBUILD_SORT
  procedure :: load     => SPBUILD_LOAD
  procedure :: save     => SPBUILD_SAVE
end type

type, extends(MATRIX_TYPE) :: BAND_TYPE
  !! Band storage example
  !!    / a_11  a_12  a_13       \
  !!  A=| a_21  a_22  a_23 a_24  |
  !!    |       a_32  a_33  a_34 |
  !!    \             a_43  a_44 /
  !!  here n=4, nupper=nlower=1
  !!
  !! LAPACK band storage
  !!
  !!          /             a_13  a_24 \    ! <- row 1
  !!          /       a_12  a_23  a_34 \    ! <- row 2
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 3
  !!          \ a_21  a_32  a_43       /    ! <- row 4
  !!
  !! This band storage (shifted rows)
  !!
  !!          /             a_13  a_24 \    ! <- row -2
  !!          /       a_12  a_23  a_34 \    ! <- row -1
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 0
  !!          \ a_21  a_32  a_43       /    ! <- row 1


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
  procedure :: scale     => BAND_SCALE
  procedure :: mul_vec   => BAND_MUL_VEC
  procedure :: mul_mat   => BAND_MUL_MAT
  procedure :: factorize => BAND_FACTORIZE
  procedure :: solve_vec => BAND_SOLVE_VEC
  procedure :: solve_mat => BAND_SOLVE_MAT
  procedure :: to_dense  => BAND_TO_DENSE
  procedure :: to_sparse => BAND_TO_SPARSE
#ifdef TCMPLX
  procedure :: to_real   => BAND_TO_REAL
#else
  procedure :: to_cmplx  => BAND_TO_CMPLX
#endif
  procedure :: BAND_ADD_BAND
  procedure :: BAND_ADD_BAND3
  generic   :: add_band  => BAND_ADD_BAND, BAND_ADD_BAND3
  procedure :: output    => BAND_OUTPUT
end type

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
  procedure :: to_dense  => HESSENBERG_TO_DENSE
  procedure :: to_sparse => HESSENBERG_TO_SPARSE
  procedure :: to_triang => HESSENBERG_TO_TRIANG
#ifdef TCMPLX
  procedure :: to_real   => HESSENBERG_TO_REAL
#else
  procedure :: to_cmplx  => HESSENBERG_TO_CMPLX
#endif
  procedure :: is_triang => HESSENBERG_IS_TRIANG
end type

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
  procedure :: to_dense  => TRIANG_TO_DENSE
  procedure :: to_sparse => TRIANG_TO_SPARSE
#ifdef TCMPLX
  procedure :: to_real   => TRIANG_TO_REAL
#else
  procedure :: to_cmplx  => TRIANG_TO_CMPLX
#endif
end type

type, extends(MATRIX_TYPE) :: BLOCK_TYPE
  ! dimensions
  integer              :: nbrows ! number of block rows
  integer              :: nbcols ! number of block cols
  integer, allocatable :: i0(:)  ! block row start indices
  integer, allocatable :: i1(:)  ! block row stop  indices
  integer, allocatable :: j0(:)  ! block col start indices
  integer, allocatable :: j1(:)  ! block col stop  indices

  ! blocks
  type(MATRIX_ALC_TYPE), allocatable :: b(:,:)

  ! flags
  logical :: diag
  logical :: tridiag
  logical :: factorized
contains
  procedure :: BLOCK_INIT
  generic   :: init           => BLOCK_INIT
  procedure :: destruct       => BLOCK_DESTRUCT
  procedure :: reset          => BLOCK_RESET
  procedure :: scale          => BLOCK_SCALE
  procedure :: mul_vec        => BLOCK_MUL_VEC
  procedure :: mul_mat        => BLOCK_MUL_MAT
  procedure :: factorize      => BLOCK_FACTORIZE
  procedure :: solve_vec      => BLOCK_SOLVE_VEC
  procedure :: solve_mat      => BLOCK_SOLVE_MAT
  procedure :: to_dense       => BLOCK_TO_DENSE
  procedure :: to_sparse      => BLOCK_TO_SPARSE
  procedure :: get_global_idx => BLOCK_GET_GLOBAL_IDX
  procedure :: get_local_idx  => BLOCK_GET_LOCAL_IDX
  procedure :: BLOCK_GET_DENSE
  procedure :: BLOCK_GET_SPARSE
  procedure :: BLOCK_GET_BAND
  procedure :: BLOCK_GET_BLOCK
  generic   :: get => BLOCK_GET_DENSE, BLOCK_GET_SPARSE, BLOCK_GET_BAND, BLOCK_GET_BLOCK
end type

interface
  subroutine MATRIX_DESTRUCT(this)
    import MATRIX_TYPE
    class(MATRIX_TYPE), intent(inout) :: this
  end subroutine

  subroutine MATRIX_RESET(this, only_factorization)
    import MATRIX_TYPE
    class(MATRIX_TYPE), intent(inout) :: this
    logical, optional,  intent(in)    :: only_factorization
  end subroutine

  subroutine MATRIX_SCALE(this, fact)
    import MATRIX_TYPE
    class(MATRIX_TYPE), intent(inout) :: this
    TT,                 intent(in)    :: fact
  end subroutine

  subroutine MATRIX_MUL_VEC(this, x, y, fact_y, trans)
    import MATRIX_TYPE
    class(MATRIX_TYPE),         intent(in)    :: this
    TT,                         intent(in)    :: x(:)
    TT,                         intent(inout) :: y(:)
    TT,               optional, intent(in)    :: fact_y
    character(len=1), optional, intent(in)    :: trans
  end subroutine

  subroutine MATRIX_MUL_MAT(this, x, y, fact_y, trans)
    import MATRIX_TYPE
    class(MATRIX_TYPE),         intent(in)    :: this
    TT,                         intent(in)    :: x(:,:)
    TT,                         intent(inout) :: y(:,:)
    TT,               optional, intent(in)    :: fact_y
    character(len=1), optional, intent(in)    :: trans
  end subroutine

  subroutine MATRIX_FACTORIZE(this)
    import MATRIX_TYPE
    class(MATRIX_TYPE), intent(inout) :: this
  end subroutine

  subroutine MATRIX_SOLVE_VEC(this, rhs, x, trans, iter_refine)
    import MATRIX_TYPE
    class(MATRIX_TYPE),         intent(in)  :: this
    TT,                         intent(in)  :: rhs(:)
    TT,                         intent(out) :: x(:)
    character(len=1), optional, intent(in)  :: trans
    logical,          optional, intent(in)  :: iter_refine
  end subroutine

  subroutine MATRIX_SOLVE_MAT(this, rhs, x, trans, iter_refine)
    import MATRIX_TYPE
    class(MATRIX_TYPE),         intent(in)  :: this
    TT,                         intent(in)  :: rhs(:,:)
    TT,                         intent(out) :: x(:,:)
    character(len=1), optional, intent(in)  :: trans
    logical,          optional, intent(in)  :: iter_refine
  end subroutine

  subroutine MATRIX_TO_DENSE(this, d, i0, j0)
    import MATRIX_TYPE, DENSE_TYPE
    class(MATRIX_TYPE), intent(in)    :: this
    type(DENSE_TYPE),   intent(inout) :: d
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
  end subroutine

  subroutine MATRIX_TO_SPARSE(this, sb, i0, j0, struct, drop_zeros)
    import MATRIX_TYPE, SPBUILD_TYPE
    class(MATRIX_TYPE), intent(in)    :: this
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
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

#undef MATRIX_TYPE
#undef MATRIX_ALC_TYPE
#undef MATRIX_INIT
#undef MATRIX_DESTRUCT
#undef MATRIX_RESET
#undef MATRIX_SCALE
#undef MATRIX_MUL_VEC
#undef MATRIX_MUL_MAT
#undef MATRIX_FACTORIZE
#undef MATRIX_SOLVE_VEC
#undef MATRIX_SOLVE_MAT
#undef MATRIX_TO_DENSE
#undef MATRIX_TO_SPARSE

#undef DENSE_TYPE
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
#undef DENSE_TO_DENSE
#undef DENSE_TO_SPARSE
#undef DENSE_TO_REAL
#undef DENSE_TO_CMPLX
#undef DENSE_TRANSPOSE
#undef DENSE_TRANSPOSE2
#undef DENSE_ADD_DENSE
#undef DENSE_ADD_DENSE3
#undef DENSE_ADD_SPARSE
#undef DENSE_ADD_SPARSE3
#undef DENSE_ADD_BAND
#undef DENSE_ADD_BAND3
#undef DENSE_MUL_DENSE
#undef DENSE_EIG
#undef DENSE_OUTPUT

#undef SPARSE_TYPE
#undef SPARSE_INIT
#undef SPARSE_DESTRUCT
#undef SPARSE_RESET
#undef SPARSE_SCALE
#undef SPARSE_MUL_VEC
#undef SPARSE_MUL_MAT
#undef SPARSE_FACTORIZE
#undef SPARSE_SOLVE_VEC
#undef SPARSE_SOLVE_MAT
#undef SPARSE_TO_DENSE
#undef SPARSE_TO_SPARSE
#undef SPARSE_TO_REAL
#undef SPARSE_TO_CMPLX
#undef SPARSE_TRANSPOSE
#undef SPARSE_TRANSPOSE2
#undef SPARSE_ADD_SPARSE
#undef SPARSE_ADD_SPARSE3
#undef SPARSE_ADD_BAND
#undef SPARSE_ADD_BAND3
#undef SPARSE_MUL_SPARSE
#undef SPARSE_OUTPUT

#undef SPBUILD_TYPE
#undef SPBUILD_INIT
#undef SPBUILD_DESTRUCT
#undef SPBUILD_RESET
#undef SPBUILD_SET
#undef SPBUILD_SET_ROW
#undef SPBUILD_ADD
#undef SPBUILD_SORT
#undef SPBUILD_LOAD
#undef SPBUILD_SAVE

#undef BAND_TYPE
#undef BAND_INIT
#undef BAND_DESTRUCT
#undef BAND_RESET
#undef BAND_SCALE
#undef BAND_MUL_VEC
#undef BAND_MUL_MAT
#undef BAND_FACTORIZE
#undef BAND_SOLVE_VEC
#undef BAND_SOLVE_MAT
#undef BAND_TO_DENSE
#undef BAND_TO_SPARSE
#undef BAND_TO_REAL
#undef BAND_TO_CMPLX
#undef BAND_ADD_BAND
#undef BAND_ADD_BAND3
#undef BAND_OUTPUT

#undef HESSENBERG_TYPE
#undef HESSENBERG_INIT
#undef HESSENBERG_DESTRUCT
#undef HESSENBERG_RESET
#undef HESSENBERG_SCALE
#undef HESSENBERG_MUL_VEC
#undef HESSENBERG_MUL_MAT
#undef HESSENBERG_FACTORIZE
#undef HESSENBERG_SOLVE_VEC
#undef HESSENBERG_SOLVE_MAT
#undef HESSENBERG_TO_DENSE
#undef HESSENBERG_TO_SPARSE
#undef HESSENBERG_TO_TRIANG
#undef HESSENBERG_IS_TRIANG
#undef HESSENBERG_TO_REAL
#undef HESSENBERG_TO_CMPLX

#undef TRIANG_TYPE
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
#undef TRIANG_TO_DENSE
#undef TRIANG_TO_SPARSE
#undef TRIANG_TO_REAL
#undef TRIANG_TO_CMPLX

#undef BLOCK_TYPE
#undef BLOCK_INIT
#undef BLOCK_DESTRUCT
#undef BLOCK_RESET
#undef BLOCK_SCALE
#undef BLOCK_MUL_VEC
#undef BLOCK_MUL_MAT
#undef BLOCK_FACTORIZE
#undef BLOCK_SOLVE_VEC
#undef BLOCK_SOLVE_MAT
#undef BLOCK_TO_DENSE
#undef BLOCK_TO_SPARSE
#undef BLOCK_GET_GLOBAL_IDX
#undef BLOCK_GET_LOCAL_IDX
#undef BLOCK_SET_BLOCK_TYPE
#undef BLOCK_GET_DENSE
#undef BLOCK_GET_SPARSE
#undef BLOCK_GET_BAND
#undef BLOCK_GET_BLOCK
