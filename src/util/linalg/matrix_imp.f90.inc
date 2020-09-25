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

subroutine MATRIX_INIT(this, dum, nrows, ncols)
  !! Initialize matrix.

  class(MATRIX_TYPE), intent(out) :: this
    !! Matrix
  character(*),       intent(in)  :: dum
    !! Dummy, unused
  integer,            intent(in)  :: nrows
    !! Number of rows
  integer, optional,  intent(in)  :: ncols
    !! Optional number of columns (default: ncols = nrows)

  IGNORE(dum)

  this%nrows = nrows
  this%ncols = nrows
  if (present(ncols)) this%ncols = ncols

  ASSERT(this%nrows > 0)
  ASSERT(this%ncols > 0)

  this%factorized = .false.
end subroutine

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
