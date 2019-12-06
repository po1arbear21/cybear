#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ARRAY_TYPE           CONCAT(array,T)
#define VECTOR_TYPE          CONCAT(vector,T)
#define DENSE_TYPE           CONCAT(dense,T)
#define SPARSE_TYPE          CONCAT(sparse,T)
#define SPBUILD_TYPE         CONCAT(spbuild,T)
#define BAND_TYPE            CONCAT(band,T)
#define HESSENBERG_TYPE      CONCAT(hessenberg,T)
#define TRIANG_TYPE          CONCAT(triang,T)
#define BLOCK_TYPE           CONCAT(block,T)

#define TRIANG_INIT      CONCAT3(triang,T,init)
#define TRIANG_DESTRUCT  CONCAT3(triang,T,destruct)
#define TRIANG_RESET     CONCAT3(triang,T,reset)
#define TRIANG_SCALE     CONCAT3(triang,T,scale)
#define TRIANG_MUL_VEC   CONCAT3(triang,T,mul_vec)
#define TRIANG_MUL_MAT   CONCAT3(triang,T,mul_mat)
#define TRIANG_FACTORIZE CONCAT3(triang,T,factorize)
#define TRIANG_SOLVE_VEC CONCAT3(triang,T,solve_vec)
#define TRIANG_SOLVE_MAT CONCAT3(triang,T,solve_mat)
#define TRIANG_TO_DENSE  CONCAT3(triang,T,to_dense)
#define TRIANG_TO_SPARSE CONCAT3(triang,T,to_sparse)
#ifdef TCMPLX
#define TRIANG_TO_REAL   CONCAT3(triang,T,to_real)
#else
#define TRIANG_TO_CMPLX  CONCAT3(triang,T,to_cmplx)
#endif

subroutine TRIANG_INIT(this, nrows, upper)
  !! Initialize triangular matrix

  class(TRIANG_TYPE), intent(out) :: this
    !! Triangular matrix
  integer,            intent(in)  :: nrows
    !! Number of rows
  logical,            intent(in)  :: upper
    !! Upper or lower triangular matrix

  ! init base
  call this%init("", nrows)

  ! set upper/lower flag
  this%upper = upper

  ! allocate memory
  allocate (this%d(nrows,nrows))
  this%d = 0
end subroutine

subroutine TRIANG_DESTRUCT(this)
  !! Destruct triangular matrix by deallocating its memory.

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix

  if (allocated(this%d)) deallocate(this%d)
end subroutine

subroutine TRIANG_RESET(this, only_factorization)
  !! Reset triangular matrix to zero

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix
  logical, optional,  intent(in)    :: only_factorization
    !! ignored

  this%d = 0
end subroutine

subroutine TRIANG_SCALE(this, fact)
  !! Scale triangular matrix by a factor.

  class(TRIANG_TYPE), intent(inout) :: this
    !! Triangular matrix
  TT,                 intent(in)    :: fact
    !! Scaling factor

  ! local variables
  integer :: i

  if (this%upper) then
    do i = 1, this%nrows
      this%d(i,i:this%nrows) = this%d(i,i:this%nrows) * fact
    end do
  else
    do i = 1, this%nrows
      this%d(i,1:i) = this%d(i,1:i) * fact
    end do
  end if
end subroutine

.....

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
#undef DENSE_TYPE
#undef SPARSE_TYPE
#undef SPBUILD_TYPE
#undef BAND_TYPE
#undef HESSENBERG_TYPE
#undef TRIANG_TYPE
#undef BLOCK_TYPE

#undef TRIANG_INIT
#undef TRIANG_DESTRUCT
#undef TRIANG_RESET
#undef TRIANG_SCALE
#undef TRIANG_MUL_VEC
#undef TRIANG_MUL_MAT
#undef TRIANG_FACTORIZE
#undef TRIANG_SOLVE_VEC
#undef TRIANG_SOLVE_MAT
#undef TRIANG_TO_DENSE
#undef TRIANG_TO_SPARSE
#undef TRIANG_TO_TRIANG
#undef TRIANG_TO_REAL
#undef TRIANG_TO_CMPLX
