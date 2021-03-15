#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)

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

#define ADD_BAND2            CONCAT(add_band2,T)
#define ADD_BAND3            CONCAT(add_band3,T)
#define ADD_BAND_DENSE       CONCAT(add_band_dense,T)
#define ADD_BAND_DENSE2      CONCAT(add_band_dense2,T)
#define ADD_BAND_SPARSE      CONCAT(add_band_sparse,T)
#define ADD_BAND_SPARSE2     CONCAT(add_band_sparse2,T)
#define ADD_DENSE2           CONCAT(add_dense2,T)
#define ADD_DENSE3           CONCAT(add_dense3,T)
#define ADD_SPARSE_DENSE     CONCAT(add_sparse_dense,T)
#define ADD_SPARSE_DENSE2    CONCAT(add_sparse_dense2,T)
#define ADD_SPARSE2          CONCAT(add_sparse2,T)
#define ADD_SPARSE3          CONCAT(add_sparse3,T)

#define DIAG_ARR_BAND        CONCAT(diag_arr_band,T)
#define DIAG_ARR_DENSE       CONCAT(diag_arr_dense,T)
#define DIAG_ARR_SPARSE      CONCAT(diag_arr_sparse,T)
#define DIAG_BAND_ARR        CONCAT(diag_band_arr,T)
#define DIAG_DENSE_ARR       CONCAT(diag_dense_arr,T)
#define DIAG_SPARSE_ARR      CONCAT(diag_sparse_arr,T)

#define APPROX_SPARSE        CONCAT(approx_sparse,T)

interface matrix_add
  module subroutine ADD_BAND2(b1, b2, fact)
    type(BAND_TYPE), intent(in)    :: b1
    type(BAND_TYPE), intent(inout) :: b2
    TT, optional,    intent(in)    :: fact
  end subroutine

  module subroutine ADD_BAND3(b1, b2, b3, fact1, fact2)
    type(BAND_TYPE), intent(in)  :: b1
    type(BAND_TYPE), intent(in)  :: b2
    type(BAND_TYPE), intent(out) :: b3
    TT, optional,    intent(in)  :: fact1
    TT, optional,    intent(in)  :: fact2
  end subroutine

  module subroutine ADD_BAND_DENSE(b, d, fact)
    type(BAND_TYPE),  intent(in)    :: b
    type(DENSE_TYPE), intent(inout) :: d
    TT, optional,     intent(in)    :: fact
  end subroutine

  module subroutine ADD_BAND_DENSE2(b1, d2, d3, fact1, fact2)
    type(BAND_TYPE),  intent(in)  :: b1
    type(DENSE_TYPE), intent(in)  :: d2
    type(DENSE_TYPE), intent(out) :: d3
    TT, optional,     intent(in)  :: fact1
    TT, optional,     intent(in)  :: fact2
  end subroutine

  module subroutine ADD_BAND_SPARSE(b, s, fact)
    type(BAND_TYPE),   intent(in)    :: b
    type(SPARSE_TYPE), intent(inout) :: s
    TT, optional,      intent(in)    :: fact
  end subroutine

  module subroutine ADD_BAND_SPARSE2(b1, s2, s3, fact1, fact2)
    type(BAND_TYPE),   intent(in)  :: b1
    type(SPARSE_TYPE), intent(in)  :: s2
    type(SPARSE_TYPE), intent(out) :: s3
    TT, optional,      intent(in)  :: fact1
    TT, optional,      intent(in)  :: fact2
  end subroutine

  module subroutine ADD_DENSE2(d1, d2, fact)
    type(DENSE_TYPE), intent(in)    :: d1
    type(DENSE_TYPE), intent(inout) :: d2
    TT, optional,     intent(in)    :: fact
  end subroutine

  module subroutine ADD_DENSE3(d1, d2, d3, fact1, fact2)
    type(DENSE_TYPE), intent(in)  :: d1
    type(DENSE_TYPE), intent(in)  :: d2
    type(DENSE_TYPE), intent(out) :: d3
    TT, optional,     intent(in)  :: fact1
    TT, optional,     intent(in)  :: fact2
  end subroutine

  module subroutine ADD_SPARSE_DENSE(s, d, fact)
    type(SPARSE_TYPE), intent(in)    :: s
    type(DENSE_TYPE),  intent(inout) :: d
    TT, optional,      intent(in)    :: fact
  end subroutine

  module subroutine ADD_SPARSE_DENSE2(s1, d2, d3, fact1, fact2)
    type(SPARSE_TYPE), intent(in)  :: s1
    type(DENSE_TYPE),  intent(in)  :: d2
    type(DENSE_TYPE),  intent(out) :: d3
    TT, optional,      intent(in)  :: fact1
    TT, optional,      intent(in)  :: fact2
  end subroutine

  module subroutine ADD_SPARSE2(s1, s2, fact)
    type(SPARSE_TYPE), intent(in)    :: s1
    type(SPARSE_TYPE), intent(inout) :: s2
    TT, optional,      intent(in)    :: fact
  end subroutine

  module subroutine ADD_SPARSE3(s1, s2, s3, fact1, fact2)
    type(SPARSE_TYPE), intent(in)  :: s1
    type(SPARSE_TYPE), intent(in)  :: s2
    type(SPARSE_TYPE), intent(out) :: s3
    TT, optional,      intent(in)  :: fact1
    TT, optional,      intent(in)  :: fact2
  end subroutine
end interface

interface matrix_diag
  module subroutine DIAG_ARR_BAND(arr_in, b_out)
    TT,              intent(in)  :: arr_in(:)
    type(BAND_TYPE), intent(out) :: b_out
  end subroutine

  module subroutine DIAG_ARR_DENSE(arr_in, d_out)
    TT,               intent(in)  :: arr_in(:)
    type(DENSE_TYPE), intent(out) :: d_out
  end subroutine

  module subroutine DIAG_ARR_SPARSE(arr_in, s_out)
    TT,                intent(in)  :: arr_in(:)
    type(SPARSE_TYPE), intent(out) :: s_out
  end subroutine

  module subroutine DIAG_BAND_ARR(b_in, arr_out)
    type(BAND_TYPE), intent(in)  :: b_in
    TT,              intent(out) :: arr_out(:)
  end subroutine

  module subroutine DIAG_DENSE_ARR(d_in, arr_out)
    type(DENSE_TYPE), intent(in)  :: d_in
    TT,               intent(out) :: arr_out(:)
  end subroutine

  module subroutine DIAG_SPARSE_ARR(s_in, arr_out)
    type(SPARSE_TYPE), intent(in)  :: s_in
    TT,                intent(out) :: arr_out(:)
  end subroutine
end interface

interface matrix_approx
  module subroutine APPROX_SPARSE(thres, s1, s2)
    real,              intent(in)          :: thres
    type(SPARSE_TYPE), intent(in),  target :: s1
    type(SPARSE_TYPE), intent(out), target :: s2
  end subroutine
end interface

#undef T
#undef TT
#undef TCMPLX

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

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

#undef ADD_BAND2
#undef ADD_BAND3
#undef ADD_BAND_DENSE
#undef ADD_BAND_DENSE2
#undef ADD_BAND_SPARSE
#undef ADD_BAND_SPARSE2
#undef ADD_DENSE2
#undef ADD_DENSE3
#undef ADD_SPARSE_DENSE
#undef ADD_SPARSE_DENSE2
#undef ADD_SPARSE2
#undef ADD_SPARSE3

#undef DIAG_ARR_BAND
#undef DIAG_ARR_DENSE
#undef DIAG_ARR_SPARSE
#undef DIAG_BAND_ARR
#undef DIAG_DENSE_ARR
#undef DIAG_SPARSE_ARR

#undef APPROX_SPARSE
