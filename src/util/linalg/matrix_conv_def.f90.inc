#define PASTE(X)              X
#define PASTE2(X)             PASTE(X)_
#define CONCATHELP(X, Y)      PASTE2(X)Y
#define CONCAT(X, Y)          CONCATHELP(X,Y)

#define BAND_TYPE             CONCAT(band,T)
#define BAND_PTR_TYPE         CONCAT(band_ptr,T)
#define BLOCK_TYPE            CONCAT(block,T)
#define BLOCK_PTR_TYPE        CONCAT(block_ptr,T)
#define DENSE_TYPE            CONCAT(dense,T)
#define DENSE_PTR_TYPE        CONCAT(dense_ptr,T)
#define HESSENBERG_TYPE       CONCAT(hessenberg,T)
#define HESSENBERG_PTR_TYPE   CONCAT(hessenberg_ptr,T)
#define MATRIX_TYPE           CONCAT(matrix,T)
#define MATRIX_PTR_TYPE       CONCAT(matrix_ptr,T)
#define SPARSE_TYPE           CONCAT(sparse,T)
#define SPARSE_PTR_TYPE       CONCAT(sparse_ptr,T)
#define SPBUILD_TYPE          CONCAT(spbuild,T)
#define TRIANG_TYPE           CONCAT(triang,T)
#define TRIANG_PTR_TYPE       CONCAT(triang_ptr,T)

#define SPARSE_TO_BAND        CONCAT(sparse_to_band,T)

#define BAND_TO_DENSE         CONCAT(band_to_dense,T)
#define BLOCK_TO_DENSE        CONCAT(block_to_dense,T)
#define DENSE_TO_DENSE        CONCAT(dense_to_dense,T)
#define HESSENBERG_TO_DENSE   CONCAT(hessenberg_to_dense,T)
#define SPARSE_TO_DENSE       CONCAT(sparse_to_dense,T)
#define TRIANG_TO_DENSE       CONCAT(triang_to_dense,T)

#define BLOCK_TO_SPARSE       CONCAT(block_to_sparse,T)

#define BAND_TO_SPBUILD       CONCAT(band_to_spbuild,T)
#define BLOCK_TO_SPBUILD      CONCAT(block_to_spbuild,T)
#define DENSE_TO_SPBUILD      CONCAT(dense_to_spbuild,T)
#define HESSENBERG_TO_SPBUILD CONCAT(hessenberg_to_spbuild,T)
#define SPARSE_TO_SPBUILD     CONCAT(sparse_to_spbuild,T)
#define TRIANG_TO_SPBUILD     CONCAT(triang_to_spbuild,T)

#define HESSENBERG_TO_TRIANG  CONCAT(hessenberg_to_triang,T)

! to_band
interface matrix_convert
  module subroutine SPARSE_TO_BAND(s, b, i0, j0)
    type(SPARSE_TYPE), intent(in)    :: s
    type(BAND_TYPE),   intent(inout) :: b
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine
end interface

! to_dense
interface matrix_convert
  module subroutine BAND_TO_DENSE(b, d, i0, j0)
    type(BAND_TYPE),   intent(in)    :: b
    type(DENSE_TYPE),  intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module recursive subroutine BLOCK_TO_DENSE(b, d, i0, j0)
    type(BLOCK_TYPE),  intent(in)    :: b
    type(DENSE_TYPE),  intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module subroutine DENSE_TO_DENSE(d1, d2, i0, j0)
    type(DENSE_TYPE),  intent(in)    :: d1
    type(DENSE_TYPE),  intent(inout) :: d2
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module subroutine HESSENBERG_TO_DENSE(h, d, i0, j0)
    type(HESSENBERG_TYPE), intent(in)    :: h
    type(DENSE_TYPE),      intent(inout) :: d
    integer, optional,     intent(in)    :: i0
    integer, optional,     intent(in)    :: j0
  end subroutine

  module subroutine SPARSE_TO_DENSE(s, d, i0, j0)
    type(SPARSE_TYPE), intent(in)    :: s
    type(DENSE_TYPE),  intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module subroutine TRIANG_TO_DENSE(t, d, i0, j0)
    type(TRIANG_TYPE), intent(in)    :: t
    type(DENSE_TYPE),  intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine
end interface

! to_sparse
interface matrix_convert
  module subroutine BLOCK_TO_SPARSE(b, s, struct, drop_zeros)
    type(BLOCK_TYPE),  intent(in)  :: b
    type(SPARSE_TYPE), intent(out) :: s
    logical, optional, intent(in)  :: struct(:,:)
    logical, optional, intent(in)  :: drop_zeros
  end subroutine
end interface

! to_spbuild
interface matrix_convert
  module subroutine BAND_TO_SPBUILD(b, sb, i0, j0, struct, drop_zeros)
    type(BAND_TYPE),    intent(in)    :: b
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
  end subroutine

  module subroutine BLOCK_TO_SPBUILD(b, sb, i0, j0, struct, drop_zeros)
    type(BLOCK_TYPE),   intent(in)    :: b
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
  end subroutine

  module subroutine DENSE_TO_SPBUILD(d, sb, i0, j0, struct, drop_zeros)
    type(DENSE_TYPE),   intent(in)    :: d
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
  end subroutine

  module subroutine HESSENBERG_TO_SPBUILD(h, sb, i0, j0, struct, drop_zeros)
    type(HESSENBERG_TYPE), intent(in)    :: h
    type(SPBUILD_TYPE),    intent(inout) :: sb
    integer, optional,     intent(in)    :: i0
    integer, optional,     intent(in)    :: j0
    logical, optional,     intent(in)    :: struct(:,:)
    logical, optional,     intent(in)    :: drop_zeros
  end subroutine

  module subroutine SPARSE_TO_SPBUILD(s, sb, i0, j0, struct, drop_zeros)
    type(SPARSE_TYPE),  intent(in)    :: s
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
  end subroutine

  module subroutine TRIANG_TO_SPBUILD(t, sb, i0, j0, struct, drop_zeros)
    type(TRIANG_TYPE),  intent(in)    :: t
    type(SPBUILD_TYPE), intent(inout) :: sb
    integer, optional,  intent(in)    :: i0
    integer, optional,  intent(in)    :: j0
    logical, optional,  intent(in)    :: struct(:,:)
    logical, optional,  intent(in)    :: drop_zeros
  end subroutine
end interface

! to_triang
interface matrix_convert
  module subroutine HESSENBERG_TO_TRIANG(h, t)
    type(HESSENBERG_TYPE), intent(in)  :: h
    type(TRIANG_TYPE),     intent(out) :: t
  end subroutine
end interface

! to_real
#ifdef TCMPLX
interface matrix_convert
  module subroutine to_real_band(bc_in, br_out)
    type(band_cmplx), intent(in)  :: bc_in
    type(band_real),  intent(out) :: br_out
  end subroutine

  module subroutine to_real_dense(dc_in, dr_out)
    type(dense_cmplx), intent(in)  :: dc_in
    type(dense_real),  intent(out) :: dr_out
  end subroutine

  module subroutine to_real_hessenberg(hc_in, hr_out)
    type(hessenberg_cmplx), intent(in)  :: hc_in
    type(hessenberg_real),  intent(out) :: hr_out
  end subroutine

  module subroutine to_real_sparse(sc_in, sr_out)
    type(sparse_cmplx), intent(in)  :: sc_in
    type(sparse_real),  intent(out) :: sr_out
  end subroutine

  module subroutine to_real_triang(tc_in, tr_out)
    type(triang_cmplx), intent(in)  :: tc_in
    type(triang_real),  intent(out) :: tr_out
  end subroutine
end interface
#endif

! to_cmplx
#ifndef TCMPLX
interface matrix_convert
  module subroutine to_cmplx_band(br_in, bc_out)
    type(band_real),  intent(in)  :: br_in
    type(band_cmplx), intent(out) :: bc_out
  end subroutine

  module subroutine to_cmplx_dense(dr_in, dc_out)
    type(dense_real),  intent(in)  :: dr_in
    type(dense_cmplx), intent(out) :: dc_out
  end subroutine

  module subroutine to_cmplx_hessenberg(hr_in, hc_out)
    type(hessenberg_real),  intent(in)  :: hr_in
    type(hessenberg_cmplx), intent(out) :: hc_out
  end subroutine

  module subroutine to_cmplx_sparse(sr_in, sc_out)
    type(sparse_real),  intent(in)  :: sr_in
    type(sparse_cmplx), intent(out) :: sc_out
  end subroutine

  module subroutine to_cmplx_triang(tr_in, tc_out)
    type(triang_real),  intent(in)  :: tr_in
    type(triang_cmplx), intent(out) :: tc_out
  end subroutine
end interface
#endif

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

#undef SPARSE_TO_BAND

#undef BAND_TO_DENSE
#undef BLOCK_TO_DENSE
#undef DENSE_TO_DENSE
#undef HESSENBERG_TO_DENSE
#undef SPARSE_TO_DENSE
#undef TRIANG_TO_DENSE

#undef BLOCK_TO_SPARSE

#undef BAND_TO_SPBUILD
#undef BLOCK_TO_SPBUILD
#undef DENSE_TO_SPBUILD
#undef HESSENBERG_TO_SPBUILD
#undef SPARSE_TO_SPBUILD
#undef TRIANG_TO_SPBUILD

#undef HESSENBERG_TO_TRIANG
