! to_band
interface matrix_convert
  module subroutine sparse_to_band_{}T{}(s, b, i0, j0)
    type(sparse_{}T),  intent(in)    :: s
    type(band_{}T),    intent(inout) :: b
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine
end interface

! to_dense
interface matrix_convert
  module subroutine band_to_dense_{}T{}(b, d, i0, j0)
    type(band_{}T),    intent(in)    :: b
    type(dense_{}T),   intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module recursive subroutine block_to_dense_{}T{}(b, d, i0, j0)
    class(block_{}T),  intent(in)    :: b
    type(dense_{}T),   intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module subroutine dense_to_dense_{}T{}(d1, d2, i0, j0)
    type(dense_{}T),   intent(in)    :: d1
    type(dense_{}T),   intent(inout) :: d2
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine

  module subroutine sparse_to_dense_{}T{}(s, d, i0, j0)
    type(sparse_{}T),  intent(in)    :: s
    type(dense_{}T),   intent(inout) :: d
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
  end subroutine
end interface

! to_sparse
interface matrix_convert
  module subroutine block_to_sparse_{}T{}(b, s, struct, drop_zeros)
    class(block_{}T),  intent(in)  :: b
    type(sparse_{}T),  intent(out) :: s
    logical, optional, intent(in)  :: struct(:,:)
    logical, optional, intent(in)  :: drop_zeros
  end subroutine
end interface

! to_spbuild
interface matrix_convert
  module subroutine band_to_spbuild_{}T{}(b, sb, i0, j0, struct, drop_zeros)
    type(band_{}T),    intent(in)    :: b
    type(spbuild_{}T), intent(inout) :: sb
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
    logical, optional, intent(in)    :: struct(:,:)
    logical, optional, intent(in)    :: drop_zeros
  end subroutine

  module subroutine block_to_spbuild_{}T{}(b, sb, i0, j0, struct, drop_zeros)
    class(block_{}T),  intent(in)    :: b
    type(spbuild_{}T), intent(inout) :: sb
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
    logical, optional, intent(in)    :: struct(:,:)
    logical, optional, intent(in)    :: drop_zeros
  end subroutine

  module subroutine dense_to_spbuild_{}T{}(d, sb, i0, j0, struct, drop_zeros)
    type(dense_{}T),   intent(in)    :: d
    type(spbuild_{}T), intent(inout) :: sb
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
    logical, optional, intent(in)    :: struct(:,:)
    logical, optional, intent(in)    :: drop_zeros
  end subroutine

  module subroutine sparse_to_spbuild_{}T{}(s, sb, i0, j0, struct, drop_zeros)
    type(sparse_{}T),  intent(in)    :: s
    type(spbuild_{}T), intent(inout) :: sb
    integer, optional, intent(in)    :: i0
    integer, optional, intent(in)    :: j0
    logical, optional, intent(in)    :: struct(:,:)
    logical, optional, intent(in)    :: drop_zeros
  end subroutine
end interface

! convert between real and complex
interface matrix_convert
m4_ifelse(T,cmplx,{
  module subroutine to_real_band(bc_in, br_out)
    type(band_cmplx), intent(in)  :: bc_in
    type(band_real),  intent(out) :: br_out
  end subroutine

  module subroutine to_real_dense(dc_in, dr_out)
    type(dense_cmplx), intent(in)  :: dc_in
    type(dense_real),  intent(out) :: dr_out
  end subroutine

  module subroutine to_real_sparse(sc_in, sr_out)
    type(sparse_cmplx), intent(in)  :: sc_in
    type(sparse_real),  intent(out) :: sr_out
  end subroutine
  },{
  module subroutine to_cmplx_band(br_in, bc_out)
    type(band_real),  intent(in)  :: br_in
    type(band_cmplx), intent(out) :: bc_out
  end subroutine

  module subroutine to_cmplx_dense(dr_in, dc_out)
    type(dense_real),  intent(in)  :: dr_in
    type(dense_cmplx), intent(out) :: dc_out
  end subroutine

  module subroutine to_cmplx_sparse(sr_in, sc_out)
    type(sparse_real),  intent(in)  :: sr_in
    type(sparse_cmplx), intent(out) :: sc_out
  end subroutine
})
end interface
