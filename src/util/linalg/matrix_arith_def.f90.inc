interface matrix_add
  module subroutine add_band2_{}T{}(b1, b2, fact)
    type(band_{}T),       intent(in)    :: b1
    type(band_{}T),       intent(inout) :: b2
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_band3_{}T{}(b1, b2, b3, fact1, fact2)
    type(band_{}T),       intent(in)  :: b1
    type(band_{}T),       intent(in)  :: b2
    type(band_{}T),       intent(out) :: b3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine

  module subroutine add_band_dense_{}T{}(b, d, fact)
    type(band_{}T),       intent(in)    :: b
    type(dense_{}T),      intent(inout) :: d
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_band_dense2_{}T{}(b1, d2, d3, fact1, fact2)
    type(band_{}T),       intent(in)  :: b1
    type(dense_{}T),      intent(in)  :: d2
    type(dense_{}T),      intent(out) :: d3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine

  module subroutine add_band_sparse_{}T{}(b, s, fact)
    type(band_{}T),       intent(in)    :: b
    type(sparse_{}T),     intent(inout) :: s
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_band_sparse2_{}T{}(b1, s2, s3, fact1, fact2)
    type(band_{}T),       intent(in)  :: b1
    type(sparse_{}T),     intent(in)  :: s2
    type(sparse_{}T),     intent(out) :: s3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine

  module subroutine add_dense2_{}T{}(d1, d2, fact)
    type(dense_{}T),      intent(in)    :: d1
    type(dense_{}T),      intent(inout) :: d2
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_dense3_{}T{}(d1, d2, d3, fact1, fact2)
    type(dense_{}T),      intent(in)  :: d1
    type(dense_{}T),      intent(in)  :: d2
    type(dense_{}T),      intent(out) :: d3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine

  module subroutine add_sparse_dense_{}T{}(s, d, fact)
    type(sparse_{}T),     intent(in)    :: s
    type(dense_{}T),      intent(inout) :: d
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_sparse_dense2_{}T{}(s1, d2, d3, fact1, fact2)
    type(sparse_{}T),     intent(in)  :: s1
    type(dense_{}T),      intent(in)  :: d2
    type(dense_{}T),      intent(out) :: d3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine

  module subroutine add_sparse2_{}T{}(s1, s2, fact)
    type(sparse_{}T),     intent(in)    :: s1
    type(sparse_{}T),     intent(inout) :: s2
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine add_sparse3_{}T{}(s1, s2, s3, fact1, fact2)
    type(sparse_{}T),     intent(in)  :: s1
    type(sparse_{}T),     intent(in)  :: s2
    type(sparse_{}T),     intent(out) :: s3
    m4_type(T), optional, intent(in)  :: fact1
    m4_type(T), optional, intent(in)  :: fact2
  end subroutine
end interface

interface matrix_diag
  module subroutine diag_arr_band_{}T{}(arr_in, b_out)
    m4_type(T),     intent(in)  :: arr_in(:)
    type(band_{}T), intent(out) :: b_out
  end subroutine

  module subroutine diag_arr_dense_{}T{}(arr_in, d_out)
    m4_type(T),      intent(in)  :: arr_in(:)
    type(dense_{}T), intent(out) :: d_out
  end subroutine

  module subroutine diag_arr_sparse_{}T{}(arr_in, s_out)
    m4_type(T),       intent(in)  :: arr_in(:)
    type(sparse_{}T), intent(out) :: s_out
  end subroutine

  module subroutine diag_band_arr_{}T{}(b_in, arr_out)
    type(band_{}T), intent(in)  :: b_in
    m4_type(T),     intent(out) :: arr_out(:)
  end subroutine

  module subroutine diag_dense_arr_{}T{}(d_in, arr_out)
    type(dense_{}T), intent(in)  :: d_in
    m4_type(T),      intent(out) :: arr_out(:)
  end subroutine

  module subroutine diag_sparse_arr_{}T{}(s_in, arr_out)
    type(sparse_{}T), intent(in)  :: s_in
    m4_type(T),       intent(out) :: arr_out(:)
  end subroutine
end interface

interface matrix_approx
  module subroutine approx_sparse_{}T{}(thres, s1, s2)
    real,             intent(in)          :: thres
    type(sparse_{}T), intent(in),  target :: s1
    type(sparse_{}T), intent(out), target :: s2
  end subroutine
end interface
