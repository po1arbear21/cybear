module jacval_snapshot_m
  !! Serialize a (state, claimed Jacobian, CSD reference Jacobian, residual)
  !! tuple to disk as a pair of files:
  !!
  !!   <stem>.json  -- manifest (schema in jacval.md section 2)
  !!   <stem>.bin   -- concatenated raw little-endian float64 / int64 arrays
  !!
  !! The split avoids adding HDF5 as a cybear dependency while remaining
  !! trivially loadable in Python (numpy.frombuffer at known offsets).
  !!
  !! Internally the module evaluates both the analytical Jacobian (via
  !! prob%eval_jac) and the CSD reference (via prob%eval_cmplx column-by-
  !! column) on a dense work array, then encodes the union sparsity pattern
  !! as CSR. Diagnostic meshes only (intended <= 500 unknowns); dense
  !! storage at that scale costs at most a few MB.

  use, intrinsic :: iso_fortran_env, only: real64, int64, output_unit
  use jacobian_validator_m, only: jacobian_problem_cmplx, validate_complex_step

  implicit none

  private
  public :: jacval_snapshot_dump
  public :: index_map_t
  public :: physics_flags_t

  type :: index_map_t
    !! Per-row metadata that resolves unknown index -> (equation kind, cell).
    !!
    !! equation(i)  is an index into labels(:), 1-based. e.g. labels=["psi","n","p"]
    !!              and equation(7)=2 means row 7 is the n-equation.
    !! cell(i)      is the grid cell associated with row i.
    integer, allocatable          :: equation(:)
    integer, allocatable          :: cell(:)
    character(len=16), allocatable :: labels(:)
  end type

  type :: physics_flags_t
    !! Which physics terms were active at snapshot time. Toggle-diff compares
    !! two snapshots whose flags differ in exactly one entry.
    logical :: thermionic    = .true.
    logical :: tunneling     = .true.
    logical :: ifbl          = .true.
    logical :: drift_contact = .true.
    logical :: srh           = .true.
    logical :: auger         = .true.
  end type

  ! Magic word written at the head of the .bin file so the Python loader
  ! can detect format / endianness mismatches before reading garbage. Eight
  ! bytes, little-endian: "JACVAL\0\1" (the trailing byte is the format version).
  character(len=8), parameter :: BIN_MAGIC = "JACVAL"//char(0)//char(1)

  integer, parameter :: FORMAT_VERSION = 1

contains

  subroutine jacval_snapshot_dump(stem, prob, x, im, flags, bias_V, build_hash, &
                                  state_psi, state_n, state_p, mesh_x,          &
                                  doping_ND, doping_NA, struct_zero)
    !! Serialize a snapshot. Required: stem, prob, x, im, flags, bias_V,
    !! build_hash. State / mesh / doping arrays are optional: tests against
    !! generic problems (e.g., the Bratu fixture) skip them; cybear DD runs
    !! pass them through.

    character(*),                  intent(in)           :: stem
    class(jacobian_problem_cmplx), intent(in)           :: prob
    real,                          intent(in)           :: x(:)
    type(index_map_t),             intent(in)           :: im
    type(physics_flags_t),         intent(in)           :: flags
    real,                          intent(in)           :: bias_V
    character(*),                  intent(in)           :: build_hash
    real,                          intent(in), optional :: state_psi(:)
    real,                          intent(in), optional :: state_n(:)
    real,                          intent(in), optional :: state_p(:)
    real,                          intent(in), optional :: mesh_x(:)
    real,                          intent(in), optional :: doping_ND(:)
    real,                          intent(in), optional :: doping_NA(:)
    real,                          intent(in), optional :: struct_zero
    !! Entries with max(|J_claim|, |J_csd|) <= struct_zero are treated as
    !! structurally zero. Default: 1e-14 * max(|J_claim|_inf, |J_csd|_inf).

    real,    allocatable :: J_claim(:,:), J_csd(:,:), F(:)
    real,    allocatable :: claim_vals(:), csd_vals_in_claim(:)
    real,    allocatable :: extra_vals(:)
    integer, allocatable :: claim_row_ptr(:), claim_col_idx(:)
    integer, allocatable :: extra_row_idx(:), extra_col_idx(:)

    real    :: max_err, struct_zero_eff, jmax
    logical :: ok
    integer :: n, nnz_claim, n_extras
    integer :: u_bin, u_json, ios
    integer(int64) :: offset
    character(len=:), allocatable :: bin_filename
    character(len=:), allocatable :: json_filename

    n = size(x)
    allocate (J_claim(n,n), J_csd(n,n), F(n))

    ! Reuse the optional out-args of validate_complex_step so the snapshot
    ! CSD reference is byte-identical to what the validator produces.
    call validate_complex_step(prob, x, ok, max_err,                         &
                               J_claim_out = J_claim, J_ref_out = J_csd)

    call prob%eval_real(x, F)

    jmax = max(maxval(abs(J_claim)), maxval(abs(J_csd)))
    if (present(struct_zero)) then
      struct_zero_eff = struct_zero
    else
      struct_zero_eff = 1.0e-14_real64 * max(jmax, 1.0_real64)
    end if

    call build_csr_from_claim(J_claim, J_csd, struct_zero_eff,               &
                              claim_vals, csd_vals_in_claim,                 &
                              claim_row_ptr, claim_col_idx, nnz_claim)

    call collect_csd_extras(J_claim, J_csd, struct_zero_eff,                 &
                            extra_vals, extra_row_idx, extra_col_idx, n_extras)

    bin_filename  = trim(stem) // ".bin"
    json_filename = trim(stem) // ".json"

    open (newunit = u_bin, file = bin_filename, access = "stream",           &
          form = "unformatted", status = "replace", iostat = ios)
    if (ios /= 0) then
      write (output_unit, "(A,A,A,I0)")                                      &
        "jacval_snapshot_dump: open failed for '", bin_filename, "' iostat=", ios
      return
    end if

    write (u_bin) BIN_MAGIC
    offset = int(len(BIN_MAGIC), kind=int64)

    open (newunit = u_json, file = json_filename, status = "replace",        &
          action = "write", iostat = ios)
    if (ios /= 0) then
      close (u_bin)
      write (output_unit, "(A,A,A,I0)")                                      &
        "jacval_snapshot_dump: open failed for '", json_filename, "' iostat=", ios
      return
    end if

    call json_header(u_json, bias_V, build_hash, n, size(im%cell),           &
                     bin_filename, im%labels, flags)

    ! Write each array to the binary file and record its offset/nbytes in the
    ! manifest. Order is irrelevant for correctness (loader uses offsets) but
    ! we keep it stable for diff-friendliness.
    write (u_json, "(A)") "  ""arrays"": {"

    call write_array_real (u_bin, u_json, "state.x", x, offset, last=.false.)
    if (present(state_psi)) call write_array_real (u_bin, u_json, "state.psi", state_psi, offset, last=.false.)
    if (present(state_n))   call write_array_real (u_bin, u_json, "state.n",   state_n,   offset, last=.false.)
    if (present(state_p))   call write_array_real (u_bin, u_json, "state.p",   state_p,   offset, last=.false.)
    if (present(mesh_x))    call write_array_real (u_bin, u_json, "mesh.x",    mesh_x,    offset, last=.false.)
    if (present(doping_ND)) call write_array_real (u_bin, u_json, "doping.ND", doping_ND, offset, last=.false.)
    if (present(doping_NA)) call write_array_real (u_bin, u_json, "doping.NA", doping_NA, offset, last=.false.)

    call write_array_real (u_bin, u_json, "J_claim.values",  claim_vals,         offset, last=.false.)
    call write_array_int  (u_bin, u_json, "J_claim.row_ptr", claim_row_ptr,      offset, last=.false.)
    call write_array_int  (u_bin, u_json, "J_claim.col_idx", claim_col_idx,      offset, last=.false.)
    call write_array_real (u_bin, u_json, "J_csd.values",    csd_vals_in_claim,  offset, last=.false.)

    ! Extras are optional; only emit if any.
    if (n_extras > 0) then
      call write_array_real (u_bin, u_json, "J_csd_extras.values",  extra_vals,    offset, last=.false.)
      call write_array_int  (u_bin, u_json, "J_csd_extras.row_idx", extra_row_idx, offset, last=.false.)
      call write_array_int  (u_bin, u_json, "J_csd_extras.col_idx", extra_col_idx, offset, last=.false.)
    end if

    call write_array_real (u_bin, u_json, "residual",            F,           offset, last=.false.)
    call write_array_int  (u_bin, u_json, "index_map.equation",  im%equation, offset, last=.false.)
    call write_array_int  (u_bin, u_json, "index_map.cell",      im%cell,     offset, last=.true.)

    write (u_json, "(A)") "  }"
    write (u_json, "(A)") "}"

    close (u_bin)
    close (u_json)

    write (output_unit, "(A,A,A,I0,A,I0,A)")                                 &
      "jacval_snapshot_dump: wrote ", trim(stem), " (n=", n,                  &
      ", nnz=", nnz_claim, ")"
  end subroutine

  subroutine build_csr_from_claim(J_claim, J_csd, struct_zero,               &
                                  claim_vals, csd_vals_in_claim,             &
                                  row_ptr, col_idx, nnz)
    !! Encode the claim nonzero entries as CSR; for each such entry, also
    !! record the matching CSD value (which may be different / zero / a bug).
    real,    intent(in)  :: J_claim(:,:), J_csd(:,:)
    real,    intent(in)  :: struct_zero
    real,    allocatable, intent(out) :: claim_vals(:), csd_vals_in_claim(:)
    integer, allocatable, intent(out) :: row_ptr(:), col_idx(:)
    integer,              intent(out) :: nnz

    integer :: n, i, j, k

    n = size(J_claim, 1)
    allocate (row_ptr(n + 1))

    nnz = 0
    do i = 1, n
      do j = 1, n
        if (abs(J_claim(i,j)) > struct_zero) nnz = nnz + 1
      end do
    end do

    allocate (claim_vals(nnz), csd_vals_in_claim(nnz), col_idx(nnz))

    k = 0
    do i = 1, n
      row_ptr(i) = k
      do j = 1, n
        if (abs(J_claim(i,j)) > struct_zero) then
          k = k + 1
          claim_vals(k)        = J_claim(i,j)
          csd_vals_in_claim(k) = J_csd(i,j)
          col_idx(k)           = j - 1     ! 0-based for Python compatibility
        end if
      end do
    end do
    row_ptr(n + 1) = k
  end subroutine

  subroutine collect_csd_extras(J_claim, J_csd, struct_zero,                 &
                                vals, row_idx, col_idx, n_extras)
    !! Entries where CSD is structurally nonzero but the claim analytical
    !! Jacobian dropped them entirely. Encoded as COO (row_idx, col_idx, vals)
    !! to keep the format trivial for the typically small extras set.
    real,    intent(in)  :: J_claim(:,:), J_csd(:,:)
    real,    intent(in)  :: struct_zero
    real,    allocatable, intent(out) :: vals(:)
    integer, allocatable, intent(out) :: row_idx(:), col_idx(:)
    integer,              intent(out) :: n_extras

    integer :: n, i, j, k

    n = size(J_claim, 1)
    n_extras = 0
    do i = 1, n
      do j = 1, n
        if ( abs(J_csd(i,j))   >  struct_zero .and.                          &
             abs(J_claim(i,j)) <= struct_zero ) then
          n_extras = n_extras + 1
        end if
      end do
    end do

    allocate (vals(n_extras), row_idx(n_extras), col_idx(n_extras))
    k = 0
    do i = 1, n
      do j = 1, n
        if ( abs(J_csd(i,j))   >  struct_zero .and.                          &
             abs(J_claim(i,j)) <= struct_zero ) then
          k = k + 1
          vals(k)    = J_csd(i,j)
          row_idx(k) = i - 1
          col_idx(k) = j - 1
        end if
      end do
    end do
  end subroutine

  subroutine json_header(unit, bias_V, build_hash, n_unknowns, n_cells,      &
                         bin_filename, labels, flags)
    integer,              intent(in) :: unit
    real,                 intent(in) :: bias_V
    character(*),         intent(in) :: build_hash
    integer,              intent(in) :: n_unknowns, n_cells
    character(*),         intent(in) :: bin_filename
    character(len=16),    intent(in) :: labels(:)
    type(physics_flags_t),intent(in) :: flags

    integer :: i

    write (unit, "(A)")               "{"
    write (unit, "(A,I0,A)")          "  ""snapshot_format_version"": ", FORMAT_VERSION, ","
    write (unit, "(A,ES23.15,A)")     "  ""bias_V"": ", bias_V, ","
    write (unit, "(A,A,A)")           "  ""cybear_build_hash"": """, trim(build_hash), ""","
    write (unit, "(A,I0,A)")          "  ""n_unknowns"": ", n_unknowns, ","
    write (unit, "(A,I0,A)")          "  ""n_cells"": ", n_cells, ","
    write (unit, "(A,A,A)")           "  ""binary_file"": """, trim(bin_filename), ""","

    ! equation_labels
    write (unit, "(A)", advance="no") "  ""equation_labels"": ["
    do i = 1, size(labels)
      if (i > 1) write (unit, "(A)", advance="no") ", "
      write (unit, "(A,A,A)", advance="no") """", trim(labels(i)), """"
    end do
    write (unit, "(A)") "],"

    ! flags
    write (unit, "(A)") "  ""flags"": {"
    write (unit, "(A,A,A)") "    ""thermionic"": ",    bool(flags%thermionic),    ","
    write (unit, "(A,A,A)") "    ""tunneling"": ",     bool(flags%tunneling),     ","
    write (unit, "(A,A,A)") "    ""ifbl"": ",          bool(flags%ifbl),          ","
    write (unit, "(A,A,A)") "    ""drift_contact"": ", bool(flags%drift_contact), ","
    write (unit, "(A,A,A)") "    ""srh"": ",           bool(flags%srh),           ","
    write (unit, "(A,A)")   "    ""auger"": ",         bool(flags%auger)
    write (unit, "(A)") "  },"
  end subroutine

  subroutine write_array_real(u_bin, u_json, name, arr, offset, last)
    !! Write a real array to the binary file at the current offset, and emit
    !! its manifest entry to the JSON file. The offset argument is advanced.
    integer,        intent(in)    :: u_bin, u_json
    character(*),   intent(in)    :: name
    real,           intent(in)    :: arr(:)
    integer(int64), intent(inout) :: offset
    logical,        intent(in)    :: last        ! last array in the manifest?

    integer(int64) :: nbytes
    real(real64), allocatable :: arr64(:)

    allocate (arr64(size(arr)))
    arr64 = real(arr, kind=real64)

    nbytes = int(size(arr), kind=int64) * 8_int64
    write (u_bin) arr64
    call emit_array_entry(u_json, name, "f8", size(arr), offset, nbytes, last)
    offset = offset + nbytes
    deallocate (arr64)
  end subroutine

  subroutine write_array_int(u_bin, u_json, name, arr, offset, last)
    integer,        intent(in)    :: u_bin, u_json
    character(*),   intent(in)    :: name
    integer,        intent(in)    :: arr(:)
    integer(int64), intent(inout) :: offset
    logical,        intent(in)    :: last

    integer(int64) :: nbytes
    integer(int64), allocatable :: arr64(:)

    allocate (arr64(size(arr)))
    arr64 = int(arr, kind=int64)

    nbytes = int(size(arr), kind=int64) * 8_int64
    write (u_bin) arr64
    call emit_array_entry(u_json, name, "i8", size(arr), offset, nbytes, last)
    offset = offset + nbytes
    deallocate (arr64)
  end subroutine

  subroutine emit_array_entry(u_json, name, dtype, shape1, offset, nbytes, last)
    integer,        intent(in) :: u_json
    character(*),   intent(in) :: name, dtype
    integer,        intent(in) :: shape1
    integer(int64), intent(in) :: offset, nbytes
    logical,        intent(in) :: last
    character(len=2) :: trail
    if (last) then
      trail = ""
    else
      trail = ","
    end if
    write (u_json, "(A,A,A,A,A,A,I0,A,A,I0,A,I0,A,A)")                       &
      "    """, trim(name), """: ",                                          &
      "{""dtype"": """, trim(dtype), """, ""shape"": [", shape1, "], ",      &
      """offset_bytes"": ", offset, ", ""nbytes"": ", nbytes, "}", trim(trail)
  end subroutine

  pure function bool(b) result(s)
    logical, intent(in) :: b
    character(len=5) :: s
    if (b) then
      s = "true "
    else
      s = "false"
    end if
  end function

end module
