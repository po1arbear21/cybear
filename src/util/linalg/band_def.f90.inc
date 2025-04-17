type, extends(matrix_{}T) :: band_{}T
  !! Band matrix.
  !! Storage example
  !!    / a_11  a_12  a_13       \
  !!  A=| a_21  a_22  a_23  a_24 |
  !!    |       a_32  a_33  a_34 |
  !!    \             a_43  a_44 /
  !!  here n=4, nlower=1, nupper=2
  !!
  !! Band storage used here (shifted rows)
  !!
  !!          /             a_13  a_24 \    ! <- row -2
  !!          |       a_12  a_23  a_34 |    ! <- row -1
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 0
  !!          \ a_21  a_32  a_43       /    ! <- row 1
  !!  fixme numbering should be reversed b.c. matlab uses row>0 for upper diagonals. has effect on set_diag routines.
  !!
  !! LAPACK's band storage
  !!
  !!          /             a_13  a_24 \    ! <- row 1
  !!          |       a_12  a_23  a_34 \    ! <- row 2
  !!  A_bs =  | a_11  a_22  a_33  a_44 |    ! <- row 3
  !!          \ a_21  a_32  a_43       /    ! <- row 4


  integer              :: nupper, nlower
    !! upper/lower number of diagonals (n>=0).
    !! nupper=nlower=0 <=> diagonal matrix

  m4_type(T), allocatable :: d(:,:)
    !! matrix entries in band storage format.
    !! index1: diagonal number (-nupper:+nlower)
    !! index2: column (1:ncols)

  ! factorization data
  m4_type(T), allocatable :: f(:,:)
  integer,    allocatable :: ipiv(:)

  m4_ifdef({m4_spike},{
  type(spike_params) :: sp
    !! spike solver params
  m4_type(T), allocatable :: work(:)
    !! spike work array
  })
contains
  procedure :: band_{}T{}_init
  generic   :: init      => band_{}T{}_init
  procedure :: destruct  => band_{}T{}_destruct
  procedure :: reset     => band_{}T{}_reset
  procedure :: reset_row => band_{}T{}_reset_row
  procedure :: scale     => band_{}T{}_scale
  procedure :: mul_vec   => band_{}T{}_mul_vec
  procedure :: mul_mat   => band_{}T{}_mul_mat
  procedure :: factorize => band_{}T{}_factorize
  procedure :: solve_vec => band_{}T{}_solve_vec
  procedure :: solve_mat => band_{}T{}_solve_mat
  procedure :: add_elem  => band_{}T{}_add_elem
  procedure :: get_elem  => band_{}T{}_get_elem
  procedure :: set_elem  => band_{}T{}_set_elem
  procedure :: band_{}T{}_set_diag_val
  procedure :: band_{}T{}_set_diag_arr
  generic   :: set_diag  => band_{}T{}_set_diag_val, band_{}T{}_set_diag_arr
  procedure :: output    => band_{}T{}_output
end type

type band_ptr_{}T
  type(band_{}T), pointer :: p => null()
end type

interface
  module subroutine band_{}T{}_init(this, nrows, nlower, nupper, d0)
    class(band_{}T),      intent(out) :: this
    integer,              intent(in)  :: nrows
    integer,              intent(in)  :: nlower
    integer,    optional, intent(in)  :: nupper
    m4_type(T), optional, intent(in)  :: d0(:,:)
  end subroutine

  module subroutine band_{}T{}_destruct(this)
    class(band_{}T), intent(inout) :: this
  end subroutine

  module subroutine band_{}T{}_reset(this, only_factorization)
    class(band_{}T),   intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine band_{}T{}_reset_row(this, i)
    class(band_{}T),  intent(inout) :: this
    integer,          intent(in)    :: i
  end subroutine

  module subroutine band_{}T{}_scale(this, fact)
    class(band_{}T), intent(inout) :: this
    m4_type(T),      intent(in)    :: fact
  end subroutine

  module subroutine band_{}T{}_mul_vec(this, x, y, fact_y, trans)
    class(band_{}T),      intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine band_{}T{}_mul_mat(this, x, y, fact_y, trans)
    class(band_{}T),      intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:,:)
    m4_type(T),           intent(inout) :: y(:,:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine band_{}T{}_factorize(this, solver)
    class(band_{}T),   intent(inout) :: this
    integer, optional, intent(in)    :: solver
  end subroutine

  module subroutine band_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
    class(band_{}T),     intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:)
    m4_type(T),          intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine band_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
    class(band_{}T),     intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:,:)
    m4_type(T),          intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine band_{}T{}_output(this, file, funit, fmt)
    class(band_{}T),        intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module subroutine band_{}T{}_add_elem(this, i, j, x)
    class(band_{}T), intent(inout) :: this
    integer,         intent(in)    :: i
    integer,         intent(in)    :: j
    m4_type(T),      intent(in)    :: x
  end subroutine

  module function band_{}T{}_get_elem(this, i, j) result(x)
    class(band_{}T), intent(in) :: this
    integer,         intent(in) :: i
    integer,         intent(in) :: j
    m4_type(T)                  :: x
  end function

  module subroutine band_{}T{}_set_elem(this, i, j, x)
    class(band_{}T), intent(inout) :: this
    integer,         intent(in)    :: i
    integer,         intent(in)    :: j
    m4_type(T),      intent(in)    :: x
  end subroutine

  module subroutine band_{}T{}_set_diag_val(this, val, k)
    class(band_{}T), intent(inout)        :: this
    m4_type(T),      intent(in)           :: val
    integer,         intent(in), optional :: k
  end subroutine

  module subroutine band_{}T{}_set_diag_arr(this, vals, k)
    class(band_{}T), intent(inout)        :: this
    m4_type(T),      intent(in)           :: vals(:)
    integer,         intent(in), optional :: k
  end subroutine

  module function band_eye_{}T{}(nrows) result(eye)
    integer, intent(in) :: nrows
    type(band_{}T)      :: eye
  end function
end interface
