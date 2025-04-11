type, extends(matrix_{}T) :: dense_{}T
  !! Dense matrix type.

  m4_type(T), allocatable :: d(:,:)
    !! Raw matrix data
  m4_type(T), allocatable :: lu(:,:)
    !! LU decomposition data
  integer,    allocatable :: ipiv(:)
    !! Permutation used in conjunction with lu
contains
  procedure :: dense_{}T{}_init_by_size, dense_{}T{}_init_by_data
  generic   :: init       => dense_{}T{}_init_by_size, dense_{}T{}_init_by_data
  procedure :: destruct   => dense_{}T{}_destruct
  procedure :: reset      => dense_{}T{}_reset
  procedure :: scale      => dense_{}T{}_scale
  procedure :: mul_vec    => dense_{}T{}_mul_vec
  procedure :: mul_mat    => dense_{}T{}_mul_mat
  procedure :: factorize  => dense_{}T{}_factorize
  procedure :: solve_vec  => dense_{}T{}_solve_vec
  procedure :: solve_mat  => dense_{}T{}_solve_mat
  procedure :: dense_{}T{}_transpose                        ! transposes this matrix inplace
  procedure :: dense_{}T{}_transpose2                       ! transpose this matrix into another
  generic   :: transpose  => dense_{}T{}_transpose, dense_{}T{}_transpose2
  procedure :: eig        => dense_{}T{}_eig
  procedure :: output     => dense_{}T{}_output
end type

type dense_ptr_{}T
  type(dense_{}T), pointer :: p => null()
end type

interface
  module subroutine dense_{}T{}_init_by_data(this, d)
    class(dense_{}T), intent(out) :: this
    m4_type(T),       intent(in)  :: d(:,:)
  end subroutine

  module subroutine dense_{}T{}_init_by_size(this, nrows, ncols)
    class(dense_{}T),  intent(out) :: this
    integer,           intent(in)  :: nrows
    integer, optional, intent(in)  :: ncols
  end subroutine

  module subroutine dense_{}T{}_destruct(this)
    class(dense_{}T), intent(inout) :: this
  end subroutine

  module subroutine dense_{}T{}_reset(this, only_factorization)
    class(dense_{}T),  intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine dense_{}T{}_scale(this, fact)
    class(dense_{}T), intent(inout) :: this
    m4_type(T),       intent(in)    :: fact
  end subroutine

  module subroutine dense_{}T{}_mul_vec(this, x, y, fact_y, trans)
    class(dense_{}T),     intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine dense_{}T{}_mul_mat(this, x, y, fact_y, trans)
    class(dense_{}T),     intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:,:)
    m4_type(T),           intent(inout) :: y(:,:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine dense_{}T{}_factorize(this, solver)
    class(dense_{}T),  intent(inout) :: this
    integer, optional, intent(in)    :: solver
  end subroutine

  module subroutine dense_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
    class(dense_{}T),    intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:)
    m4_type(T),          intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine dense_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
    class(dense_{}T),    intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:,:)
    m4_type(T),          intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine dense_{}T{}_transpose(this)
    class(dense_{}T), intent(inout) :: this
  end subroutine

  module subroutine dense_{}T{}_transpose2(this, d2)
    class(dense_{}T), intent(in) :: this
    type(dense_{}T),  intent(out) :: d2
  end subroutine

  module subroutine dense_{}T{}_eig(this, e, R, L, sort)
    class(dense_{}T),  intent(in)  :: this
    complex,           intent(out) :: e(:)
    complex, optional, intent(out) :: R(:,:)
    complex, optional, intent(out) :: L(:,:)
    logical, optional, intent(in)  :: sort
  end subroutine

  module subroutine dense_{}T{}_output(this, file, funit, fmt)
    class(dense_{}T),       intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module function dense_eye_{}T{}(nrows, ncols) result(eye)
    integer,           intent(in)  :: nrows
    integer, optional, intent(in)  :: ncols
    type(dense_{}T)                :: eye
  end function
end interface
