type, extends(matrix_{}T) :: sparse_{}T
  !! Sparse matrix type in CSR format.

  integer(SPARSE_IDX), allocatable :: ia(:)
  integer,             allocatable :: ja(:)
  m4_type(T),                  allocatable :: a(:)

  integer :: solver = 0
    !! SOLVER_PARDISO, SOLVER_MUMPS or SOLVER_ILUPACK
  integer :: solver_handle = 0
    !! pardiso, mumps or ilupack solver handle
contains
  procedure :: sparse_{}T{}_init
  generic   :: init          => sparse_{}T{}_init
  procedure :: init_solver   => sparse_{}T{}_init_solver
  procedure :: destruct      => sparse_{}T{}_destruct
  procedure :: reset         => sparse_{}T{}_reset
  procedure :: scale         => sparse_{}T{}_scale
  procedure :: mul_vec       => sparse_{}T{}_mul_vec
  procedure :: mul_vec_slice => sparse_{}T{}_mul_vec_slice
  procedure :: mul_mat       => sparse_{}T{}_mul_mat
  procedure :: factorize     => sparse_{}T{}_factorize
  procedure :: solve_vec     => sparse_{}T{}_solve_vec
  procedure :: solve_mat     => sparse_{}T{}_solve_mat
  procedure :: sparse_{}T{}_transpose
  procedure :: sparse_{}T{}_transpose2
  generic   :: transpose     => sparse_{}T{}_transpose, &
    &                           sparse_{}T{}_transpose2
  procedure :: mul_sparse    => sparse_{}T{}_mul_sparse
  procedure :: input         => sparse_{}T{}_input
  procedure :: output        => sparse_{}T{}_output
  procedure :: nnz           => sparse_{}T{}_nnz
  procedure :: is_empty      => sparse_{}T{}_is_empty
end type

type sparse_ptr_{}T
  type(sparse_{}T), pointer :: p => null()
end type

type spbuild_{}T
  !! Sparse matrix builder which is more flexible than simple CSR.
  type(sparse_{}T), pointer     :: sp
    !! Non-owning pointer to sparse matrix this builder is tied to
  type(vector_int), allocatable :: cols(:)
    !! Column indices of non-zero elements in each row
  type(vector_{}T), allocatable :: vals(:)
    !! Values of non-zero elements in each row
  type(vector_log), allocatable :: filled(:)
    !! Indicates whether a value is actually used (keep structure after reset)
  logical,          allocatable :: sorted(:)
    !! Indicates whether a row is sorted or not
  logical,          allocatable :: keep_struct(:)
    !! Indicates whether a row should be deleted completely on reset, or structure should be kept
  logical                       :: assume_unique
    !! user ensures no duplicate entries
contains
  procedure :: init      => spbuild_{}T{}_init
  procedure :: destruct  => spbuild_{}T{}_destruct
  procedure :: reset     => spbuild_{}T{}_reset
  procedure :: reset_row => spbuild_{}T{}_reset_row
  procedure :: set       => spbuild_{}T{}_set
  procedure :: set_row   => spbuild_{}T{}_set_row
  procedure :: add       => spbuild_{}T{}_add
  procedure :: load      => spbuild_{}T{}_load
  procedure :: save      => spbuild_{}T{}_save

  procedure, private :: sort => spbuild_{}T{}_sort
end type

interface
  module subroutine sparse_{}T{}_init(this, nrows, ncols)
    class(sparse_{}T), intent(out) :: this
    integer,           intent(in)  :: nrows
    integer, optional, intent(in)  :: ncols
  end subroutine

  module subroutine sparse_{}T{}_destruct(this)
    class(sparse_{}T), intent(inout) :: this
  end subroutine

  module subroutine sparse_{}T{}_reset(this, only_factorization)
    class(sparse_{}T), intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine sparse_{}T{}_scale(this, fact)
    class(sparse_{}T), intent(inout) :: this
    m4_type(T),        intent(in)    :: fact
  end subroutine

  module subroutine sparse_{}T{}_mul_vec_slice(this, x, y, fact_y, i0, i1)
    class(sparse_{}T),    intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    integer,    optional, intent(in)    :: i0
    integer,    optional, intent(in)    :: i1
  end subroutine

  module subroutine sparse_{}T{}_mul_vec(this, x, y, fact_y, trans)
    class(sparse_{}T),    intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine sparse_{}T{}_mul_mat(this, x, y, fact_y, trans)
    class(sparse_{}T),    intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:,:)
    m4_type(T),           intent(inout) :: y(:,:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine sparse_{}T{}_init_solver(this)
    class(sparse_{}T), intent(inout) :: this
  end subroutine

  module subroutine sparse_{}T{}_factorize(this)
    class(sparse_{}T), intent(inout) :: this
  end subroutine

  module subroutine sparse_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
    class(sparse_{}T),   intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:)
    m4_type(T),          intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine sparse_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
    class(sparse_{}T),   intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:,:)
    m4_type(T),          intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine sparse_{}T{}_transpose(this)
    class(sparse_{}T), intent(inout) :: this
  end subroutine

  module subroutine sparse_{}T{}_transpose2(this, s2)
    class(sparse_{}T), intent(in)  :: this
    type(sparse_{}T),  intent(out) :: s2
  end subroutine

  module subroutine sparse_{}T{}_mul_sparse(this, s2, s3, sort)
    class(sparse_{}T),  intent(in)  :: this
    type(sparse_{}T),   intent(in)  :: s2
    type(sparse_{}T),   intent(out) :: s3
    logical, optional,  intent(in)  :: sort
  end subroutine

  module subroutine sparse_{}T{}_input(this, file, funit)
    class(sparse_{}T),      intent(out) :: this
    character(*), optional, intent(in)  :: file
    integer,      optional, intent(in)  :: funit
  end subroutine

  module subroutine sparse_{}T{}_output(this, file, funit, fmt)
    class(sparse_{}T),      intent(in) :: this
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: funit
    character(*), optional, intent(in) :: fmt
  end subroutine

  module function sparse_{}T{}_nnz(this, only_nonzeros) result(nnz)
    class(sparse_{}T), intent(in) :: this
    logical, optional, intent(in) :: only_nonzeros
    integer(SPARSE_IDX)           :: nnz
  end function

  module function sparse_{}T{}_is_empty(this) result(empty)
    class(sparse_{}T), intent(in) :: this
    logical                       :: empty
  end function

  module function sparse_eye_{}T{}(nrows, ncols) result(eye)
    integer,           intent(in)  :: nrows
    integer, optional, intent(in)  :: ncols
    type(sparse_{}T)               :: eye
  end function

  module function sparse_zero_{}T{}(nrows, ncols) result(zero)
    integer,           intent(in) :: nrows
    integer, optional, intent(in) :: ncols
    type(sparse_{}T)             :: zero
  end function

  module subroutine spbuild_{}T{}_init(this, s, n_per_row, assume_unique)
    class(spbuild_{}T),       intent(out) :: this
    type(sparse_{}T), target, intent(in)  :: s
    integer, optional,        intent(in)  :: n_per_row
    logical, optional,        intent(in)  :: assume_unique
  end subroutine

  module subroutine spbuild_{}T{}_destruct(this)
    class(spbuild_{}T), intent(inout) :: this
  end subroutine

  module subroutine spbuild_{}T{}_reset(this)
    class(spbuild_{}T), intent(inout) :: this
  end subroutine

  module subroutine spbuild_{}T{}_reset_row(this, row)
    class(spbuild_{}T), intent(inout) :: this
    integer,            intent(in)    :: row
  end subroutine

  module subroutine spbuild_{}T{}_set(this, row, col, val, fact, search)
    class(spbuild_{}T),   intent(inout) :: this
    integer,              intent(in)    :: row
    integer,              intent(in)    :: col
    m4_type(T),           intent(in)    :: val
    m4_type(T), optional, intent(in)    :: fact
    logical,    optional, intent(in)    :: search
  end subroutine

  module subroutine spbuild_{}T{}_set_row(this, row, vals, j0, j1)
    class(spbuild_{}T), intent(inout) :: this
    integer,            intent(in)    :: row
    m4_type(T),         intent(in)    :: vals(:)
    integer, optional,  intent(in)    :: j0
    integer, optional,  intent(in)    :: j1
  end subroutine

  module subroutine spbuild_{}T{}_add(this, row, col, val)
    class(spbuild_{}T), intent(inout) :: this
    integer,            intent(in)    :: row
    integer,            intent(in)    :: col
    m4_type(T),         intent(in)    :: val
  end subroutine

  module subroutine spbuild_{}T{}_sort(this)
    class(spbuild_{}T), intent(inout) :: this
  end subroutine

  module subroutine spbuild_{}T{}_load(this)
    class(spbuild_{}T), intent(inout) :: this
  end subroutine

  module subroutine spbuild_{}T{}_save(this)
    class(spbuild_{}T), intent(inout) :: this
  end subroutine
end interface
