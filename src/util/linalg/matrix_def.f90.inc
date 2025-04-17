type, abstract :: matrix_{}T
  !! Base type for matrices.
  !! Can be extended if more specialized matrices are needed besides the built-in types.
  integer :: nrows
    !! Number of rows
  integer :: ncols
    !! Number of columns
  logical :: factorized
    !! Is matrix factorized and ready to solve?
  integer, private :: solver = 0
    !! Solver used when computing the stored factorization
contains
  procedure :: matrix_{}T{}_init
  generic   :: init => matrix_{}T{}_init

  procedure(matrix_{}T{}_destruct),  deferred :: destruct
  procedure(matrix_{}T{}_reset),     deferred :: reset
  procedure(matrix_{}T{}_scale),     deferred :: scale
  procedure(matrix_{}T{}_mul_vec),   deferred :: mul_vec
  procedure(matrix_{}T{}_mul_mat),   deferred :: mul_mat
  procedure(matrix_{}T{}_factorize), deferred :: factorize
  procedure(matrix_{}T{}_solve_vec), deferred :: solve_vec
  procedure(matrix_{}T{}_solve_mat), deferred :: solve_mat
end type

type matrix_ptr_{}T
  class(matrix_{}T), pointer :: p => null()
end type

interface
  module subroutine matrix_{}T{}_destruct(this)
    class(matrix_{}T), intent(inout) :: this
  end subroutine

  module subroutine matrix_{}T{}_reset(this, only_factorization)
    class(matrix_{}T), intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine matrix_{}T{}_scale(this, fact)
    class(matrix_{}T), intent(inout) :: this
    m4_type(T),        intent(in)    :: fact
  end subroutine

  module subroutine matrix_{}T{}_mul_vec(this, x, y, fact_y, trans)
    class(matrix_{}T),    intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine matrix_{}T{}_mul_mat(this, x, y, fact_y, trans)
    class(matrix_{}T),    intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:,:)
    m4_type(T),           intent(inout) :: y(:,:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine matrix_{}T{}_factorize(this, solver)
    class(matrix_{}T), intent(inout) :: this
    integer, optional, intent(in)    :: solver
  end subroutine

  module subroutine matrix_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
    class(matrix_{}T),   intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:)
    m4_type(T),          intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine matrix_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
    class(matrix_{}T),   intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:,:)
    m4_type(T),          intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine
end interface
