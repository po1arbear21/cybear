type, extends(matrix_{}T) :: block_{}T
  !! Block matrix
  !!
  !! Example
  !!        / A B \
  !!    M = | C D | for A=(2x3), B=(2x4), C=(3x3), D=(3x4), E=(1x3), F=(1x4)-matrix
  !!        \ E F /
  !!
  !!    => nbrows=3, nbcols=2, i0=[1,3,6], i1=[2,5,6], j0=[1,4], j1=[3,7]

  ! dimensions
  integer              :: nbrows
    !! number of block rows
  integer              :: nbcols
    !! number of block cols
  integer, allocatable :: i0(:)
    !! block row start indices
  integer, allocatable :: i1(:)
    !! block row stop  indices
  integer, allocatable :: j0(:)
    !! block col start indices
  integer, allocatable :: j1(:)
    !! block col stop  indices

  ! blocks
  type(matrix_ptr_{}T), allocatable :: b(:,:)
  logical,              allocatable :: alc(:,:)
    !! internal flags to keep track of allocations of b(i,j)%p (delete b(i,j)%p in destruct if alc(i,j) == true)

  ! flags
  logical :: diag
  logical :: tridiag

  ! factorization data
  type(sparse_{}T{}) :: sparse_conversion
  type(dense_{}T{})  :: dense_conversion
contains
  procedure :: block_{}T{}_init
  generic   :: init           => block_{}T{}_init
  procedure :: destruct       => block_{}T{}_destruct
  procedure :: reset          => block_{}T{}_reset
  procedure :: delete_block   => block_{}T{}_delete_block
  procedure :: scale          => block_{}T{}_scale
  procedure :: mul_vec        => block_{}T{}_mul_vec
  procedure :: mul_mat        => block_{}T{}_mul_mat
  procedure :: factorize      => block_{}T{}_factorize
  procedure :: solve_vec      => block_{}T{}_solve_vec
  procedure :: solve_mat      => block_{}T{}_solve_mat
  procedure :: solve_lsqr     => block_{}T{}_solve_lsqr
  procedure :: get_global_idx => block_{}T{}_get_global_idx
  procedure :: get_local_idx  => block_{}T{}_get_local_idx
  procedure :: block_{}T{}_get_dense
  procedure :: block_{}T{}_get_sparse
  procedure :: block_{}T{}_get_band
  procedure :: block_{}T{}_get_block
  generic   :: get            => block_{}T{}_get_dense, &
    &                            block_{}T{}_get_sparse, &
    &                            block_{}T{}_get_band, &
    &                            block_{}T{}_get_block
  procedure :: set            => block_{}T{}_set
  procedure :: set_ptr        => block_{}T{}_set_ptr
  procedure :: copy           => block_{}T{}_copy
  procedure :: copy_deep      => block_{}T{}_copy_deep

  procedure, private :: block_allowed => block_{}T{}_block_allowed
  procedure, private :: check_sizes   => block_{}T{}_check_sizes
end type

type block_ptr_{}T
  type(block_{}T), pointer :: p => null()
end type

interface
  module subroutine block_{}T{}_init(this, row_dim, col_dim, diag, tridiag)
    class(block_{}T),  intent(out) :: this
    integer,           intent(in)  :: row_dim(:)
    integer, optional, intent(in)  :: col_dim(:)
    logical, optional, intent(in)  :: diag
    logical, optional, intent(in)  :: tridiag
  end subroutine

  module subroutine block_{}T{}_destruct(this)
    class(block_{}T), intent(inout) :: this
  end subroutine

  module subroutine block_{}T{}_delete_block(this, ib, jb)
    class(block_{}T), intent(inout) :: this
    integer,          intent(in)    :: ib, jb
  end subroutine

  module subroutine block_{}T{}_reset(this, only_factorization)
    class(block_{}T),  intent(inout) :: this
    logical, optional, intent(in)    :: only_factorization
  end subroutine

  module subroutine block_{}T{}_scale(this, fact)
    class(block_{}T), intent(inout) :: this
    m4_type(T),       intent(in)    :: fact
  end subroutine

  module subroutine block_{}T{}_mul_vec(this, x, y, fact_y, trans)
    class(block_{}T),     intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:)
    m4_type(T),           intent(inout) :: y(:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine block_{}T{}_mul_mat(this, x, y, fact_y, trans)
    class(block_{}T),     intent(in)    :: this
    m4_type(T),           intent(in)    :: x(:,:)
    m4_type(T),           intent(inout) :: y(:,:)
    m4_type(T), optional, intent(in)    :: fact_y
    character,  optional, intent(in)    :: trans
  end subroutine

  module subroutine block_{}T{}_factorize(this, solver)
    class(block_{}T),  intent(inout) :: this
    integer, optional, intent(in)    :: solver
  end subroutine

  module subroutine block_{}T{}_solve_vec(this, rhs, x, trans, iter_refine)
    class(block_{}T),    intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:)
    m4_type(T),          intent(out) :: x(:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine block_{}T{}_solve_mat(this, rhs, x, trans, iter_refine)
    class(block_{}T),    intent(in)  :: this
    m4_type(T),          intent(in)  :: rhs(:,:)
    m4_type(T),          intent(out) :: x(:,:)
    character, optional, intent(in)  :: trans
    logical,   optional, intent(in)  :: iter_refine
  end subroutine

  module subroutine block_{}T{}_solve_lsqr(this, rhs, damp, x, istop, atol, btol)
    class(block_{}T), intent(in)  :: this
    m4_type(T),       intent(in)  :: rhs(:)
    real,             intent(in)  :: damp
    m4_type(T),       intent(out) :: x(:)
    integer,          intent(out) :: istop
    real, optional,   intent(in)  :: atol
    real, optional,   intent(in)  :: btol
  end subroutine

  module subroutine block_{}T{}_get_global_idx(this, ib, jb, il, jl, ig, jg)
    class(block_{}T), intent(in)  :: this
    integer,          intent(in)  :: ib
    integer,          intent(in)  :: jb
    integer,          intent(in)  :: il
    integer,          intent(in)  :: jl
    integer,          intent(out) :: ig
    integer,          intent(out) :: jg
  end subroutine

  module subroutine block_{}T{}_get_local_idx(this, ig, jg, ib, jb, il, jl)
    class(block_{}T), intent(in)  :: this
    integer,          intent(in)  :: ig
    integer,          intent(in)  :: jg
    integer,          intent(out) :: ib
    integer,          intent(out) :: jb
    integer,          intent(out) :: il
    integer,          intent(out) :: jl
  end subroutine

  module subroutine block_{}T{}_get_dense(this, ib, jb, p)
    class(block_{}T),         intent(inout) :: this
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: jb
    type(dense_{}T), pointer, intent(out)   :: p
  end subroutine

  module subroutine block_{}T{}_get_sparse(this, ib, jb, p)
    class(block_{}T),          intent(inout) :: this
    integer,                   intent(in)    :: ib
    integer,                   intent(in)    :: jb
    type(sparse_{}T), pointer, intent(out)   :: p
  end subroutine

  module subroutine block_{}T{}_get_band(this, ib, jb, p)
    class(block_{}T),         intent(inout) :: this
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: jb
    type(band_{}T), pointer,  intent(out)   :: p
  end subroutine

  module subroutine block_{}T{}_get_block(this, ib, jb, p)
    class(block_{}T),         intent(inout) :: this
    integer,                  intent(in)    :: ib
    integer,                  intent(in)    :: jb
    type(block_{}T), pointer, intent(out)   :: p
  end subroutine

  module subroutine block_{}T{}_set(this, ib, jb, b, fact)
    class(block_{}T),     intent(inout) :: this
    integer,              intent(in)    :: ib
    integer,              intent(in)    :: jb
    class(matrix_{}T),    intent(in)    :: b
    m4_type(T), optional, intent(in)    :: fact
  end subroutine

  module subroutine block_{}T{}_set_ptr(this, ib, jb, b)
    class(block_{}T),          intent(inout) :: this
    integer,                   intent(in)    :: ib
    integer,                   intent(in)    :: jb
    class(matrix_{}T), target, intent(in)    :: b
  end subroutine

  module subroutine block_{}T{}_copy(this, copy)
    class(block_{}T), intent(in)  :: this
    type(block_{}T),  intent(out) :: copy
  end subroutine

  module subroutine block_{}T{}_copy_deep(this, copy)
    class(block_{}T), intent(in)  :: this
    type(block_{}T),  intent(out) :: copy
  end subroutine

  module function block_{}T{}_block_allowed(this, ib, jb) result(tf)
    class(block_{}T), intent(in) :: this
    integer,          intent(in) :: ib
    integer,          intent(in) :: jb
    logical                      :: tf
  end function

  module function block_{}T{}_check_sizes(this) result(tf)
    class(block_{}T), intent(in) :: this
    logical                      :: tf
  end function
end interface
