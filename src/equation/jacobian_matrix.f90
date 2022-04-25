m4_include(../util/macro.f90.inc)

module jacobian_matrix_m
  use matrix_m,    only: matrix_add, sparse_real, block_real, spbuild_real, sparse_ptr_real
  use vselector_m, only: vselector

  implicit none

  private
  public jacobian_matrix, jacobian_matrix_ptr, vector_jacobian_matrix_ptr

  type, extends(block_real) :: jacobian_matrix
    !! jacobian in block matrix representation

    type(vselector), pointer :: v1 => null()
      !! result variable selector
    type(vselector), pointer :: v2 => null()
      !! dependency variable selector

    type(sparse_ptr_real), allocatable :: s(:,:)
      !! size: (v1%ntab, v2%ntab)
      !! size of matrix s(itab1,itab2): (v1%nvals(itab1), v2%nvals(itab2))
    type(spbuild_real),    allocatable :: sb(:,:)
      !! size: (v1%ntab, v2%ntab)
      !! size of matrix sb(itab1,itab2): (v1%nvals(itab1), v2%nvals(itab2))

    logical, allocatable :: const(:,:)
      !! const flags for each block (true means linear equation)
    logical, allocatable :: zero(:,:)
      !! zero flags for each block
    logical, allocatable :: valmsk(:,:,:,:)
      !! value masks for each block, false means entry is always zero, for all grid points (v1%nval x v2%nval x size(v1%tab) x size(v2%tab))
  contains
    procedure :: jacobian_matrix_init
    generic   :: init     => jacobian_matrix_init
    procedure :: reset    => jacobian_matrix_reset
    procedure :: add_jaco => jacobian_matrix_add_jaco
    procedure :: mul_jaco => jacobian_matrix_mul_jaco
  end type

  type jacobian_matrix_ptr
    !! pointer to jacobian matrix
    type(jacobian_matrix), pointer :: p => null()
  end type

  m4_define({T},{jacobian_matrix_ptr})
  m4_include(../util/vector_def.f90.inc)

contains

  m4_define({T},{jacobian_matrix_ptr})
  m4_include(../util/vector_imp.f90.inc)

  subroutine jacobian_matrix_init(this, v1, v2, const, zero, valmsk)
    !! initialize jacobian matrix
    class(jacobian_matrix),  intent(out) :: this
    type(vselector), target, intent(in)  :: v1
      !! result variable selector
    type(vselector), target, intent(in)  :: v2
      !! dependency variable selector
    logical,                 intent(in)  :: const(:,:)
      !! const flags (v1%ntab x v2%ntab)
    logical,                 intent(in)  :: zero(:,:)
      !! zero flags (v1%ntab x v2%ntab)
    logical,                 intent(in)  :: valmsk(:,:,:,:)
      !! value mask (v1%nval x v2%nval x v1%ntab x v2%ntab)

    integer :: itab1, itab2

    ! init base
    call this%block_real_init(v1%nvals, col_dim = v2%nvals)

    ! save var selectors
    this%v1 => v1
    this%v2 => v2

    ! allocate flags
    this%const  = const
    this%zero   = zero
    this%valmsk = valmsk

    ! init blocks as empty sparse matrices
    allocate (this%s( v1%ntab,v2%ntab))
    allocate (this%sb(v1%ntab,v2%ntab))
    do itab1 = 1, v1%ntab; do itab2 = 1, v2%ntab
      call this%get(itab1, itab2, this%s(itab1,itab2)%p)
      call this%s( itab1,itab2)%p%init(v1%nvals(itab1), v2%nvals(itab2))
      call this%sb(itab1,itab2)%init(this%s(itab1,itab2)%p, assume_unique=.true.)
    end do; end do
  end subroutine

  subroutine jacobian_matrix_reset(this, only_factorization)
    !! reset jacobian matrix
    class(jacobian_matrix), intent(inout) :: this
    logical, optional,      intent(in)    :: only_factorization

    integer :: itab1, itab2

    ! reset base
    call this%block_real%reset(only_factorization=only_factorization)

    ! reset spbuild
    if (allocated(this%sb)) then
      do itab1 = 1, this%v1%ntab; do itab2 = 1, this%v2%ntab
        call this%sb(itab1,itab2)%reset()
      end do; end do
    end if
  end subroutine

  subroutine jacobian_matrix_add_jaco(this, jaco, reset, itab1, itab2)
    !! add jacobian matrices to this matrix
    class(jacobian_matrix),    intent(inout) :: this
      !! result
    type(jacobian_matrix_ptr), intent(in)    :: jaco(:)
      !! matrices to add
    logical, optional,         intent(in)    :: reset
      !! reset result before adding (default: true)
    integer, optional,         intent(in)    :: itab1
      !! 1st table/block index
    integer, optional,         intent(in)    :: itab2
      !! 2nd table/block index

    integer :: i, j, k, i0, i1, j0, j1
    logical :: reset_

    ! process optional arguments
    reset_ = .true.
    if (present(reset)) reset_ = reset
    if (present(itab1)) then
      i0 = itab1
      i1 = itab1
    else
      i0 = 1
      i1 = this%v1%ntab
    end if
    if (present(itab2)) then
      j0 = itab2
      j1 = itab2
    else
      j0 = 1
      j1 = this%v2%ntab
    end if

    ! add matrices
    do i = i0, i1; do j = j0, j1
      ! reset matrix
      if (reset_) call this%s(i,j)%p%reset()

      ! add matrices
      do k = 1, size(jaco)
        call matrix_add(jaco(k)%p%s(i,j)%p, this%s(i,j)%p)
      end do
    end do; end do
  end subroutine

  subroutine jacobian_matrix_mul_jaco(this, jaco1, jaco2, reset, itab1, itab2)
    !! multiply two jacobians and add result to this
    class(jacobian_matrix), target, intent(inout) :: this
      !! result matrix
    class(jacobian_matrix), target, intent(in)    :: jaco1
      !! first jacobian
    class(jacobian_matrix), target, intent(in)    :: jaco2
      !! second jacobian
    logical, optional,              intent(in)    :: reset
      !! reset result first (default: true)
    integer, optional,              intent(in)    :: itab1
      !! 1st table/block index
    integer, optional,              intent(in)    :: itab2
      !! 2nd table/block index

    integer           :: i, j, k, i0, i1, j0, j1
    logical           :: reset_
    type(sparse_real) :: tmp

    ! process optional arguments
    reset_ = .true.
    if (present(reset)) reset_ = reset
    if (present(itab1)) then
      i0 = itab1
      i1 = itab1
    else
      i0 = 1
      i1 = this%v1%ntab
    end if
    if (present(itab2)) then
      j0 = itab2
      j1 = itab2
    else
      j0 = 1
      j1 = this%v2%ntab
    end if

    ! multiply matrices
    do i = i0, i1; do j = j0, j1
      ! reset result
      if (reset_) call this%s(i,j)%p%reset()

      ! loop over jaco1%ncols = jaco2%nrows
      do k = 1, jaco1%v2%ntab
        call jaco1%s(i,k)%p%mul_sparse(jaco2%s(k,j)%p, tmp)
        call matrix_add(tmp, this%s(i,j)%p)
      end do
    end do; end do
  end subroutine
end module
