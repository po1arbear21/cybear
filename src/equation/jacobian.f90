#include "../util/macro.f90.inc"

module jacobian_m

  use array_m,           only: array4_real
  use error_m
  use jacobian_matrix_m, only: jacobian_matrix
  use stencil_m,         only: stencil_ptr, static_stencil, dynamic_stencil
  use vselector_m,       only: vselector

  implicit none

  private
  public jacobian, jacobian_ptr

  type jacobian
    !! jacobian structure to store derivatives

    type(jacobian_matrix) :: matr
      !! jacobian in blockmatrix representation

    type(stencil_ptr), allocatable :: st(:)
      !! stencils (v1%ntab)
    type(array4_real), allocatable :: d(:)
      !! derivatives (v1%nval x v2%nval x st(itab1)%nmax x v1%tab(itab1)%n) x (v1%ntab)
  contains
    procedure :: init     => jacobian_init
    procedure :: destruct => jacobian_destruct
    procedure :: reset    => jacobian_reset
    procedure :: set_matr => jacobian_set_matr
    generic   :: set      => jacobian_set_multi_stat,  &
      &                      jacobian_set_multi_dyn,   &
      &                      jacobian_set_single_stat, &
      &                      jacobian_set_single_dyn,  &
      &                      jacobian_set_scalar_stat, &
      &                      jacobian_set_scalar_dyn
    generic   :: update   => jacobian_update_multi_stat,  &
      &                      jacobian_update_multi_dyn,   &
      &                      jacobian_update_single_stat, &
      &                      jacobian_update_single_dyn,  &
      &                      jacobian_update_scalar_stat, &
      &                      jacobian_update_scalar_dyn

    procedure, private :: jacobian_set_multi_stat
    procedure, private :: jacobian_set_multi_dyn
    procedure, private :: jacobian_set_single_stat
    procedure, private :: jacobian_set_single_dyn
    procedure, private :: jacobian_set_scalar_stat
    procedure, private :: jacobian_set_scalar_dyn

    procedure, private :: jacobian_update_multi_stat
    procedure, private :: jacobian_update_multi_dyn
    procedure, private :: jacobian_update_single_stat
    procedure, private :: jacobian_update_single_dyn
    procedure, private :: jacobian_update_scalar_stat
    procedure, private :: jacobian_update_scalar_dyn

    procedure, private :: jacobian_edit_multi_stat
    procedure, private :: jacobian_edit_multi_dyn
    procedure, private :: jacobian_edit_single_stat
    procedure, private :: jacobian_edit_single_dyn
  end type

  type jacobian_ptr
    type(jacobian), pointer :: p => null()
  end type

contains

  subroutine jacobian_init(this, v1, v2, st, const, zero, valmsk)
    !! initialize jacobian
    class(jacobian), target, intent(out) :: this
    type(vselector), target, intent(in)  :: v1
      !! result variable selector
    type(vselector), target, intent(in)  :: v2
      !! dependency variable selector
    type(stencil_ptr),       intent(in)  :: st(:)
      !! stencils (v1%ntab)
    logical, optional,       intent(in)  :: const
      !! const flag; default = false; the same for all blocks
    logical, optional,       intent(in)  :: zero(:,:)
      !! zero flags (v1%ntab x v2%ntab); default: set automatically by checking stencils
    logical, optional,       intent(in)  :: valmsk(:,:)
      !! value mask (v1%nval x v2%nval); default = true; the same for all blocks

    integer :: i, j, itab1, itab2
    logical :: const_(v1%ntab,v2%ntab), zero_(v1%ntab,v2%ntab), valmsk_(v1%nval,v2%nval,v1%ntab,v2%ntab), status

    ASSERT(size(st) == v1%ntab)

    ! optional arguments
    const_ = .false.
    if (present(const)) const_ = const
    if (present(zero)) then
      zero_ = zero
    else
      ! set zero flags automatically by checking stencils
      zero_ = .true.
      do itab1 = 1, v1%ntab
        select type (st_ptr => st(itab1)%p)
          class is (static_stencil)
            ! do nothing for empty stencil
            if (.not. associated(st_ptr)) cycle

            block
              integer :: idx1(v1%g%idx_dim), idx2(v2%g%idx_dim)

              ! loop over result points
              do i = 1, v1%tab(itab1)%p%n
                ! get result grid indices
                idx1 = v1%tab(itab1)%p%get_idx(i)

                ! loop over dependency grid indices
                do j = 1, st_ptr%nmax
                  ! get dependency grid indices
                  call st_ptr%get(idx1, j, idx2, status)
                  if (.not. status) exit

                  ! stencil connects itab1 to itab2 => not zero
                  itab2 = v2%itab%get(idx2)
                  zero_(itab1,itab2) = .false.
                end do
              end do
            end block

          class is (dynamic_stencil)
            zero_ = .false.

          class default
            call program_error("stencil type not supported; derive from static or dynamic stencil!")
        end select
      end do
    end if
    const_  = (const_ .or. zero_) ! zero blocks are also constant
    valmsk_ = .true.
    if (present(valmsk)) then
      do itab2 = 1, v2%ntab; do itab1 = 1, v1%ntab
        valmsk_(:,:,itab1,itab2) = valmsk
      end do; end do
    end if

    ! init matrix
    call this%matr%init(v1, v2, const_, zero_, valmsk_)

    ! save stencil pointers
    this%st = st

    ! allocate memory for derivatives of static stencil data
    allocate (this%d(v1%ntab))
    do itab1 = 1, v1%ntab
      select type (st_ptr => st(itab1)%p)
        class is (static_stencil)
          allocate (this%d(itab1)%d(v1%nval, v2%nval, st_ptr%nmax, v1%tab(itab1)%p%n), source = 0.0)
        class is (dynamic_stencil)
          ASSERT(.not. all(const_))
      end select
    end do
  end subroutine

  subroutine jacobian_destruct(this)
    !! destruct jacobian
    class(jacobian), intent(inout) :: this

    call this%matr%destruct()
    if (allocated(this%st)) deallocate (this%st)
    if (allocated(this%d )) deallocate (this%d )
  end subroutine

  subroutine jacobian_reset(this, const, nonconst)
    !! reset jacobian to zeros
    class(jacobian),   intent(inout) :: this
    logical, optional, intent(in)    :: const
      !! reset constant blocks (default: true)
    logical, optional, intent(in)    :: nonconst
      !! reset non-constant blocks (default: true)

    integer :: i, j, itab1, itab2
    logical :: const_, nonconst_, reset(this%matr%v1%ntab,this%matr%v2%ntab)

    ! optional arguments
    const_ = .true.
    if (present(const)) const_ = const
    nonconst_ = .true.
    if (present(nonconst)) nonconst_ = nonconst

    ! decide which blocks to reset
    reset = (this%matr%const .and. const_) .or. (.not. this%matr%const .and. nonconst_)

    associate (v1 => this%matr%v1, v2 => this%matr%v2)
      ! loop over result tables
      do itab1 = 1, v1%ntab
        if (.not. associated(this%st(itab1)%p)) cycle
        if (all(.not. reset(itab1,:))) cycle

        select type (st => this%st(itab1)%p)
          class is (static_stencil)
            if (all(reset(itab1,:))) then
              this%d(itab1)%d = 0
              cycle
            end if

            block
              integer :: idx1(v1%g%idx_dim), idx2(v2%g%idx_dim)
              logical :: status

              ! loop over result points
              do i = 1, v1%tab(itab1)%p%n
                ! get result grid indices
                idx1 = v1%tab(itab1)%p%get_idx(i)

                ! loop over dependency grid indices
                do j = 1, st%nmax
                  ! get dependency grid indices
                  call st%get(idx1, j, idx2, status)
                  if (.not. status) exit

                  ! get dependency table index
                  itab2 = v2%itab%get(idx2)

                  ! reset only if flag is set
                  if (reset(itab1,itab2)) this%d(itab1)%d(:,:,j,i) = 0
                end do
              end do
            end block

          class is (dynamic_stencil)
            do itab2 = 1, v2%ntab
              if (reset(itab1,itab2)) then
                call this%matr%s( itab1,itab2)%p%reset()
                call this%matr%sb(itab1,itab2)%reset()
              end if
            end do
        end select
      end do
    end associate
  end subroutine

  subroutine jacobian_set_matr(this, const, nonconst)
    !! generate sparse matrix representation
    class(jacobian),   intent(inout) :: this
    logical, optional, intent(in)    :: const
      !! enable processing of const blocks (default: true)
    logical, optional, intent(in)    :: nonconst
      !! enable processing of non-const blocks (default: true)

    integer :: i, j, itab1, itab2, ival1, ival2, row, row0, col, col0
    logical :: const_, nonconst_, set(this%matr%v1%ntab,this%matr%v2%ntab)

    ! optional arguments
    const_ = .true.
    if (present(const)) const_ = const
    nonconst_ = .true.
    if (present(nonconst)) nonconst_ = nonconst

    ! decide which blocks to set
    set = (this%matr%const .and. const_) .or. (.not. this%matr%const .and. nonconst_)

    associate (v1 => this%matr%v1, v2 => this%matr%v2)
      ! loop over all blocks
      do itab1 = 1, v1%ntab
        ! decide whether to process any blocks of this block row
        if (.not. any(set(itab1,:))) cycle

        ! insert data into spbuilders for static stencils
        select type (st => this%st(itab1)%p)
          class is (static_stencil)
            block
              integer :: idx1(v1%g%idx_dim), idx2(v2%g%idx_dim)
              logical :: status

              ! reset blocks before saving values
              do itab2 = 1, v2%ntab
                if (set(itab1,itab2)) then
                  call this%matr%s( itab1,itab2)%p%reset()
                  call this%matr%sb(itab1,itab2)%reset()
                end if
              end do

              if (.not. all(this%matr%zero(itab1,:))) then
                ! loop over points in result table
                do i = 1, v1%tab(itab1)%p%n
                  ! get result grid indices
                  idx1 = v1%tab(itab1)%p%get_idx(i)

                  ! get row base index
                  row0 = (i - 1) * v1%nval

                  ! loop over dependency points
                  do j = 1, st%nmax
                    ! get dependency points
                    call st%get(idx1, j, idx2, status)
                    if (.not. status) exit

                    ! get dependency table index
                    itab2 = v2%itab%get(idx2)

                    ! do nothing if dependency matrix is not used (const or zero flag prevents it)
                    if ((.not. set(itab1,itab2)) .or. (this%matr%zero(itab1,itab2))) cycle

                    ! get column base index
                    col0 = (v2%tab(itab2)%p%get_flat(idx2) - 1) * v2%nval

                    ! set derivatives
                    do ival1 = 1, v1%nval
                      row = row0 + ival1
                      do ival2 = 1, v2%nval
                        if (.not. this%matr%valmsk(ival1,ival2,itab1,itab2)) cycle
                        col = col0 + ival2
                        call this%matr%sb(itab1,itab2)%set(row, col, this%d(itab1)%d(ival1,ival2,j,i), search = .false.)
                      end do
                    end do
                  end do
                end do
              end if
            end block
        end select

        ! save data in sparse_matrix
        do itab2 = 1, v2%ntab
          if (set(itab1,itab2)) call this%matr%sb(itab1,itab2)%save()
        end do
      end do
    end associate
  end subroutine

  !
  ! stationary stencil procedures
  !
  subroutine jacobian_set_multi_stat(this, itab1, idx1, j, d)
    !! set one element of multidimensional derivative for static stencil.
    !! see arguments' documentation in jacobian_edit_multi_stat.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: j
    real,            intent(in)    :: d(:,:)

    call this%jacobian_edit_multi_stat(itab1, idx1, j, d, set=.true.)
  end subroutine

  subroutine jacobian_update_multi_stat(this, itab1, idx1, j, d)
    !! update one element of multidimensional derivative for static stencil.
    !! see arguments' documentation in jacobian_edit_multi_stat.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: j
    real,            intent(in)    :: d(:,:)

    call this%jacobian_edit_multi_stat(itab1, idx1, j, d, update=.true.)
  end subroutine

  subroutine jacobian_edit_multi_stat(this, itab1, idx1, j, d, set, update)
    !! set/update full (multidimensional) derivative for static stencil
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! grid table index
    integer,           intent(in)    :: idx1(:)
      !! result index
    integer,           intent(in)    :: j
      !! stencil dependency index
    real,            intent(in)      :: d(:,:)
      !! derivatives. size: (v1%nval, v2%nval)
    logical, optional, intent(in)    :: set, update
      !! flags to either set or update the derivative. (default: false)
      !! exactly one must be true.

    integer :: i
    logical :: set_, update_

    ASSERT(size(d, dim=1) == this%matr%v1%nval)
    ASSERT(size(d, dim=2) == this%matr%v2%nval)

    ! optional flags
    set_ = .false.
    if (present(set)) set_ = set
    update_ = .false.
    if (present(update)) update_ = update
    ASSERT(set_ .neqv. update_)

    ! result flat index
    i = this%matr%v1%tab(itab1)%p%get_flat(idx1)

    ! edit derivative
    if      (set_   ) then
      this%d(itab1)%d(:,:,j,i) = d
    else if (update_) then
      this%d(itab1)%d(:,:,j,i) = d + this%d(itab1)%d(:,:,j,i)
    end if
  end subroutine


  subroutine jacobian_set_single_stat(this, itab1, idx1, j, ival1, ival2, d)
    !! set one element of multidimensional derivative for static stencil.
    !! see arguments' documentation in jacobian_edit_single_stat.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: j
    integer,         intent(in)    :: ival1
    integer,         intent(in)    :: ival2
    real,            intent(in)    :: d

    call this%jacobian_edit_single_stat(itab1, idx1, j, ival1, ival2, d, set=.true.)
  end subroutine

  subroutine jacobian_update_single_stat(this, itab1, idx1, j, ival1, ival2, d)
    !! update one element of multidimensional derivative for static stencil.
    !! see arguments' documentation in jacobian_edit_single_stat.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: j
    integer,         intent(in)    :: ival1
    integer,         intent(in)    :: ival2
    real,            intent(in)    :: d

    call this%jacobian_edit_single_stat(itab1, idx1, j, ival1, ival2, d, update=.true.)
  end subroutine

  subroutine jacobian_edit_single_stat(this, itab1, idx1, j, ival1, ival2, d, set, update)
    !! set/update one element of multidimensional derivative for static stencil
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! grid table index
    integer,           intent(in)    :: idx1(:)
      !! result index
    integer,           intent(in)    :: j
      !! stencil dependency index
    integer,           intent(in)    :: ival1
      !! result dimension
    integer,           intent(in)    :: ival2
      !! dependency dimension
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: set, update
    !! flags to either set or update the derivative. (default: false)
    !! exactly one must be true.

    integer :: i
    logical :: set_, update_

    ! optional flags
    set_ = .false.
    if (present(set)) set_ = set
    update_ = .false.
    if (present(update)) update_ = update
    ASSERT(set_ .neqv. update_)

    ! result flat index
    i = this%matr%v1%tab(itab1)%p%get_flat(idx1)

    ! edit derivative
    if      (set_   ) then
      this%d(itab1)%d(ival1,ival2,j,i) = d
    else if (update_) then
      this%d(itab1)%d(ival1,ival2,j,i) = d + this%d(itab1)%d(ival1,ival2,j,i)
    end if
  end subroutine


  subroutine jacobian_set_scalar_stat(this, itab1, idx1, j, d)
    !! set scalar derivative for static stencil of one-dimensional variables.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
      !! grid table index
    integer,         intent(in)    :: idx1(:)
      !! result index
    integer,         intent(in)    :: j
      !! stencil dependency index
    real,            intent(in)    :: d
      !! derivative

    ASSERT(this%matr%v1%nval == 1)
    ASSERT(this%matr%v2%nval == 1)

    call this%set(itab1, idx1, j, 1, 1, d)
  end subroutine

  subroutine jacobian_update_scalar_stat(this, itab1, idx1, j, d)
    !! update scalar derivative for static stencil of one-dimensional variables.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
      !! grid table index
    integer,         intent(in)    :: idx1(:)
      !! result index
    integer,         intent(in)    :: j
      !! stencil dependency index
    real,            intent(in)    :: d
      !! derivative

    ASSERT(this%matr%v1%nval == 1)
    ASSERT(this%matr%v2%nval == 1)

    call this%update(itab1, idx1, j, 1, 1, d)
  end subroutine


  !
  ! dynamical stencil procedures
  !
  subroutine jacobian_set_multi_dyn(this, itab1, idx1, idx2, d)
    !! set full (multidimensional) derivative for dynamic stencil
    !! see arguments' documentation in jacobian_edit_multi_dyn.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: idx2(:)
    real,            intent(in)    :: d(:,:)

    call this%jacobian_edit_multi_dyn(itab1, idx1, idx2, d, set=.true.)
  end subroutine

  subroutine jacobian_update_multi_dyn(this, itab1, idx1, idx2, d)
    !! update full (multidimensional) derivative for dynamic stencil
    !! see arguments' documentation in jacobian_edit_multi_dyn.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: idx2(:)
    real,            intent(in)    :: d(:,:)

    call this%jacobian_edit_multi_dyn(itab1, idx1, idx2, d, update=.true.)
  end subroutine

  subroutine jacobian_edit_multi_dyn(this, itab1, idx1, idx2, d, set, update)
    !! set/update full (multidimensional) derivative for dynamic stencil
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! grid table index
    integer,           intent(in)    :: idx1(:)
      !! result index
    integer,           intent(in)    :: idx2(:)
      !! dependency index
    real,              intent(in)    :: d(:,:)
      !! derivatives. size: (v1%nval, v2%nval)
    logical, optional, intent(in)    :: set, update
      !! flags to either set or update the derivative. (default: false)
      !! exactly one must be true.

    integer :: itab2, ival1, ival2, col, col0, row, row0
    logical :: set_, update_

    ASSERT(size(d, dim=1) == this%matr%v1%nval)
    ASSERT(size(d, dim=2) == this%matr%v2%nval)

    ! optional flags
    set_ = .false.
    if (present(set)) set_ = set
    update_ = .false.
    if (present(update)) update_ = update
    ASSERT(set_ .neqv. update_)

    ! get row base index
    row0 = (this%matr%v1%tab(itab1)%p%get_flat(idx1) - 1) * this%matr%v1%nval

    ! get column base index
    itab2 = this%matr%v2%itab%get(idx2)
    col0  = (this%matr%v2%tab(itab2)%p%get_flat(idx2) - 1) * this%matr%v2%nval

    ! edit derivatives
    do ival1 = 1, this%matr%v1%nval
      row = row0 + ival1
      do ival2 = 1, this%matr%v2%nval
        col = col0 + ival2
        if      (set_   ) then
          call this%matr%sb(itab1,itab2)%set(row, col, d(ival1,ival2), search = .false.)
        else if (update_) then
          call this%matr%sb(itab1,itab2)%add(row, col, d(ival1,ival2))
        end if
      end do
    end do
  end subroutine


  subroutine jacobian_set_single_dyn(this, itab1, idx1, idx2, ival1, ival2, d)
    !! set one element of multidimensional derivative for dynamic stencil.
    !! see arguments' documentation in jacobian_edit_single_dyn.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: idx2(:)
    integer,         intent(in)    :: ival1
    integer,         intent(in)    :: ival2
    real,            intent(in)    :: d

    call this%jacobian_edit_single_dyn(itab1, idx1, idx2, ival1, ival2, d, set=.true.)
  end subroutine

  subroutine jacobian_update_single_dyn(this, itab1, idx1, idx2, ival1, ival2, d)
    !! update one element of multidimensional derivative for dynamic stencil.
    !! see arguments' documentation in jacobian_edit_single_dyn.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
    integer,         intent(in)    :: idx1(:)
    integer,         intent(in)    :: idx2(:)
    integer,         intent(in)    :: ival1
    integer,         intent(in)    :: ival2
    real,            intent(in)    :: d

    call this%jacobian_edit_single_dyn(itab1, idx1, idx2, ival1, ival2, d, update=.true.)
  end subroutine

  subroutine jacobian_edit_single_dyn(this, itab1, idx1, idx2, ival1, ival2, d, set, update)
    !! set/update one element of multidimensional derivative for dynamic stencil
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! grid table index
    integer,           intent(in)    :: idx1(:)
      !! result index
    integer,           intent(in)    :: idx2(:)
      !! dependency index
    integer,           intent(in)    :: ival1
      !! result dimension
    integer,           intent(in)    :: ival2
      !! dependency dimension
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: set, update
      !! flags to either set or update the derivative. (default: false)
      !! exactly one must be true.

    integer :: itab2, col, col0, row, row0
    logical :: set_, update_

    ! optional flags
    set_ = .false.
    if (present(set)) set_ = set
    update_ = .false.
    if (present(update)) update_ = update
    ASSERT(set_ .neqv. update_)

    ! get row base index
    row0 = (this%matr%v1%tab(itab1)%p%get_flat(idx1) - 1) * this%matr%v1%nval

    ! get column base index
    itab2 = this%matr%v2%itab%get(idx2)
    col0  = (this%matr%v2%tab(itab2)%p%get_flat(idx2) - 1) * this%matr%v2%nval

    ! set derivatives
    row = row0 + ival1
    col = col0 + ival2
    if      (set_   ) then
      call this%matr%sb(itab1,itab2)%set(row, col, d, search = .false.)
    else if (update_) then
      call this%matr%sb(itab1,itab2)%add(row, col, d)
    end if
  end subroutine


  subroutine jacobian_set_scalar_dyn(this, itab1, idx1, idx2, d)
    !! set scalar derivative for dynamic stencil of one-dimensional variables.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
      !! grid table index
    integer,         intent(in)    :: idx1(:)
      !! result index
    integer,         intent(in)    :: idx2(:)
      !! dependency index
    real,            intent(in)    :: d
      !! derivative

    ASSERT(this%matr%v1%nval == 1)
    ASSERT(this%matr%v2%nval == 1)

    call this%set(itab1, idx1, idx2, 1, 1, d)
  end subroutine

  subroutine jacobian_update_scalar_dyn(this, itab1, idx1, idx2, d)
    !! update scalar derivative for dynamic stencil of one-dimensional variables.
    class(jacobian), intent(inout) :: this
    integer,         intent(in)    :: itab1
      !! grid table index
    integer,         intent(in)    :: idx1(:)
      !! result index
    integer,         intent(in)    :: idx2(:)
      !! dependency index
    real,            intent(in)    :: d
      !! derivative

    ASSERT(this%matr%v1%nval == 1)
    ASSERT(this%matr%v2%nval == 1)

    call this%update(itab1, idx1, idx2, 1, 1, d)
  end subroutine

end module
