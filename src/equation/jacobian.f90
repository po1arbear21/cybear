#include "../util/macro.f90.inc"

module jacobian_m
  use error_m
  use hashmap_m,         only: hashmap_int
  use jacobian_matrix_m, only: jacobian_matrix
  use stencil_m,         only: stencil_ptr, static_stencil, dynamic_stencil
  use vselector_m,       only: vselector

  implicit none

  private
  public jacobian, jacobian_ptr

  type static_data
    !! derivative data for static stencils
    real, allocatable :: d(:,:,:,:)
      !! derivatives (v1%nval x v2%nval x st%nmax x tab1%n)
    type(hashmap_int) :: hmap
      !! (i1, itab2, i2) -> j (= 3rd index in d)
  contains
    procedure :: init     => static_data_init
    procedure :: destruct => static_data_destruct
  end type

  type jacobian
    !! jacobian structure to store derivatives

    type(jacobian_matrix) :: matr
      !! jacobian in blockmatrix representation

    type(stencil_ptr), allocatable :: st(:)
      !! stencils (v1%ntab)
    type(static_data), allocatable :: sd(:)
      !! derivatives for static stencils (v1%ntab)
  contains
    procedure :: init     => jacobian_init
    procedure :: destruct => jacobian_destruct
    procedure :: reset    => jacobian_reset
    procedure :: set_matr => jacobian_set_matr
    generic   :: set      => jacobian_set_itab_nval, &
      &                      jacobian_set_itab_ival, &
      &                      jacobian_set_itab,      &
      &                      jacobian_set_nval,      &
      &                      jacobian_set_ival,      &
      &                      jacobian_set
    generic   :: add      => jacobian_add_itab_nval, &
      &                      jacobian_add_itab_ival, &
      &                      jacobian_add_itab,      &
      &                      jacobian_add_nval,      &
      &                      jacobian_add_ival,      &
      &                      jacobian_add

    procedure, private :: jacobian_set_itab_nval, &
      &                   jacobian_set_itab_ival, &
      &                   jacobian_set_itab,      &
      &                   jacobian_set_nval,      &
      &                   jacobian_set_ival,      &
      &                   jacobian_set
    procedure, private :: jacobian_add_itab_nval, &
      &                   jacobian_add_itab_ival, &
      &                   jacobian_add_itab,      &
      &                   jacobian_add_nval,      &
      &                   jacobian_add_ival,      &
      &                   jacobian_add
  end type

  type jacobian_ptr
    type(jacobian), pointer :: p => null()
  end type

contains

  subroutine static_data_init(this, v1, v2, itab1, st)
    class(static_data),    intent(out) :: this
    type(vselector),       intent(in)  :: v1
      !! result variable selector
    type(vselector),       intent(in)  :: v2
      !! dependency variable selector
    integer,               intent(in)  :: itab1
      !! result table index
    class(static_stencil), intent(in)  :: st
      !! stencil

    integer :: i1, itab2, i2, j, idx1(v1%g%idx_dim), idx2(v2%g%idx_dim)
    logical :: status

    allocate (this%d(v1%nval, v2%nval, st%nmax, v1%tab(itab1)%p%n), source = 0.0)

    ! init hashmap
    call this%hmap%init(c = v1%tab(itab1)%p%n * st%nmax)
    do i1 = 1, v1%tab(itab1)%p%n
      ! get result grid indices
      idx1 = v1%tab(itab1)%p%get_idx(i1)

      ! loop over dependency grid indices
      do j = 1, st%nmax
        ! get dependency grid indices
        call st%get(idx1, j, idx2, status)
        if (.not. status) exit

        ! lookup idx2
        itab2 = v2%itab%get(idx2)
        i2    = v2%tab(itab2)%p%get_flat(idx2)

        ! save stencil dependency index (j) in hashmap
        call this%hmap%set([i1, itab2, i2], j)
      end do
    end do
  end subroutine

  subroutine static_data_destruct(this)
    class(static_data), intent(inout) :: this

    if (allocated(this%d)) deallocate(this%d)
    call this%hmap%destruct()
  end subroutine

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
      ASSERT(all(shape(zero) == [v1%ntab, v2%ntab]))
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
            zero_(itab1,:) = .false.

          class default
            call program_error("stencil type not supported; derive from static or dynamic stencil!")
        end select
      end do
    end if
    const_  = (const_ .or. zero_) ! zero blocks are also constant
    valmsk_ = .true.
    if (present(valmsk)) then
      ASSERT(all(shape(valmsk) == [v1%ntab, v2%ntab]))
      do itab2 = 1, v2%ntab; do itab1 = 1, v1%ntab
        valmsk_(:,:,itab1,itab2) = valmsk
      end do; end do
    end if

    ! init matrix
    call this%matr%init(v1, v2, const_, zero_, valmsk_)

    ! save stencil pointers
    this%st = st

    ! allocate memory for derivatives of static stencil data
    allocate (this%sd(v1%ntab))
    do itab1 = 1, v1%ntab
      select type (st_ptr => st(itab1)%p)
        class is (static_stencil)
          call this%sd(itab1)%init(v1, v2, itab1, st_ptr)
        class is (dynamic_stencil)
          ASSERT(.not. all(const_))
      end select
    end do
  end subroutine

  subroutine jacobian_destruct(this)
    !! destruct jacobian
    class(jacobian), intent(inout) :: this

    integer :: i

    call this%matr%destruct()
    if (allocated(this%st)) deallocate (this%st)
    if (allocated(this%sd)) then
      do i = 1, size(this%sd)
        call this%sd(i)%destruct()
      end do
    end if
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
              this%sd(itab1)%d = 0
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
                  if (reset(itab1,itab2)) this%sd(itab1)%d(:,:,j,i) = 0
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
                        call this%matr%sb(itab1,itab2)%set(row, col, this%sd(itab1)%d(ival1,ival2,j,i), search = .false.)
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

  subroutine jacobian_set_itab_nval(this, itab1, i1, idx2, d, add)
    !! set/update derivatives (select result grid_table, non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d(:,:)
      !! derivatives (nval1 x nval2)
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    integer :: itab2, i2, j, row, row0, col, col0, ival1, ival2
    logical :: add_

    ASSERT(size(d, dim=1) == this%matr%v1%nval)
    ASSERT(size(d, dim=2) == this%matr%v2%nval)

    ! optional addition flag
    add_ = .false.
    if (present(add)) add_ = add

    ! lookup idx2
    itab2 = this%matr%v2%itab%get(idx2)
    i2    = this%matr%v2%tab(itab2)%p%get_flat(idx2)

    select type (st => this%st(itab1)%p)
      class is (static_stencil)
        ! get stencil dependency index
        call this%sd(itab1)%hmap%get([i1, itab2, i2], j)

        ! edit derivatives
        if (add_) then
          this%sd(itab1)%d(:,:,j,i1) = this%sd(itab1)%d(:,:,j,i1) + d
        else
          this%sd(itab1)%d(:,:,j,i1) = d
        end if

      class is (dynamic_stencil)
        ! get row base index
        row0 = (i1 - 1) * this%matr%v1%nval

        ! get column base index
        col0 = (i2 - 1) * this%matr%v2%nval

        ! edit derivatives
        do ival1 = 1, this%matr%v1%nval
          row = row0 + ival1
          do ival2 = 1, this%matr%v2%nval
            col = col0 + ival2
            if (add_) then
              call this%matr%sb(itab1,itab2)%add(row, col, d(ival1, ival2))
            else
              call this%matr%sb(itab1,itab2)%set(row, col, d(ival1, ival2))
            end if
          end do
        end do
    end select
  end subroutine

  subroutine jacobian_set_itab_ival(this, itab1, i1, idx2, ival1, ival2, d, add)
    !! set/update one element of derivatives (select result grid_table, non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    integer,           intent(in)    :: ival1
      !! result value index
    integer,           intent(in)    :: ival2
      !! dependency value index
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    integer :: itab2, i2, j, row, col
    logical :: add_

    ! optional addition flag
    add_ = .false.
    if (present(add)) add_ = add

    ! lookup idx2
    itab2 = this%matr%v2%itab%get(idx2)
    i2    = this%matr%v2%tab(itab2)%p%get_flat(idx2)

    select type (st => this%st(itab1)%p)
      class is (static_stencil)
        ! get stencil dependency index
        call this%sd(itab1)%hmap%get([i1, itab2, i2], j)

        ! edit derivatives
        if (add_) then
          this%sd(itab1)%d(ival1,ival2,j,i1) = this%sd(itab1)%d(ival1,ival2,j,i1) + d
        else
          this%sd(itab1)%d(ival1,ival2,j,i1) = d
        end if

      class is (dynamic_stencil)
        ! get row index
        row = (i1 - 1) * this%matr%v1%nval + ival1

        ! get column index
        col = (i2 - 1) * this%matr%v2%nval + ival2

        ! edit derivative
        if (add_) then
          call this%matr%sb(itab1,itab2)%add(row, col, d)
        else
          call this%matr%sb(itab1,itab2)%set(row, col, d)
        end if
    end select
  end subroutine

  subroutine jacobian_set_itab(this, itab1, i1, idx2, d, add)
    !! set/update derivative (select result grid_table, scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    ASSERT(this%matr%v1%nval == 1)
    ASSERT(this%matr%v2%nval == 1)

    call this%jacobian_set_itab_ival(itab1, i1, idx2, 1, 1, d, add = add)
  end subroutine

  subroutine jacobian_set_nval(this, idx1, idx2, d, add)
    !! set/update derivatives (non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d(:,:)
      !! derivatives (nval1 x nval2)
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    integer :: itab1, i1

    ! lookup idx1
    itab1 = this%matr%v1%itab%get(idx1)
    i1    = this%matr%v1%tab(itab1)%p%get_flat(idx1)

    call this%jacobian_set_itab_nval(itab1, i1, idx2, d, add = add)
  end subroutine

  subroutine jacobian_set_ival(this, idx1, idx2, ival1, ival2, d, add)
    !! set/update one element of derivatives (non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    integer,           intent(in)    :: ival1
      !! result value index
    integer,           intent(in)    :: ival2
      !! dependency value index
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    integer :: itab1, i1

    ! lookup idx1
    itab1 = this%matr%v1%itab%get(idx1)
    i1    = this%matr%v1%tab(itab1)%p%get_flat(idx1)

    call this%jacobian_set_itab_ival(itab1, i1, idx2, ival1, ival2, d, add = add)
  end subroutine

  subroutine jacobian_set(this, idx1, idx2, d, add)
    !! set/update derivative (scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d
      !! derivative
    logical, optional, intent(in)    :: add
      !! addition flag (default: false)

    integer :: itab1, i1

    ! lookup idx1
    itab1 = this%matr%v1%itab%get(idx1)
    i1    = this%matr%v1%tab(itab1)%p%get_flat(idx1)

    call this%jacobian_set_itab(itab1, i1, idx2, d, add = add)
  end subroutine


  subroutine jacobian_add_itab_nval(this, itab1, i1, idx2, d)
    !! update derivatives (select result grid_table, non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d(:,:)
      !! derivatives (nval1 x nval2)

    call this%jacobian_set_itab_nval(itab1, i1, idx2, d, add = .true.)
  end subroutine

  subroutine jacobian_add_itab_ival(this, itab1, i1, idx2, ival1, ival2, d)
    !! update one element of derivatives (select result grid_table, non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    integer,           intent(in)    :: ival1
      !! result value index
    integer,           intent(in)    :: ival2
      !! dependency value index
    real,              intent(in)    :: d
      !! derivative

    call this%jacobian_set_itab_ival(itab1, i1, idx2, ival1, ival2, d, add = .true.)
  end subroutine

  subroutine jacobian_add_itab(this, itab1, i1, idx2, d)
    !! update derivative (select result grid_table, scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: itab1
      !! result grid_table index
    integer,           intent(in)    :: i1
      !! result grid_table entry index
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d
      !! derivative

    call this%jacobian_set_itab(itab1, i1, idx2, d, add = .true.)
  end subroutine

  subroutine jacobian_add_nval(this, idx1, idx2, d)
    !! update derivatives (non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d(:,:)
      !! derivatives (nval1 x nval2)

    call this%jacobian_set_nval(idx1, idx2, d, add = .true.)
  end subroutine

  subroutine jacobian_add_ival(this, idx1, idx2, ival1, ival2, d)
    !! update one element of derivatives (non-scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    integer,           intent(in)    :: ival1
      !! result value index
    integer,           intent(in)    :: ival2
      !! dependency value index
    real,              intent(in)    :: d
      !! derivative

    call this%jacobian_set_ival(idx1, idx2, ival1, ival2, d, add = .true.)
  end subroutine

  subroutine jacobian_add(this, idx1, idx2, d)
    !! update derivative (scalar variables)
    class(jacobian),   intent(inout) :: this
    integer,           intent(in)    :: idx1(:)
      !! result grid indices
    integer,           intent(in)    :: idx2(:)
      !! dependency grid indices
    real,              intent(in)    :: d
      !! derivative

    call this%jacobian_set(idx1, idx2, d, add = .true.)
  end subroutine
end module
