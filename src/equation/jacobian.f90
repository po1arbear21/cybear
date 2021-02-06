module jacobian_m
  use array_m
  use jacobian_matrix_m
  use stencil_m
  implicit none

  type jacobian
    !! jacobian structure to store derivatives

    type(jacobian_matrix) :: matr
      !! jacobian in blockmatrix representation

    type(stencil_ptr), allocatable :: st(:)
      !! stencils (v1%ntab)
    type(array4_real), allocatable :: d(:)
      !! derivatives (v1%nval x v2%nval x st(itab1)%max_ndep x v1%tab(itab1)%n) x (v1%ntab)
  contains
    procedure :: init     => jacobian_init
    procedure :: destruct => jacobian_destruct
    procedure :: reset    => jacobian_reset
    procedure :: set_matr => jacobian_set_matr
  end type

  type jacobian_ptr
    type(jacobian), pointer :: p => null()
  end type

contains

  subroutine jacobian_init(this, v1, v2, st, const, zero, valmsk)
    !! initialize jacobian
    class(jacobian),         intent(out) :: this
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

    integer :: i, j, itab1, itab2, ndep
    logical :: const_(v1%ntab,v2%ntab), zero_(v1%ntab,v2%ntab), valmsk_(v1%nval,v2%nval,v1%ntab,v2%ntab)

    ! optional arguments
    const_ = .false.
    if (present(const)) const_ = const
    if (present(zero)) then
      zero_ = zero
    else
      ! set zero flags automatically by checking stencils
      zero_ = .true.
      do itab1 = 1, v1%ntab
        ! do nothing for empty stencil
        if (.not. associated(st(itab1)%p)) cycle

        block
          integer :: idx1(v1%g%idx_dim), idx2(v2%g%idx_dim, st(itab1)%p%max_ndep)

          ! loop over result points
          do i = 1, v1%tab(itab1)%p%n
            ! get result grid indices
            idx1 = v1%tab(itab1)%p%get_idx(i)

            ! get dependency grid indices
            call st(itab1)%p%get_dep(idx1, idx2, ndep)

            ! loop over dependency grid indices
            do j = 1, ndep
              ! stencil connects itab1 to itab2 => not zero
              itab2 = v2%itab%get(idx2(:,j))
              zero_(itab1,itab2) = .false.
            end do
          end do
        end block
      end do
    end if
    const_ = (const_ .or. zero_) ! zero blocks are also constant
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

    ! allocate memory for derivatives
    allocate (this%d(v1%ntab))
    do itab1 = 1, v1%ntab
      if (.not. associated(st(itab1)%p)) cycle
      allocate (this%d(itab1)%d(v1%nval, v2%nval, st(itab1)%p%max_ndep, v1%tab(itab1)%p%n), source = 0.0)
    end do
  end subroutine

  subroutine jacobian_destruct(this)
    !! destruct jacobian
    class(jacobian), intent(inout) :: this

    call this%matr%destruct()
    if (allocated(this%st)) deallocate (this%st)
    if (allocated(this%d )) deallocate (this%d)
  end subroutine

  subroutine jacobian_reset(this, const, nonconst)
    !! reset jacobian to zeros
    class(jacobian),   intent(inout) :: this
    logical, optional, intent(in)    :: const
      !! reset constant blocks (default: true)
    logical, optional, intent(in)    :: nonconst
      !! reset non-constant blocks (default: true)

    integer :: i, j, itab1, itab2, ndep
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
        if (all(reset(itab1,:))) then
          this%d(itab1)%d = 0
          cycle
        end if

        associate (st => this%st(itab1)%p)
          block
            integer :: idx1(v1%g%idx_dim), idx2(v2%g%idx_dim, st%max_ndep)

            ! loop over result points
            do i = 1, v1%tab(itab1)%p%n
              ! get result grid indices
              idx1 = v1%tab(itab1)%p%get_idx(i)

              ! get dependency grid indices
              call st%get_dep(idx1, idx2, ndep)

              ! loop over dependency grid indices
              do j = 1, ndep
                ! get dependency table index
                itab2 = v2%itab%get(idx2(:,j))

                ! reset only if flag is set
                if (reset(itab1,itab2)) this%d(itab1)%d(:,:,j,i) = 0
              end do
            end do
          end block
        end associate
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

    integer               :: i, j, itab1, itab2, ival1, ival2, ndep, row, row0, col, col0
    logical               :: const_, nonconst_
    type(sparse_ptr_real) :: s(this%matr%v2%ntab)
    type(spbuild_real)    :: sb(this%matr%v2%ntab)

    ! optional arguments
    const_ = .true.
    if (present(const)) const_ = const
    nonconst_ = .true.
    if (present(nonconst)) nonconst_ = nonconst

    associate (v1 => this%matr%v1, v2 => this%matr%v2)
      ! loop over all blocks
      do itab1 = 1, v1%ntab
        ! decide whether to process any blocks of this block row
        if (all(      this%matr%const(itab1,:)) .and. .not.    const_) cycle
        if (all(.not. this%matr%const(itab1,:)) .and. .not. nonconst_) cycle

        ! get sparse matrices and builders
        do itab2 = 1, v2%ntab
          nullify (s(itab2)%p)

          ! decide whether to process block
          if (this%matr%const(itab1, itab2) .and. .not. const_) cycle
          if (.not. this%matr%const(itab1, itab2) .and. .not. nonconst_) cycle

          ! get sparse matrix block
          call this%matr%get(itab1, itab2, s(itab2)%p)

          ! init sparse builder
          call sb(itab2)%init(s(itab2)%p)

          if ((.not. associated(this%st(itab1)%p)) .or. this%matr%zero(itab1, itab2)) then
            ! finish this block early
            call sb(itab2)%save()
            nullify (s(itab2)%p)
          end if
        end do

        ! do nothing for empty stencil or if all blocks in this block row are zero
        if ((.not. associated(this%st(itab1)%p)) .or. all(this%matr%zero(itab1,:))) cycle

        associate (st => this%st(itab1)%p)
          block
            integer :: idx1(v1%g%idx_dim)
            integer :: idx2(v2%g%idx_dim,st%max_ndep)

            ! loop over points in result table
            do i = 1, v1%tab(itab1)%p%n
              ! get result grid indices
              idx1 = v1%tab(itab1)%p%get_idx(i)

              ! get row base index
              row0 = (i - 1) * v1%nval

              ! get dependency points
              call st%get_dep(idx1, idx2, ndep)

              ! loop over dependency points
              do j = 1, ndep
                ! get dependency table index
                itab2 = v2%itab%get(idx2(:,j))

                ! do nothing if dependency matrix is not used (const or zero flag prevents it)
                if (.not. associated(s(itab2)%p)) cycle

                ! get column base index
                col0 = (v2%tab(itab2)%p%get_flat(idx2(:,j)) - 1) * v2%nval

                ! set derivatives
                do ival1 = 1, v1%nval
                  row = row0 + ival1
                  do ival2 = 1, v2%nval
                    if (.not. this%matr%valmsk(ival1,ival2,itab1,itab2)) cycle
                    col = col0 + ival2
                    call sb(itab2)%set(row, col, this%d(itab1)%d(ival1,ival2,j,i), search = .false.)
                  end do
                end do
              end do
            end do
          end block
        end associate

        ! save values in sparse matrices
        do itab2 = 1, v2%ntab
          if (associated(s(itab2)%p)) call sb(itab2)%save()
        end do
      end do
    end associate
  end subroutine

end module
