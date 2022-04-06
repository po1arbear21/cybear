#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define GRID_DATA_TYPE         CONCAT(grid_data,T)
#define GRID_DATA0_TYPE        CONCAT(grid_data0,T)
#define GRID_DATA1_TYPE        CONCAT(grid_data1,T)
#define GRID_DATA2_TYPE        CONCAT(grid_data2,T)
#define GRID_DATA3_TYPE        CONCAT(grid_data3,T)
#define GRID_DATA4_TYPE        CONCAT(grid_data4,T)
#define GRID_DATA5_TYPE        CONCAT(grid_data5,T)
#define GRID_DATA6_TYPE        CONCAT(grid_data6,T)
#define GRID_DATA7_TYPE        CONCAT(grid_data7,T)
#define GRID_DATA8_TYPE        CONCAT(grid_data8,T)

#define GRID_DATA_INIT         CONCAT3(grid_data,T,init)
#define GRID_DATA_DESTRUCT     CONCAT3(grid_data,T,destruct)
#define GRID_DATA_RESET        CONCAT3(grid_data,T,reset)
#define GRID_DATA_GET_PTR0     CONCAT3(grid_data,T,get_ptr0)
#define GRID_DATA_GET_PTR1     CONCAT3(grid_data,T,get_ptr1)
#define GRID_DATA_GET_PTR2     CONCAT3(grid_data,T,get_ptr2)
#define GRID_DATA_GET_PTR3     CONCAT3(grid_data,T,get_ptr3)
#define GRID_DATA_GET_PTR4     CONCAT3(grid_data,T,get_ptr4)
#define GRID_DATA_GET_PTR5     CONCAT3(grid_data,T,get_ptr5)
#define GRID_DATA_GET_PTR6     CONCAT3(grid_data,T,get_ptr6)
#define GRID_DATA_GET_PTR7     CONCAT3(grid_data,T,get_ptr7)
#define GRID_DATA_GET_PTR8     CONCAT3(grid_data,T,get_ptr8)
#define GRID_DATA_GET_POINT    CONCAT3(grid_data,T,get_point)
#define GRID_DATA_GET_ALL      CONCAT3(grid_data,T,get_all)
#define GRID_DATA_SET_POINT    CONCAT3(grid_data,T,set_point)
#define GRID_DATA_SET_ALL      CONCAT3(grid_data,T,set_all)
#define GRID_DATA_UPDATE_POINT CONCAT3(grid_data,T,update_point)
#define GRID_DATA_UPDATE_ALL   CONCAT3(grid_data,T,update_all)

#define ALLOCATE_GRID_DATA0    CONCAT(allocate_grid_data0,T)
#define ALLOCATE_GRID_DATA1    CONCAT(allocate_grid_data1,T)
#define ALLOCATE_GRID_DATA2    CONCAT(allocate_grid_data2,T)
#define ALLOCATE_GRID_DATA3    CONCAT(allocate_grid_data3,T)
#define ALLOCATE_GRID_DATA4    CONCAT(allocate_grid_data4,T)

module subroutine ALLOCATE_GRID_DATA0(gdata, dim)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata
  integer,                            intent(in)  :: dim
    !! grid data dimension

  select case (dim)
  case (0)
    allocate (GRID_DATA0_TYPE :: gdata)
  case (1)
    allocate (GRID_DATA1_TYPE :: gdata)
  case (2)
    allocate (GRID_DATA2_TYPE :: gdata)
  case (3)
    allocate (GRID_DATA3_TYPE :: gdata)
  case (4)
    allocate (GRID_DATA4_TYPE :: gdata)
  case (5)
    allocate (GRID_DATA5_TYPE :: gdata)
  case (6)
    allocate (GRID_DATA6_TYPE :: gdata)
  case (7)
    allocate (GRID_DATA7_TYPE :: gdata)
  case (8)
    allocate (GRID_DATA8_TYPE :: gdata)
  case default
    call program_error("grid data dimension must be in range 0:8")
  end select
end subroutine

module subroutine ALLOCATE_GRID_DATA1(gdata, dim, i0, i1)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata(:)
  integer,                            intent(in)  :: dim
    !! grid data dimension
  integer,                            intent(in)  :: i0
    !! lower array bound
  integer,                            intent(in)  :: i1
    !! upper array bound

  select case (dim)
  case (0)
    allocate (GRID_DATA0_TYPE :: gdata(i0:i1))
  case (1)
    allocate (GRID_DATA1_TYPE :: gdata(i0:i1))
  case (2)
    allocate (GRID_DATA2_TYPE :: gdata(i0:i1))
  case (3)
    allocate (GRID_DATA3_TYPE :: gdata(i0:i1))
  case (4)
    allocate (GRID_DATA4_TYPE :: gdata(i0:i1))
  case (5)
    allocate (GRID_DATA5_TYPE :: gdata(i0:i1))
  case (6)
    allocate (GRID_DATA6_TYPE :: gdata(i0:i1))
  case (7)
    allocate (GRID_DATA7_TYPE :: gdata(i0:i1))
  case (8)
    allocate (GRID_DATA8_TYPE :: gdata(i0:i1))
  case default
    call program_error("grid data dimension must be in range 0:8")
  end select
end subroutine

module subroutine ALLOCATE_GRID_DATA2(gdata, dim, i0, i1)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata(:,:)
  integer,                            intent(in)  :: dim
    !! grid data dimension
  integer,                            intent(in)  :: i0(2)
    !! lower array bounds
  integer,                            intent(in)  :: i1(2)
    !! upper array bounds

  select case (dim)
  case (0)
    allocate (GRID_DATA0_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (1)
    allocate (GRID_DATA1_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (2)
    allocate (GRID_DATA2_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (3)
    allocate (GRID_DATA3_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (4)
    allocate (GRID_DATA4_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (5)
    allocate (GRID_DATA5_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (6)
    allocate (GRID_DATA6_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (7)
    allocate (GRID_DATA7_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case (8)
    allocate (GRID_DATA8_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2)))
  case default
    call program_error("grid data dimension must be in range 0:8")
  end select
end subroutine

module subroutine ALLOCATE_GRID_DATA3(gdata, dim, i0, i1)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata(:,:,:)
  integer,                            intent(in)  :: dim
    !! grid data dimension
  integer,                            intent(in)  :: i0(3)
    !! lower array bounds
  integer,                            intent(in)  :: i1(3)
    !! upper array bounds

  select case (dim)
  case (0)
    allocate (GRID_DATA0_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (1)
    allocate (GRID_DATA1_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (2)
    allocate (GRID_DATA2_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (3)
    allocate (GRID_DATA3_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (4)
    allocate (GRID_DATA4_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (5)
    allocate (GRID_DATA5_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (6)
    allocate (GRID_DATA6_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (7)
    allocate (GRID_DATA7_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case (8)
    allocate (GRID_DATA8_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  case default
    call program_error("grid data dimension must be in range 0:8")
  end select
end subroutine

module subroutine ALLOCATE_GRID_DATA4(gdata, dim, i0, i1)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata(:,:,:,:)
  integer,                            intent(in)  :: dim
    !! grid data dimension
  integer,                            intent(in)  :: i0(4)
    !! lower array bounds
  integer,                            intent(in)  :: i1(4)
    !! upper array bounds

  select case (dim)
  case (0)
    allocate (GRID_DATA0_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (1)
    allocate (GRID_DATA1_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (2)
    allocate (GRID_DATA2_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (3)
    allocate (GRID_DATA3_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (4)
    allocate (GRID_DATA4_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (5)
    allocate (GRID_DATA5_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (6)
    allocate (GRID_DATA6_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (7)
    allocate (GRID_DATA7_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case (8)
    allocate (GRID_DATA8_TYPE :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  case default
    call program_error("grid data dimension must be in range 0:8")
  end select
end subroutine

module subroutine GRID_DATA_INIT(this, g, idx_type, idx_dir, d0)
  !! initialize and allocate grid data
  class(GRID_DATA_TYPE), intent(out) :: this
  class(grid),           intent(in)  :: g
  integer,               intent(in)  :: idx_type
    !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
  integer,               intent(in)  :: idx_dir
    !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
  TT, optional,          intent(in)  :: d0
    !! initial value

  TT :: d0_

  ! get index bounds
  allocate (this%idx_bnd(g%idx_dim))
  call g%get_idx_bnd(idx_type, idx_dir, this%idx_bnd)

  ! set total number of values
  this%n = product(this%idx_bnd) ! = 1 for idx_dim == 0

  ! get initial value
  if (present(d0)) then
    d0_ = d0
  else
#ifdef TLOG
    d0_ = .false.
#else
    d0_ = 0
#endif
  end if

  ! allocate data array
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = d0_
  type is (GRID_DATA1_TYPE)
    allocate (this%data(this%idx_bnd(1)), source = d0_)
  type is (GRID_DATA2_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2)), source = d0_)
  type is (GRID_DATA3_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3)), source = d0_)
  type is (GRID_DATA4_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3), &
      &                 this%idx_bnd(4)), source = d0_)
  type is (GRID_DATA5_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3), &
      &                 this%idx_bnd(4), &
      &                 this%idx_bnd(5)), source = d0_)
  type is (GRID_DATA6_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3), &
      &                 this%idx_bnd(4), &
      &                 this%idx_bnd(5), &
      &                 this%idx_bnd(6)), source = d0_)
  type is (GRID_DATA7_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3), &
      &                 this%idx_bnd(4), &
      &                 this%idx_bnd(5), &
      &                 this%idx_bnd(6), &
      &                 this%idx_bnd(7)), source = d0_)
  type is (GRID_DATA8_TYPE)
    allocate (this%data(this%idx_bnd(1), &
      &                 this%idx_bnd(2), &
      &                 this%idx_bnd(3), &
      &                 this%idx_bnd(4), &
      &                 this%idx_bnd(5), &
      &                 this%idx_bnd(6), &
      &                 this%idx_bnd(7), &
      &                 this%idx_bnd(8)), source = d0_)
  end select
end subroutine

module function GRID_DATA_GET_PTR0(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA0_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA0_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR1(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA1_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA1_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR2(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA2_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA2_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR3(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA3_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA3_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR4(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA4_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA4_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR5(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA5_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA5_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR6(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA6_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA6_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR7(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA7_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA7_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module function GRID_DATA_GET_PTR8(this) result(p)
  class(GRID_DATA_TYPE), target, intent(in) :: this
  type(GRID_DATA8_TYPE), pointer            :: p

  select type (this)
  type is (GRID_DATA8_TYPE)
    p => this
  class default
    p => null()
  end select
end function

module subroutine GRID_DATA_DESTRUCT(this)
  !! destruct grid data
  class(GRID_DATA_TYPE), intent(inout) :: this

  if (allocated(this%idx_bnd)) deallocate (this%idx_bnd)
  select type (this)
  type is (GRID_DATA1_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA2_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA3_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA4_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA5_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA6_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA7_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  type is (GRID_DATA8_TYPE)
    if (allocated(this%data)) deallocate (this%data)
  end select
end subroutine

module subroutine GRID_DATA_RESET(this, d0)
  !! reset data for all points
  class(GRID_DATA_TYPE), intent(inout) :: this
  TT, optional,          intent(in)    :: d0
    !! initial value

  TT :: d0_

  ! get initial value
  if (present(d0)) then
    d0_ = d0
  else
#ifdef TLOG
    d0_ = .false.
#else
    d0_ = 0
#endif
  end if

  ! set data array
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = d0_
  type is (GRID_DATA1_TYPE)
    this%data = d0_
  type is (GRID_DATA2_TYPE)
    this%data = d0_
  type is (GRID_DATA3_TYPE)
    this%data = d0_
  type is (GRID_DATA4_TYPE)
    this%data = d0_
  type is (GRID_DATA5_TYPE)
    this%data = d0_
  type is (GRID_DATA6_TYPE)
    this%data = d0_
  type is (GRID_DATA7_TYPE)
    this%data = d0_
  type is (GRID_DATA8_TYPE)
    this%data = d0_
  end select
end subroutine

module function GRID_DATA_GET_POINT(this, idx) result(d)
  !! get data for single point with bounds check (out of bounds: return default value)
  class(GRID_DATA_TYPE), intent(in) :: this
  integer,               intent(in) :: idx(:)
    !! grid indices (idx_dim)
  TT                                :: d
    !! return data

  ASSERT(size(idx) == size(this%idx_bnd)             )
  ASSERT(all(idx >= 1) .and. all(idx <= this%idx_bnd))

  ! default value
#ifdef TLOG
  d = .false.
#else
  d = 0
#endif

  ! get data
  select type (this)
  type is (GRID_DATA0_TYPE)
    d = this%data
  type is (GRID_DATA1_TYPE)
    d = this%data(idx(1))
  type is (GRID_DATA2_TYPE)
    d = this%data(idx(1),idx(2))
  type is (GRID_DATA3_TYPE)
    d = this%data(idx(1),idx(2),idx(3))
  type is (GRID_DATA4_TYPE)
    d = this%data(idx(1),idx(2),idx(3),idx(4))
  type is (GRID_DATA5_TYPE)
    d = this%data(idx(1),idx(2),idx(3),idx(4),idx(5))
  type is (GRID_DATA6_TYPE)
    d = this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
  type is (GRID_DATA7_TYPE)
    d = this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7))
  type is (GRID_DATA8_TYPE)
    d = this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8))
  end select
end function

module function GRID_DATA_GET_ALL(this) result(d)
  !! get data for all points in flat array
  class(GRID_DATA_TYPE), intent(in) :: this
  TT                                :: d(this%n)
    !! return all data

  select type (this)
  type is (GRID_DATA0_TYPE)
    d = this%data
  type is (GRID_DATA1_TYPE)
    d = this%data
  type is (GRID_DATA2_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA3_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA4_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA5_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA6_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA7_TYPE)
    d = reshape(this%data, [this%n])
  type is (GRID_DATA8_TYPE)
    d = reshape(this%data, [this%n])
  end select
end function

module subroutine GRID_DATA_SET_POINT(this, idx, d)
  !! set data for single point with bounds check (do nothing if out of bounds)
  class(GRID_DATA_TYPE), intent(inout) :: this
  integer,               intent(in)    :: idx(:)
    !! grid indices (idx_dim)
  TT,                    intent(in)    :: d
    !! new value

  ASSERT(size(idx) == size(this%idx_bnd)             )
  ASSERT(all(idx >= 1) .and. all(idx <= this%idx_bnd))

  ! set data
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = d
  type is (GRID_DATA1_TYPE)
    this%data(idx(1)) = d
  type is (GRID_DATA2_TYPE)
    this%data(idx(1),idx(2)) = d
  type is (GRID_DATA3_TYPE)
    this%data(idx(1),idx(2),idx(3)) = d
  type is (GRID_DATA4_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4)) = d
  type is (GRID_DATA5_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5)) = d
  type is (GRID_DATA6_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6)) = d
  type is (GRID_DATA7_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7)) = d
  type is (GRID_DATA8_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8)) = d
  end select
end subroutine

module subroutine GRID_DATA_SET_ALL(this, d)
  !! set data for all points
  class(GRID_DATA_TYPE), intent(inout) :: this
  TT,                    intent(in)    :: d(:)
    !! new values

  ASSERT(size(d) == this%n)

  ! set data
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = d(1)
  type is (GRID_DATA1_TYPE)
    this%data = d
  type is (GRID_DATA2_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2)])
  type is (GRID_DATA3_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3)])
  type is (GRID_DATA4_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4)])
  type is (GRID_DATA5_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5)])
  type is (GRID_DATA6_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6)])
  type is (GRID_DATA7_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6), this%idx_bnd(7)])
  type is (GRID_DATA8_TYPE)
    this%data = reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6), this%idx_bnd(7), this%idx_bnd(8)])
  end select
end subroutine

#ifndef TLOG
module subroutine GRID_DATA_UPDATE_POINT(this, idx, d)
  !! update data for single point with bounds check (do nothing if out of bounds)
  class(GRID_DATA_TYPE), intent(inout) :: this
  integer,               intent(in)    :: idx(:)
    !! grid indices (idx_dim)
  TT,                    intent(in)    :: d
    !! delta value

  ASSERT(size(idx) == size(this%idx_bnd)             )
  ASSERT(all(idx >= 1) .and. all(idx <= this%idx_bnd))

  ! set data
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = this%data + d
  type is (GRID_DATA1_TYPE)
    this%data(idx(1)) = &
      & this%data(idx(1)) + d
  type is (GRID_DATA2_TYPE)
    this%data(idx(1),idx(2)) = &
      & this%data(idx(1),idx(2)) + d
  type is (GRID_DATA3_TYPE)
    this%data(idx(1),idx(2),idx(3)) = &
      & this%data(idx(1),idx(2),idx(3)) + d
  type is (GRID_DATA4_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4)) = &
      & this%data(idx(1),idx(2),idx(3),idx(4)) + d
  type is (GRID_DATA5_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5)) = &
      & this%data(idx(1),idx(2),idx(3),idx(4),idx(5)) + d
  type is (GRID_DATA6_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6)) = &
      & this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6)) + d
  type is (GRID_DATA7_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7)) = &
      & this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7)) + d
  type is (GRID_DATA8_TYPE)
    this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8)) = &
      & this%data(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6),idx(7),idx(8)) + d
  end select
end subroutine

module subroutine GRID_DATA_UPDATE_ALL(this, d)
  !! update data for all points
  class(GRID_DATA_TYPE), intent(inout) :: this
  TT,                    intent(in)    :: d(:)
    !! delta values

  ASSERT(size(d) == this%n)

  ! set data
  select type (this)
  type is (GRID_DATA0_TYPE)
    this%data = this%data + d(1)
  type is (GRID_DATA1_TYPE)
    this%data = this%data + d
  type is (GRID_DATA2_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2)])
  type is (GRID_DATA3_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3)])
  type is (GRID_DATA4_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4)])
  type is (GRID_DATA5_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5)])
  type is (GRID_DATA6_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6)])
  type is (GRID_DATA7_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6), this%idx_bnd(7)])
  type is (GRID_DATA8_TYPE)
    this%data = this%data + reshape(d, [this%idx_bnd(1), this%idx_bnd(2), this%idx_bnd(3), this%idx_bnd(4), this%idx_bnd(5), this%idx_bnd(6), this%idx_bnd(7), this%idx_bnd(8)])
  end select
end subroutine
#endif

#undef T
#undef TT
#undef TLOG

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
#undef CONCATHELP3
#undef CONCAT3

#undef GRID_DATA_TYPE
#undef GRID_DATA0_TYPE
#undef GRID_DATA1_TYPE
#undef GRID_DATA2_TYPE
#undef GRID_DATA3_TYPE
#undef GRID_DATA4_TYPE
#undef GRID_DATA5_TYPE
#undef GRID_DATA6_TYPE
#undef GRID_DATA7_TYPE
#undef GRID_DATA8_TYPE

#undef GRID_DATA_INIT
#undef GRID_DATA_DESTRUCT
#undef GRID_DATA_RESET
#undef GRID_DATA_GET_PTR0
#undef GRID_DATA_GET_PTR1
#undef GRID_DATA_GET_PTR2
#undef GRID_DATA_GET_PTR3
#undef GRID_DATA_GET_PTR4
#undef GRID_DATA_GET_PTR5
#undef GRID_DATA_GET_PTR6
#undef GRID_DATA_GET_PTR7
#undef GRID_DATA_GET_PTR8
#undef GRID_DATA_GET_POINT
#undef GRID_DATA_GET_ALL
#undef GRID_DATA_SET_POINT
#undef GRID_DATA_SET_ALL
#undef GRID_DATA_UPDATE_POINT
#undef GRID_DATA_UPDATE_ALL

#undef ALLOCATE_GRID_DATA0
#undef ALLOCATE_GRID_DATA1
#undef ALLOCATE_GRID_DATA2
#undef ALLOCATE_GRID_DATA3
#undef ALLOCATE_GRID_DATA4
