#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define ALLOCATE_GRID_DATA     CONCAT(allocate_grid_data,T)

#define GRID_DATA_TYPE         CONCAT(grid_data,T)
#define GRID_DATA_INIT         CONCAT3(grid_data,T,init)
#define GRID_DATA_DESTRUCT     CONCAT3(grid_data,T,destruct)
#define GRID_DATA_GET_POINT    CONCAT3(grid_data,T,get_point)
#define GRID_DATA_GET_ALL      CONCAT3(grid_data,T,get_all)
#define GRID_DATA_SET_POINT    CONCAT3(grid_data,T,set_point)
#define GRID_DATA_SET_ALL      CONCAT3(grid_data,T,set_all)
#define GRID_DATA_UPDATE_POINT CONCAT3(grid_data,T,update_point)
#define GRID_DATA_UPDATE_ALL   CONCAT3(grid_data,T,update_all)

#define GRID_DATA1_TYPE        CONCAT(grid_data1,T)
#define GRID_DATA2_TYPE        CONCAT(grid_data2,T)
#define GRID_DATA3_TYPE        CONCAT(grid_data3,T)
#define GRID_DATA4_TYPE        CONCAT(grid_data4,T)
#define GRID_DATA5_TYPE        CONCAT(grid_data5,T)
#define GRID_DATA6_TYPE        CONCAT(grid_data6,T)
#define GRID_DATA7_TYPE        CONCAT(grid_data7,T)
#define GRID_DATA8_TYPE        CONCAT(grid_data8,T)

subroutine ALLOCATE_GRID_DATA(gdata, dim)
  !! allocate grid data
  class(GRID_DATA_TYPE), allocatable, intent(out) :: gdata
  integer,                            intent(in)  :: dim
    !! grid data dimension

  select case (dim)
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
      call program_error("grid data dimension must be in range 1:8")
  end select
end subroutine

subroutine GRID_DATA_INIT(this, g, idx_type, idx_dir)
  !! initialize and allocate grid data
  class(GRID_DATA_TYPE), intent(out) :: this
  class(grid),           intent(in)  :: g
  integer,               intent(in)  :: idx_type
    !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
  integer,               intent(in)  :: idx_dir
    !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)

  TT :: d0

  ! get index bounds
  allocate (this%idx_bnd(g%idx_dim))
  call g%get_idx_bnd(idx_type, idx_dir, this%idx_bnd)

  ! set total number of values
  this%n = product(this%idx_bnd)

  ! get initial value
#ifdef TLOG
  d0 = .false.
#else
  d0 = 0
#endif

  ! allocate data array
  select type (this)
    type is (GRID_DATA1_TYPE)
      allocate (this%data(this%idx_bnd(1)), source = d0)
    type is (GRID_DATA2_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2)), source = d0)
    type is (GRID_DATA3_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3)), source = d0)
    type is (GRID_DATA4_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3), &
        &                 this%idx_bnd(4)), source = d0)
    type is (GRID_DATA5_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3), &
        &                 this%idx_bnd(4), &
        &                 this%idx_bnd(5)), source = d0)
    type is (GRID_DATA6_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3), &
        &                 this%idx_bnd(4), &
        &                 this%idx_bnd(5), &
        &                 this%idx_bnd(6)), source = d0)
    type is (GRID_DATA7_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3), &
        &                 this%idx_bnd(4), &
        &                 this%idx_bnd(5), &
        &                 this%idx_bnd(6), &
        &                 this%idx_bnd(7)), source = d0)
    type is (GRID_DATA8_TYPE)
      allocate (this%data(this%idx_bnd(1), &
        &                 this%idx_bnd(2), &
        &                 this%idx_bnd(3), &
        &                 this%idx_bnd(4), &
        &                 this%idx_bnd(5), &
        &                 this%idx_bnd(6), &
        &                 this%idx_bnd(7), &
        &                 this%idx_bnd(8)), source = d0)
  end select
end subroutine

subroutine GRID_DATA_DESTRUCT(this)
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

function GRID_DATA_GET_POINT(this, idx) result(d)
  !! get data for single point with bounds check (out of bounds: return default value)
  class(GRID_DATA_TYPE), intent(in) :: this
  integer,               intent(in) :: idx(:)
    !! grid indices (idx_dim)
  TT                                :: d
    !! return data

  ! default value
#ifdef TLOG
  d = .false.
#else
  d = 0
#endif

  ! bounds check
  if (any(idx < 1) .or. any(idx > this%idx_bnd)) return

  ! get data
  select type (this)
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

function GRID_DATA_GET_ALL(this) result(d)
  !! get data for all points in flat array
  class(GRID_DATA_TYPE), intent(in) :: this
  TT                                :: d(this%n)
    !! return all data

  select type (this)
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

subroutine GRID_DATA_SET_POINT(this, idx, d)
  !! set data for single point with bounds check (do nothing if out of bounds)
  class(GRID_DATA_TYPE), intent(inout) :: this
  integer,               intent(in)    :: idx(:)
    !! grid indices (idx_dim)
  TT,                    intent(in)    :: d
    !! new value

  ! bounds check
  if (any(idx < 1) .or. any(idx > this%idx_bnd)) return

  ! set data
  select type (this)
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

subroutine GRID_DATA_SET_ALL(this, d)
  !! set data for all points
  class(GRID_DATA_TYPE), intent(inout) :: this
  TT,                    intent(in)    :: d(this%n)
    !! new values

  ! set data
  select type (this)
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
subroutine GRID_DATA_UPDATE_POINT(this, idx, d)
  !! update data for single point with bounds check (do nothing if out of bounds)
  class(GRID_DATA_TYPE), intent(inout) :: this
  integer,               intent(in)    :: idx(:)
    !! grid indices (idx_dim)
  TT,                    intent(in)    :: d
    !! delta value

  ! bounds check
  if (any(idx < 1) .or. any(idx > this%idx_bnd)) return

  ! set data
  select type (this)
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

subroutine GRID_DATA_UPDATE_ALL(this, d)
  !! update data for all points
  class(GRID_DATA_TYPE), intent(inout) :: this
  TT,                    intent(in)    :: d(this%n)
    !! delta values

  ! set data
  select type (this)
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

#undef ALLOCATE_GRID_DATA

#undef GRID_DATA_TYPE
#undef GRID_DATA_INIT
#undef GRID_DATA_DESTRUCT
#undef GRID_DATA_GET_POINT
#undef GRID_DATA_GET_ALL
#undef GRID_DATA_SET_POINT
#undef GRID_DATA_SET_ALL
#undef GRID_DATA_UPDATE_POINT
#undef GRID_DATA_UPDATE_ALL

#undef GRID_DATA1_TYPE
#undef GRID_DATA2_TYPE
#undef GRID_DATA3_TYPE
#undef GRID_DATA4_TYPE
#undef GRID_DATA5_TYPE
#undef GRID_DATA6_TYPE
#undef GRID_DATA7_TYPE
#undef GRID_DATA8_TYPE
