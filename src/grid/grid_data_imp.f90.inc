#define PASTE(X)             X
#define PASTE2(X)            PASTE(X)_
#define CONCATHELP(X, Y)     PASTE2(X)Y
#define CONCAT(X, Y)         CONCATHELP(X,Y)
#define CONCATHELP3(X, Y, Z) PASTE2(CONCATHELP(X,Y))Z
#define CONCAT3(X, Y, Z)     CONCATHELP3(X,Y,Z)

#define GRID_DATA_TYPE       CONCAT(grid_data,T)
#define GRID_DATA_INIT       CONCAT3(grid_data,T,init)
#define GRID_DATA_DESTRUCT   CONCAT3(grid_data,T,destruct)
#define GRID_DATA_GET        CONCAT3(grid_data,T,get)
#define GRID_DATA_SET        CONCAT3(grid_data,T,set)

subroutine GRID_DATA_INIT(this, g, idx_type, dir)
  !! initialize and allocate grid data
  class(GRID_DATA_TYPE), intent(out) :: this
  class(grid),           intent(in)  :: g
  integer,               intent(in)  :: idx_type
    !! grid index type (e.g. IDX_VERTEX)
  integer,               intent(in)  :: dir
    !! index of direction (only used for IDX_EDGE and IDX_FACE; range = 1:idx_dim)

  TT :: d0

  ! get index bounds
  allocate (this%idx_bnd(g%idx_dim))
  call g%get_idx_bnd(idx_type, dir, this%idx_bnd)

#ifdef TLOG
  d0 = .false.
#else
  d0 = 0
#endif

  ! allocate data array
  allocate (this%data(product(this%idx_bnd)), source = d0)
end subroutine

subroutine GRID_DATA_DESTRUCT(this)
  !! destruct/deallocate grid data
  class(GRID_DATA_TYPE), intent(inout) :: this

  if (allocated(this%idx_bnd)) deallocate (this%idx_bnd)
  if (allocated(this%data   )) deallocate (this%data   )
end subroutine

function GRID_DATA_GET(this, idx) result(d)
  !! get data with bounds check (out of bounds: return default value)
  class(GRID_DATA_TYPE), intent(in) :: this
  integer,               intent(in) :: idx(:)
    !! grid indices (idx_dim)
  TT                                :: d
    !! return data

  integer :: i, j

  ! default value
#ifdef TLOG
  d = .false.
#else
  d = 0
#endif

  ! bounds check
  if (any(idx < 1) .or. any(idx > this%idx_bnd)) return

  ! get data
  i = idx(1)
  do j = 2, size(idx)
    i = i + idx(j) * this%idx_bnd(j-1)
  end do
  d = this%data(i)
end function

subroutine GRID_DATA_SET(this, idx, d)
  !! set data with bounds check (do nothing if out of bounds)
  class(GRID_DATA_TYPE), intent(inout) :: this
  integer,               intent(in)    :: idx(:)
    !! grid indices (idx_dim)
  TT,                    intent(in)    :: d
    !! new value

  integer :: i, j

  ! bounds check
  if (any(idx < 1) .or. any(idx > this%idx_bnd)) return

  ! set data
  i = idx(1)
  do j = 2, size(idx)
    i = i + idx(j) * this%idx_bnd(j-1)
  end do
  this%data(i) = d
end subroutine

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
#undef GRID_DATA_INIT
#undef GRID_DATA_DESTRUCT
#undef GRID_DATA_GET
#undef GRID_DATA_SET