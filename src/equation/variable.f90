module variable_m

  use grid_m, only: grid, grid_data_real, allocate_grid_data

  implicit none

  private
  public variable, variable_ptr

  type, abstract :: variable
    !! base variable
    !! derive for specific variable, e.g. potential...

    character(:), allocatable :: name
      !! name of variable
    character(:), allocatable :: unit
      !! physical unit token

    class(grid), pointer :: g => null()
      !! pointer to grid

    integer :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer :: idx_dir
      !! index direction for edges and faces

    class(grid_data_real), allocatable :: data
      !! values
  contains
    procedure :: variable_init
    procedure :: get_ptr => variable_get_ptr

    procedure, private :: variable_get_point
    procedure, private :: variable_get_all
    generic            :: get => variable_get_point, variable_get_all

    procedure, private :: variable_set_point
    procedure, private :: variable_set_all
    generic            :: set => variable_set_point, variable_set_all
  end type

  type variable_ptr
    class(variable), pointer :: p => null()
  end type

contains

  subroutine variable_init(this, name, unit, g, idx_type, idx_dir)
    !! initialize base variable
    class(variable),     intent(out) :: this
    character(*),        intent(in)  :: name
      !! variable name
    character(*),        intent(in)  :: unit
      !! physical unit token
    class(grid), target, intent(in)  :: g
      !! grid this variable is defined on
    integer,             intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,             intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)

    this%name     =  name
    this%unit     =  unit
    this%g        => g
    this%idx_type =  idx_type
    this%idx_dir  =  idx_dir

    ! allocate data
    call allocate_grid_data(this%data, g%idx_dim)
    call this%data%init(g, idx_type, idx_dir)
  end subroutine

  function variable_get_ptr(this) result(ptr)
    !! returns pointer type to this variable
    class(variable), target, intent(in) :: this
    type(variable_ptr)                  :: ptr

    ptr%p => this
  end function

  function variable_get_point(this, idx) result(d)
    !! get data for single point with bounds check (out of bounds: return default value)
    class(variable), intent(in) :: this
    integer,         intent(in) :: idx(:)
      !! grid indices
    real                        :: d
      !! return data

    d = this%data%get(idx)
  end function

  function variable_get_all(this) result(d)
    !! get data for all points in flat array
    class(variable), intent(in) :: this
    real                        :: d(this%data%n)
      !! return all data

    d = this%data%get()
  end function

  subroutine variable_set_point(this, idx, d)
    !! set data for single point with bounds check (do nothing if out of bounds)
    class(variable), intent(inout) :: this
    integer,         intent(in)    :: idx(:)
      !! grid indices (idx_dim)
    real,            intent(in)    :: d
      !! new value

    call this%data%set(idx, d)
  end subroutine

  subroutine variable_set_all(this, d)
    !! set data for all points
    class(variable), intent(inout) :: this
    real,            intent(in)    :: d(:)
      !! new values

    call this%data%set(d)
  end subroutine

end module
