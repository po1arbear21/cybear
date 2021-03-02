module variable_m

  use grid_m,      only: grid
  use grid_data_m, only: grid_data_real, allocate_grid_data

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
    procedure :: get     => variable_get
    procedure :: set     => variable_set
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

  function variable_get(this, idx) result(d)
    !! returns data for given grid indices
    class(variable), intent(in) :: this
    integer,         intent(in) :: idx(:)
      !! grid indices
    real                        :: d
      !! value

    d = this%data%get(idx)
  end function

  subroutine variable_set(this, idx, d)
    !! set data for given grid indices
    class(variable), intent(inout) :: this
    integer,         intent(in)    :: idx(:)
      !! grid indices
    real,            intent(in)    :: d
      !! value

    call this%data%set(idx, d)
  end subroutine

end module
