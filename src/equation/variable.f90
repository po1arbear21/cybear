module variable_m
  use grid_m
  use grid_data_m
  use grid_table_m
  implicit none

  type, abstract :: variable
    !! abstract base variable

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

end module