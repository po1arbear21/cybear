#include "../util/macro.f90.inc"

module variable_m

  use error_m,         only: assert_failed
  use grid_m,          only: allocate_grid_data, grid, grid_data_real, IDX_VERTEX
  use grid0D_m,        only: get_dummy_grid
  use normalization_m, only: denorm, norm

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
    procedure :: get_ptr   => variable_get_ptr
    generic   :: get       => variable_get_point, variable_get_all
    generic   :: set       => variable_set_point, variable_set_all
    procedure :: load_data => variable_load_data
    procedure :: save_data => variable_save_data

    procedure, private :: variable_get_point, variable_get_all
    procedure, private :: variable_set_point, variable_set_all
  end type

  type variable_ptr
    class(variable), pointer :: p => null()
  end type

contains

  subroutine variable_init(this, name, unit, g, idx_type, idx_dir)
    !! initialize base variable
    class(variable),               intent(out) :: this
    character(*),                  intent(in)  :: name
      !! variable name
    character(*),                  intent(in)  :: unit
      !! physical unit token
    class(grid), target, optional, intent(in)  :: g
      !! grid this variable is defined on
      !! default: type(grid0D) dummy_grid
    integer,             optional, intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,             optional, intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)

    this%name = name
    this%unit = unit

    ! optional arguments
    if (present(g)) then
      ASSERT(present(idx_type))
      ASSERT(present(idx_dir))
      this%g        => g
      this%idx_type =  idx_type
      this%idx_dir  =  idx_dir
    else
      ASSERT(.not. present(idx_type))
      ASSERT(.not. present(idx_dir))
      this%g        => get_dummy_grid()
      this%idx_type =  IDX_VERTEX
      this%idx_dir  =  0
    end if

    ! allocate data
    call allocate_grid_data(this%data, this%g%idx_dim)
    call this%data%init(this%g, this%idx_type, this%idx_dir)
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

  subroutine variable_save_data(this, fname)
    !! save variable's data into file as denormalized values.
    class(variable), intent(in) :: this
    character(*),    intent(in) :: fname
      !! file name, e.g. "output/tmp/pot.csv"

    integer           :: iounit, i
    real, allocatable :: d(:)

    d = denorm(this%get(), this%unit)

    open (newunit = iounit, file = fname, action = 'WRITE')
    do i = 1, size(d)
      write (iounit, '(E32.24)') d(i)
    end do
    close (unit = iounit)
  end subroutine

  subroutine variable_load_data(this, fname)
    !! load variable's data from file.
    !! assume denormalized data in file and perform normalization.
    class(variable), intent(inout) :: this
    character(*),    intent(in)    :: fname
      !! file name, e.g. "output/tmp/pot.csv"

    integer           :: iounit, n, ios
    real              :: tmp
    real, allocatable :: d(:)

    open (newunit = iounit, file = fname, action = 'READ')

    ! count lines in file
    n   = -1
    ios = 0
    do while (ios == 0)
      n = n + 1
      read (iounit, *, iostat = ios) tmp
    end do

    allocate (d(n))
    rewind (iounit)
    read (iounit, *) d
    close (unit = iounit)

    call this%set(norm(d, this%unit))
  end subroutine

end module
