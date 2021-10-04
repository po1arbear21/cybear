#include "../util/macro.f90.inc"

module variable_m

  use error_m,         only: assert_failed
  use grid_m,          only: allocate_grid_data, grid, grid_data_real, grid_table, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
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
    procedure :: get_ptr     => variable_get_ptr
    generic   :: get         => variable_get_point, variable_get_all
    generic   :: set         => variable_set_point, variable_set_all
    procedure :: output_data => variable_output_data

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

  subroutine variable_output_data(this, fname, grd_unit, tab)
    !! save variable's data into file as denormalized values.
    class(variable), target,            intent(in) :: this
    character(*),                       intent(in) :: fname
      !! file name, e.g. "output/tmp/pot.csv"
    character(*),             optional, intent(in) :: grd_unit
      !! physical unit token (default: "nm")
    type(grid_table), target, optional, intent(in) :: tab
      !! optional grid table for certain data (default: tab_all)

    integer                   :: iounit, idx(this%g%idx_dim), i, j
    real                      :: coord(this%g%dim)
    type(grid_table), pointer :: tab_
    character(:), allocatable :: grd_unit_

    ! optional arguments
    allocate (character(2) :: grd_unit_)
    grd_unit_ = "nm"
    if (present(grd_unit)) grd_unit_ = grd_unit
    tab_ => this%g%tab_all(this%idx_type, this%idx_dir)
    if (present(tab)) tab_ => tab

    open (newunit = iounit, file = fname, action = 'write')
    do i = 1, tab_%n
      idx = tab_%get_idx(i)

      select case (this%idx_type)
        case(IDX_CELL)
          block
            real :: coord_cl(this%g%dim, this%g%cell_dim)
            call this%g%get_cell(idx, coord_cl)
            coord = sum(coord_cl, dim = 2) / size(coord_cl, dim = 2)
          end block
        case(IDX_EDGE)
          block
            real :: coord_ed(this%g%dim ,2)
            call this%g%get_edge(idx, this%idx_dir, coord_ed)
            coord = sum(coord_ed, dim = 2) / size(coord_ed, dim = 2)
          end block
        case(IDX_FACE)
          block
            real :: coord_fc(this%g%dim, this%g%face_dim(this%idx_dir))
            call this%g%get_face(idx, this%idx_dir, coord_fc)
            coord = sum(coord_fc, dim = 2) / size(coord_fc, dim = 2)
          end block
        case(IDX_VERTEX)
          call this%g%get_vertex(idx, coord)
      end select

      coord = denorm(coord, grd_unit_)

      do j = 1 , size(coord)
        write (iounit, '(ES24.16)', advance = "no") coord(j)
      end do
      write (iounit, '(ES24.16)') denorm(this%get(idx), this%unit)
    end do
    close (unit = iounit)
  end subroutine

end module
