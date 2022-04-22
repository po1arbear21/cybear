m4_include(../util/macro.f90.inc)

module variable_m

  use error_m,         only: assert_failed
  use grid_m,          only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use grid_data_m,     only: allocate_grid_data, grid_data_real, grid_data_cmplx
  use grid_table_m,    only: grid_table
  use grid0D_m,        only: get_dummy_grid
  use iso_fortran_env, only: int64, int32
  use json_m,          only: json_object
  use normalization_m, only: denorm, norm
  use output_file_m,   only: output_file

  implicit none

  private
  public variable, variable_ptr, vector_variable_ptr, variable_real, variable_cmplx

  type, abstract :: variable
    !! base variable

    character(:), allocatable :: name
      !! name of variable
    character(:), allocatable :: unit
      !! physical unit token

    class(grid), pointer :: g => null()
      !! pointer to grid

    type(grid_table) :: tab_all
      !! table that selects all possible grid indices

    integer :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer :: idx_dir
      !! index direction for edges and faces
  contains
    procedure :: variable_init
    procedure :: reset       => variable_reset
    procedure :: get_ptr     => variable_get_ptr
    procedure :: hashkey     => variable_hashkey
    procedure :: output_info => variable_output_info
    procedure :: output_data => variable_output_data
  end type

  type variable_ptr
    class(variable), pointer :: p => null()
  end type

  m4_define({T},{variable_ptr})
  m4_include(../util/vector_def.f90.inc)

  type, abstract, extends(variable) :: variable_real
    !! real valued variable

    class(grid_data_real), allocatable :: data
      !! real values
  contains
    generic :: get => variable_real_get_point, variable_real_get_all
    generic :: set => variable_real_set_point, variable_real_set_all

    procedure, private :: variable_real_get_point, variable_real_get_all
    procedure, private :: variable_real_set_point, variable_real_set_all
  end type

  type, abstract, extends(variable) :: variable_cmplx
    !! complex valued variable

    class(grid_data_cmplx), allocatable :: data
      !! complex values
  contains
    generic :: get => variable_cmplx_get_point, variable_cmplx_get_all
    generic :: set => variable_cmplx_set_point, variable_cmplx_set_all

    procedure, private :: variable_cmplx_get_point, variable_cmplx_get_all
    procedure, private :: variable_cmplx_set_point, variable_cmplx_set_all
  end type

contains

  m4_define({T},{variable_ptr})
  m4_include(../util/vector_imp.f90.inc)

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

    character(32) :: tab_name

    this%name = name
    this%unit = unit

    ! optional arguments
    if (present(g)) then
      m4_assert(present(idx_type))
      m4_assert(present(idx_dir))
      this%g        => g
      this%idx_type =  idx_type
      this%idx_dir  =  idx_dir
    else
      m4_assert(.not. present(idx_type))
      m4_assert(.not. present(idx_dir))
      this%g        => get_dummy_grid()
      this%idx_type =  IDX_VERTEX
      this%idx_dir  =  0
    end if

    ! allocate data
    select type (this)
    class is (variable_real)
      call allocate_grid_data(this%data, this%g%idx_dim)
      call this%data%init(this%g, this%idx_type, this%idx_dir)
    class is (variable_cmplx)
      call allocate_grid_data(this%data, this%g%idx_dim)
      call this%data%init(this%g, this%idx_type, this%idx_dir)
    end select

    ! init tab_all
    if (this%idx_dir > 0) then
      write (tab_name, "(A,A,I0)") "All", trim(IDX_NAME(this%idx_type)), this%idx_dir
    else
      write (tab_name, "(A,A)") "All", trim(IDX_NAME(this%idx_type))
    end if
    call this%tab_all%init(trim(tab_name), this%g, this%idx_type, this%idx_dir, initial_flags = .true.)
    call this%tab_all%init_final()
  end subroutine

  subroutine variable_reset(this)
    !! reset variable data
    class(variable), intent(inout) :: this

    select type (this)
    class is (variable_real)
      call this%data%reset()
    class is (variable_cmplx)
      call this%data%reset()
    end select
  end subroutine

  function variable_get_ptr(this) result(ptr)
    !! returns pointer type to this variable
    class(variable), target, intent(in) :: this
    type(variable_ptr)                  :: ptr

    ptr%p => this
  end function

  function variable_hashkey(this) result(hkey)
    !! get key for hashmap (based on address)
    class(variable), target, intent(in) :: this
    integer                             :: hkey(m4_ptrsize)

    integer(int64)           :: iptr
    class(variable), pointer :: vptr

    vptr => this
    iptr =  loc(vptr)

    m4_ifelse(m4_intsize,32,{
      hkey(1) = int(iptr, kind = int32)
      hkey(2) = int(ishft(iptr, -32), kind = int32)
    },{
      hkey(1) = iptr
    })
  end function

  subroutine variable_output_info(this, of)
    !! output variable info
    class(variable),   intent(in)    :: this
    type(output_file), intent(inout) :: of
      !! output file handle

    type(json_object), pointer :: obj

    obj => of%new_object("Variables")
    call obj%add_string("Name", this%name)
    select type (this)
    class is (variable_real)
      call obj%add_string("Type", "Real")
    class is (variable_cmplx)
      call obj%add_string("Type", "Complex")
    end select
    call obj%add_string("Grid", this%g%name)
    call obj%add_string("IdxType", trim(IDX_NAME(this%idx_type)))
    call obj%add_int("IdxDir", this%idx_dir)
  end subroutine

  subroutine variable_output_data(this, of, obj)
    !! output variable data
    class(variable),            intent(in)    :: this
    type(output_file),          intent(inout) :: of
      !! output file handle
    type(json_object), pointer, intent(inout) :: obj
      !! parent object in output file

    select type (this)
    class is (variable_real)
      call this%data%output(of, obj, this%name, this%unit)
    class is (variable_cmplx)
      call this%data%output(of, obj, this%name, this%unit)
    end select
  end subroutine

  function variable_real_get_point(this, idx) result(d)
    !! get data for single point with bounds check (out of bounds: return default value)
    class(variable_real), intent(in) :: this
    integer,              intent(in) :: idx(:)
      !! grid indices
    real                             :: d
      !! return data

    d = this%data%get(idx)
  end function

  function variable_real_get_all(this) result(d)
    !! get data for all points in flat array
    class(variable_real), intent(in) :: this
    real                             :: d(this%data%n)
      !! return all data

    d = this%data%get()
  end function

  subroutine variable_real_set_point(this, idx, d)
    !! set data for single point with bounds check (do nothing if out of bounds)
    class(variable_real), intent(inout) :: this
    integer,              intent(in)    :: idx(:)
      !! grid indices (idx_dim)
    real,                 intent(in)    :: d
      !! new value

    call this%data%set(idx, d)
  end subroutine

  subroutine variable_real_set_all(this, d)
    !! set data for all points
    class(variable_real), intent(inout) :: this
    real,                 intent(in)    :: d(:)
      !! new values

    call this%data%set(d)
  end subroutine

  function variable_cmplx_get_point(this, idx) result(d)
    !! get data for single point with bounds check (out of bounds: return default value)
    class(variable_cmplx), intent(in) :: this
    integer,               intent(in) :: idx(:)
      !! grid indices
    complex                           :: d
      !! return data

    d = this%data%get(idx)
  end function

  function variable_cmplx_get_all(this) result(d)
    !! get data for all points in flat array
    class(variable_cmplx), intent(in) :: this
    complex                           :: d(this%data%n)
      !! return all data

    d = this%data%get()
  end function

  subroutine variable_cmplx_set_point(this, idx, d)
    !! set data for single point with bounds check (do nothing if out of bounds)
    class(variable_cmplx), intent(inout) :: this
    integer,               intent(in)    :: idx(:)
      !! grid indices (idx_dim)
    complex,               intent(in)    :: d
      !! new value

    call this%data%set(idx, d)
  end subroutine

  subroutine variable_cmplx_set_all(this, d)
    !! set data for all points
    class(variable_cmplx), intent(inout) :: this
    complex,               intent(in)    :: d(:)
      !! new values

    call this%data%set(d)
  end subroutine

end module
