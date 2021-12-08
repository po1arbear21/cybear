#include "../util/macro.f90.inc"

module vselector_m

  use error_m,         only: assert_failed, program_error
  use grid_m,          only: grid, grid_data_int, grid_table, grid_table_ptr, allocate_grid_data
  use iso_fortran_env, only: int64, int32
  use util_m,          only: hash
  use variable_m,      only: variable, variable_ptr, variable_real, variable_cmplx

  implicit none

  private
  public vselector, vselector_ptr, vector_vselector_ptr

  type vselector
    !! variable selector: select data from one or more variables using one or more grid tables

    character(:), allocatable :: name
      !! variable selector name

    class(grid), pointer :: g => null()
      !! pointer to grid

    integer :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer :: idx_dir
      !! index direction for edges and faces

    type(variable_ptr),   allocatable :: v(:)
      !! pointers to selected variables
    type(grid_table_ptr), allocatable :: tab(:)
      !! pointers to selected grid tables

    integer              :: nvar
      !! number of variables: nvar = size(v)
    integer              :: nval
      !! number of values per grid point
    integer              :: ntab
      !! number of tables/blocks: ntab = size(tab)
    integer, allocatable :: nvals(:)
      !! number of values per block: nvals(itab) = tab(itab)%n * nval
    integer              :: n
      !! total number of values selected: n = sum(nvals)

    class(grid_data_int), allocatable :: itab
      !! get table index from grid indices
  contains
    generic   :: init         => vselector_init_nvar_ntab, &
      &                          vselector_init_var_ntab,  &
      &                          vselector_init_nvar_tab,  &
      &                          vselector_init_var_tab
    procedure :: get_ptr      => vselector_get_ptr
    procedure :: reset        => vselector_reset
    procedure :: compare      => vselector_compare
    procedure :: hashkey_size => vselector_hashkey_size
    procedure :: hashkey      => vselector_hashkey
    generic   :: get          => vselector_get_single,    vselector_get_block,    vselector_get_all
    generic   :: set          => vselector_set_single,    vselector_set_block,    vselector_set_all
    generic   :: update       => vselector_update_single, vselector_update_block, vselector_update_all
    procedure :: print        => vselector_print

    procedure, private :: vselector_init_nvar_ntab, vselector_init_var_ntab, vselector_init_nvar_tab, vselector_init_var_tab
    procedure, private :: vselector_get_single,    vselector_get_block,    vselector_get_all
    procedure, private :: vselector_set_single,    vselector_set_block,    vselector_set_all
    procedure, private :: vselector_update_single, vselector_update_block, vselector_update_all
  end type

  type vselector_ptr
    type(vselector), pointer :: p => null()
  end type

#define T vselector_ptr
#define TT type(vselector_ptr)
#include "../util/vector_def.f90.inc"

contains

#define T vselector_ptr
#define TT type(vselector_ptr)
#include "../util/vector_imp.f90.inc"

  subroutine vselector_init_nvar_ntab(this, v, tab, name)
    !! initialize variable selector given multiple variables and multiple tables.
    class(vselector),     intent(out) :: this
    type(variable_ptr),   intent(in)  :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)  :: tab(:)
      !! grid table pointers
    character(*),         intent(in)  :: name
      !! selector name

    integer :: i, itab, idx(v(1)%p%g%idx_dim)

    this%name     =  name
    this%g        => v(1)%p%g
    this%idx_type =  v(1)%p%idx_type
    this%idx_dir  =  v(1)%p%idx_dir
    this%v        =  v
    this%tab      =  tab

#ifdef DEBUG
    do i = 2, size(v)
      if (.not. associated(v(i)%p%g, this%g)) call program_error("variables defined on different grids")
      if (v(i)%p%idx_type /= this%idx_type  ) call program_error("different variable index types"      )
      if (v(i)%p%idx_dir  /= this%idx_dir   ) call program_error("different variable index directions" )
    end do
    do i = 1, size(tab)
      if (.not. associated(tab(i)%p%g, this%g)) call program_error("grid table defined on wrong grid"     )
      if (tab(i)%p%idx_type /= this%idx_type  ) call program_error("grid table has wrong index type"      )
      if (tab(i)%p%idx_dir  /= this%idx_dir   ) call program_error("grid table has wrong index direction" )
    end do
#endif

    ! number of variables
    this%nvar = size(v)

    ! number of values per grid point
    this%nval = 0
    do i = 1, this%nvar
      select type (v => v(i)%p)
      class is (variable_real)
        this%nval = this%nval + 1
      class is (variable_cmplx)
        this%nval = this%nval + 2 ! real + imag part
      end select
    end do

    ! number of tables
    this%ntab = size(tab)

    ! number of values per table
    allocate (this%nvals(size(tab)), source = 0)
    do itab = 1, this%ntab
      this%nvals(itab) = tab(itab)%p%n * this%nval
    end do

    ! total number of values
    this%n = sum(this%nvals)

    ! initialize table index data
    call allocate_grid_data(this%itab, this%g%idx_dim)
    call this%itab%init(this%g, this%idx_type, this%idx_dir)
    do itab = 1, this%ntab
      do i = 1, tab(itab)%p%n
        idx = tab(itab)%p%get_idx(i)
        ASSERT(this%itab%get(idx) == 0)
        call this%itab%set(idx, itab)
      end do
    end do
  end subroutine

  subroutine vselector_init_var_ntab(this, v, tab, name)
    !! initialize variable selector given one variable and multiple tables.
    class(vselector),       intent(out) :: this
    class(variable),        intent(in)  :: v
      !! variable
    type(grid_table_ptr),   intent(in)  :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)  :: name
      !! selector name (default: variable%name)

    character(:), allocatable :: name_

    ! optional name parsing
    if (present(name)) then
      allocate (name_, source = name)
    else
      allocate (name_, source = v%name)
    end if

    call this%init([v%get_ptr()], tab, name_)
  end subroutine

  subroutine vselector_init_nvar_tab(this, v, name, tab)
    !! initialize variable selector given multiple variables and one table.
    class(vselector),           intent(out) :: this
    type(variable_ptr),         intent(in)  :: v(:)
      !! variable pointers
    character(*),               intent(in)  :: name
      !! selector name
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variables' whole grids via v(1)%g%tab_all)

    if (present(tab)) then
      call this%init(v, [tab%get_ptr()], name)
    else
      call this%init(v, [v(1)%p%g%tab_all(v(1)%p%idx_type,v(1)%p%idx_dir)%get_ptr()], name)
    end if
  end subroutine

  subroutine vselector_init_var_tab(this, v, tab, name)
    !! initialize variable selector given one variable and one table.
    class(vselector),           intent(out) :: this
    class(variable),            intent(in)  :: v
      !! variable
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    character(*),     optional, intent(in)  :: name
      !! selector name (default: variable%name)

    character(:), allocatable :: name_

    ! optional name parsing
    if (present(name)) then
      allocate (name_, source = name)
    else
      allocate (name_, source = v%name)
    end if

    if (present(tab)) then
      call this%init([v%get_ptr()], [tab%get_ptr()], name_)
    else
      call this%init([v%get_ptr()], [v%g%tab_all(v%idx_type,v%idx_dir)%get_ptr()], name_)
    end if
  end subroutine

  function vselector_get_ptr(this) result(ptr)
    !! returns pointer type to this vselector
    class(vselector), target, intent(in) :: this
    type(vselector_ptr)                  :: ptr

    ptr%p => this
  end function

  subroutine vselector_reset(this, itab)
    !! reset all data selected to zero
    class(vselector),  intent(inout) :: this
    integer, optional, intent(in)    :: itab
      !! select table to reset (default: reset all tables)

    integer :: i, j, itab0, itab1, itab_, idx(this%g%idx_dim)

    ! optional table index
    itab0 = 1
    itab1 = this%ntab
    if (present(itab)) then
      itab0 = itab
      itab1 = itab
    end if

    do itab_ = itab0, itab1
      associate (tab => this%tab(itab_)%p)
        do i = 1, tab%n
          idx = tab%get_idx(i)
          do j = 1, this%nvar
            select type (v => this%v(j)%p)
            class is (variable_real)
              call v%data%set(idx, 0.0)
            class is (variable_cmplx)
              call v%data%set(idx, (0.0, 0.0))
            end select
          end do
        end do
      end associate
    end do
  end subroutine

  function vselector_compare(this, other) result(c)
    !! compare two variable selectors (names are not checked)
    class(vselector), intent(in) :: this
    type(vselector),  intent(in) :: other
    logical                      :: c
      !! return this == other

    integer :: i

    c = .false.

    if (.not. associated(this%g, other%g)) return
    if (this%idx_type /= other%idx_type) return
    if (this%idx_dir /= other%idx_dir) return
    if (this%nvar /= other%nvar) return
    if (this%nval /= other%nval) return
    do i = 1, this%nvar
      if (.not. associated(this%v(i)%p, other%v(i)%p)) return
    end do
    if (this%ntab /= other%ntab) return
    do i = 1, this%ntab
      if (.not. associated(this%tab(i)%p, other%tab(i)%p)) return
    end do

    c = .true.
  end function

  pure function vselector_hashkey_size(this) result(n)
    !! get key size
    class(vselector), intent(in) :: this
    integer                      :: n

#ifdef INTSIZE32
    n = 2 + 2*(1 + size(this%v) + size(this%tab))
#endif
#ifdef INTSIZE64
    n = 3 + size(this%v) + size(this%tab)
#endif
  end function

  function vselector_hashkey(this) result(hkey)
    !! get key for hashmap (name is not included)
    class(vselector), intent(in) :: this
    integer                      :: hkey(this%hashkey_size())
      !! return hashmap key

    integer                    :: i, n
    integer(int64)             :: iptr
    class(grid),       pointer :: gptr ! use local pointers to avoid ifort bug
    class(variable),   pointer :: vptr ! otherwise loc does not return address of target
    class(grid_table), pointer :: tptr

    ! convert to integer array (convert pointers to integers using loc function)
    hkey(1) = this%idx_type
    hkey(2) = this%idx_dir
    gptr => this%g
    iptr =  loc(gptr)
#ifdef INTSIZE32
    hkey(3) = int(iptr, kind = int32)
    hkey(4) = int(ishft(iptr, -32), kind = int32)
    n = 4
#endif
#ifdef INTSIZE64
    hkey(3) = iptr
    n = 3
#endif
    do i = 1, size(this%v)
      vptr => this%v(i)%p
      iptr =  loc(vptr)
#ifdef INTSIZE32
      hkey(n+1) = int(iptr, kind = int32)
      hkey(n+2) = int(ishft(iptr, -32), kind = int32)
      n = n + 2
#endif
#ifdef INTSIZE64
      hkey(n+1) = iptr
      n = n + 1
#endif
    end do
    do i = 1, size(this%tab)
      tptr => this%tab(i)%p
      iptr =  loc(tptr)
#ifdef INTSIZE32
      hkey(n+1) = int(iptr, kind = int32)
      hkey(n+2) = int(ishft(iptr, -32), kind = int32)
      n = n + 2
#endif
#ifdef INTSIZE64
      hkey(n+1) = iptr
      n = n + 1
#endif
    end do
  end function

  function vselector_get_single(this, idx) result(x)
    !! get data for single grid index tuple
    class(vselector), intent(in) :: this
    integer,          intent(in) :: idx(:)
      !! grid indices
    real                         :: x(this%nval)
      !! return variable data at idx

    integer :: i, j
    complex :: c

    j = 1
    do i = 1, this%nvar
      select type (v => this%v(i)%p)
      class is (variable_real)
        x(j) = v%get(idx)
        j = j + 1
      class is (variable_cmplx)
        c      = v%get(idx)
        x(j  ) = real(c)
        x(j+1) = aimag(c)
        j = j + 2
      end select
    end do
  end function

  function vselector_get_block(this, itab) result(x)
    !! get data for whole block
    class(vselector), intent(in) :: this
    integer,          intent(in) :: itab
      !! grid table index
    real                         :: x(this%nvals(itab))

    integer :: i, i0, i1, idx(this%g%idx_dim)

    associate (tab => this%tab(itab)%p)
      i1 = 0
      do i = 1, tab%n
        i0 = i1 + 1
        i1 = i1 + this%nval

        idx = tab%get_idx(i)
        x(i0:i1) = this%vselector_get_single(idx)
      end do
    end associate
  end function

  function vselector_get_all(this) result(x)
    !! get data for all blocks
    class(vselector), intent(in) :: this
    real                         :: x(this%n)

    integer :: itab, i0, i1

    i1 = 0
    do itab = 1, this%ntab
      i0 = i1 + 1
      i1 = i1 + this%nvals(itab)

      x(i0:i1) = this%vselector_get_block(itab)
    end do
  end function

  subroutine vselector_set_single(this, idx, x)
    !! set data for single grid index tuple
    class(vselector), intent(inout) :: this
    integer,          intent(in)    :: idx(:)
      !! grid indices
    real,             intent(in)    :: x(:)
      !! new variable data (nval)

    integer :: i, j

    ASSERT(size(x) == this%nval)

    j = 1
    do i = 1, this%nvar
      select type (v => this%v(i)%p)
      class is (variable_real)
        call v%set(idx, x(j))
        j = j + 1
      class is (variable_cmplx)
        call v%set(idx, complex(x(j), x(j+1)))
        j = j + 2
      end select
    end do
  end subroutine

  subroutine vselector_set_block(this, itab, x)
    !! set data for whole block
    class(vselector), intent(inout) :: this
    integer,          intent(in)    :: itab
      !! table index
    real,             intent(in)    :: x(:)
      !! new variable data (nvals(itab))

    integer :: i, i0, i1, idx(this%g%idx_dim)

    ASSERT(size(x) == this%nvals(itab))

    associate (tab => this%tab(itab)%p)
      i1 = 0
      do i = 1, tab%n
        i0 = i1 + 1
        i1 = i1 + this%nval

        idx = tab%get_idx(i)
        call this%vselector_set_single(idx, x(i0:i1))
      end do
    end associate
  end subroutine

  subroutine vselector_set_all(this, x)
    !! set data for all blocks
    class(vselector), intent(inout) :: this
    real,             intent(in)    :: x(:)
      !! new variable data (n)

    integer :: itab, i0, i1

    ASSERT(size(x) == this%n)

    i1 = 0
    do itab = 1, this%ntab
      i0 = i1 + 1
      i1 = i1 + this%nvals(itab)

      call this%vselector_set_block(itab, x(i0:i1))
    end do
  end subroutine

  subroutine vselector_update_single(this, idx, dx)
    !! update data for single grid index tuple
    class(vselector), intent(inout) :: this
    integer,          intent(in)    :: idx(:)
      !! grid indices
    real,             intent(in)    :: dx(:)
      !! delta variable data (nval)

    integer :: i, j

    ASSERT(size(dx) == this%nval)

    j = 1
    do i = 1, this%nval
      select type (v => this%v(i)%p)
      class is (variable_real)
        call v%data%update(idx, dx(j))
        j = j + 1
      class is (variable_cmplx)
        call v%data%update(idx, complex(dx(j), dx(j+1)))
        j = j + 2
      end select
    end do
  end subroutine

  subroutine vselector_update_block(this, itab, dx)
    !! update data for whole block
    class(vselector), intent(inout) :: this
    integer,          intent(in)    :: itab
      !! table index
    real,             intent(in)    :: dx(:)
      !! delta variable data (nvals(itab))

    integer :: i, i0, i1, idx(this%g%idx_dim)

    ASSERT(size(dx) == this%nvals(itab))

    associate (tab => this%tab(itab)%p)
      i1 = 0
      do i = 1, tab%n
        i0 = i1 + 1
        i1 = i1 + this%nval

        idx = tab%get_idx(i)
        call this%vselector_update_single(idx, dx(i0:i1))
      end do
    end associate
  end subroutine

  subroutine vselector_update_all(this, dx)
    !! update data for all blocks
    class(vselector), intent(inout) :: this
    real,             intent(in)    :: dx(:)
      !! delta variable data (n)

    integer :: itab, i0, i1

    ASSERT(size(dx) == this%n)

    i1 = 0
    do itab = 1, this%ntab
      i0 = i1 + 1
      i1 = i1 + this%nvals(itab)

      call this%vselector_update_block(itab, dx(i0:i1))
    end do
  end subroutine

  subroutine vselector_print(this)
    !! print information to console
    class(vselector), intent(in) :: this

    integer :: i

    print "(A)", "  name: "//this%name
    print "(A)", "  vars:"
    do i = 1, this%nval
      print "(A)", "    "//this%v(i)%p%name
    end do
    print "(A)", "  tab:"
    do i = 1, this%ntab
      print "(A)", "    "//this%tab(i)%p%name
    end do
  end subroutine

end module
