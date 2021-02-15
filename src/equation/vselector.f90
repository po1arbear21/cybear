#include "../util/macro.f90.inc"

module vselector_m

  use error_m
  use grid_m,       only: grid
  use grid_table_m, only: grid_table_ptr
  use grid_data_m,  only: grid_data_int, allocate_grid_data
  use variable_m,   only: variable_ptr

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

    integer              :: nval
      !! number of values per grid point: nval = size(v)
    integer              :: ntab
      !! number of tables/blocks: ntab = size(tab)
    integer, allocatable :: nvals(:)
      !! number of values per block: nvals(itab) = tab(itab)%n * nval
    integer              :: n
      !! total number of values selected: n = sum(nvals)

    class(grid_data_int), allocatable :: itab
      !! get table index from grid indices
  contains
    procedure :: init    => vselector_init
    procedure :: reset   => vselector_reset
    procedure :: compare => vselector_compare

    procedure :: vselector_get_single
    procedure :: vselector_get_block
    procedure :: vselector_get_all
    generic   :: get => vselector_get_single, vselector_get_block, vselector_get_all

    procedure :: vselector_set_single
    procedure :: vselector_set_block
    procedure :: vselector_set_all
    generic   :: set => vselector_set_single, vselector_set_block, vselector_set_all

    procedure :: vselector_update_single
    procedure :: vselector_update_block
    procedure :: vselector_update_all
    generic   :: update => vselector_update_single, vselector_update_block, vselector_update_all

    procedure :: print => vselector_print
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

  subroutine vselector_init(this, name, v, tab)
    !! initialize variable selector
    class(vselector),     intent(out) :: this
    character(*),         intent(in)  :: name
      !! selector name
    type(variable_ptr),   intent(in)  :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)  :: tab(:)
      !! grid table pointers

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

    ! get number of selected values
    this%nval = size(v)
    this%ntab = size(tab)
    allocate (this%nvals(size(tab)), source = 0)
    do itab = 1, this%ntab
      this%nvals(itab) = tab(itab)%p%n * this%nval
    end do
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

  subroutine vselector_reset(this)
    !! reset all data selected to zero
    class(vselector), intent(inout) :: this

    integer :: i, j, itab, idx(this%g%idx_dim)

    do itab = 1, this%ntab
      associate (tab => this%tab(itab)%p)
        do i = 1, tab%n
          idx = tab%get_idx(i)
          do j = 1, this%nval
            call this%v(j)%p%data%set(idx, 0.0)
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
    if (this%nval /= other%nval) return
    do i = 1, this%nval
      if (.not. associated(this%v(i)%p, other%v(i)%p)) return
    end do
    if (this%ntab /= other%ntab) return
    do i = 1, this%ntab
      if (.not. associated(this%tab(i)%p, other%tab(i)%p)) return
    end do

    c = .true.
  end function

  function vselector_get_single(this, idx) result(x)
    !! get data for single grid index tuple
    class(vselector), intent(in) :: this
    integer,          intent(in) :: idx(:)
      !! grid indices
    real                         :: x(this%nval)
      !! return variable data at idx

    integer :: i

    do i = 1, this%nval
      x(i) = this%v(i)%p%data%get(idx)
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

    integer :: i

    ASSERT(size(x) == this%nval)

    do i = 1, this%nval
      call this%v(i)%p%data%set(idx, x(i))
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

    integer :: i

    ASSERT(size(dx) == this%nval)

    do i = 1, this%nval
      call this%v(i)%p%data%update(idx, dx(i))
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
