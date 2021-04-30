#include "../util/macro.f90.inc"

submodule (grid_m) grid_table_m

contains

  module subroutine grid_table_init(this, name, g, idx_type, idx_dir, initial_flags)
    !! initialize grid table
    class(grid_table),   intent(out) :: this
    character(*),        intent(in)  :: name
      !! grid table name
    class(grid), target, intent(in)  :: g
    integer,             intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,             intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    logical, optional,   intent(in)  :: initial_flags
      !! initialize flags to this value (default: false)

    ! set members
    this%name     =  name
    this%g        => g
    this%idx_type =  idx_type
    this%idx_dir  =  idx_dir

    ! initialize flags data
    call allocate_grid_data(this%flags, g%idx_dim)
    call this%flags%init(g, idx_type, idx_dir, d0 = initial_flags)
  end subroutine

  module subroutine grid_table_init_final(this)
    !! initialize internal tables (call this after flags have been set)
    class(grid_table), intent(inout) :: this

    integer :: i, j, idx(this%g%idx_dim)

    associate (idx_dim => this%g%idx_dim, idx_bnd => this%flags%idx_bnd)
      ! count number of entries
      this%n = count(this%flags%get())

      ! allocate internal tables
      allocate (this%flat2idx(idx_dim,this%n), source = 0)
      call allocate_grid_data(this%idx2flat, this%g%idx_dim)
      call this%idx2flat%init(this%g, this%idx_type, this%idx_dir)

      ! collect indices for which flag is set
      if (idx_dim == 0) then
        call this%idx2flat%set(idx, 1)
      else
        idx = 1
        i   = 0
        do while (idx(idx_dim) <= idx_bnd(idx_dim))
          if (this%flags%get(idx)) then
            i = i + 1
            this%flat2idx(:,i) = idx
            call this%idx2flat%set(idx, i)
          end if

          ! go to next index
          idx(1) = idx(1) + 1
          do j = 1, idx_dim-1
            if (idx(j) <= idx_bnd(j)) exit
            idx(j  ) = 1
            idx(j+1) = idx(j+1) + 1
          end do
        end do
      end if
    end associate
  end subroutine

  module function grid_table_get_ptr(this) result(ptr)
    !! returns pointer type to this grid_table
    class(grid_table), target, intent(in) :: this
    type(grid_table_ptr)                  :: ptr

    ptr%p => this
  end function

  module function grid_table_get_idx(this, i) result(idx)
    !! convert flat index to grid indices
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: i
      !! flat index
    integer                        :: idx(this%g%idx_dim)
      !! return grid indices

    ASSERT((i > 0 ) .and. (i <= this%n))

    idx = this%flat2idx(:,i)
  end function

  module function grid_table_get_flat(this, idx) result(i)
    !! convert grid indices to flat index
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: idx(:)
      !! grid indices (idx_dim)
    integer                        :: i
      !! return flat index (0 if grid point is not part of table)

    ASSERT(this%g%idx_allowed(this%idx_type, this%idx_dir, idx=idx))

    i = this%idx2flat%get(idx)
  end function

end submodule
