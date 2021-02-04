module grid_table_m
  use grid_data_m
  implicit none

  type grid_table
    !! table to select points from grid

    class(grid), pointer :: g => null()
      !! pointer to grid

    integer             :: idx_type
      !! index type
    integer             :: dir
      !! direction index
    type(grid_data_log) :: flags
      !! include/exclude point (product idx_bnd)

    integer              :: n
      !! number of entries
    integer, allocatable :: flat2idx(:,:)
      !! flat index to grid indices (idx_dim x n)
    type(grid_data_int)  :: idx2flat
      !! grid indices to flat index
  contains
    procedure :: init       => grid_table_init
    procedure :: init_final => grid_table_init_final
    procedure :: get_idx    => grid_table_get_idx
    procedure :: get_flat   => grid_table_get_flat
  end type

contains

  subroutine grid_table_init(this, g, idx_type, dir)
    !! initialize grid table
    class(grid_table),   intent(out) :: this
    class(grid), target, intent(in)  :: g
    integer,             intent(in)  :: idx_type
      !! grid index type
    integer,             intent(in)  :: dir
      !! index of direction (only used for IDX_EDGE and IDX_FACE; range = 1:idx_dim)

    ! save members
    this%g        => g
    this%idx_type = idx_type
    this%dir      = dir

    ! initialize flags data
    call this%flags%init(g, idx_type, dir)
  end subroutine

  subroutine grid_table_init_final(this)
    !! initialize internal tables (call this after flags have been set)
    class(grid_table), intent(inout) :: this

    integer :: i, j, idx(this%g%idx_dim)

    associate (idx_dim => this%g%idx_dim, idx_bnd => this%flags%idx_bnd)
      ! count number of entries
      this%n = count(this%flags%data)

      ! allocate internal tables
      allocate (this%flat2idx(idx_dim, this%n), source = 0)
      call this%idx2flat%init(this%g, this%idx_type, this%dir)

      ! collect indices for which flag is set
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
    end associate
  end subroutine

  subroutine grid_table_get_idx(this, i, idx)
    !! convert flat index to grid indices
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: i
      !! flat index
    integer,           intent(out) :: idx(:)
      !! return grid indices

    idx = this%flat2idx(:,i)
  end subroutine

  subroutine grid_table_get_flat(this, idx, i)
    !! convert grid indices to flat index
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: idx(:)
      !! grid indices (idx_dim)
    integer,           intent(out) :: i
      !! return flat index (0 if grid point is not part of table)

    i = this%idx2flat%get(idx)
  end subroutine

end module