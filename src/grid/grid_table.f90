module grid_table_m
  use grid_data_m
  implicit none

  type grid_table
    !! table to select points from grid

    character(:), allocatable :: name
      !! grid table name

    class(grid), pointer :: g => null()
      !! pointer to grid

    integer             :: idx_type
      !! index type
    integer             :: idx_dir
      !! index direction for edges and faces
    class(grid_data_log), allocatable :: flags
      !! include/exclude point (product idx_bnd)

    integer                           :: n
      !! number of entries
    integer,              allocatable :: flat2idx(:,:)
      !! flat index to grid indices (idx_dim x n)
    class(grid_data_int), allocatable :: idx2flat
      !! grid indices to flat index
  contains
    procedure :: init       => grid_table_init
    procedure :: init_final => grid_table_init_final
    procedure :: get_idx    => grid_table_get_idx
    procedure :: get_flat   => grid_table_get_flat
  end type

  type grid_table_ptr
    type(grid_table), pointer :: p => null()
  end type

contains

  subroutine grid_table_init(this, name, g, idx_type, idx_dir)
    !! initialize grid table
    class(grid_table),   intent(out) :: this
    character(*),        intent(in)  :: name
      !! grid table name
    class(grid), target, intent(in)  :: g
    integer,             intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,             intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)

    ! set members
    this%name     =  name
    this%g        => g
    this%idx_type =  idx_type
    this%idx_dir  =  idx_dir

    ! initialize flags data
    call allocate_grid_data(this%flags, g%idx_dim)
    call this%flags%init(g, idx_type, idx_dir)
  end subroutine

  subroutine grid_table_init_final(this)
    !! initialize internal tables (call this after flags have been set)
    class(grid_table), intent(inout) :: this

    integer :: i, j, idx(this%g%idx_dim)

    associate (idx_dim => this%g%idx_dim, idx_bnd => this%flags%idx_bnd)
      ! count number of entries
      this%n = count(this%flags%get())

      ! allocate internal tables
      allocate (this%flat2idx(idx_dim, this%n), source = 0)
      call allocate_grid_data(this%idx2flat, this%g%idx_dim)
      call this%idx2flat%init(this%g, this%idx_type, this%idx_dir)

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

  function grid_table_get_idx(this, i) result(idx)
    !! convert flat index to grid indices
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: i
      !! flat index
    integer                        :: idx(this%g%idx_dim)
      !! return grid indices

    idx = this%flat2idx(:,i)
  end function

  function grid_table_get_flat(this, idx) result(i)
    !! convert grid indices to flat index
    class(grid_table), intent(in)  :: this
    integer,           intent(in)  :: idx(:)
      !! grid indices (idx_dim)
    integer                        :: i
      !! return flat index (0 if grid point is not part of table)

    i = this%idx2flat%get(idx)
  end function

end module
