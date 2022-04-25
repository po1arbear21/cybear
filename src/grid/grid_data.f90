m4_include(../util/macro.f90.inc)

module grid_data_m

  use error_m,       only: assert_failed, program_error
  use grid_m,        only: grid
  use json_m,        only: json_object
  use output_file_m, only: output_file

  implicit none

  ! supported grid_data types
  m4_define({m4_typelist},{
    m4_X(int)
    m4_X(log)
    m4_X(real)
    m4_X(cmplx)
  })

  ! physical unit?
  m4_define({m4_unit},{m4_divert(m4_ifelse($1,real,0,m4_ifelse($1,cmplx,0,-1)))$2{}m4_divert(0)})

  ! grid_data dimensions (0..max_dim)
  m4_define({m4_max_dim},{8})
  m4_define({m4_dimlist},{
    m4_ifelse($1,0,,{m4_dimlist(m4_decr($1),$2)})
    m4_Y($1,$2)
  })

  ! combine type-list with dimension-list to create full list
  m4_define({m4_X},{m4_dimlist(m4_max_dim,$1)})
  m4_define({m4_list},m4_typelist)

  ! make base grid_data types public
  m4_define({m4_X},{public grid_data_$1})
  m4_typelist

  ! make derived grid_data types public
  m4_define({m4_Y},{public grid_data$1_$2})
  m4_list

  ! idx(1), idx(2), ....
  m4_define({m4_idx},{m4_ifelse($1,1,,{m4_idx(m4_decr($1)),}){idx($1)}})

  ! (idx(1), idx(2), ...)
  m4_define({m4_pidx},{m4_ifelse($1,0,,{(m4_idx($1))})})

  ! this%idx_bnd(1), this%idx_bnd(2), ...
  m4_define({m4_this_idx_bnd},{m4_ifelse($1,1,,{m4_this_idx_bnd(m4_decr($1)),}){this%idx_bnd($1)}})

  ! abstract grid_data type definitions
  m4_define({m4_Y},{procedure :: get_ptr$1 => grid_data_$2_get_ptr$1})
  m4_define({m4_X},{
    type, abstract :: grid_data_$1
      !! grid based array with bounds checks
      integer, allocatable :: idx_bnd(:)
        !! index bounds
      integer              :: n
        !! total number of values
    contains
      procedure :: init     => grid_data_$1_init
      procedure :: destruct => grid_data_$1_destruct
      procedure :: reset    => grid_data_$1_reset

      m4_dimlist(m4_max_dim,$1)

      generic   :: get      => grid_data_$1_get_point, grid_data_$1_get_all
      generic   :: set      => grid_data_$1_set_point, grid_data_$1_set_all
      m4_ifelse($1,log,,{generic   :: update   => grid_data_$1_update_point, grid_data_$1_update_all})

      procedure :: output   => grid_data_$1_output

      procedure, private :: grid_data_$1_get_point, grid_data_$1_get_all
      procedure, private :: grid_data_$1_set_point, grid_data_$1_set_all
      m4_ifelse($1,log,,{procedure, private :: grid_data_$1_update_point, grid_data_$1_update_all})
    end type
  })
  m4_typelist

  ! derived grid_data type definitions
  m4_define({m4_Y}, {
    type, extends(grid_data_$2) :: grid_data$1_$2
      m4_type($2){}m4_ifelse($1,0,,{, allocatable}) :: data{}m4_pshape($1)
        !! data
    end type
  })
  m4_list

  interface allocate_grid_data
    m4_define({m4_X},{
      module procedure :: allocate_grid_data0_$1
      module procedure :: allocate_grid_data1_$1
      module procedure :: allocate_grid_data2_$1
      module procedure :: allocate_grid_data3_$1
      module procedure :: allocate_grid_data4_$1
    })
    m4_typelist
  end interface

contains

  m4_define({m4_Y},{
    case ($1)
      allocate (grid_data$1_$2 :: gdata)
  })
  m4_define({m4_X},{subroutine allocate_grid_data0_$1(gdata, dim)
    !! allocate grid data
    class(grid_data_$1), allocatable, intent(out) :: gdata
    integer,                          intent(in)  :: dim
      !! grid data dimension

    select case (dim)
      m4_dimlist(m4_max_dim,$1)
    case default
      call program_error("grid data dimension must be in range 0:{}m4_max_dim{}")
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    case ($1)
      allocate (grid_data$1_$2 :: gdata(i0:i1))
  })
  m4_define({m4_X},{subroutine allocate_grid_data1_$1(gdata, dim, i0, i1)
    !! allocate grid data
    class(grid_data_$1), allocatable, intent(out) :: gdata(:)
    integer,                          intent(in)  :: dim
      !! grid data dimension
    integer,                          intent(in)  :: i0
      !! lower array bound
    integer,                          intent(in)  :: i1
      !! upper array bound

    select case (dim)
      m4_dimlist(m4_max_dim,$1)
    case default
      call program_error("grid data dimension must be in range 0:{}m4_max_dim{}")
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    case ($1)
      allocate (grid_data$1_$2 :: gdata(i0(1):i1(1),i0(2):i1(2)))
  })
  m4_define({m4_X},{subroutine allocate_grid_data2_$1(gdata, dim, i0, i1)
    !! allocate grid data
    class(grid_data_$1), allocatable, intent(out) :: gdata(:,:)
    integer,                          intent(in)  :: dim
      !! grid data dimension
    integer,                          intent(in)  :: i0(2)
      !! lower array bound
    integer,                          intent(in)  :: i1(2)
      !! upper array bound

    select case (dim)
      m4_dimlist(m4_max_dim,$1)
    case default
      call program_error("grid data dimension must be in range 0:{}m4_max_dim{}")
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    case ($1)
      allocate (grid_data$1_$2 :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3)))
  })
  m4_define({m4_X},{subroutine allocate_grid_data3_$1(gdata, dim, i0, i1)
    !! allocate grid data
    class(grid_data_$1), allocatable, intent(out) :: gdata(:,:,:)
    integer,                          intent(in)  :: dim
      !! grid data dimension
    integer,                          intent(in)  :: i0(3)
      !! lower array bound
    integer,                          intent(in)  :: i1(3)
      !! upper array bound

    select case (dim)
      m4_dimlist(m4_max_dim,$1)
    case default
      call program_error("grid data dimension must be in range 0:{}m4_max_dim{}")
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    case ($1)
      allocate (grid_data$1_$2 :: gdata(i0(1):i1(1),i0(2):i1(2),i0(3):i1(3),i0(4):i1(4)))
  })
  m4_define({m4_X},{subroutine allocate_grid_data4_$1(gdata, dim, i0, i1)
    !! allocate grid data
    class(grid_data_$1), allocatable, intent(out) :: gdata(:,:,:,:)
    integer,                          intent(in)  :: dim
      !! grid data dimension
    integer,                          intent(in)  :: i0(4)
      !! lower array bound
    integer,                          intent(in)  :: i1(4)
      !! upper array bound

    select case (dim)
      m4_dimlist(m4_max_dim,$1)
    case default
      call program_error("grid data dimension must be in range 0:{}m4_max_dim{}")
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      m4_ifelse($1,0,{this%data = d0_},{allocate (this%data(m4_this_idx_bnd($1)), source = d0_)})
  })
  m4_define({m4_X},{subroutine grid_data_$1_init(this, g, idx_type, idx_dir, d0)
    !! initialize and allocate grid data
    class(grid_data_$1),   intent(out) :: this
    class(grid),           intent(in)  :: g
    integer,               intent(in)  :: idx_type
      !! grid index type (IDX_VERTEX, IDX_EDGE, IDX_FACE or IDX_CELL)
    integer,               intent(in)  :: idx_dir
      !! index direction for edges and faces (must be 0 for IDX_VERTEX and IDX_CELL)
    m4_type($1), optional, intent(in)  :: d0
      !! initial value

    m4_type($1) :: d0_

    ! get index bounds
    allocate (this%idx_bnd(g%idx_dim))
    call g%get_idx_bnd(idx_type, idx_dir, this%idx_bnd)

    ! set total number of values
    this%n = product(this%idx_bnd) ! = 1 for idx_dim == 0

    ! get initial value
    d0_ = m4_ifelse($1,log,{.false.},0)
    if (present(d0)) d0_ = d0

    ! allocate data array
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{m4_ifelse($1,0,,{
    type is (grid_data$1_$2)
      if (allocated(this%data)) deallocate (this%data)
  })})
  m4_define({m4_X},{subroutine grid_data_$1_destruct(this)
    class(grid_data_$1), intent(inout) :: this

    if (allocated(this%idx_bnd)) deallocate (this%idx_bnd)
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      this%data = d0_
  })
  m4_define({m4_X},{subroutine grid_data_$1_reset(this, d0)
    !! reset data for all points
    class(grid_data_$1),   intent(inout) :: this
    m4_type($1), optional, intent(in)    :: d0
      !! optional: reset to this value

    m4_type($1) :: d0_

    ! get initial value
    d0_ = m4_ifelse($1,log,{.false.},0)
    if (present(d0)) d0_ = d0

    ! set data array
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{function grid_data_$2_get_ptr$1(this) result(p)
    class(grid_data_$2), target, intent(in) :: this
    type(grid_data$1_$2), pointer           :: p

    select type (this)
    type is (grid_data$1_$2)
      p => this
    class default
      p => null()
    end select
  end function})
  m4_list

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      d = this%data{}m4_pidx($1)
  })
  m4_define({m4_X},{function grid_data_$1_get_point(this, idx) result(d)
    !! get data for single point with bounds check (out of bounds: return default value)
    class(grid_data_$1), intent(in) :: this
    integer,             intent(in) :: idx(:)
      !! grid indices (idx_dim)
    m4_type($1)                     :: d
      !! return data

    m4_assert(size(idx) == size(this%idx_bnd)             )
    m4_assert(all(idx >= 1) .and. all(idx <= this%idx_bnd))

    ! default value
    d = m4_ifelse($1,log,{.false.},0)

    ! get data
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end function})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      d = m4_ifelse($1,0,{this%data},{reshape(this%data, [this%n])})
  })
  m4_define({m4_X},{function grid_data_$1_get_all(this) result(d)
    !! get data for all points in flat array
    class(grid_data_$1), intent(in) :: this
    m4_type($1)                     :: d(this%n)
      !! return all data

    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end function})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      this%data{}m4_pidx($1) = d
  })
  m4_define({m4_X},{subroutine grid_data_$1_set_point(this, idx, d)
    !! set data for single point with bounds check (do nothing if out of bounds)
    class(grid_data_$1), intent(inout) :: this
    integer,             intent(in)    :: idx(:)
      !! grid indices (idx_dim)
    m4_type($1),         intent(in)    :: d
      !! new value

    m4_assert(size(idx) == size(this%idx_bnd)             )
    m4_assert(all(idx >= 1) .and. all(idx <= this%idx_bnd))

    ! set data
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      this%data = m4_ifelse($1,0,{d(1)},{m4_ifelse($1,1,{d},{reshape(d, [m4_this_idx_bnd($1)])})})
  })
  m4_define({m4_X},{subroutine grid_data_$1_set_all(this, d)
    !! set data for all points
    class(grid_data_$1), intent(inout) :: this
    m4_type($1),         intent(in)    :: d(:)
      !! new values

    m4_assert(size(d) == this%n)

    ! set data
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      this%data{}m4_pidx($1) = this%data{}m4_pidx($1) + d
  })
  m4_define({m4_X},{m4_divert(m4_ifelse($1,log,-1,0))subroutine grid_data_$1_update_point(this, idx, d)
    !! update data for single point with bounds check (do nothing if out of bounds)
    class(grid_data_$1), intent(inout) :: this
    integer,             intent(in)    :: idx(:)
      !! grid indices (idx_dim)
    m4_type($1),         intent(in)    :: d
      !! delta value

    m4_assert(size(idx) == size(this%idx_bnd)             )
    m4_assert(all(idx >= 1) .and. all(idx <= this%idx_bnd))

    ! update data
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine{}m4_divert(0)})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      this%data = this%data + m4_ifelse($1,0,{d(1)},{m4_ifelse($1,1,{d},{reshape(d, [m4_this_idx_bnd($1)])})})
  })
  m4_define({m4_X},{m4_divert(m4_ifelse($1,log,-1,0))subroutine grid_data_$1_update_all(this, d)
    !! update data for all points
    class(grid_data_$1), intent(inout) :: this
    m4_type($1),         intent(in)    :: d(:)
      !! delta values

    m4_assert(size(d) == this%n)

    ! update data
    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine{}m4_divert(0)})
  m4_typelist

  m4_define({m4_Y},{
    type is (grid_data$1_$2)
      call of%write(obj, name, this%data{}m4_unit($2,{, unit = unit}))
  })
  m4_define({m4_X},{subroutine grid_data_$1_output(this, of, obj, name{}m4_unit($1,{, unit}))
    !! output grid data
    class(grid_data_$1),        intent(in)    :: this
    type(output_file),          intent(inout) :: of
      !! output file handle
    type(json_object), pointer, intent(inout) :: obj
      !! parent object in output file
    character(*),               intent(in)    :: name
      !! data name
    m4_unit($1,{
    character(*),               intent(in)    :: unit
      !! physical unit token
    })

    select type (this)
      m4_dimlist(m4_max_dim,$1)
    end select
  end subroutine})
  m4_typelist

end module
