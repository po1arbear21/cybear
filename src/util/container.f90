m4_include(macro.f90.inc)

module container_m

  use grid_m,          only: grid, grid_ptr, IDX_CELL, IDX_EDGE, IDX_FACE, IDX_VERTEX
  use grid0D_m,        only: grid0D
  use grid1D_m,        only: grid1D
  use tensor_grid_m,   only: tensor_grid
  use triang_grid_m,   only: triang_grid
  use grid_table_m,    only: grid_table
  use grid_data_m
  use normalization_m, only: norm, denorm
  use storage_m,       only: storage, variable, STORAGE_READ, STORAGE_WRITE, COMPR_DEFAULT, DYNAMIC_NO, DYNAMIC_EXT, DYNAMIC_APP, ERR_ALREADY_EXISTS, ERR_FILE_OP, ERR_INTERNAL, ERR_INVALID_ARGUMENTS, ERR_NOT_FOUND
  use string_m
  use variable_m,      only: variable_real
  use vector_m,        only: vector_string
  use vselector_m,     only: vselector

  implicit none

  m4_define({m4_typelist},{
    m4_X(int)
    m4_X(log)
    m4_X(real)
    m4_X(cmplx)
  })
  m4_define({m4_max_dim},{8})
  m4_define({m4_dimlist},{
    m4_ifelse($1,0,,{m4_dimlist(m4_decr($1),$2)})
    m4_Y($1,$2)
  })
  m4_define({m4_X},{m4_dimlist(m4_max_dim,$1)})
  m4_define({m4_list},m4_typelist)

  type, extends(storage) :: container
    type(vector_string) :: structs
  contains
    procedure :: open => container_open
    procedure :: has_struct => container_has_struct

    procedure :: lgrid => container_load_grid

    generic :: save => container_save_grid, &
                     & container_save_grid_table, &
                     & container_save_variable, &
                     & container_save_vselector
    
    generic :: load => container_load_grid_table, &
                     & container_load_grid0D, &
                     & container_load_grid1D, &
                     & container_load_tensor_grid, &
                     & container_load_triang_grid, &
                     & container_load_variable, &
                     & container_load_vselector
            
    m4_define({m4_X},{generic :: save => container_save_grid_data_$1})
    m4_typelist

    m4_define({m4_X},{generic :: load => container_load_grid_data_$1})
    m4_typelist

    procedure, private :: container_save_grid, container_save_grid_table, container_save_variable, container_save_vselector
    m4_define({m4_X},{procedure, private :: container_save_grid_data_$1})
    m4_typelist

    procedure, private :: container_load_grid0D, container_load_grid1D, container_load_tensor_grid, container_load_triang_grid, container_load_grid_table, container_load_variable, container_load_vselector
    m4_define({m4_X},{procedure, private :: container_load_grid_data_$1})
    m4_typelist
  end type

  m4_define({T},{string})
  m4_define({U},{grid_ptr})
  m4_include(../util/map_def.f90.inc)

contains

  subroutine container_open(this, file, flag, stat, err_msg)
    class(container),          intent(inout) :: this
    character(*),              intent(in)    :: file
    integer,         optional, intent(in)    :: flag
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    integer :: i, flag_
    type(variable), allocatable :: vars(:)

    flag_ = STORAGE_WRITE
    if (present(flag)) flag_ = flag

    call this%storage_open(file, flag_, stat, err_msg)
    allocate(vars(this%variables%n))
    
    call this%structs%init(0)
    ! Load existing structs
    call this%variables%to_array(values=vars)
    do i = 1, this%variables%n
      if (0 /= index(vars(i)%name%s, "sys/variables") .or. 0 /= index(vars(i)%name%s, "sys/vselectors")) then
        call this%structs%push(vars(i)%name)
      end if
    end do
  end subroutine

  function container_has_struct(this, name) result(found)
    class(container), intent(inout) :: this
    character(*),     intent(in)    :: name

    logical :: found
    integer :: i
    type(string), allocatable :: list(:)

    found = .false.

    list = this%structs%to_array()
    do i = 1, size(list)
      if (list(i)%s == name) then 
        found = .true.
        return
      end if
    end do

  end function

  recursive subroutine container_save_grid(this, g, stat, err_msg)
    class(container),          intent(inout) :: this
    class(grid), target,       intent(in)    :: g
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg
 
    select type (g)
      class is (grid0D)
        call this%write("sys/grids/" // g%name, new_string("grid0D"), stat=stat, err_msg=err_msg)
        if (present(stat)) then
          if (stat /= 0) return
        end if

      class is (grid1D)
        call this%write("sys/grids/" // g%name, new_string("grid1D"), stat=stat, err_msg=err_msg)
        if (present(stat)) then
          if (stat /= 0) return
        end if

        call this%write("sys/grids/" // g%name // "/i0",       g%i0)
        call this%write("sys/grids/" // g%name // "/unit",     g%unit(1))
        call this%write("sys/grids/" // g%name // "/vertices", denorm(g%x, g%unit(1)%s), compression=COMPR_DEFAULT)

      class is (tensor_grid)
      block
        integer      :: i
        type(string) :: grids(size(g%g))

        call this%write("sys/grids/" // g%name, new_string("tensor_grid"), stat=stat, err_msg=err_msg)
        if (present(stat)) then
          if (stat /= 0) return
        end if
        
        do i = 1, size(g%g)
          grids(i)%s = g%g(i)%p%name
          if (.not. this%contains(new_string("sys/grids/") // grids(i))) then
            call this%container_save_grid(g%g(i)%p)
          end if
        end do
        
        call this%write("sys/grids/" // g%name // "/grids", grids)
      end block

      class is (triang_grid)
        call this%write("sys/grids/" // g%name, new_string("triang_grid"), stat=stat, err_msg=err_msg)
        if (present(stat)) then
          if (stat /= 0) return
        end if

        call this%write("sys/grids/" // g%name // "/units", g%unit)

        block
          real, allocatable :: vertices(:,:)
          allocate(vertices(2, g%nvert))
          vertices = g%vert
          vertices(1,:) = denorm(vertices(1,:), g%unit(1)%s)
          vertices(2,:) = denorm(vertices(2,:), g%unit(2)%s)

          call this%write("sys/grids/" // g%name // "/vertices", vertices)
          call this%write("sys/grids/" // g%name // "/cells",    g%cell2vert)
        end block


      class default
        ! Generic writing of all the data of the grid, disable for now
        call program_error("Cannot write grid of this kind")
    end select
  end subroutine

  m4_define({m4_Y},{
  type is (grid_data$1_$2)
    call this%write("sys/data/" // name // "/data",   data%data, unit=unit, compression=COMPR_DEFAULT, dynamic=dflag)
})
m4_define({m4_X},{subroutine container_save_grid_data_$1(this, name, data, unit, dynamic, stat, err_msg)
  class(container),          intent(inout) :: this
  class(grid_data_$1),       intent(in)    :: data
  character(*),              intent(in)    :: name
  character(*),    optional, intent(in)    :: unit
  logical,         optional, intent(in)    :: dynamic
  integer(kind=4), optional, intent(out)   :: stat
  type(string),    optional, intent(out)   :: err_msg

  integer :: dflag

  dflag = DYNAMIC_NO
  if (present(dynamic)) then
    if (dynamic) dflag = DYNAMIC_EXT
  end if

  if (.not. this%contains(new_string("sys/data/" // name))) then
    call this%write("sys/data/" // name, new_string("$2_$1"), stat=stat, err_msg=err_msg)
    if (present(stat)) then
      if (stat /= 0) return
    end if

    call this%write("sys/data/" // name // "/bounds", data%idx_bnd)
  end if
  select type (data)
    m4_dimlist(m4_max_dim,$1)
  end select
end subroutine})
m4_typelist

  subroutine container_save_grid_table(this, tab, stat, err_msg)
    class(container),          intent(inout) :: this
    type(grid_table), target,  intent(in)    :: tab
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    ! Write all data in a tree under grid_tables/<name>
    call this%write("sys/tables/" // tab%name, new_string("grid_table"), stat=stat, err_msg=err_msg)
    if (present(stat)) then
          if (stat /= 0) return
        end if

    call this%write("sys/tables/" // tab%name // "/grid",     new_string(tab%g%name))
    call this%write("sys/tables/" // tab%name // "/type",     tab%idx_type)
    call this%write("sys/tables/" // tab%name // "/dir",      tab%idx_dir)
    call this%write("sys/tables/" // tab%name // "/flat2idx", tab%flat2idx, compression=COMPR_DEFAULT)
  end subroutine

  subroutine container_save_variable(this, var, parent, dynamic, stat, err_msg)
    class(container),             intent(inout) :: this
    class(variable_real), target, intent(in)    :: var
    character(*),                 intent(in)    :: parent
    logical,            optional, intent(in)    :: dynamic
    integer(kind=4),     optional, intent(out)   :: stat
    type(string),       optional, intent(out)   :: err_msg

    integer :: dflag
    
    dflag = DYNAMIC_NO
    if (present(dynamic)) then
      if (dynamic) dflag = DYNAMIC_EXT
    end if

    ! If not done already, the struct will be saved
    if (.not. this%has_struct("sys/variables/" // var%name)) then
      call this%write("sys/variables/" // var%name, new_string("variable_real"))
      call this%write("sys/variables/" // var%name // "/grid",  new_string(var%g%name))
      call this%write("sys/variables/" // var%name // "/unit",  new_string(var%unit))
      call this%write("sys/variables/" // var%name // "/type",  var%idx_type)
      call this%write("sys/variables/" // var%name // "/dir",   var%idx_dir)
      call this%structs%push(new_string("sys/variables/" // var%name))
    end if

    call this%write(parent // "/" // var%name, var%get(), unit=var%unit, compression=COMPR_DEFAULT, dynamic=dflag, stat=stat, err_msg=err_msg)

    ! Maybe we can shove in the array data%data directly without calling %get().
    ! However, the array bounds are not read correctly
    ! select type(data => var%data)
    !   type is (grid_data0_real)
    !     call this%write(parent // "/" // var%name, data%data)
    !   type is (grid_data1_real) ...
  end subroutine

  subroutine container_save_vselector(this, vsel, parent, dynamic, stat, err_msg)
    class(container),          intent(inout) :: this
    type(vselector), target,   intent(in)    :: vsel
    character(*),              intent(in)    :: parent
    logical,         optional, intent(in)    :: dynamic
    integer(kind=4), optional, intent(out)   :: stat
    type(string),    optional, intent(out)   :: err_msg

    integer      :: i, j, k, lbnd, ubnd
    real         :: values(vsel%n)
    type(string) :: var_names(vsel%nvar), tab_names(vsel%ntab)
    integer :: dflag = DYNAMIC_NO
    if (present(dynamic)) then
      if (dynamic) dflag = DYNAMIC_EXT
    end if

    ! If not done already, the struct will be saved
    if (.not. this%has_struct("sys/vselectors/" // vsel%name)) then
      do i = 1, vsel%nvar
        var_names(i) = new_string(vsel%v(i)%p%name)
      end do
      do i = 1, vsel%ntab
        tab_names(i) = new_string(vsel%tab(i)%p%name)
      end do

      call this%write("sys/vselectors/" // vsel%name, new_string("vselector"))
      call this%write("sys/vselectors/" // vsel%name // "/grid",  new_string(vsel%g%name))
      call this%write("sys/vselectors/" // vsel%name // "/variables", var_names)
      call this%write("sys/vselectors/" // vsel%name // "/tables", tab_names)
      call this%structs%push(new_string("sys/vselectors/" // vsel%name))
    end if

    ! Denormalize every variable on each block/tab. CAUTION: Assumed only variable_real
    k = 0
    values = vsel%get()
    do i = 1, vsel%ntab
      do j = 1, vsel%nvar
        lbnd = k + j
        ubnd = k + vsel%nvals(i) + j - vsel%nvar
        values(lbnd:ubnd:vsel%nvar) = denorm(values(lbnd:ubnd:vsel%nvar), vsel%v(j)%p%unit) ! PLEASE LET ME DO THIS
      end do
      k = k + vsel%nvals(i)
    end do
    call this%write(new_string(parent // "/" // vsel%name), values, compression=COMPR_DEFAULT, dynamic=dflag, stat=stat, err_msg=err_msg)
  end subroutine

  subroutine container_load_grid0D(this, name, g)
    class(container), intent(inout) :: this
    character(*),     intent(in)    :: name
    type(grid0D),     intent(inout) :: g

    integer      :: i
    type(string) :: var_name, gtype
    
    i = index(name, "/")
    if (i /= 0) call program_error("Grid name with '/' is not allowed")
    var_name = new_string("sys/grids/") // name
    
    ! Do we need to check the allocation status of the grid?

    call this%read(var_name, gtype)
    if (gtype%s /= "grid0D") call program_error("Grid variable (grid0D) does not match the stored grid class: " // gtype%s)
    call g%init(name)
  end subroutine

  subroutine container_load_grid1D(this, name, g)
    class(container), intent(inout) :: this
    character(*),     intent(in)    :: name
    type(grid1D),     intent(inout) :: g

    integer           :: i
    real, allocatable :: x(:)
    type(string)      :: var_name, gtype, unit
    
    ! Do we need to check the allocation status of the grid?
    
    i = index(name, "/")
    if (i /= 0) call program_error("Grid name with '/' is not allowed")
    var_name = new_string("sys/grids/") // name

    call this%read(var_name, gtype)
    if (gtype%s /= "grid1D") call program_error("Grid variable does not match the stored grid class")
    call this%read(var_name%s // "/i0",       i)
    call this%read(var_name%s // "/unit",     unit)
    call this%read(var_name%s // "/vertices", x)

    x = norm(x, unit%s)

    call g%init(name, x, i, unit=unit%s)
  end subroutine

  subroutine container_load_grid(this, name, g)
    class(container), intent(inout) :: this
    character(*),     intent(in)    :: name
    class(grid),      intent(inout) :: g

    select type (g)
      class is (grid0D)
        call this%container_load_grid0D(name, g)
      class is (grid1D)
        call this%container_load_grid1D(name, g)
      class is (tensor_grid)
        call this%container_load_tensor_grid(name, g)
      class is (triang_grid)
        call this%container_load_triang_grid(name, g)
    end select
  end subroutine

  subroutine container_load_tensor_grid(this, name, g, map)
    class(container),          intent(inout) :: this
    character(*),              intent(in)    :: name
    type(tensor_grid),         intent(inout) :: g
    type(map_string_grid_ptr), optional, target, intent(inout) :: map

    type(map_string_grid_ptr), pointer :: map_
    integer                   :: i
    type(string), allocatable :: grids(:)
    type(string)              :: var_name, gtype
    
    type(grid_ptr)              :: gptr
    type(grid_ptr), allocatable :: gs(:)
    type(mapnode_string_grid_ptr), pointer :: p => null()

    if (present(map)) map_ => map
    
    i = index(name, "/")
    if (i /= 0) call program_error("Grid name with '/' is not allowed")
    var_name = new_string("sys/grids/") // name
    
    ! Do we need to check the allocation status of the grid?

    call this%read(var_name, gtype)
    if (gtype%s /= "tensor_grid") call program_error("Grid variable does not match the stored grid class")
    call this%read(var_name%s // "/grids", grids)
    
    allocate(gs(size(grids)))
    do i = 1, size(grids)
      ! Look for the grid in the map
      if (associated(map_)) p => map_%find(grids(i))
      if (associated(p)) then
        gs(i) = p%value
        continue
      end if

      gptr%p => null()
      call this%read("sys/grids/" // grids(i)%s, gtype)
      select case (gtype%s)
        case ("grid0D")
          allocate(grid0D :: gptr%p)
        case ("grid1D")
          allocate(grid1D :: gptr%p)
        case ("tensor_grid")
          call program_error("Inception not allowed")
        case ("triang_grid")
          allocate(triang_grid :: gptr%p)
      end select
      call this%lgrid(grids(i)%s, gptr%p)

      if (associated(map_)) call map_%set(grids(i), gptr)
      gs(i) = gptr
    end do

    call g%init(name, gs)
  end subroutine

  subroutine container_load_triang_grid(this, name, g)
    class(container),  intent(inout) :: this
    character(*),      intent(in)    :: name
    type(triang_grid), intent(inout) :: g
    
    integer                   :: i
    real,         allocatable :: vert(:,:)
    integer,      allocatable :: cells(:,:)
    type(string)              :: var_name, gtype
    type(string), allocatable :: units(:)
    
    i = index(name, "/")
    if (i /= 0) call program_error("Grid name with '/' is not allowed")
    var_name = new_string("sys/grids/") // name
    
    ! Do we need to check the allocation status of the grid?

    call this%read(var_name, gtype)
    if (gtype%s /= "triang_grid") call program_error("Grid variable does not match the stored grid class")
    call this%read(var_name%s // "/units", units)
    call this%read(var_name%s // "/vertices", vert)
    call this%read(var_name%s // "/cells", cells)

    if (size(units, dim=1) /= 2) then
      call program_error("Got invalid amount of units back")
    end if

    vert(1,:) = norm(vert(1,:), units(1)%s)
    vert(2,:) = norm(vert(2,:), units(2)%s)

    call g%init(name, vert, cells, unit=units)
  end subroutine

  subroutine container_load_grid_table(this, name, tab)
    class(container), intent(inout) :: this
    character(*),     intent(in)    :: name
    type(grid_table), intent(inout) :: tab

    integer :: i, idx(tab%g%idx_dim)
    type(string) :: s, var_name

    i = index(name, "/")
    if (i /= 0) call program_error("Grid table name with '/' is not allowed")
    var_name = new_string("sys/tables/") // name

    call this%read(var_name, s)
    if (s%s /= "grid_table") call program_error("Variable not of type grid_table")
 
    ! Check the initialization status, we want to be called after init() and instead of init_final()
    if (tab%finished_init)       call program_error("The grid table is already fully initialized")
    if (allocated(tab%flat2idx)) call program_error("The flat indices are already allocated")

    ! Compare the defined grids
    call this%read(var_name // new_string("/grid"), s)
    if (s%s /= tab%g%name) call program_error("Loading a grid table requires the same grids")
    call this%read(var_name // new_string("/type"), i)
    if (i /= tab%idx_type) call program_error("The index types do not match")
    call this%read(var_name // new_string("/dir"), i)
    if (i /= tab%idx_dir) call program_error("The index directions do not match")

    ! Set the data, we make the assumption that tab%flags is initialized and tab%idx2flat is not!!!
    call allocate_grid_data(tab%idx2flat, tab%g%idx_dim)
    call tab%idx2flat%init(tab%g, tab%idx_type, tab%idx_dir, 0)
    call allocate_grid_data(tab%flags, tab%g%idx_dim)
    call tab%flags%init(tab%g, tab%idx_type, tab%idx_dir, .false.)
    call this%read(var_name // new_string("/flat2idx"), tab%flat2idx)
    
    do i = 1, size(tab%flat2idx, 2)
      idx = tab%flat2idx(:,i)
      call tab%flags%set(idx, .true.)
      call tab%idx2flat%set(idx, i)
    end do
    tab%finished_init = .true. 
  end subroutine

  subroutine container_load_variable(this, name, var)
    class(container),             intent(inout) :: this
    character(*),                 intent(in)    :: name
    class(variable_real), target, intent(inout) :: var

    integer           :: i
    type(string)      :: s
    real, allocatable :: values(:)
    
    call this%read("sys/variables/" // var%name, s)
    if (s%s /= "variable_real") call program_error("The stored entry at '"// name //"' is not a variable_real")
    
    ! Compare stored metadata
    call this%read("sys/variables/" // var%name // "/grid", s)
    if (s%s /= var%g%name) call program_error("The grids do not match")
    call this%read("sys/variables/" // var%name // "/unit", s)
    if (s%s /= var%unit) call program_error("The units do not match")
    call this%read("sys/variables/" // var%name // "/type", i)
    if (i /= var%idx_type) call program_error("The index type does not match")
    call this%read("sys/variables/" // var%name // "/dir", i)
    if (i /= var%idx_dir) call program_error("The index direction does not match")

    ! Ideally no extra allocation here and set the variable grid data directly
    call this%read(name, values)
    call var%set(values)
  end subroutine

  m4_define({m4_Y},{
  type is (grid_data$1_$2)
    call this%read(new_string("grid_data/" // name // "/data"), data%data)
})
m4_define({m4_X},{subroutine container_load_grid_data_$1(this, name, data)
  !! get data for all points in flat array
  class(container),    intent(inout) :: this
  class(grid_data_$1), intent(inout) :: data
  character(*),        intent(in)    :: name

  call this%read(new_string("grid_data/" // name // "/bounds"), data%idx_bnd)
  select type (data)
    m4_dimlist(m4_max_dim,$1)
  end select
end subroutine})
m4_typelist

  subroutine container_load_vselector(this, name, vsel)
    class(container),        intent(inout) :: this
    type(string),            intent(in)    :: name
    type(vselector), target, intent(inout) :: vsel

    integer           :: i, j, k, lbnd, ubnd
    real, allocatable :: values(:)
    
    ! TODO: Is vsel initialized and matches the struct definition?

    k = 0
    call this%read(name, values)
    do i = 1, vsel%ntab
      do j = 1, vsel%nvar
        lbnd = k + j
        ubnd = k + vsel%nvals(i) + j - vsel%nvar
        values(lbnd:ubnd:vsel%nvar) = norm(values(lbnd:ubnd:vsel%nvar), vsel%v(j)%p%unit)
      end do
      k = k + vsel%nvals(i)
    end do
    call vsel%set(values)
    
  end subroutine
  
  m4_define({T},{string})
  m4_define({U},{grid_ptr})
  m4_include(../util/map_imp.f90.inc)

end module