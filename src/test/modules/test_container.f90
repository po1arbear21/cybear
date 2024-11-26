module test_container_m

  use test_case_m,     only: test_case
  use container_m,     only: container
  use error_m,         only: program_error
  use storage_m,       only: STORAGE_NEW, STORAGE_UNPACK
  use string_m,        only: new_string
  use grid_m,          only: grid, IDX_VERTEX, IDX_CELL, IDX_EDGE, IDX_FACE
  use grid0D_m,        only: grid0D
  use grid1D_m,        only: grid1D
  use grid_data_m,     only: grid_data0_real, grid_data2_real
  use grid_table_m,    only: grid_table
  use normalization_m, only: norm, denorm
  use tensor_grid_m,   only: tensor_grid
  use triang_grid_m,   only: triang_grid
  use variable_m,      only: variable_real
  use vselector_m,     only: vselector
  use math_m,          only: linspace
  
  implicit none

  private
  public test_container

  type, extends(variable_real) :: var
    type(grid_data0_real), pointer :: d0 => null()
    type(grid_data2_real), pointer :: d2 => null()
  contains
    procedure :: init => var_init
  end type

contains

  subroutine var_init(this, name, unit, g, idx_type, idx_dir)
    class(var),   intent(out) :: this
    character(*), intent(in)  :: name
    class(grid), target, intent(in) :: g
    character(*), intent(in) :: unit
    integer, intent(in) :: idx_type
    integer, intent(in) :: idx_dir

    call this%variable_init(name, unit, g, idx_type, idx_dir)

    this%d0 => this%data%get_ptr0()
    this%d2 => this%data%get_ptr2()
  end subroutine

  subroutine test_container()
    type(test_case) :: tc
    type(container) :: depot, topedo
    
    integer :: i, j, idx(2)
    real, allocatable :: x(:), y(:)
    type(grid0D) :: g0, g0_readback
    type(grid1D) :: gx, gy, gx_read, gy_read
    type(var)    :: var0, vart, vart2, vart3, vart4
    type(vselector) :: vsel, vsel_edge
    type(tensor_grid) :: gt, gt_read
    type(grid_table) :: tab, tab_read, tab_edge1, tab_edge2

    call tc%init("container")

    ! Save to file
    call depot%open("container.fbs", STORAGE_NEW)

    call g0%init("0D-test-grid")
    call gx%init("x", linspace(-1.0, 2.0, 31), -10)
    call gy%init("y", linspace(0.0, 10.0, 101))
    call gt%init("test-tensor", [gx%get_ptr(), gy%get_ptr()])

    call tab%init("tab", gt, IDX_CELL, 0, initial_flags = .true.)
    call tab%init_final()

    call depot%save(g0)
    call depot%save(gx)
    call depot%save(gy)
    call depot%save(gt)

    call var0%init("0var", "1", g0, IDX_VERTEX, 0)
    call var0%set([1.0])

    call vart%init("tensor_var", "1", gt, IDX_CELL, 0)
    do i = 1, tab%n
      idx = tab%get_idx(i)
      call vart%set(idx, 9.9)
    end do
    call vart%set([1,1], 7.7)

    call vart2%init("vart_edge2", "V", gt, IDX_EDGE, 1)
    call vart3%init("vart_edge3", "1/cm^3", gt, IDX_EDGE, 1)
    call vart4%init("vart_edge4", "1", gt, IDX_EDGE, 1)
    call tab_edge1%init("TABE1", gt, IDX_EDGE, 1, initial_flags=.true.)
    call tab_edge2%init("TABE2", gt, IDX_EDGE, 1)
    
    do i = -10, 19
      do j = 1, 101
        if (i < 1 .and. j < 21) call tab_edge1%flags%set([i, j], .false.)
        if (i < 1 .and. j < 10) call tab_edge2%flags%set([i, j], .true.)

        call vart2%set([i, j], norm(1.0, "V"))
        call vart3%set([i, j], norm(3.3e17, "1/cm^3"))
        call vart4%set([i, j], norm(-0.752, "1"))
      end do
    end do
    call tab_edge1%init_final()
    call tab_edge2%init_final()
    call vsel_edge%init([vart2%get_ptr(), vart3%get_ptr(), vart4%get_ptr()], [tab_edge1%get_ptr(), tab_edge2%get_ptr()], "edgy")

    call depot%save(var0, "vars")
    call depot%save(vart, "vars")
    call depot%save(tab)

    call vsel%init(vart, name="vselll", tab=tab)
    call vart%set([1,1], 9.9)
    call depot%save(vsel, "sim")
    call depot%save(vsel_edge, "full")

    call depot%close()
    
    call var0%set([3.3])
    do i = 1, tab%n
      idx = tab%get_idx(i)
      call vart%set(idx, 0.0)
    end do

    allocate(x(vsel_edge%n), source=0.0)
    y = vsel_edge%get()
    call vsel_edge%set(x)

    ! Load from file
    call topedo%open("container.fbs", STORAGE_UNPACK)

    call topedo%load("0D-test-grid", g0_readback)
    call topedo%load("x", gx_read)
    call topedo%load("y", gy_read)
    call topedo%load("test-tensor", gt_read)
    
    call assert_eq_grid(g0, g0_readback)
    call assert_eq_grid(gx, gx_read)
    call assert_eq_grid(gy, gy_read)
    call assert_eq_grid(gt, gt_read)
    
    call tab_read%init("tab", gt_read, IDX_CELL, 0)
    call topedo%load("tab", tab_read)
    call tc%assert(all(tab%flags%get() .eqv. tab_read%flags%get()), "Grid tables do not match")

    call topedo%load("vars/0var", var0)
    call tc%assert_eq([1.0], var0%get(), 1e-16, 1e-16, "Variable values do not match")
    call topedo%load("vars/tensor_var", vart)
    call tc%assert_eq(7.7, vart%get([1, 1]), 1e-16, 1e-16, "Variable values do not match")
    call vart%set([1,1], 9.9)
    call tc%assert(all(vart%get() == 9.9), "Variable tensor_var values do not match")
    
    call topedo%load(new_string("sim/vselll"), vsel)
    call tc%assert(all(vsel%get() == 9.9), "Vselector values do not match")

    deallocate(x)
    call topedo%read(new_string("full/edgy"), x)
    call tc%assert_eq(x(1), 1.0, 1e-15, 1e-15, "vart2 in vsel not as expected")
    call tc%assert_eq(x(2), 3.3e17, 1e-15, 1e-15, "vart3 in vsel not as expected")
    call tc%assert_eq(x(3), -0.752, 1e-15, 1e-15, "vart4 in vsel not as expected")

    call topedo%load(new_string("full/edgy"), vsel_edge)
    x = vsel_edge%get()
    call tc%assert_eq(y, x, 1e-15, 1e-15, "Some edgy values do not match expectations")

    call topedo%close()
    call tc%finish()

  contains
    subroutine assert_eq_grid(expected, value)
      class(grid),  intent(in) :: expected
      class(grid),  intent(in) :: value

      integer :: idx_bnd(2,expected%idx_dim), idx_bnd_(2,value%idx_dim)
      integer :: idx(expected%idx_dim), i, j
      real :: evertex(expected%dim), vvertex(expected%dim)

      call tc%assert_eq(new_string(expected%name), new_string(value%name), "The grid names do not match")
      call tc%assert_eq(expected%idx_dim, value%idx_dim, "The grid index dimensions do not match")
      call tc%assert_eq(expected%dim,     value%dim, "The grid dimensions do not match")
      
      call expected%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
      call value%get_idx_bnd(   IDX_VERTEX, 0, idx_bnd_)
      call tc%assert_eq(idx_bnd, idx_bnd_, "The index bounds do not match")

      ! collect indices for which flag is set
      if (expected%idx_dim > 0) then
        idx = idx_bnd(1,:)
        i   = 0
        do while (idx(expected%idx_dim) <= idx_bnd(2,expected%idx_dim))
          
          call expected%get_vertex(idx, evertex) 
          call value%get_vertex(   idx, vvertex) 
          call tc%assert_eq(evertex, vvertex, 1e-16, 1e-16, "Grid vertices do not match")

          ! go to next index
          idx(1) = idx(1) + 1
          do j = 1, expected%idx_dim-1
            if (idx(j) <= idx_bnd(2,j)) exit
            idx(j  ) = idx_bnd(1,j)
            idx(j+1) = idx(j+1) + 1
          end do
        end do
      end if

    end subroutine
  end subroutine

end module