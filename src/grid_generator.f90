module grid_generator_m

  use error_m,       only: program_error
  use grid_m,        only: grid, grid_ptr
  use grid1D_m,      only: grid1D
  use input_m,       only: input_file
  use tensor_grid_m, only: tensor_grid
  use triang_grid_m, only: triang_grid
  use triangle_m,    only: triangulation
  use region_m,      only: region, region_ptr, region_poisson, region_transport, region_doping, region_contact
  use vector_m,      only: vector_real
  use math_m,        only: linspace
  use qsort_m,       only: qsort

  implicit none

  private
  public  :: DIR_NAME, generate_1D_grid, generate_triangle_grid

   ! parameters
  character(*), parameter :: DIR_NAME(3)  = [ "x", "y", "z"]

contains

  subroutine generate_1D_grid(file, load, dir, reg, g1D, gptr, ngptr)
    !! create 1D grid x, y, or z
    type(input_file),     intent(in)    :: file
      !! device file
    logical,              intent(in)    :: load
      !! load or generate?
    integer,              intent(in)    :: dir
      !! axis direction
    type(region_ptr),     intent(in)    :: reg(:)
      !! pointers to all regions
    type(grid1D), target, intent(out)   :: g1D
      !! output newly generated 1D grid
    type(grid_ptr),       intent(inout) :: gptr(:)
      !! save pointer to g1D in gptr (used to create tensor grid afterwards)
    integer,              intent(inout) :: ngptr
      !! number of already used gptr entries (gets increased by 1)

    real              :: max_dx
    real, allocatable :: x0(:), x(:)

    ! load or generate axis
    if (load) then
      call file%get("grid", DIR_NAME(dir), x)
    else
      call file%get("grid", "max_d"//DIR_NAME(dir), max_dx)
      call get_bounds(dir, reg, x0)
      call generate_axis(x0, max_dx, x)
    end if

    ! initialize grid
    call g1D%init(DIR_NAME(dir), x)

    ! update pointer list
    ngptr         =  ngptr + 1
    gptr(ngptr)%p => g1D
  end subroutine

  subroutine generate_triangle_grid(file, load, reg, tr, gtr, gptr, ngptr)
    !! create triangle grid
    type(input_file),          intent(in)    :: file
      !! device file
    logical,                   intent(in)    :: load
      !! load or generate? (load not implemented yet)
    type(region_ptr),          intent(in)    :: reg(:)
      !! pointers to all regions
    type(triangulation),       intent(out)   :: tr
      !! output triangulation
    type(triang_grid), target, intent(out)   :: gtr
      !! output corresponding triangle grid
    type(grid_ptr),            intent(inout) :: gptr(:)
      !! save pointer to g1D in gptr (used to create tensor grid afterwards)
    integer,                   intent(inout) :: ngptr
      !! number of already used gptr entries (gets increased by 1)

    integer :: i
    real    :: max_area

    ! load or generate triangle grid
    if (load) then
      call program_error("loading triangle grid not implemented yet")
    else
      call tr%init()
      do i = 1, size(reg)
        call tr%add_polygon(reg(i)%p%xyz(1,:), reg(i)%p%xyz(2,:), closed = .true.)
      end do
      call file%get("grid", "max_area", max_area)
      call tr%triangulate(max_area = max_area)
      call tr%get_grid("tr", gtr)
    end if

    ! update pointer list
    ngptr = ngptr + 1
    gptr(ngptr)%p => gtr
  end subroutine

  subroutine get_bounds(dim, reg, x0)
    integer,           intent(in)  :: dim
    type(region_ptr),  intent(in)  :: reg(:)
    real, allocatable, intent(out) :: x0(:)

    integer :: i, j

    ! collect bounds of all regions in array
    allocate (x0(2*size(reg)))
    j = 0
    do i = 1, size(reg)
      x0(j+1:j+2) = reg(i)%p%xyz(dim,1:2)
      j = j + 2
    end do
  end subroutine

  subroutine generate_axis(x0, dx, x)
    real,              intent(in)  :: x0(:)
    real,              intent(in)  :: dx
    real, allocatable, intent(out) :: x(:)

    real, parameter :: TOL = 1e-10

    real, allocatable :: xtmp(:), xx(:)
    integer           :: n, m, j, i
    type(vector_real) :: xv

    ! sort region bounds
    allocate (xtmp(size(x0)))
    xtmp = x0
    call qsort(xtmp)
    n = size(xtmp)

    ! remove duplicates
    j = 1
    do i = 2 , n
      if (abs(xtmp(i)-xtmp(j)) > TOL) then
        j = j + 1
        xtmp(j) = xtmp(i)
      end if
    end do
    n = j

    ! concat linearly spaced arrays together
    call xv%init(0, c = 2 * ceiling((xtmp(n) - xtmp(1)) / dx))
    call xv%push(xtmp(1))
    do i = 1, n-1
      m  = ceiling((xtmp(i+1)-xtmp(i))/dx) + 1
      if (allocated(xx)) deallocate (xx)
      xx = linspace(xtmp(i), xtmp(i+1), m)
      call xv%push(xx(2:m))
    end do
    x= xv%to_array()
  end subroutine

end module
