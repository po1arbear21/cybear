module grid_generator_m

  use grid_m,        only: grid, grid_ptr
  use grid1D_m,      only: grid1D
  use tensor_grid_m, only: tensor_grid
  use triang_grid_m, only: triang_grid
  use triangle_m,    only: triangulation
  use region_m,      only: region, region_ptr, region_poisson, region_transport, region_doping, region_contact
  use vector_m,      only: vector_real
  use math_m,        only: linspace
  use qsort_m,       only: qsort

  implicit none

  private
  public  :: DIR_NAME, generate_cartesian_grid, generate_triangle_grid

   ! parameters
  character(*), parameter :: DIR_NAME(3)  = [ "x", "y", "z"]

contains

  subroutine generate_cartesian_grid(dim, reg, max_dxyz, g1D, tg, g)
    !! create grids "x", "xy", "xyz"
    integer,                    intent(in)  :: dim
      !! grid dimension
    type(region_ptr),           intent(in)  :: reg(:)
      !! pointers to all regions
    real,                       intent(in)  :: max_dxyz(:)
      !! [max_dx, max_dy, max_dz]
    type(grid1D),      target,  intent(out) :: g1D(:)
      !! output x, y, z grid
    type(tensor_grid), target,  intent(out) :: tg
      !! output tensor grid
    class(grid),       pointer, intent(out) :: g
      !! output pointer to tg

    integer           :: i
    real, allocatable :: x0(:), xyz(:)
    type(grid_ptr)    :: gptr(3)

    do i = 1, dim
      ! get bounds
      call get_bounds(i, reg, x0)
      if (allocated(xyz)) deallocate(xyz)
      call generate_axis(x0, max_dxyz(i), xyz)

      call g1D(i)%init(DIR_NAME(i), xyz)
      gptr(i) = g1D(i)%get_ptr()
    end do
    if (dim > 1) then
      call tg%init("grid", gptr(1:dim))
      g => tg
    else
      g => g1D(1)
    end if
  end subroutine

  subroutine generate_triangle_grid(dim, reg, max_areadz, tr, gtr, g1D, tg, g)
    !! create triangle grid "tr_xy", "tr_xyz"
    integer,                    intent(in)  :: dim
      !! grid dimension
    type(region_ptr),           intent(in)  :: reg(:)
      !! region pointer
    real,                       intent(in)  :: max_areadz(:)
      !! [max_area, max_dz]
    type(triangulation),        intent(out) :: tr
      !! output triangulation
    type(triang_grid), target,  intent(out) :: gtr
      !! output corresponding triangle grid
    type(grid1D),      target,  intent(out) :: g1D(:)
      !! output z grid
    type(tensor_grid), target,  intent(out) :: tg
      !! output tensor grid
    class(grid),       pointer, intent(out) :: g
      !! output pointer to gtr or tg

    integer           :: i
    real, allocatable :: z0(:), z(:)
    type(grid_ptr)    :: gptr(2)

    ! triangle grid
    call tr%init()
    do i = 1, size(reg)
      call tr%add_polygon(reg(i)%p%xyz(1,:), reg(i)%p%xyz(2,:), closed = .true.)
    end do
    call tr%triangulate(max_area = max_areadz(1))
    call tr%get_grid("tr", gtr)
    g => gtr
    gptr(1) = gtr%get_ptr()

    ! add z axis if gtype == "tr_xyz"
    if (dim == 3) then
      call get_bounds(dim, reg, z0)
      call generate_axis(z0, max_areadz(2), z)
      call g1D(3)%init(DIR_NAME(3), z)
      gptr(2) = g1D(3)%get_ptr()

      ! tensor grid
      call tg%init("grid", gptr(1:2))
      g => tg
    end if
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
