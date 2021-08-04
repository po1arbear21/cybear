module example_device_m

  use bin_search_m,    only: bin_search
  use grid_m,          only: IDX_VERTEX, IDX_CELL, grid_data1_real
  use grid1D_m,        only: grid1D
  use input_m,         only: input_file
  use math_m,          only: linspace
  use normalization_m, only: init_normconst

  implicit none

  private
  public init_device
  public adj_v, n_intrin, grd, eps, dop, dop_v

  real                  :: n_intrin
  type(grid1D)          :: grd
  type(grid_data1_real) :: adj_v, dop, dop_v, eps


contains
  subroutine init_device(f)
    type(input_file), intent(in) :: f
      !! input file for initialisation

    call init_globals(f)
    call init_grid(f)
    call init_adj_v()
    call init_permittivity(f)
    call init_doping(f)

  end subroutine

  subroutine init_globals(f)
    type(input_file), intent(in) :: f
    !! input file for initialisation

    real :: T

    call f%get("", "temperature", T, normalize = .false.)
    call init_normconst(T)
    call f%get("", "n_intrin", n_intrin)
  end subroutine

  subroutine init_grid(f)
    type(input_file), intent(in) :: f
    !! input file for initialisation

    integer           :: nx
    real              :: x0, x1
    real, allocatable :: x(:)

    ! getting input for the grid
    call f%get("grid", "Nx", nx)
    call f%get("grid", "x0", x0)
    call f%get("grid", "x1", x1)

    x = linspace(x0, x1, nx)
    call grd%init(x)
  end subroutine

  subroutine init_adj_v()
    integer :: i

    ! initialise the grid_data for the adj_v
    call adj_v%init(grd, IDX_VERTEX, 0)

    ! set values in the given range
    do i = 1, size(grd%x)-1
      call adj_v%update([i],   grd%get_vol([i])/2)
      call adj_v%update([i+1], grd%get_vol([i])/2)
    end do
  end subroutine

  subroutine init_permittivity(f)
    type(input_file), intent(in) :: f
    !! input file for initialisation

    integer :: i, ind0, ind1
    real :: x0, x1, value

    ! initialise the grid_data for the permittivity
    call eps%init(grd, IDX_CELL, 0)

    ! getting input for the permittivity
    call f%get("permittivity", "eps", value)
    call f%get("permittivity", "x0",  x0)
    call f%get("permittivity", "x1",  x1)


    ! set values in the given range
    ind0 = bin_search(grd%x, x0)
    ind1 = bin_search(grd%x, x1) - 1
    do i = ind0, ind1
      call eps%set([i], value)
    end do
  end subroutine

  subroutine init_doping(f)
    type(input_file), intent(in) :: f
    !! input file for initialisation

    integer              :: ind0, ind1, i, j
    integer, allocatable :: sid(:)
    real                 :: x0, x1, value

    ! initialise the grid_data for the doping and dop_v
    call dop%init(grd, IDX_CELL, 0)
    call dop_v%init(grd, IDX_VERTEX, 0)

    ! get different sections of doping
    call f%get_sections("doping", sid)

    ! getting input for the doping at each section
    do i = 1, size(sid)
      call f%get(sid(i), "dcon", value)
      call f%get(sid(i), "x0",   x0)
      call f%get(sid(i), "x1",   x1)

      ! set values in the given range
      ind0 = bin_search(grd%x, x0)
      ind1 = bin_search(grd%x, x1) - 1
      do j = ind0, ind1
        call dop%set([j], value)
      end do
    end do

    ! transfering the doping from the cells to the vertices
    do i = 1, dop%n
      call dop_v%update([i],   0.5*grd%get_vol([i])*dop%get([i])/adj_v%get([i]))
      call dop_v%update([i+1], 0.5*grd%get_vol([i])*dop%get([i])/adj_v%get([i+1]))
    end do
  end subroutine

end module
