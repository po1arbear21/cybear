module beam_generation_m
  !! External beam generation for STEM-EBIC simulation
  !!
  !! This module provides carrier generation due to electron beam excitation.
  !! The beam generates electron-hole pairs at a specified location (point or line).

  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use semiconductor_m, only: CR_NAME, CR_ELEC, CR_HOLE
  use variable_m,      only: variable_real

  implicit none

  private
  public beam_generation, calc_beam_generation

  type, extends(variable_real) :: beam_generation
    !! External beam generation rate variable
    !!
    !! Unit: 1/cm^3/s (generation rate per unit volume)
    !! Applied equally to both electrons and holes

    integer :: ci
      !! carrier index (CR_ELEC or CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for 1D grids
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for 2D grids
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for 3D grids
  contains
    procedure :: init => beam_generation_init
  end type

  type, extends(equation) :: calc_beam_generation
    !! Calculate external beam generation rate
    !!
    !! Sets generation rate based on beam position and shape.
    !! Currently supports line generation (constant G along a vertical line).

    type(device_params), pointer :: par => null()
      !! device parameters

    type(beam_generation), pointer :: bgen => null()
      !! beam generation variable to update

    integer :: ci
      !! carrier index

  contains
    procedure :: init => calc_beam_generation_init
    procedure :: eval => calc_beam_generation_eval
  end type

contains

  subroutine beam_generation_init(this, par, ci)
    !! Initialize beam generation variable
    class(beam_generation), intent(out) :: this
    type(device_params),    intent(in)  :: par
      !! device parameters
    integer,                intent(in)  :: ci
      !! carrier index (CR_ELEC or CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init variable base
    call this%variable_init("bgen_"//CR_NAME(ci), "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

    ! get pointer to data for direct access
    select case (par%g%idx_dim)
    case (1)
      p1 => this%data%get_ptr1()
      this%x1 => p1%data
    case (2)
      p2 => this%data%get_ptr2()
      this%x2 => p2%data
    case (3)
      p3 => this%data%get_ptr3()
      this%x3 => p3%data
    end select
  end subroutine

  subroutine calc_beam_generation_init(this, par, bgen)
    !! Initialize beam generation calculation equation
    class(calc_beam_generation),    intent(out) :: this
    type(device_params),    target, intent(in)  :: par
      !! device parameters
    type(beam_generation),  target, intent(in)  :: bgen
      !! beam generation variable

    integer :: iprov

    print "(A)", "calc_beam_generation_init"

    ! init equation
    call this%equation_init("calc_"//bgen%name)
    this%par  => par
    this%bgen => bgen
    this%ci   = bgen%ci

    ! provides beam generation on transport vertices
    iprov = this%provide(bgen, par%transport(IDX_VERTEX, 0))

    ! No dependencies - beam generation is externally specified
    ! (no jacobians needed since G doesn't depend on solution)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_beam_generation_eval(this)
    !! Evaluate beam generation - copy from device_params to variable
    class(calc_beam_generation), intent(inout) :: this

    integer              :: i
    integer, allocatable :: idx(:)
    real                 :: G

    allocate (idx(this%par%g%idx_dim))

    ! Copy beam generation rate from device_params to variable
    ! The generation rate is pre-computed in device_params based on beam region
    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
      G = this%par%beam_rate%get(idx)
      call this%bgen%set(idx, G)
    end do
  end subroutine

end module
