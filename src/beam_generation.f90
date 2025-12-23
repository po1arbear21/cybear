module beam_generation_m
  !! External beam generation for STEM-EBIC simulation
  !!
  !! This module provides carrier generation due to electron beam excitation.
  !! The beam generates electron-hole pairs at a specified location (point or line).

  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data0_real, grid_data1_real, grid_data2_real, grid_data3_real
  use semiconductor_m, only: CR_NAME, CR_ELEC, CR_HOLE
  use variable_m,      only: variable
  use normalization_m, only: denorm

  implicit none

  private
  public beam_generation, beam_position, calc_beam_generation

  type, extends(variable) :: beam_generation
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

  type, extends(variable) :: beam_position
    !! Beam position variable for STEM-EBIC sweep
    !!
    !! Represents the current beam position (y-coordinate) during a linescan.
    !! Unit: cm (normalized internally)

    real, pointer :: x => null()
      !! pointer to position value (normalized to cm)
  contains
    procedure :: init => beam_position_init
  end type

  type, extends(equation) :: calc_beam_generation
    !! Calculate external beam generation rate
    !!
    !! Sets generation rate based on beam position and shape.
    !! Supports dynamic beam position for EBIC linescan sweeps.

    type(device_params), pointer :: par => null()
      !! device parameters

    type(beam_generation), pointer :: bgen => null()
      !! beam generation variable to update

    type(beam_position), pointer :: beam_pos => null()
      !! beam position variable (for sweep mode)

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

  subroutine beam_position_init(this, name)
    !! Initialize beam position variable for STEM-EBIC sweep
    class(beam_position), intent(out) :: this
    character(*),         intent(in)  :: name
      !! variable name (e.g., "Y_BEAM")

    type(grid_data0_real), pointer :: p

    ! init base with length unit (matches input format: Y_BEAM = ... : nm)
    call this%variable_init(name, "nm")

    ! set pointer to value for direct access
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

  subroutine calc_beam_generation_init(this, par, bgen, beam_pos)
    !! Initialize beam generation calculation equation
    class(calc_beam_generation),    intent(out) :: this
    type(device_params),    target, intent(in)  :: par
      !! device parameters
    type(beam_generation),  target, intent(in)  :: bgen
      !! beam generation variable
    type(beam_position), optional, target, intent(in) :: beam_pos
      !! beam position variable for sweep mode (optional)

    integer :: iprov

    print "(A)", "calc_beam_generation_init"

    ! init equation
    call this%equation_init("calc_"//bgen%name)
    this%par  => par
    this%bgen => bgen
    this%ci   = bgen%ci

    ! store beam position pointer if provided (for sweep mode)
    if (present(beam_pos)) then
      this%beam_pos => beam_pos
    end if

    ! provides beam generation on transport vertices
    iprov = this%provide(bgen, par%transport(IDX_VERTEX, 0))

    ! No dependencies - beam generation is externally specified
    ! (no jacobians needed since G doesn't depend on solution)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_beam_generation_eval(this)
    !! Evaluate beam generation
    !!
    !! If beam_pos is set (sweep mode): compute generation dynamically
    !! based on current beam position.
    !! Otherwise (static mode): copy pre-computed rates from device_params.
    class(calc_beam_generation), intent(inout) :: this

    integer              :: i, n_gen
    integer, allocatable :: idx(:)
    real                 :: G, point(3)
    real                 :: beam_x, beam_y, y_min, y_max, G0
    real                 :: mesh_tol_x, mesh_tol_y

    allocate (idx(this%par%g%idx_dim))

    ! Static mode: copy from device_params (original behavior)
    if (.not. associated(this%beam_pos)) then
      do i = 1, this%par%transport(IDX_VERTEX, 0)%n
        idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
        G = this%par%beam_rate%get(idx)
        call this%bgen%set(idx, G)
      end do
      return
    end if

    ! Dynamic mode: compute generation based on current beam position
    ! Get current beam y-position from sweep input
    beam_y = this%beam_pos%x

    ! Get fixed beam x-position from device file (center of x-range)
    beam_x = (this%par%reg_beam(1)%xyz(1,1) + this%par%reg_beam(1)%xyz(1,2)) / 2.0

    ! Get y-bounds from device file
    y_min = this%par%reg_beam(1)%xyz(2,1)
    y_max = this%par%reg_beam(1)%xyz(2,2)

    ! Get generation rate from device file
    G0 = this%par%reg_beam(1)%G0

    ! Use mesh spacing for tolerance (half grid spacing)
    ! Fallback to 1e-5 (~10nm normalized) if max_dx/max_dy not set
    mesh_tol_x = this%par%max_dx / 2.0
    mesh_tol_y = this%par%max_dy / 2.0
    if (mesh_tol_x < 1e-10) mesh_tol_x = 1e-5
    if (mesh_tol_y < 1e-10) mesh_tol_y = 1e-5

    ! Debug: print tolerance values once
    if (this%ci == 1) then
      print "(A,ES12.4,A,ES12.4)", "  max_dx=", this%par%max_dx, " max_dy=", this%par%max_dy
      print "(A,ES12.4,A,ES12.4)", "  tol_x=", mesh_tol_x, " tol_y=", mesh_tol_y
    end if

    ! Compute generation for each transport vertex
    n_gen = 0
    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
      call this%par%g%get_vertex(idx, point(1:this%par%g%dim))

      G = 0.0
      ! Check if vertex is at beam position (within mesh tolerance)
      if (abs(point(1) - beam_x) < mesh_tol_x) then
        if (abs(point(2) - beam_y) < mesh_tol_y) then
          if (point(2) >= y_min .and. point(2) <= y_max) then
            G = G0
            n_gen = n_gen + 1
          end if
        end if
      end if

      call this%bgen%set(idx, G)
    end do

    ! Debug: print once per eval (first carrier only)
    if (this%ci == 1) then
      print "(A,ES12.4,A,ES12.4,A,I0,A)", "  beam @ (", denorm(beam_x,'nm'), ",", denorm(beam_y,'nm'), ") -> ", n_gen, " vertices"
    end if
  end subroutine

end module
