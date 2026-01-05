module beam_generation_m
  !! Point beam generation for STEM-EBIC simulation
  !!
  !! Simple model: G = I_beam / tr_vol at the beam vertex
  !! The beam is a point that can sweep along y.

  use bin_search_m,    only: bin_search
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
    !! Beam generation rate variable
    !! Unit: 1/cm^3/s

    integer :: ci
      !! carrier index (CR_ELEC or CR_HOLE)

    real, pointer :: x1(:)     => null()
    real, pointer :: x2(:,:)   => null()
    real, pointer :: x3(:,:,:) => null()
  contains
    procedure :: init => beam_generation_init
  end type

  type, extends(variable) :: beam_position
    !! Beam y-position for STEM-EBIC sweep
    !! Unit: cm (normalized)

    real, pointer :: x => null()
  contains
    procedure :: init => beam_position_init
  end type

  type, extends(equation) :: calc_beam_generation
    !! Calculate G = I_beam / tr_vol at beam position

    type(device_params), pointer :: par => null()
    type(beam_generation), pointer :: bgen => null()
    type(beam_position), pointer :: beam_pos => null()

    integer :: ci
      !! carrier index
    integer :: ix_beam
      !! beam x-index (fixed)
    integer :: iy_min
      !! beam y-range minimum index
    integer :: iy_max
      !! beam y-range maximum index
  contains
    procedure :: init => calc_beam_generation_init
    procedure :: eval => calc_beam_generation_eval
  end type

contains

  subroutine beam_generation_init(this, par, ci)
    class(beam_generation), intent(out) :: this
    type(device_params),    intent(in)  :: par
    integer,                intent(in)  :: ci

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    call this%variable_init("bgen_"//CR_NAME(ci), "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

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
    class(beam_position), intent(out) :: this
    character(*),         intent(in)  :: name

    type(grid_data0_real), pointer :: p

    call this%variable_init(name, "nm")
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

  subroutine calc_beam_generation_init(this, par, bgen, beam_pos)
    class(calc_beam_generation),    intent(out) :: this
    type(device_params),    target, intent(in)  :: par
    type(beam_generation),  target, intent(in)  :: bgen
    type(beam_position),    target, intent(in)  :: beam_pos

    integer :: iprov
    real    :: x_actual, y_min_actual, y_max_actual

    print "(A)", "calc_beam_generation_init"

    call this%equation_init("calc_"//bgen%name)
    this%par      => par
    this%bgen     => bgen
    this%beam_pos => beam_pos
    this%ci       = bgen%ci

    ! Find beam x-index (fixed for all y positions)
    this%ix_beam = bin_search(par%g1D(1)%x, par%reg_beam(1)%beam_x)
    this%iy_min  = bin_search(par%g1D(2)%x, par%reg_beam(1)%beam_y_min)
    this%iy_max  = bin_search(par%g1D(2)%x, par%reg_beam(1)%beam_y_max)

    ! Get actual grid positions
    x_actual     = par%g1D(1)%x(this%ix_beam)
    y_min_actual = par%g1D(2)%x(this%iy_min)
    y_max_actual = par%g1D(2)%x(this%iy_max)

    ! Print alignment info (first carrier only)
    if (this%ci == CR_ELEC) then
      print "(A,I0,A,ES12.4,A)", "  ix_beam = ", this%ix_beam, " (x = ", denorm(x_actual, 'nm'), " nm)"
      print "(A,I0,A,I0)", "  iy range: ", this%iy_min, " to ", this%iy_max
    end if

    ! Provide beam generation on transport vertices
    iprov = this%provide(bgen, par%transport(IDX_VERTEX, 0))

    call this%init_final()
  end subroutine

  subroutine calc_beam_generation_eval(this)
    !! Evaluate G = I_beam / (tr_vol * curr_fact) at beam position
    class(calc_beam_generation), intent(inout) :: this

    integer              :: i, iy_beam
    integer, allocatable :: idx(:)
    real                 :: I_beam, tr_vol, G, beam_y, y_actual

    allocate(idx(this%par%g%idx_dim))

    ! Get beam current and y-position
    I_beam = this%par%reg_beam(1)%I_beam
    beam_y = this%beam_pos%x

    ! Find nearest y-index
    iy_beam  = bin_search(this%par%g1D(2)%x, beam_y)
    y_actual = this%par%g1D(2)%x(iy_beam)

    ! Check bounds
    if (iy_beam < this%iy_min .or. iy_beam > this%iy_max) then
      if (this%ci == CR_ELEC) then
        print "(A,I0,A,I0,A,I0)", "  WARNING: beam iy=", iy_beam, " outside range [", this%iy_min, ",", this%iy_max, "]"
      end if
    end if

    ! Get tr_vol at beam vertex and calculate G
    idx    = [this%ix_beam, iy_beam]
    tr_vol = this%par%tr_vol%get(idx)
    G      = I_beam / (tr_vol * this%par%curr_fact)

    ! Debug output (first carrier only)
    if (this%ci == CR_ELEC) then
      print "(A,ES12.4,A,ES12.4,A,I0,A)", "  beam_y = ", denorm(beam_y, 'nm'), &
        " -> ", denorm(y_actual, 'nm'), " nm (iy=", iy_beam, ")"
      print "(A,ES12.4,A,ES12.4)", "  tr_vol = ", denorm(tr_vol, 'nm^2'), &
        ", G = ", denorm(G, '1/cm^3/s')
    end if

    ! Set G at beam vertex, 0 elsewhere
    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)
      if (idx(1) == this%ix_beam .and. idx(2) == iy_beam) then
        call this%bgen%set(idx, G)
      else
        call this%bgen%set(idx, 0.0)
      end if
    end do
  end subroutine

end module
