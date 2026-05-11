module electric_field_m

  use device_params_m,  only: device_params
  use equation_m,       only: equation
  use error_m,          only: program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m, only: DIR_NAME
  use jacobian_m,       only: jacobian, jacobian_ptr
  use potential_m,      only: potential
  use stencil_m,        only: near_neighb_stencil, dirichlet_stencil
  use variable_m,       only: variable
  use vselector_m,      only: vselector

  implicit none

  private
  public electric_field, electric_field_abs, calc_efield, calc_efield_abs

  type, extends(variable) :: electric_field
    !! directional electric field strength on vertices

    integer       :: dir
      !! direction

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => electric_field_init
  end type

  type, extends(variable) :: electric_field_abs
    !! absolute electrical field strength

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => electric_field_abs_init
  end type

  type, extends(equation) :: calc_efield
    !! calculate electric field component on vertices

    type(device_params), pointer :: par => null()
      !! device parameters
    type(vselector)              :: efield
      !! directional electric field
    type(vselector)              :: pot
      !! potential

    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot => null()

    real, allocatable :: work(:)
      !! pre-allocated work array for evaluation
  contains
    procedure :: init => calc_efield_init
    procedure :: eval => calc_efield_eval
  end type

  type, extends(equation) :: calc_efield_abs
    !! calculate electric field on vertices

    type(device_params), pointer :: par => null()
      !! device parameters
    type(vselector), allocatable :: efield(:)
      !! directional electric fields
    type(vselector)              :: efield_abs
      !! absolute electric field

    type(dirichlet_stencil) :: st_dir

    ! TODO: extend jacobian support so it can carry jacobians for the
    ! magnitude as well as for individual directions (currently only
    ! d|E|/dE_dir per direction).
    type(jacobian_ptr), allocatable :: jaco_efield(:)
  contains
    procedure :: init => calc_efield_abs_init
    procedure :: eval => calc_efield_abs_eval
  end type
contains

  subroutine electric_field_init(this, par, dir)
    !! initialize electric field component
    class(electric_field), intent(out) :: this
    type(device_params),   intent(in)  :: par
      !! device parameters
    integer,               intent(in)  :: dir
      !! direction

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("electric_field_"//DIR_NAME(dir), "V/m", g = par%g, idx_type = IDX_VERTEX, idx_dir=0)
    this%dir = dir

    ! get pointer to data
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
    case default
      call program_error("Maximal 3 dimensions allowed")
    end select
  end subroutine

  subroutine electric_field_abs_init(this, par)
    !! initialize absolute electric field
    class(electric_field_abs), intent(out) :: this
    type(device_params),       intent(in)  :: par
      !! device parameters

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("electric_field_abs", "V/m", g = par%g, idx_type = IDX_VERTEX, idx_dir=0)

    ! get pointer to data
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
    case default
      call program_error("Maximal 3 dimensions allowed")
    end select
  end subroutine

  subroutine calc_efield_init(this, par, efield, pot)
    class(calc_efield),           intent(out) :: this
    type(device_params),  target, intent(in)  :: par
      !! device parameters
    type(electric_field),         intent(in)  :: efield
      !! electric field component variable
    type(potential),              intent(in)  :: pot
      !! potential variable

    integer              :: dim, dir, i, idx_dim, idx_dir, iprov
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    logical              :: status
    real                 :: dx, edge_len, dEFdPot, surf
    real,    allocatable :: p1(:), p2(:)

    print "(A)", "calc_efield_init"

    ! init base
    call this%equation_init("calc_"//efield%name)
    this%par => par

    ! create variable selectors
    call this%efield%init(efield, par%transport(IDX_VERTEX, 0))
    call this%pot%init(pot, par%transport(IDX_VERTEX, 0))

    ! init stencil
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)

    ! provide electric field
    iprov = this%provide(efield, par%transport(IDX_VERTEX, 0))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()], const = .true.)

    ! calculate jacobian
    idx_dim = par%g%idx_dim
    dir = efield%dir
    dim = par%g%dim

    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), p1(dim), p2(dim))

    if ((this%par%gtype%s == "x") .or. (this%par%gtype%s == "xy") .or. (this%par%gtype%s == "xyz")) then
      idx_dir = dir
      ! loop over edges
      do i = 1, this%par%transport(IDX_EDGE, dir)%n
        idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)

        surf = this%par%tr_surf(idx_dir)%get(idx)
        edge_len = this%par%g%get_len(idx, idx_dir)

        ! get vertices and distance
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)
        call this%par%g%get_vertex(idx1, p1)
        call this%par%g%get_vertex(idx2, p2)

        dx = p2(dir) - p1(dir)

        ! add values to jacobian
        dEFdPot =  (surf * dx) / (2 * this%par%tr_vol%get(idx1) * edge_len)
        call this%jaco_pot%add(idx1, idx1, dEFdPot)
        call this%jaco_pot%add(idx1, idx2, -dEFdPot)

        dEFdPot = (surf * dx) / (2 * this%par%tr_vol%get(idx2) * edge_len)
        call this%jaco_pot%add(idx2, idx2, -dEFdPot)
        call this%jaco_pot%add(idx2, idx1, dEFdPot)
      end do
    else
      ! TODO: need to loop differently if dir does not match idx_dir
      call program_error("electric field not implemented for triangular grids")
    end if

    ! pre-allocate work array for evaluation
    allocate(this%work(this%efield%n))

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_efield_abs_init(this, par, efield, efield_abs)
    class(calc_efield_abs),      intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(electric_field),        intent(in)  :: efield(:)
      !! electric field component variable
    type(electric_field_abs),    intent(in)  :: efield_abs
      !! electric field variable

    integer :: iprov, dir

    print "(A)", "calc_efield_abs_init"

    ! init base
    call this%equation_init("calc_"//efield_abs%name)
    this%par => par

    ! init variable selectors
    allocate (this%efield(par%g%dim))
    call this%efield_abs%init(efield_abs, par%transport(IDX_VERTEX, 0))
    do dir = 1, par%g%dim
      call this%efield(dir)%init(efield(dir), par%transport(IDX_VERTEX, 0))
    end do

    ! init stencils
    call this%st_dir%init(par%g)

    ! provide electric field
    iprov = this%provide(efield_abs, par%transport(IDX_VERTEX, 0))

    ! depend on electric field components
    allocate (this%jaco_efield(par%g%dim))
    do dir = 1, par%g%dim
      this%jaco_efield(dir)%p => this%init_jaco(iprov, this%depend(efield(dir), par%transport(IDX_VERTEX, 0)), st = [this%st_dir%get_ptr()])
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_efield_eval(this)
    class(calc_efield), intent(inout) :: this

    call this%jaco_pot%matr%mul_vec(this%pot%get(), this%work)
    call this%efield%set(this%work)
  end subroutine

  subroutine calc_efield_abs_eval(this)
    class(calc_efield_abs), intent(inout) :: this

    integer              :: i, dir
    integer, allocatable :: idx(:)
    real                 :: fieldsum(1), efield(1)

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)

      fieldsum(1) = 0.0
      ! sum all directional field squares
      do dir = 1, this%par%g%dim
        fieldsum = fieldsum + this%efield(dir)%get(idx)**2
      end do

      ! add minimum for numeric stability
      fieldsum = fieldsum + this%par%smc%ii_ef_min**2

      ! get total field strength
      fieldsum = sqrt(fieldsum)
      call this%efield_abs%set(idx, fieldsum)

      ! set jacobians
      do dir = 1, this%par%g%dim
        efield = this%efield(dir)%get(idx)
        call this%jaco_efield(dir)%p%add(idx, idx, efield(1)/fieldsum(1))
      end do
    end do
  end subroutine

end module
