module imref_m

  use density_m,       only: density
  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use semiconductor_m, only: CR_NAME, CR_CHARGE
  use variable_m,      only: variable
  use error_m,         only: program_error

  implicit none

  private
  public :: chemical_pot, imref, calc_chemical_pot_dens, calc_chemical_pot_imref, calc_imref, calc_density

  type, extends(variable) :: chemical_pot
    !! chemical potential

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => chemical_pot_init
  end type

  type, extends(variable) :: imref
    !! quasi-fermi potential

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => imref_init
  end type

  type, extends(equation) :: calc_chemical_pot_dens
    !! calculate chemical potential from density: eta = F_inv(dens / edos)

    type(device_params), pointer :: par  => null()
    type(density),       pointer :: dens => null()
    type(chemical_pot),  pointer :: eta  => null()

    type(jacobian), pointer :: jaco_dens => null()
  contains
    procedure :: init => calc_chemical_pot_dens_init
    procedure :: eval => calc_chemical_pot_dens_eval
  end type

  type, extends(equation) :: calc_imref
    !! calculate imref from chemical potential: iref = pot - band_edge/e +/- eta

    type(device_params), pointer :: par  => null()
    type(chemical_pot),  pointer :: eta  => null()
    type(potential),     pointer :: pot  => null()
    type(imref),         pointer :: iref => null()

    type(jacobian), pointer :: jaco_eta => null()
    type(jacobian), pointer :: jaco_pot => null()
  contains
    procedure :: init => calc_imref_init
    procedure :: eval => calc_imref_eval
  end type

  type, extends(equation) :: calc_chemical_pot_imref
    !! calculate chemical potential from imref: eta = +/- (pot - iref - band_edge)

    type(device_params), pointer :: par  => null()
    type(potential),     pointer :: pot  => null()
    type(imref),         pointer :: iref => null()
    type(chemical_pot),  pointer :: eta  => null()

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_iref => null()
  contains
    procedure :: init => calc_chemical_pot_imref_init
    procedure :: eval => calc_chemical_pot_imref_eval
  end type

  type, extends(equation) :: calc_density
    !! calculate density from chemical potential: dens = edos * F(eta)

    type(device_params), pointer :: par  => null()
    type(chemical_pot),  pointer :: eta  => null()
    type(density),       pointer :: dens => null()

    type(jacobian), pointer :: jaco_eta  => null()
  contains
    procedure :: init => calc_density_init
    procedure :: eval => calc_density_eval
  end type

contains

  subroutine chemical_pot_init(this, par, ci)
    !! initialize eta
    class(chemical_pot), intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("eta_"//CR_NAME(ci), "1", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

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

  subroutine imref_init(this, par, ci)
    !! initialize imref
    class(imref),        intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init(CR_NAME(ci)//"imref", "V", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

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

  subroutine calc_chemical_pot_dens_init(this, par, dens, eta)
    !! initialize eta calculation equation
    class(calc_chemical_pot_dens), intent(out) :: this
    type(device_params), target,   intent(in)  :: par
      !! device parameters
    type(density),       target,   intent(in)  :: dens
      !! density variable
    type(chemical_pot),  target,   intent(in)  :: eta
      !! chemical potential variable

    integer :: iprov

    print "(A)", "calc_chemical_pot_from_dens_init"

    ! init equation
    call this%equation_init("calc_"//eta%name//"_from_"//dens%name)
    this%par  => par
    this%dens => dens
    this%eta  => eta

    ! provides imref
    iprov = this%provide(eta, par%transport(IDX_VERTEX,0))

    ! depend on density
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX,0)), const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_chemical_pot_dens_eval(this)
    !! evaluate eta calculation equation
    class(calc_chemical_pot_dens), intent(inout) :: this

    integer              :: i
    integer, allocatable :: idx(:)
    real                 :: dens, eta, detadF, edos_loc

    allocate (idx(this%par%g%idx_dim))
    !$omp parallel do default(none) schedule(dynamic) private(i,idx,dens,eta,detadF,edos_loc) shared(this)
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate eta with per-vertex edos so n = dens/edos is consistent across heterojunctions
      dens     = this%dens%get(idx)
      edos_loc = this%par%edos_v(this%eta%ci)%get(idx)
      call this%par%smc(this%par%smc_default)%get_inv_dist(dens / edos_loc, eta, detadF)
      call this%eta%set(idx, eta)
      call this%jaco_dens%add(idx, idx, detadF / edos_loc)
    end do
    !$omp end parallel do
  end subroutine

  subroutine calc_imref_init(this, par, eta, pot, iref)
    !! initialize imref calculation equation
    class(calc_imref),           intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(chemical_pot),  target, intent(in)  :: eta
      !! chemical potential variable
    type(potential),     target, intent(in)  :: pot
      !! potential variable
    type(imref),         target, intent(in)  :: iref
      !! imref variable

    integer              :: i, iprov
    integer, allocatable :: idx(:)

    print "(A)", "calc_imref_init"

    ! init equation
    call this%equation_init("calc_"//iref%name)
    this%par  => par
    this%eta  => eta
    this%pot  => pot
    this%iref => iref

    ! provides imref
    iprov = this%provide(iref, par%transport(IDX_VERTEX,0))

    ! depend on chemical potential
    this%jaco_eta => this%init_jaco(iprov, this%depend(eta, par%transport(IDX_VERTEX,0)), const = .true.)

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX,0)), const = .true.)

    ! set jacobian entries (constant)
    allocate (idx(par%g%idx_dim))
    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx = par%transport(IDX_VERTEX,0)%get_idx(i)
      call this%jaco_pot%add(idx, idx, 1.0)
      call this%jaco_eta%add(idx, idx, CR_CHARGE(this%iref%ci))
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_imref_eval(this)
    !! evaluate imref calculation equation
    class(calc_imref), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ch

    ci = this%iref%ci
    ch = CR_CHARGE(this%iref%ci)

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate imref using per-vertex band edge
      call this%iref%set(idx, this%pot%get(idx) - this%par%band_edge_v(ci)%get(idx) + ch * this%eta%get(idx))
    end do
  end subroutine

  subroutine calc_chemical_pot_imref_init(this, par, pot, iref, eta)
    !! initialize eta calculation equation
    class(calc_chemical_pot_imref), intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    type(potential),     target,    intent(in)  :: pot
      !! potential variable
    type(imref),         target,    intent(in)  :: iref
      !! imref variable
    type(chemical_pot),  target,    intent(in)  :: eta
      !! chemical potential variable

    integer              :: i, iprov
    integer, allocatable :: idx(:)

    print "(A)", "calc_chemical_pot_from_imref_init"

    ! init equation
    call this%equation_init("calc_"//eta%name//"_from_"//iref%name)
    this%par  => par
    this%pot  => pot
    this%iref => iref
    this%eta  => eta

    ! provides imref
    iprov = this%provide(eta, par%transport(IDX_VERTEX,0))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX,0)), const = .true.)

    ! depend on imref
    this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%transport(IDX_VERTEX,0)), const = .true.)

    ! set jacobian entries (constant)
    allocate (idx(par%g%idx_dim))
    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx = par%transport(IDX_VERTEX,0)%get_idx(i)
      call this%jaco_pot%add( idx, idx, -CR_CHARGE(this%iref%ci))
      call this%jaco_iref%add(idx, idx,  CR_CHARGE(this%iref%ci))
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_chemical_pot_imref_eval(this)
    !! evaluate eta calculation equation
    class(calc_chemical_pot_imref), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ch

    ci = this%iref%ci
    ch = CR_CHARGE(ci)

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate eta from imref and potential using per-vertex band edge
      call this%eta%set(idx, -ch * (this%pot%get(idx) - this%iref%get(idx) - this%par%band_edge_v(ci)%get(idx)))
    end do
  end subroutine

  subroutine calc_density_init(this, par, eta, dens)
    !! initialize density calculation equation
    class(calc_density),         intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(chemical_pot),  target, intent(in)  :: eta
      !! chemical potential variable
    type(density),       target, intent(in)  :: dens
      !! density variable

    integer :: iprov

    print "(A)", "calc_density_init"

    ! init equation
    call this%equation_init("calc_"//dens%name)
    this%par  => par
    this%eta  => eta
    this%dens => dens

    ! provides density
    iprov = this%provide(dens, par%transport(IDX_VERTEX,0))

    ! depend on chemical potential
    this%jaco_eta => this%init_jaco(iprov, this%depend(eta, par%transport(IDX_VERTEX,0)), const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_density_eval(this)
    !! evaluate density calculation equation
    class(calc_density), intent(inout) :: this

    integer              :: i
    real                 :: F, dFdeta, edos_loc
    integer, allocatable :: idx(:)

    allocate (idx(this%par%g%idx_dim))
    !$omp parallel do default(none) schedule(dynamic) private(i,idx,F,dFdeta,edos_loc) shared(this)
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate density with per-vertex edos
      edos_loc = this%par%edos_v(this%eta%ci)%get(idx)
      call this%par%smc(this%par%smc_default)%get_dist(this%eta%get(idx), 0, F, dFdeta)
      call this%dens%set(idx, edos_loc * F)
      call this%jaco_eta%add(idx, idx, edos_loc * dFdeta)
    end do
    !$omp end parallel do
  end subroutine

end module
