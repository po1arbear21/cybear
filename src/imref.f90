module imref_m

  use density_m,       only: density
  use device_params_m, only: device_params
  use distributions_m, only: fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use semiconductor_m, only: CR_NAME, CR_CHARGE
  use variable_m,      only: variable_real
  use error_m,         only: assert_failed, program_error

  implicit none

  private
  public :: imref, calc_imref, calc_density

  type, extends(variable_real) :: imref
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

  type, extends(equation) :: calc_imref
    !! iref = pot +/- V_T * log(dens / n_intrin) or equivalent formula for degenerate case

    type(device_params), pointer :: par  => null()
    type(potential),     pointer :: pot  => null()
    type(density),       pointer :: dens => null()
    type(imref),         pointer :: iref => null()

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_dens => null()
  contains
    procedure :: init => calc_imref_init
    procedure :: eval => calc_imref_eval
  end type

  type, extends(equation) :: calc_density
    !! dens = n_intrin * exp(+/- (iref - pot) / V_T) or equivalent formula for degenerate case

    type(device_params), pointer :: par  => null()
    type(potential),     pointer :: pot  => null()
    type(density),       pointer :: dens => null()
    type(imref),         pointer :: iref => null()

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_iref => null()
  contains
    procedure :: init => calc_density_init
    procedure :: eval => calc_density_eval
  end type

contains

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

  subroutine calc_imref_init(this, par, pot, dens, iref)
    !! initialize imref calculation equation
    class(calc_imref),           intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(potential),     target, intent(in)  :: pot
      !! potential variable
    type(density),       target, intent(in)  :: dens
      !! density variable
    type(imref),         target, intent(in)  :: iref
      !! imref variable

    integer               :: i, iprov
    integer, allocatable  :: idx(:)

    ! init equation
    call this%equation_init("calc_"//iref%name)
    this%par  => par
    this%pot  => pot
    this%dens => dens
    this%iref => iref

    ! provides imref
    iprov = this%provide(iref, par%transport(IDX_VERTEX,0))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX,0)), const = .true.)

    ! depend on density
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX,0)), const = .false.)

    ! set jaco_pot entries (constant)
    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx = par%transport(IDX_VERTEX,0)%get_idx(i)
      call this%jaco_pot%set(idx, idx, 1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_imref_eval(this)
    !! evaluate imref calculation equation
    class(calc_imref), intent(inout) :: this

    integer              :: i, ci
    real                 :: ch, pot, dens, iref, IF12, dIF12
    integer, allocatable :: idx(:)

    ci = this%iref%ci
    ch = CR_CHARGE(this%iref%ci)

    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate imref
      pot  = this%pot%get(idx)
      dens = this%dens%get(idx)
      if (this%par%smc%degen) then
        call inv_fermi_dirac_integral_1h(dens / this%par%smc%edos(ci), IF12, dIF12)
        iref = pot - this%par%smc%band_edge(ci) + ch * IF12

        call this%iref%set(idx, iref)
        call this%jaco_dens%set(idx, idx, ch * dIF12 / this%par%smc%edos(ci))
      else
        iref = pot + ch * log(dens / this%par%smc%n_intrin)

        call this%iref%set(idx, iref)
        call this%jaco_dens%set(idx, idx, ch / dens)
      end if
    end do
  end subroutine

  subroutine calc_density_init(this, par, pot, dens, iref)
    !! initialize density calculation equation
    class(calc_density),         intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(potential),     target, intent(in)  :: pot
      !! potential variable
    type(density),       target, intent(in)  :: dens
      !! density variable
    type(imref),         target, intent(in)  :: iref
      !! imref variable

    integer :: iprov

    ! init equation
    call this%equation_init("calc_"//dens%name)
    this%par  => par
    this%pot  => pot
    this%dens => dens
    this%iref => iref

    ! provides density
    iprov = this%provide(dens, par%transport(IDX_VERTEX,0))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX,0)), const = .false.)

    ! depend on imref
    this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%transport(IDX_VERTEX,0)), const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_density_eval(this)
    !! evaluate density calculation equation
    class(calc_density), intent(inout) :: this

    integer               :: i, ci
    real                  :: ch, pot, dens, iref, F12, dF12
    integer, allocatable  :: idx(:)

    ci = this%iref%ci
    ch = CR_CHARGE(ci)
    allocate (idx(this%par%g%idx_dim))

    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate density
      pot  = this%pot%get(idx)
      iref = this%iref%get(idx)
      if (this%par%smc%degen) then
        call fermi_dirac_integral_1h(ch * (iref - pot + this%par%smc%band_edge(ci)), F12, dF12)
        dens = this%par%smc%edos(ci) * F12

        call this%dens%set(idx, dens)
        call this%jaco_pot%set( idx, idx, - ch * this%par%smc%edos(ci) * dF12)
        call this%jaco_iref%set(idx, idx,   ch * this%par%smc%edos(ci) * dF12)
      else
        dens = this%par%smc%n_intrin * exp(ch * (iref - pot))

        call this%dens%set(idx, dens)
        call this%jaco_pot%set( idx, idx, - ch*dens)
        call this%jaco_iref%set(idx, idx,   ch*dens)
      end if
    end do
  end subroutine

end module
