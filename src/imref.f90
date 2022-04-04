module imref_m

  use density_m,       only: density
  use device_params_m, only: device_params, CR_NAME, CR_CHARGE
  use equation_m,      only: equation
  use grid_m,          only: grid_data2_real, IDX_VERTEX
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use variable_m,      only: variable_real

  implicit none

  private
  public :: imref, calc_imref, calc_density

  type, extends(variable_real) :: imref
    !! quasi-fermi potential
    integer       :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => imref_init
  end type

  type, extends(equation) :: calc_imref
    !! iref = pot +/- V_T * log(dens / n_intrin)

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
    !! dens = n_intrin * exp(+/- (iref - pot) / V_T)

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

    type(grid_data2_real), pointer :: p

    ! init base
    call this%variable_init(CR_NAME(ci)//"imref", "V", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
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

    integer :: i, idx(2), iprov

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
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX,0)), const = .true.)

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

    integer :: i, idx(2)
    real    :: ch, pot, dens, iref

    ch = CR_CHARGE(this%iref%ci)

    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate imref
      pot  = this%pot%get(idx)
      dens = this%dens%get(idx)
      iref = pot + ch * log(dens / this%par%n_intrin)
      call this%iref%set(idx, iref)

      ! set jaco_dens entries
      call this%jaco_dens%set(idx, idx, ch / dens)
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

    integer :: i, idx(2)
    real    :: ch, pot, dens, iref

    ch = CR_CHARGE(this%iref%ci)

    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! calculate density
      pot  = this%pot%get(idx)
      iref = this%iref%get(idx)
      dens = this%par%n_intrin * exp(ch * (iref - pot))
      call this%dens%set(idx, dens)

      ! set jacobian entries
      call this%jaco_pot%set( idx, idx, -ch*dens)
      call this%jaco_iref%set(idx, idx,  ch*dens)
    end do
  end subroutine

end module
