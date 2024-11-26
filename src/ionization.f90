module ionization_m

  use density_m,       only: density
  use device_params_m, only: device_params
  use fermi_m,         only: fermi_dirac_integral_1h, fermi_dirac_generation_reg
  use equation_m,      only: equation
  use error_m,         only: program_error
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use high_precision_m
  use imref_m,         only: imref
  use ieee_arithmetic, only: ieee_is_finite
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use res_equation_m,  only: res_equation
  use semiconductor_m, only: DOP_NAME, CR_CHARGE, CR_ELEC, CR_HOLE
  use stencil_m,       only: dirichlet_stencil, empty_stencil
  use variable_m,      only: variable_real
  use vselector_m,     only: vselector

  implicit none

  private
  public ionization, calc_ionization, ion_continuity
  public generation_recombination, calc_generation_recombination

  type, extends(variable_real) :: ionization
    !! ionization concentration

    integer :: ci
      !! dopant index (DOP_DCON, DOP_ACON)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => ionization_init
  end type

  type, extends(equation) :: calc_ionization
    !! calculate the stationary dopant ionization concentration

    type(device_params), pointer :: par  => null()
    type(ionization),    pointer :: ion  => null()
    type(potential),     pointer :: pot  => null()
    type(imref),         pointer :: iref => null()

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_iref => null()
  contains
    procedure :: init => calc_ionization_init
    procedure :: eval => calc_ionization_eval
  end type

  type, extends(res_equation) :: ion_continuity
    !! continuity equation for dopant ionization concentration: d/dt ion - (G - R) = 0

    type(device_params), pointer :: par => null()

    integer :: ci

    type(vselector) :: ion
    type(vselector) :: genrec

    type(dirichlet_stencil) :: st_dir

    type(jacobian), pointer :: jaco_t      => null()
    type(jacobian), pointer :: jaco_genrec => null()
  contains
    procedure :: init => ion_continuity_init
    procedure :: eval => ion_continuity_eval
  end type

  type, extends(variable_real) :: generation_recombination
    !! generation - recombination rate

    integer      :: ci
      !! carrier index
    character(2) :: kind
      !! name of rate: "DC" or "AV"
  contains
    procedure :: init => generation_recombination_init
  end type

  type, extends(equation) :: calc_generation_recombination
    !! calculate G - R

    real :: tau
      !! time constant of the generation/recombination process

    type(device_params), pointer :: par => null()

    type(generation_recombination), pointer :: genrec => null()
      !! main variable recombination

    type(potential),  pointer :: pot  => null()
      !! potential variable
    type(imref),      pointer :: iref => null()
      !! imref variable
    type(ionization), pointer :: ion  => null()
      !! donor/acceptor ionization ratio

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_iref => null()
    type(jacobian), pointer :: jaco_ion  => null()
  contains
    procedure :: init => calc_generation_recombination_init
    procedure :: eval => calc_generation_recombination_eval
  end type

contains

  subroutine ionization_init(this, par, ci)
    !! initialize ionization ratio
    class(ionization),   intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: ci
      !! carrier/dopant index (DOP_ELEC, DOP_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("ion"//DOP_NAME(ci), "1/cm^3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
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

  subroutine calc_ionization_init(this, par, pot, ion, iref)
    !! initialize ionization ratio calculation equation
    class(calc_ionization),      intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(ionization),    target, intent(in)  :: ion
      !! donor/acceptor ionization ratio variable
    type(potential),     target, intent(in)  :: pot
      !! potential variable
    type(imref),         target, intent(in)  :: iref
      !! imref variable

    integer :: iprov

    ! init equation
    call this%equation_init("calc_"//ion%name)
    this%par  => par
    this%ion  => ion
    this%pot  => pot
    this%iref => iref

    ! provides ionization
    iprov = this%provide(ion, par%ionvert(ion%ci))

    ! depends on potential and imref if incomplete ionization is enabled
    if (par%smc%incomp_ion) then
      this%jaco_pot  => this%init_jaco(iprov, this%depend(pot,  par%ionvert(ion%ci)), const = .false.)
      this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%ionvert(ion%ci)), const = .false.)
    end if

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_ionization_eval(this)
    !! evaluate ionization calculation equation
    class(calc_ionization), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ch, dop, Edop, ii_g, iref, pot, eta, tmp, ion, dion

    allocate (idx(this%par%g%idx_dim))

    ci   = this%ion%ci
    ch   = CR_CHARGE(ci)
    ii_g = this%par%smc%ii_g(ci) ! degeneracy factor (2 or 4)

    do i = 1, this%par%ionvert(ci)%n
      idx = this%par%ionvert(ci)%get_idx(i)

      dop  = this%par%dop(IDX_VERTEX,0,ci)%get(idx)
      Edop = this%par%ii_E_dop(ci)%get(idx) ! E_C - E_D or E_A - E_V , should be > 0
      iref = this%iref%get(idx)
      pot  = this%pot%get(idx)
      eta  = - ch * (pot - iref - this%par%smc%band_edge(ci))

      if (eta < 300.0) then
        tmp  = ii_g * exp(Edop) * exp(eta)
        ion  = dop / (1.0 + tmp)
        dion = - dop * tmp / (1.0 + tmp)**2
      else
        ion  = 0
        dion = 0
      end if
      call this%ion%set(idx, ion)
      call this%jaco_pot%set( idx, idx, - ch * dion)
      call this%jaco_iref%set(idx, idx,   ch * dion)
    end do
  end subroutine

  subroutine ion_continuity_init(this, par, ion, genrec)
    !! initialize dopant residual equation
    class(ion_continuity),                  intent(out) :: this
    type(device_params),            target, intent(in)  :: par
      !! device parameters
    type(ionization),               target, intent(in)  :: ion
      !! dopant ionization density variable
    type(generation_recombination), target, intent(in)  :: genrec
      !! recombination variable

    integer              :: ci, i
    integer, allocatable :: idx(:)

    if (.not. par%smc%incomp_ion) call program_error("use only if incomplete ionization is turned on")

    ci = ion%ci

    ! init base
    call this%equation_init(ion%name//"_contin")
    this%par => par
    this%ci = ci

    ! init variable selectors
    call this%ion%init(ion, par%ionvert(ci))
    call this%genrec%init(genrec, par%ionvert(ci))

    ! init residuals using this%ion as main variable
    call this%init_f(this%ion)

    ! init dirichlet stencil
    call this%st_dir%init(par%g)

    ! init jacobians
    this%jaco_t      => this%init_jaco_f(this%depend(this%ion),    st = [this%st_dir%get_ptr()], const = .true., dtime = .true. )
    this%jaco_genrec => this%init_jaco_f(this%depend(this%genrec), st = [this%st_dir%get_ptr()], const = .true., dtime = .false.)

    ! set jacobian entries
    allocate (idx(par%g%idx_dim))
    do i = 1, par%ionvert(ci)%n
      idx = par%ionvert(ci)%get_idx(i)
      call this%jaco_t%set(     idx, idx,  1.0)
      call this%jaco_genrec%set(idx, idx, -1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine ion_continuity_eval(this)
    !! evaluate ionization continuity equation
    class(ion_continuity), intent(inout) :: this

    call this%f%set(-1.0*this%genrec%get())
  end subroutine

  subroutine generation_recombination_init(this, par, ci)
    !! initialize generation-recombination variable
    class(generation_recombination), intent(inout) :: this
    type(device_params),             intent(in)    :: par
      !! device parameters
    integer,                         intent(in)    :: ci
      !! carrier index

    character(2) :: kind

    ! kind string
    select case (ci)
    case (CR_ELEC)
      kind = "DC"
    case (CR_HOLE)
      kind = "AV"
    end select

    ! init variable base
    call this%variable_init("genrec_"//kind, "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci   = ci
    this%kind = kind
  end subroutine

  subroutine calc_generation_recombination_init(this, par, tau, genrec, pot, iref, ion)
    !! initialize generation-recombination equation
    class(calc_generation_recombination),   intent(inout) :: this
    type(device_params),            target, intent(in)    :: par
      !! device parameters
    real,                                   intent(in)    :: tau
      !! time constant
    type(generation_recombination), target, intent(in)    :: genrec
      !! recombination variable
    type(potential),                target, intent(in)    :: pot
      !! potential variable
    type(imref),                    target, intent(in)    :: iref
      !! imref variable
    type(ionization),               target, intent(in)    :: ion
      !! ionization ratio

    integer :: ci, iprov

    ! init equation and variable pointers
    call this%equation_init("calc_"//genrec%name)
    this%tau    =  tau
    this%par    => par
    this%genrec => genrec
    this%pot    => pot
    this%iref   => iref
    this%ion    => ion

    ci    = genrec%ci
    iprov = this%provide(genrec, par%ionvert(ci))

    ! dependencies
    this%jaco_pot  => this%init_jaco(iprov, this%depend(pot,  par%ionvert(ci)), const = .false.)
    this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%ionvert(ci)), const = .false.)
    this%jaco_ion  => this%init_jaco(iprov, this%depend(ion,  par%ionvert(ci)), const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_generation_recombination_eval(this)
    !! evaluate generation-recombination equation
    class(calc_generation_recombination), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ch, dop, Edop, ii_g, iref, pot, ion, eta, f, df, g, dg
    real                 :: gen, dgen, rec, drec, genrec, dgenrecdeta, dgenrecdion

    ci    = this%genrec%ci
    ch    = CR_CHARGE(ci)
    ii_g  = this%par%smc%ii_g(ci) ! degeneracy factor (2 or 4)

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%ionvert(ci)%n
      idx = this%par%ionvert(ci)%get_idx(i)

      dop  = this%par%dop(IDX_VERTEX,0,ci)%get(idx)
      Edop = this%par%ii_E_dop(ci)%get(idx) ! E_C - E_D or E_A - E_V , should be > 0
      iref = this%iref%get(idx)
      pot  = this%pot%get(idx)
      ion  = this%ion%get(idx)
      eta  = - ch * (pot - iref - this%par%smc%band_edge(ci))

      ! generation
      if (this%par%smc%degen) then
        call fermi_dirac_integral_1h(eta, f, df) ! use non-regularized value in rec for correct stationary ionization rate
        call fermi_dirac_generation_reg(eta, g, dg)
        gen  = exp(- Edop) * g
        dgen = exp(- Edop) * dg
      else
        f    = exp(eta)
        df   = f
        gen  = exp(- Edop)
        dgen = 0
      end if

      ! recombination
      rec  = ii_g * f
      drec = ii_g * df

      ! combined rate
      genrec      = hp_to_real((gen * dop - TwoSum(gen, rec) * ion) / this%tau)
      dgenrecdeta = (dgen * dop - (dgen + drec) * ion) / this%tau
      dgenrecdion = - (gen + rec) / this%tau

      ! save
      call this%genrec%set(idx, genrec)
      call this%jaco_pot%set( idx, idx, -ch * dgenrecdeta)
      call this%jaco_iref%set(idx, idx,  ch * dgenrecdeta)
      call this%jaco_ion%set( idx, idx,       dgenrecdion)
    end do
  end subroutine

end module
