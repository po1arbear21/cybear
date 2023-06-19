module ionization_m

  use density_m,       only: density
  use device_params_m, only: device_params
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
    !! ionization ratio

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
    !! calculate the stationary dopant ionization ratio

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
    !! continuity equation for dopant ionization ratio: d/dt ion + (G - R) / Ndop = 0

    type(device_params), pointer :: par => null()

    integer :: ci

    type(vselector) :: ion
    type(vselector) :: genrec

    type(dirichlet_stencil) :: st_dir

    type(jacobian), pointer :: jaco_t      => null()
    type(jacobian), pointer :: jaco_ion    => null()
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

    real :: tau_gen
      !! Time constant of the generation process
    real :: tau_rec
      !! Time constant of the recombination process

    type(device_params), pointer :: par => null()

    type(generation_recombination), pointer :: genrec => null()
      !! main variable recombination

    type(potential),  pointer :: pot  => null()
      !! potential variable
    type(imref),      pointer :: iref => null()
      !! imref variable
    type(density),    pointer :: dens => null()
      !! eletron/hole density variable
    type(ionization), pointer :: ion  => null()
      !! donor/acceptor ionization ratio

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_iref => null()
    type(jacobian), pointer :: jaco_dens => null()
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
    call this%variable_init("ion"//DOP_NAME(ci), "1", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
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
    iprov = this%provide(ion, par%dopvert(ion%ci))

    ! depends on potential and imref if incomplete ionization is enabled
    if (par%smc%incomp_ion) then
      this%jaco_pot  => this%init_jaco(iprov, this%depend(pot,  par%dopvert(ion%ci)), const = .false.)
      this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%dopvert(ion%ci)), const = .false.)
    end if

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_ionization_eval(this)
    !! evaluate ionization calculation equation
    class(calc_ionization), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: b, g, ch, edop, e, f, ion, dion, iref, pot

    ci = this%ion%ci
    ch = CR_CHARGE(ci)
    allocate (idx(this%par%g%idx_dim))

    ! full ionization
    if (.not. this%par%smc%incomp_ion) then
      do i = 1, this%par%dopvert(ci)%n
        idx = this%par%dopvert(ci)%get_idx(i)
        call this%ion%set(idx, 1.0)
      end do
      return
    end if

    g  = this%par%smc%g_dop(ci)
    do i = 1, this%par%dopvert(ci)%n
      idx  = this%par%dopvert(ci)%get_idx(i)

      b    = this%par%asb(ci)%get(idx)
      edop = this%par%smc%band_edge(ci) + ch * this%par%edop(ci)%get(idx)
      iref = this%iref%get(idx)
      pot  = this%pot%get(idx)

      ! exponent
      e = - ch * (iref - pot + edop)

      ! factor (exp(e) can overflow)
      f = 1.0 / (1.0 + g * exp(e))

      ! ionization and derivative
      ion  = 1.0 - b * f
      dion = b * f * (1.0 - f)

      ! save
      call this%ion%set(idx, ion)
      call this%jaco_pot%set(idx, idx,    ch * dion)
      call this%jaco_iref%set(idx, idx, - ch * dion)
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

    ci = ion%ci

    ! init base
    call this%equation_init(ion%name//"_contin")
    this%par => par
    this%ci = ci

    ! init variable selectors
    call this%ion%init(ion, par%dopvert(ci))
    call this%genrec%init(genrec, par%dopvert(ci))

    ! init residuals using this%ion as main variable
    call this%init_f(this%ion)

    ! init dirichlet stencil
    call this%st_dir%init(par%g)

    ! init jacobians
    if (par%smc%incomp_ion) then
      this%jaco_t      => this%init_jaco_f(this%depend(this%ion),    st = [this%st_dir%get_ptr()], const = .true., dtime = .true. )
      this%jaco_genrec => this%init_jaco_f(this%depend(this%genrec), st = [this%st_dir%get_ptr()], const = .true., dtime = .false.)
    else
      this%jaco_ion => this%init_jaco_f(this%depend(this%ion), st = [this%st_dir%get_ptr()], const = .true., dtime = .false.)
    end if

    ! set jacobian entries
    allocate (idx(par%g%idx_dim))
    do i = 1, par%dopvert(ci)%n
      idx = par%dopvert(ci)%get_idx(i)
      if (par%smc%incomp_ion) then
        call this%jaco_t%set(idx, idx, par%dop(IDX_VERTEX,0,ci)%get(idx))
        call this%jaco_genrec%set(idx, idx, -1.0)
      else
        call this%jaco_ion%set(idx, idx, 1.0)
      end if
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine ion_continuity_eval(this)
    !! evaluate ionization continuity equation
    class(ion_continuity), intent(inout) :: this

    real, allocatable :: tmp(:)

    allocate (tmp(this%f%n))
    if (this%par%smc%incomp_ion) then
      call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp)
    else
      call this%jaco_ion%matr%mul_vec(this%ion%get(), tmp)
      tmp = tmp - 1.0
    end if
    call this%f%set(tmp)
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

  subroutine calc_generation_recombination_init(this, par, tau, genrec, pot, iref, dens, ion)
    !! initialize generation-recombination equation
    class(calc_generation_recombination),   intent(inout) :: this
    type(device_params),            target, intent(in)    :: par
      !! device parameters
    real,                                   intent(in)    :: tau(2)
      !! time constant for generation (1) and recombination (2)
    type(generation_recombination), target, intent(in)    :: genrec
      !! recombination variable
    type(potential),                target, intent(in)    :: pot
      !! potential variable
    type(imref),                    target, intent(in)    :: iref
      !! imref variable
    type(density),                  target, intent(in)    :: dens
      !! eletron/hole density variable
    type(ionization),               target, intent(in)    :: ion
      !! ionization ratio

    integer :: ci, iprov

    ! init equation and variable pointers
    call this%equation_init("calc_"//genrec%name)
    this%tau_gen =  tau(1)
    this%tau_rec =  tau(2)
    this%par     => par
    this%genrec  => genrec
    this%pot     => pot
    this%iref    => iref
    this%dens    => dens
    this%ion     => ion

    ci    = genrec%ci
    iprov = this%provide(genrec, par%dopvert(ci))

    ! dependencies
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%dopvert(ci)), const = .false.)
    this%jaco_pot  => this%init_jaco(iprov, this%depend(pot,  par%dopvert(ci)), const = .false.)
    this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%dopvert(ci)), const = .false.)
    this%jaco_ion  => this%init_jaco(iprov, this%depend(ion,  par%dopvert(ci)), const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_generation_recombination_eval(this)
    !! evaluate generation-recombination equation
    class(calc_generation_recombination), intent(inout) :: this

    real, parameter :: alpha = sqrt(0.125)

    integer              :: ci, i
    integer, allocatable :: idx(:)
    logical              :: degen
    real                 :: ch, dens, dop, edop, edos, g, ion, iref, pot, rec, t
    real                 :: dgenddens, dgendpot, dgendiref, dgendion
    type(hp_real)        :: e, eta, fac, gen

    ci   = this%genrec%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)
    g    = this%par%smc%g_dop(ci)
    degen = this%par%smc%degen

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%dopvert(ci)%n
      idx = this%par%dopvert(ci)%get_idx(i)

      ! doping and doping energy level
      dop  = this%par%dop(IDX_VERTEX,0,ci)%get(idx)
      edop = this%par%smc%band_edge(ci) + ch * this%par%edop(ci)%get(idx)

      ! ionization ratio
      ion  = this%ion%get( idx)

      ! generation
      dgenddens = 0
      dgendpot  = 0
      dgendiref = 0
      if (degen) then
        dens = this%dens%get(idx)
        pot  = this%pot%get( idx)
        iref = this%iref%get(idx)

        e   = TwoSum(pot, -iref)
        eta = - ch * (e - this%par%smc%band_edge(ci))
        t   = hp_to_real(eta)
        if (t >= -17) then
          fac       = exp(ch * (e - edop))
          gen       = fac * dens / edos
          dgenddens = hp_to_real(fac / edos)
          dgendpot  = hp_to_real(  ch * gen)
          dgendiref = hp_to_real(- ch * gen)
        else
          gen = real_to_hp(exp(- ch * (edop - this%par%smc%band_edge(ci))))
          if (t >= -36) then
            ! correction for -36 < eta < -17
            gen       =                  gen / (1.0 + exp( eta) * alpha)
            dgendpot  = hp_to_real( ch * gen / (1.0 + exp(-eta) / alpha))
            dgendiref = hp_to_real(-ch * gen / (1.0 + exp(-eta) / alpha))
          end if
        end if
      else
        gen = real_to_hp(exp(- ch * (edop - this%par%smc%band_edge(ci))))
      end if

      dgendion  = - hp_to_real(gen)       * dop / this%tau_gen
      gen       = gen * TwoSum(1.0, -ion) * dop / this%tau_gen
      dgenddens = dgenddens * (1.0 - ion) * dop / this%tau_gen
      dgendpot  = dgendpot  * (1.0 - ion) * dop / this%tau_gen
      dgendiref = dgendiref * (1.0 - ion) * dop / this%tau_gen

      ! recombination
      rec = dens * ion * dop / (edos * g * this%tau_rec)

      ! combined rate
      gen       = gen - rec
      dgenddens = dgenddens -  ion * dop / (edos * g * this%tau_rec)
      dgendion  = dgendion  - dens * dop / (edos * g * this%tau_rec)

      ! save
      call this%genrec%set(idx, hp_to_real(gen))
      call this%jaco_dens%set(idx, idx, dgenddens)
      call this%jaco_pot%set( idx, idx, dgendpot)
      call this%jaco_iref%set(idx, idx, dgendiref)
      call this%jaco_ion%set( idx, idx, dgendion)
    end do
  end subroutine

end module
