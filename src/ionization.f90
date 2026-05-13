module ionization_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  use device_params_m,  only: device_params
  use fukushima_m,      only: fd1h, fdm1h
  use electric_field_m, only: electric_field_abs
  use equation_m,       only: equation
  use error_m,          only: program_error
  use grid_m,           only: IDX_VERTEX
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use imref_m,          only: chemical_pot
  use jacobian_m,       only: jacobian
  use math_m,           only: PI
  use res_equation_m,   only: res_equation
  use semiconductor_m,  only: DOP_NAME, RATE_NAME, CR_ELEC, CR_HOLE, DOS_PARABOLIC, DIST_MAXWELL, DIST_FERMI
  use stencil_m,        only: dirichlet_stencil
  use variable_m,       only: variable
  use vselector_m,      only: vselector

  implicit none

  private
  public ionization, calc_ionization, ion_continuity
  public generation_recombination, calc_generation_recombination

  type, extends(variable) :: ionization
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

    type(device_params), pointer :: par => null()
    type(chemical_pot),  pointer :: eta => null()
    type(ionization),    pointer :: ion => null()

    type(jacobian), pointer :: jaco_eta => null()
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

  type, extends(variable) :: generation_recombination
    !! generation - recombination rate

    integer      :: ci
      !! carrier index
  contains
    procedure :: init => generation_recombination_init
  end type

  type, extends(equation) :: calc_generation_recombination
    !! calculate G - R

    type(device_params),            pointer :: par    => null()
      !! device parameters
    type(generation_recombination), pointer :: genrec => null()
      !! net generation rate (G-R)
    type(chemical_pot),             pointer :: eta    => null()
      !! chemical potential
    type(ionization),               pointer :: ion    => null()
      !! donor/acceptor ionization ratio
    type(electric_field_abs),       pointer :: efield_abs => null()
      !! electric field strength

    type(jacobian), pointer :: jaco_eta => null()
    type(jacobian), pointer :: jaco_ion => null()
    type(jacobian), pointer :: jaco_efield_abs => null()
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

  subroutine calc_ionization_init(this, par, eta, ion)
    !! initialize ionization ratio calculation equation
    class(calc_ionization),      intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(chemical_pot),  target, intent(in)  :: eta
      !! chemical potential variable
    type(ionization),    target, intent(in)  :: ion
      !! donor/acceptor ionization variable

    integer :: iprov

    ! init equation
    call this%equation_init("calc_"//ion%name)
    this%par => par
    this%eta => eta
    this%ion => ion

    ! provides ionization
    iprov = this%provide(ion, par%ionvert(ion%ci))

    ! depends on eta if incomplete ionization is enabled
    if (par%smc(par%smc_default)%incomp_ion) then
      this%jaco_eta => this%init_jaco(iprov, this%depend(eta, par%ionvert(ion%ci)), const = .false.)
    end if

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_ionization_eval(this)
    !! evaluate ionization calculation equation
    class(calc_ionization), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ii_g, dop, Edop, eta, ion, diondeta

    allocate (idx(this%par%g%idx_dim))

    ci   = this%ion%ci
    ii_g = this%par%smc(this%par%smc_default)%ii_g(ci) ! degeneracy factor (2 or 4)

    do i = 1, this%par%ionvert(ci)%n
      idx = this%par%ionvert(ci)%get_idx(i)

      dop  = this%par%dop(IDX_VERTEX,0,ci)%get(idx)
      Edop = this%par%ii_E_dop(ci)%get(idx) ! E_C - E_D or E_A - E_V
      eta  = this%eta%get(idx)

      ion = dop / (1 + ii_g*exp(eta+Edop))
      diondeta = -dop / (exp(-eta-Edop)/ii_g + 2 + ii_g*exp(eta+Edop))

      call this%ion%set(idx, ion)
      call this%jaco_eta%add(idx, idx, diondeta)
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

    if (.not. par%smc(par%smc_default)%incomp_ion) call program_error("use only if incomplete ionization is turned on")

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
      call this%jaco_t%add(     idx, idx,  1.0)
      call this%jaco_genrec%add(idx, idx, -1.0)
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

    ! init variable base
    call this%variable_init("genrec_"//RATE_NAME(ci), "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci   = ci
  end subroutine

  subroutine calc_generation_recombination_init(this, par, genrec, eta, ion, efield_abs)
    !! initialize generation-recombination equation
    class(calc_generation_recombination),   intent(inout) :: this
    type(device_params),            target, intent(in)    :: par
      !! device parameters
    type(generation_recombination), target, intent(in)    :: genrec
      !! recombination variable
    type(chemical_pot),             target, intent(in)    :: eta
      !! chemical potential variable
    type(ionization),               target, intent(in)    :: ion
      !! delta ionization ratio
    type(electric_field_abs),       target, intent(in)    :: efield_abs
      !! electric field strength

    integer :: ci, iprov

    ! init equation and variable pointers
    call this%equation_init("calc_"//genrec%name)
    ci          = genrec%ci
    this%par    => par
    this%genrec => genrec
    this%eta    => eta
    this%ion    => ion
    this%efield_abs => efield_abs

    iprov = this%provide(genrec, par%ionvert(ci))

    ! dependencies
    this%jaco_eta => this%init_jaco(iprov, this%depend(eta, par%ionvert(ci)), const = .false.)
    this%jaco_ion => this%init_jaco(iprov, this%depend(ion, par%ionvert(ci)), const = .false.)
    if (this%par%smc(this%par%smc_default)%ii_pf .or. this%par%smc(this%par%smc_default)%ii_tun) then
      this%jaco_efield_abs => this%init_jaco(iprov, this%depend(efield_abs, par%ionvert(ci)), const = .false.)
    end if

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_generation_recombination_eval(this)
    !! evaluate generation-recombination equation
    class(calc_generation_recombination), intent(inout) :: this

    integer              :: i, ci
    integer, allocatable :: idx(:)
    real                 :: ii_g, dop, Edop, eta, ion, F, dF, G, dG, ef, eps, a, dE0, ddE0def, dE, ddEdef
    real                 :: prefac, dprefacdeta, dprefacdef, WKB, genrec, dgenrecdeta, dgenrecdef, dgenrecdion

    ci   = this%genrec%ci
    ii_g = this%par%smc(this%par%smc_default)%ii_g(ci) ! degeneracy factor (2 or 4)
    a    = this%par%smc(this%par%smc_default)%ii_pf_a

    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%ionvert(ci)%n
      idx = this%par%ionvert(ci)%get_idx(i)

      dop  = this%par%dop(IDX_VERTEX,0,ci)%get(idx) ! N_D or N_A
      Edop = this%par%ii_E_dop(ci)%get(idx) ! E_C - E_D or E_A - E_V
      eta  = this%eta%get(idx)
      ion  = this%ion%get(idx)
      if (this%par%smc(this%par%smc_default)%ii_pf .or. this%par%smc(this%par%smc_default)%ii_tun) then
        ef   = this%efield_abs%get(idx)
        eps  = this%par%eps(IDX_VERTEX,0)%get(idx)
      end if

      ! Poole-Frenkel barrier lowering if Edop > 0
      dE = 0.0
      ddEdef = 0.0
      if (this%par%smc(this%par%smc_default)%ii_pf .and. Edop > 0.0) then
        dE0 = sqrt(ef / (PI * eps))
        ddE0def = 0.5 / sqrt(PI * eps * ef)
        ! dE = min(dE0, Edop) with a differentiable "softmin" function
        ! the maximum difference from min(dE0, Edop) is at dE0 = Edop where dE = Edop - a*ln(2)
        dE = -a * log(exp(-dE0/a) + exp(-Edop/a))
        ddEdef = exp(-dE0/a) / (exp(-dE0/a) + exp(-Edop/a)) * ddE0def
      end if

      ! calculate distribution function F(eta+dE) and G = exp(-eta-dE)*F(eta+dE) and derivatives
      if ((this%par%smc(this%par%smc_default)%dos == DOS_PARABOLIC) .and. (this%par%smc(this%par%smc_default)%dist == DIST_MAXWELL)) then
        F    = exp(eta+dE)
        dF   = F
        G  = 1.0
        dG = 0.0
      elseif ((this%par%smc(this%par%smc_default)%dos == DOS_PARABOLIC) .and. (this%par%smc(this%par%smc_default)%dist == DIST_FERMI)) then
        ! use non-regularized functions in ionization calculation
        F  =  fd1h(eta+dE) / gamma(1.5)
        dF = fdm1h(eta+dE) / gamma(0.5)
        ! G = exp(-eta-dE)*F(eta+dE) ~ 1.0 for small eta
        if (eta > -40) then
          G = exp(-eta-dE) * F
          dG = dF * exp(-eta-dE) - F * exp(-eta-dE)
        else
          G  = 1.0
          dG = 0.0
        end if
      else
        ! FIXME: find and implement more general formula
        call program_error("Not implemented for general case")
      end if

      ! SRH prefactor with or without Poole-Frenkel effect
      prefac = (exp(-Edop+dE)*G + ii_g*F) / this%par%smc(this%par%smc_default)%ii_tau(ci)
      dprefacdeta = (exp(-Edop+dE)*dG + ii_g*dF) / this%par%smc(this%par%smc_default)%ii_tau(ci)
      dprefacdef = (exp(-Edop+dE)*G + exp(-Edop+dE)*dG + ii_g*dF) * ddEdef / this%par%smc(this%par%smc_default)%ii_tau(ci)

      ! add tunneling terms if Edop > 0
      if (this%par%smc(this%par%smc_default)%ii_tun .and. Edop > 0.0) then
        ! WKB transparency factor for a triangular barrier
        WKB = 4.0/3.0 * sqrt(2*this%par%smc(this%par%smc_default)%ii_m_tun(ci)) * Edop**1.5
        prefac = prefac + exp(-WKB/ef) / this%par%smc(this%par%smc_default)%ii_tau_tun(ci) * (1 + (ii_g-1)/(1+exp(-eta-Edop)))
        dprefacdeta = dprefacdeta + exp(-WKB/ef) / this%par%smc(this%par%smc_default)%ii_tau_tun(ci) * (ii_g-1) / (exp(eta+Edop) + 2 + exp(-eta-Edop))
        dprefacdef = dprefacdef + WKB/ef**2 * exp(-WKB/ef) / this%par%smc(this%par%smc_default)%ii_tau_tun(ci) * (1 + (ii_g-1)/(1+exp(-eta-Edop)))
      end if

      ! total rate G-R
      genrec = prefac * (dop/(1+ii_g*exp(eta+Edop)) - ion)
      dgenrecdeta = dprefacdeta * (dop/(1+ii_g*exp(eta+Edop)) - ion) - prefac * dop/(exp(-eta-Edop)/ii_g + 2 + exp(eta+Edop)*ii_g)
      dgenrecdion = -prefac
      dgenrecdef = dprefacdef * (dop/(1+ii_g*exp(eta+Edop)) - ion)

      ! save
      call this%genrec%set(idx, genrec)
      call this%jaco_eta%add(idx, idx, dgenrecdeta)
      call this%jaco_ion%add(idx, idx, dgenrecdion)
      if (this%par%smc(this%par%smc_default)%ii_pf .or. this%par%smc(this%par%smc_default)%ii_tun) call this%jaco_efield_abs%add(idx, idx, dgenrecdef)
    end do
  end subroutine

end module
