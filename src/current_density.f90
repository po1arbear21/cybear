m4_include(util/macro.f90.inc)

module current_density_m

  use density_m,        only: density
  use degen_table_m,    only: degen_table
  use device_params_m,  only: device_params
  use dual_m
  use equation_m,       only: equation
  use error_m,          only: assert_failed, program_error
  use fermi_m,          only: inv_fermi_dirac_integral_1h_reg
  use grid_m,           only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m, only: DIR_NAME
  use high_precision_m
  use ieee_arithmetic,  only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_value
  use jacobian_m,       only: jacobian
  use math_m,           only: ber, dberdx, log1p
  use mobility_m,       only: mobility
  use potential_m,      only: potential
  use radau5_m,         only: ode_options, ode_result, radau5
  use semiconductor_m,  only: CR_NAME, CR_CHARGE
  use stencil_m,        only: dirichlet_stencil, near_neighb_stencil
  use variable_m,       only: variable_real

  implicit none

  private
  public current_density, calc_current_density

  type, extends(variable_real) :: current_density
    !! electron/hole current density

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_current_density
    !! calculate current density by drift-diffusion model

    logical :: stat

    type(device_params),   pointer :: par   => null()
    type(potential),       pointer :: pot   => null()
    type(density),         pointer :: dens  => null()
    type(current_density), pointer :: cdens => null()
    type(mobility),        pointer :: mob   => null()

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_dens => null()
    type(jacobian), pointer :: jaco_mob  => null()
  contains
    procedure :: init => calc_current_density_init
    procedure :: eval => calc_current_density_eval

    procedure, private :: get_curr    => calc_current_density_get_curr
    procedure, private :: get_curr_sg => calc_current_density_get_curr_sg
  end type

contains

  subroutine current_density_init(this, par, ci, idx_dir)
    !! initialize current density
    class(current_density), intent(out) :: this
    type(device_params),    intent(in)  :: par
      !! device parameters
    integer,                intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    integer,                intent(in)  :: idx_dir
      !! edge direction

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
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

  subroutine calc_current_density_init(this, par, pot, dens, cdens, mob)
    !! initialize drift-diffusion equation
    class(calc_current_density),   intent(out) :: this
    type(device_params),   target, intent(in)  :: par
      !! device parameters
    type(potential),       target, intent(in)  :: pot
      !! potential variable
    type(density),         target, intent(in)  :: dens
      !! density variable
    type(current_density), target, intent(in)  :: cdens
      !! current density variable
    type(mobility),        target, intent(in)  :: mob
      !! mobility variable

    integer :: idx_dir, ci, iprov

    print "(A)", "calc_current_density_init"

    idx_dir = cdens%idx_dir
    ci      = cdens%ci

    ! init base
    call this%equation_init("calc_"//cdens%name)
    this%par   => par
    this%pot   => pot
    this%dens  => dens
    this%cdens => cdens
    this%mob   => mob

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_nn%init(par%g, IDX_EDGE, idx_dir, IDX_VERTEX, 0)

    ! provide current density
    iprov = this%provide(cdens, par%transport(IDX_EDGE, idx_dir))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on density
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on mobility
    this%jaco_mob => this%init_jaco(iprov, this%depend(mob, par%transport(IDX_EDGE, idx_dir)), st = [this%st_dir%get_ptr()])

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    !! evaluate drift-diffusion equation
    class(calc_current_density), intent(inout) :: this

    integer               :: i, idx_dir, idx_dim
    real                  :: pot(2), dens(2), j, djdpot(2), djddens(2), djdmob, len, mob
    integer, allocatable  :: idx(:), idx1(:), idx2(:)
    logical               :: status

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim))

    ! loop over transport edges in parallel
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(i,pot,dens,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    !$omp shared(this,idx_dir,idx_dim)
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

      ! parameters
      len     = this%par%g%get_len(idx1, idx_dir)
      pot( 1) = this%pot%get( idx1)
      pot( 2) = this%pot%get( idx2)
      dens(1) = this%dens%get(idx1)
      dens(2) = this%dens%get(idx2)
      mob     = this%mob%get( idx )

      ! get current along edge
      call this%get_curr(this%par%smc%degen, len, pot, dens, mob, j, djdpot, djddens, djdmob)

      ! set current density + derivatives
      call this%cdens%set(idx, j)
      call this%jaco_pot%set( idx, idx1, djdpot(1))
      call this%jaco_pot%set( idx, idx2, djdpot(2))
      call this%jaco_dens%set(idx, idx1, djddens(1))
      call this%jaco_dens%set(idx, idx2, djddens(2))
      call this%jaco_mob%set( idx, idx,  djdmob)
    end do
    !$omp end parallel do
  end subroutine

  subroutine calc_current_density_get_curr(this, degen, len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! get current along edge
    class(calc_current_density), intent(in)  :: this
    logical,                     intent(in)  :: degen
      !! degeneracy flag (true: FD statistics, false: MB statistics)
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    integer :: ci
    real    :: ch, dpot, edos, n(2), eta(2), detadn(2), djdeta(2), djdn(2), djddpot

    ! abbreviations
    ci   = this%cdens%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)

    ! normalize potential drop and density
    dpot = - ch * (pot(2) - pot(1))
    n    = dens / edos

    if (degen) then
      ! generalized Scharfetter-Gummel
      call inv_fermi_dirac_integral_1h_reg(n(1), eta(1), detadn(1))
      call inv_fermi_dirac_integral_1h_reg(n(2), eta(2), detadn(2))
      call this%par%degen_tab%get(eta, dpot, j, djdeta, djddpot)
      djdn = djdeta * detadn
    else
      ! Scharfetter-Gummel
      call this%get_curr_sg(n, dpot, j, djdn, djddpot)
    end if

    ! denormalize
    j       = j * mob * edos / len
    djdpot  = ch * [1.0, -1.0] * djddpot * mob * edos / len
    djddens = djdn * mob / len
    djdmob  = j / mob
  end subroutine

  subroutine calc_current_density_get_curr_sg(this, n, dpot, j, djdn, djddpot)
    use math_m, only: ber

    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real :: B1, B2, dB1, dB2

    m4_ignore(this)

    B1  = ber(    dpot)
    B2  = ber(   -dpot)
    dB1 = dberdx( dpot)
    dB2 = dberdx(-dpot)

    j       = B2 * n(1) - B1 * n(2)
    djdn    = [B2, -B1]
    djddpot = - dB2 * n(1) - dB1 * n(2)
  end subroutine

end module
