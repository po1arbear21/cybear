module current_density_m

  use density_m,       only: density
  use device_params_m, only: device_params, CR_NAME, CR_CHARGE, DIR_NAME
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,     only: grid_data2_real
  use jacobian_m,      only: jacobian
  use math_m,          only: ber, dberdx
  use mobility_m,      only: mobility
  use potential_m,     only: potential
  use stencil_m,       only: dirichlet_stencil, near_neighb_stencil
  use variable_m,      only: variable_real

  implicit none

  private
  public current_density, calc_current_density

  type, extends(variable_real) :: current_density
    !! electron/hole current density
    integer       :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_current_density
    !! calculate current density by drift-diffusion model

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

    type(grid_data2_real), pointer :: p => null()

    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
    this%ci = ci

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
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

    integer, parameter :: eye(2,2) = reshape([1, 0, 0, 1], [2, 2])

    integer :: i, idx1(2), idx2(2), idx_dir, ci
    real    :: ber1, ber2, dber1, dber2, pot1, pot2, dens1, dens2, len, mob, ch

    idx_dir = this%cdens%idx_dir
    ci      = this%cdens%ci
    ch      = CR_CHARGE(ci)

    ! loop over transport edges
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx1 = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      idx2 = idx1 + eye(:,idx_dir)

      ! parameters
      len   = this%par%g%get_len(idx1, idx_dir)
      pot1  = this%pot%get(idx1)
      pot2  = this%pot%get(idx2)
      dens1 = this%dens%get(idx1)
      dens2 = this%dens%get(idx2)
      mob   = this%mob%get(idx1)

      ! Bernoulli function
      ber1  = ber(ch * (pot1 - pot2))
      ber2  = ber(ch * (pot2 - pot1))
      dber1 = ch * dberdx(ch * (pot1 - pot2))
      dber2 = ch * dberdx(ch * (pot2 - pot1))

      ! set current density
      call this%cdens%set(idx1, - mob * (ber1 * dens2 - ber2 * dens1) / len)

      ! set jaco_pot entries
      call this%jaco_pot%set(idx1, idx1, - mob * (dber1 * dens2 + dber2 * dens1) / len)
      call this%jaco_pot%set(idx1, idx2,   mob * (dber1 * dens2 + dber2 * dens1) / len)

      ! set jaco_dens entries
      call this%jaco_dens%set(idx1, idx1,   mob * ber2 / len)
      call this%jaco_dens%set(idx1, idx2, - mob * ber1 / len)

      ! set jaco_mob entries
      call this%jaco_mob%set(idx1, idx1, -(ber1 * dens2 - ber2 * dens1) / len)
    end do
  end subroutine

end module
