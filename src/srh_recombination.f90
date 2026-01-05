module srh_recombination_m
  !! Shockley-Read-Hall recombination (simplified form)
  !!
  !! R_SRH = (n*p - ni²) / [τp*(n + ni) + τn*(p + ni)]
  !!
  !! This is the simplified SRH formula with trap at midgap (n1 = p1 = ni).
  !! Positive R means net recombination (removes electron-hole pairs).
  !! Negative R means net generation.

  use density_m,       only: density
  use device_params_m, only: device_params
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use jacobian_m,      only: jacobian
  use semiconductor_m, only: CR_ELEC, CR_HOLE
  use variable_m,      only: variable
  use normalization_m, only: denorm

  implicit none

  private
  public srh_recombination, calc_srh_recombination

  type, extends(variable) :: srh_recombination
    !! SRH recombination rate variable
    !!
    !! Unit: 1/cm^3/s (recombination rate per unit volume)
    !! Positive = net recombination, Negative = net generation

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for 1D grids
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for 2D grids
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for 3D grids
  contains
    procedure :: init => srh_recombination_init
  end type

  type, extends(equation) :: calc_srh_recombination
    !! Calculate SRH recombination rate
    !!
    !! R = (n*p - ni²) / [τp*(n + ni) + τn*(p + ni)]

    type(device_params), pointer :: par => null()
      !! device parameters

    type(srh_recombination), pointer :: srh => null()
      !! SRH recombination rate variable

    type(density), pointer :: dens_n => null()
      !! electron density variable
    type(density), pointer :: dens_p => null()
      !! hole density variable

    real :: tau_n
      !! electron lifetime (normalized)
    real :: tau_p
      !! hole lifetime (normalized)
    real :: ni_sq
      !! intrinsic concentration squared (normalized)
    real :: ni
      !! intrinsic concentration (normalized)

    type(jacobian), pointer :: jaco_dens_n => null()
      !! Jacobian dR/dn
    type(jacobian), pointer :: jaco_dens_p => null()
      !! Jacobian dR/dp
  contains
    procedure :: init => calc_srh_recombination_init
    procedure :: eval => calc_srh_recombination_eval
  end type

contains

  subroutine srh_recombination_init(this, par)
    !! Initialize SRH recombination rate variable
    class(srh_recombination), intent(out) :: this
    type(device_params),      intent(in)  :: par
      !! device parameters

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init variable base
    call this%variable_init("R_SRH", "1/cm^3/s", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data for direct access
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

  subroutine calc_srh_recombination_init(this, par, dens_n, dens_p, srh)
    !! Initialize SRH recombination calculation equation
    class(calc_srh_recombination),  intent(out)   :: this
    type(device_params),    target, intent(in)    :: par
      !! device parameters
    type(density),          target, intent(in)    :: dens_n
      !! electron density
    type(density),          target, intent(in)    :: dens_p
      !! hole density
    type(srh_recombination), target, intent(in)   :: srh
      !! SRH recombination variable

    integer :: iprov, idep_n, idep_p
    real    :: Eg

    print "(A)", "calc_srh_recombination_init"

    ! init equation
    call this%equation_init("calc_R_SRH")
    this%par    => par
    this%srh    => srh
    this%dens_n => dens_n
    this%dens_p => dens_p

    ! store lifetime parameters
    this%tau_n = par%smc%srh_tau_n
    this%tau_p = par%smc%srh_tau_p

    ! calculate ni² = Nc * Nv * exp(-Eg/kT)
    Eg = par%smc%band_gap  ! already normalized to kT
    this%ni_sq = par%smc%edos(CR_ELEC) * par%smc%edos(CR_HOLE) * exp(-Eg)
    this%ni    = sqrt(this%ni_sq)

    print "(A,ES12.4)", "  ni  = ", denorm(this%ni, 'cm^-3')
    print "(A,ES12.4,A,ES12.4)", "  tau_n  = ", denorm(this%tau_n, 's'), ", tau_p  = ", denorm(this%tau_p, 's')

    ! provide SRH rate on transport vertices
    iprov = this%provide(srh, par%transport(IDX_VERTEX, 0))

    ! dependencies on both carrier densities
    idep_n = this%depend(dens_n, par%transport(IDX_VERTEX, 0))
    idep_p = this%depend(dens_p, par%transport(IDX_VERTEX, 0))

    ! initialize Jacobians
    this%jaco_dens_n => this%init_jaco(iprov, idep_n, const = .false.)
    this%jaco_dens_p => this%init_jaco(iprov, idep_p, const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_srh_recombination_eval(this)
    !! Evaluate SRH recombination rate
    !!
    !! R = (n*p - ni²) / [τp*(n + ni) + τn*(p + ni)]
    !!
    !! Jacobians:
    !! dR/dn = [p*D - (np - ni²)*τp] / D²
    !! dR/dp = [n*D - (np - ni²)*τn] / D²
    !! where D = τp*(n + ni) + τn*(p + ni)
    class(calc_srh_recombination), intent(inout) :: this

    integer              :: i, i_max
    integer, allocatable :: idx(:)
    real                 :: n, p, np, numer, denom, R, dRdn, dRdp
    real                 :: R_max, n_max, p_max

    allocate (idx(this%par%g%idx_dim))

    R_max = 0.0
    i_max = 0

    do i = 1, this%par%transport(IDX_VERTEX, 0)%n
      idx = this%par%transport(IDX_VERTEX, 0)%get_idx(i)

      ! get carrier densities
      n = this%dens_n%get(idx)
      p = this%dens_p%get(idx)

      ! SRH recombination rate
      np    = n * p
      numer = np - this%ni_sq
      denom = this%tau_p * (n + this%ni) + this%tau_n * (p + this%ni)

      if (denom > 1e-30) then
        R = numer / denom

        ! Jacobians: dR/dn and dR/dp
        dRdn = (p * denom - numer * this%tau_p) / (denom * denom)
        dRdp = (n * denom - numer * this%tau_n) / (denom * denom)
      else
        R    = 0.0
        dRdn = 0.0
        dRdp = 0.0
      end if

      call this%srh%set(idx, R)
      call this%jaco_dens_n%set(idx, idx, dRdn)
      call this%jaco_dens_p%set(idx, idx, dRdp)

      ! Track max |R| to find where recombination is strongest
      if (abs(R) > abs(R_max)) then
        R_max = R
        n_max = n
        p_max = p
        i_max = i
      end if
    end do

    ! Debug: print max R (where beam effect should be visible)
    if (i_max > 0) then
      print "(A,I6,A,ES10.3,A,ES10.3,A,ES10.3)", "  SRH max @ i=", i_max, &
        & ": n=", denorm(n_max,'cm^-3'), " p=", denorm(p_max,'cm^-3'), " R=", denorm(R_max,'cm^-3/s')
    end if
  end subroutine

end module
