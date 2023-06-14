module dopant_m

  use device_params_m, only: device_params
  use dual_m
  use equation_m,      only: equation
  use error_m,         only: program_error
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use imref_m,         only: imref
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use semiconductor_m, only: DOP_NAME, CR_CHARGE
  use variable_m,      only: variable_real

  implicit none

  type, extends(variable_real) :: ionization
    !! ionized dopant density

    integer :: di
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
    !! Calculate the stationary dopant ionization

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

contains

  subroutine ionization_init(this, par, di)
    !! initialize density
    class(ionization),   intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: di
      !! dopant index (DOP_ELEC, DOP_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("ion"//DOP_NAME(di), "1", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%di = di

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
    !! initialize iondens calculation equation
    class(calc_ionization),      intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(ionization),    target, intent(in)  :: ion
      !! donor/acceptor density variable
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
    iprov = this%provide(ion, par%transport(IDX_VERTEX,0))

    ! Avoid jacobians if incomplete ionization is turned off
    if (par%smc%ii) then
      ! Dependent on potential
      this%jaco_pot  => this%init_jaco(iprov, this%depend(pot,  par%transport(IDX_VERTEX,0)), const = .false.)
      this%jaco_iref => this%init_jaco(iprov, this%depend(iref, par%transport(IDX_VERTEX,0)), const = .false.)
    end if


    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_ionization_eval(this)
    !! evaluate iondens calculation equation
    class(calc_ionization), intent(inout) :: this

    integer      :: i, idx(this%par%g%idx_dim), di
    real         :: b, g, ch, edop, iref
    type(dual_1) :: pot, e, ion

    ! Complete ionization
    if (.not. this%par%smc%ii) then
      do i = 1, this%par%transport(IDX_VERTEX,0)%n
        call this%ion%set(this%par%transport(IDX_VERTEX,0)%get_idx(i), 1.0)
      end do
      return
    end if

    di = this%ion%di
    ch = CR_CHARGE(di)
    g  = this%par%smc%g_dop(di)

    call pot%init(0.0, i = 1)

    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx  = this%par%transport(IDX_VERTEX,0)%get_idx(i)
      b    = this%par%asb(di)%get(idx)
      edop = this%par%smc%band_edge(di) + ch * this%par%edop(di)%get(idx)
      iref = this%iref%get(idx)

      pot%x = this%pot%get(idx)

      ! Calculate
      e   = exp(- ch * (iref - pot + edop))
      ion = 1 - b / (1 + g * e)

      call this%ion%set(idx, ion%x)
      call this%jaco_pot%set( idx, idx,  ion%dx(1))
      call this%jaco_iref%set(idx, idx, -ion%dx(1))
    end do
  end subroutine

end module
