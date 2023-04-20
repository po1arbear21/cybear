module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params, CR_CHARGE, CR_NAME
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use jacobian_m,        only: jacobian, jacobian_ptr
  use res_equation_m,    only: res_equation
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,       only: vselector
  use error_m,           only: assert_failed, program_error

  implicit none

  private
  public continuity

  type, extends(res_equation) :: continuity
    !! dn/dt + div j = 0
    integer :: ci
      !! carrier index

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    type(vselector)              :: dens
      !! main variable: density
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities in 2 directions

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, dens, cdens)
    !! initialize continuity equation
    class(continuity),           intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(density),               intent(in)  :: dens
      !! electron/hole density
    type(current_density),       intent(in)  :: cdens(:)
      !! electron/hole current density

    integer              :: i, ict, idx_dir, idens, idx_dim
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status
    real                 :: surf

    ! init base
    call this%equation_init(CR_NAME(dens%ci)//"continuity")
    this%par => par
    this%ci  = dens%ci

    idx_dim = par%g%idx_dim
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), icdens(idx_dim))
    allocate (this%cdens(idx_dim), this%st_nn(idx_dim), this%jaco_cdens(idx_dim))

    ! init variable selectors
    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, size(par%contacts))])
    do idx_dir = 1, idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
    end do

    ! init residuals using this%dens as main variable
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(par%g)
    do idx_dir = 1, idx_dim
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
    end do
    call this%st_em%init()

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do

    ! init jacobians
    this%jaco_dens   => this%init_jaco_f(idens, st = [this%st_em%get_ptr(), (this%st_dir%get_ptr(), ict = 1, size(par%contacts))], const = .true., dtime = .false.)
    this%jaco_dens_t => this%init_jaco_f(idens, st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, size(par%contacts))], const = .true., dtime = .true. )
    do idx_dir = 1, idx_dim
      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), st = [this%st_nn(idx_dir)%get_ptr(), (this%st_em%get_ptr(), ict = 1, size(par%contacts))], const = .true., dtime = .false.)
    end do

    ! loop over transport edges
    do idx_dir = 1, idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        ! get indices of edge, first and second vertex
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get edge surface
        surf = par%tr_surf(idx_dir)%get(idx)

        ! set values if uncontacted
        if (par%ict%get(idx1) == 0) call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
        if (par%ict%get(idx2) == 0) call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
      end do
    end do

    ! loop over transport vertices
    do i = 1, par%transport(IDX_VERTEX,0)%n
      idx1 = par%transport(IDX_VERTEX,0)%get_idx(i)
      ict = par%ict%get(idx1)
      if (ict == 0) then
        ! set adjoint volume for time derivative factor
        call this%jaco_dens_t%set(idx1, idx1, par%tr_vol%get(idx1))
      else
        ! dirichlet condition
        call this%jaco_dens%set(idx1, idx1, 1.0)
      end if
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    class(continuity), intent(inout) :: this

    integer              :: i, ict, idx_dir, idx_dim
    integer, allocatable :: idx(:)
    real,    allocatable :: tmp(:)

    idx_dim = this%par%g%idx_dim
    allocate(idx(idx_dim), tmp(this%dens%n))

    ! calculate residuals excluding boundary conditions
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    call this%f%set(tmp)

    ! apply ideal ohmic boundary conditions
    do ict = 1, size(this%par%contacts)
      do i = 1, this%par%transport_vct(ict)%n
        idx = this%par%transport_vct(ict)%get_idx(i)
        call this%f%update(idx, [-this%par%n_intrin * exp(- CR_CHARGE(this%ci) * this%par%contacts(ict)%phims)])
      end do
    end do
  end subroutine

end module
