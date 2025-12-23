module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use electric_field_m,  only: electric_field
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use beam_generation_m, only: beam_generation
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,       only: vselector

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
    type(vselector)              :: genrec
      !! generation - recombination
    type(vselector)              :: bgen
      !! external beam generation (STEM-EBIC)

    ! Schottky support
    type(vselector), allocatable :: efield(:)
      !! electric field components (for Schottky)
    integer, allocatable :: schottky_normal(:)
      !! normal direction for each contact (0 if not Schottky)
    real, allocatable :: accum_TE(:)
      !! accumulated TE current per contact (set during eval)
    real, allocatable :: accum_TN(:)
      !! accumulated TN current per contact (set during eval)
    real, allocatable :: v_surf(:)
      !! surface recombination velocity per contact (for Robin BC)

    real, allocatable :: b(:)
      !! right-hand side

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
    type(jacobian),     pointer     :: jaco_bgen     => null()
      !! Jacobian for beam generation (constant, no dependency on solution)
    type(jacobian),     pointer     :: jaco_iref     => null()
      !! Direct Jacobian for iref (Schottky contacts only, bypasses chain rule)
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, dens, cdens, genrec, efield, bgen)
    !! initialize continuity equation
    class(continuity),              intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    type(density),                  intent(in)  :: dens
      !! carrier density
    type(current_density),          intent(in)  :: cdens(:)
      !! electron/hole current density
    type(generation_recombination), intent(in)  :: genrec
      !! generation-recombination rate
    type(electric_field), optional, intent(in)  :: efield(:)
      !! electric field components (for Schottky, optional)
    type(beam_generation), optional, intent(in) :: bgen
      !! external beam generation (STEM-EBIC, optional)

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, ibgen, j, dir
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status, has_schottky, has_beam
    real                 :: surf, F, dF

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    call this%equation_init(CR_NAME(dens%ci)//"continuity")
    this%par => par
    this%ci  = ci

    idx_dim = par%g%idx_dim
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), icdens(idx_dim))
    allocate (this%cdens(idx_dim), this%st_nn(idx_dim), this%jaco_cdens(idx_dim))

    ! init variable selectors
    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    do idx_dir = 1, idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
    end do
    if (par%smc%incomp_ion) call this%genrec%init(genrec, par%ionvert(ci))

    ! init beam generation selector if enabled
    has_beam = present(bgen) .and. par%has_beam_gen
    if (has_beam) call this%bgen%init(bgen, par%transport(IDX_VERTEX, 0))

    ! init residuals using this%dens or this%iref as main variable
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
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (this%st_dir%get_ptr(), ict = 1, par%nct)], &
      & const = .true., dtime = .false.)
    this%jaco_dens_t => this%init_jaco_f(idens, &
      & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
      & const = .true., dtime = .true. )
    do idx_dir = 1, idx_dim
      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = [this%st_nn(idx_dir)%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! beam generation jacobian
    if (has_beam) then
      ibgen = this%depend(this%bgen)
      this%jaco_bgen => this%init_jaco_f(ibgen, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! set current density jacobian entries
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

    ! density time derivative factor
    do i = 1, par%transport_vct(0)%n
      idx1 = par%transport_vct(0)%get_idx(i)
      call this%jaco_dens_t%set(idx1, idx1, par%tr_vol%get(idx1))
    end do

    ! generation-recombination
    if (par%smc%incomp_ion) then
      do i = 1, par%ionvert(ci)%n
        idx1 = par%ionvert(ci)%get_idx(i)
        call this%jaco_genrec%set(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! beam generation (STEM-EBIC)
    ! F = dn/dt*V + div(j)*V - G_beam*V = 0  =>  dF/dG = -V
    if (has_beam) then
      do i = 1, par%transport_vct(0)%n
        idx1 = par%transport_vct(0)%get_idx(i)
        call this%jaco_bgen%set(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! boundary conditions: Dirichlet for Ohmic/Gate, handled in eval for Schottky
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      do i = 1, par%transport_vct(ict)%n
        j = j + 1
        idx1 = par%transport_vct(ict)%get_idx(i)
        call this%jaco_dens%set(idx1, idx1, 1.0)
        if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
          this%b(j) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
        else
          call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
          this%b(j) = par%smc%edos(ci) * F
        end if
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer           :: idx_dir
    real, allocatable :: tmp(:)

    allocate (tmp(this%f%n))
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)
    if (associated(this%jaco_bgen)) call this%jaco_bgen%matr%mul_vec(this%bgen%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp - this%b)
  end subroutine

end module
