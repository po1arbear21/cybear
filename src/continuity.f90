module continuity_m

  use contact_m,         only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use res_equation_m,    only: res_equation
  use schottky_m,        only: schottky_bc
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil, stencil_ptr
  use vselector_m,       only: vselector

  implicit none

  private
  public continuity

  type, extends(res_equation) :: continuity
    !! dn/dt + div(J) - (G - R) + sum_ict A_ct * j_bc(ict) = 0
    !!
    !! - ohmic/gate vertices: Dirichlet BC (dens = b)
    !! - Schottky boundary current j_bc enters via dependency on schottky_bc(ict)

    integer :: ci
      !! carrier index

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    type(vselector)              :: dens
      !! main variable: density
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities per spatial direction
    type(vselector)              :: genrec
      !! generation - recombination

    type(vselector), allocatable :: schottky_bc_sel(:)
      !! per-contact wrapper for schottky_bc(ict) on transport_vct(ict);
      !! initialized only for CT_SCHOTTKY contacts (others are unused defaults)

    real, allocatable :: b(:)
      !! right-hand side

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
    type(jacobian_ptr), allocatable :: jaco_schottky_bc(:)
      !! per-contact Jacobian (size par%nct); %p associated only for CT_SCHOTTKY contacts
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, dens, cdens, genrec, sbc)
    !! initialize continuity equation
    class(continuity),              intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    type(density),                  intent(in)  :: dens
      !! electron/hole density
    type(current_density),          intent(in)  :: cdens(:)
      !! electron/hole current density
    type(generation_recombination), intent(in)  :: genrec
      !! generation-recombination rate
    type(schottky_bc),              intent(in)  :: sbc(:)
      !! Schottky boundary current variables, indexed by contact (sized par%nct);
      !! entries for non-Schottky contacts are unused

    integer              :: ci, i, ict, ict2, idx_dir, idens, igenrec, idep_sbc, j
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status
    real                 :: surf, F, dF
    type(stencil_ptr), allocatable :: st_dens_ct(:), st_dens_t_ct(:), st_cdens_ct(:), st_sbc(:)

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    call this%equation_init(CR_NAME(dens%ci)//"continuity")
    this%par => par
    this%ci  = ci

    allocate (idx(par%g%idx_dim), idx1(par%g%idx_dim), idx2(par%g%idx_dim), icdens(par%g%idx_dim))
    allocate (this%cdens(par%g%idx_dim), this%st_nn(par%g%idx_dim), this%jaco_cdens(par%g%idx_dim))

    ! init variable selectors
    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    do idx_dir = 1, par%g%idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
    end do
    if (par%smc%incomp_ion) call this%genrec%init(genrec, par%ionvert(ci))

    ! per-contact schottky_bc selectors (only Schottky entries init'd)
    allocate(this%schottky_bc_sel(par%nct), this%jaco_schottky_bc(par%nct))
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        call this%schottky_bc_sel(ict)%init(sbc(ict), par%transport_vct(ict))
      end if
    end do

    ! init residual variable (= dens)
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(par%g)
    do idx_dir = 1, par%g%idx_dim
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
    end do
    call this%st_em%init()

    ! per-contact stencil tables for dens / dens_t / cdens
    allocate(st_dens_ct(par%nct), st_dens_t_ct(par%nct), st_cdens_ct(par%nct))
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        st_dens_ct(ict)   = this%st_em%get_ptr()
        st_dens_t_ct(ict) = this%st_dir%get_ptr()
      else
        st_dens_ct(ict)   = this%st_dir%get_ptr()
        st_dens_t_ct(ict) = this%st_em%get_ptr()
      end if
    end do

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, par%g%idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians (all const after refactor; no per-iteration rebuilds)
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (st_dens_ct(ict), ict = 1, par%nct)], &
      & const = .true., dtime = .false.)
    this%jaco_dens_t => this%init_jaco_f(idens, &
      & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
      & const = .true., dtime = .true. )
    do idx_dir = 1, par%g%idx_dim
      ! cdens at contact: Schottky has flux-divergence contribution; ohmic/gate have none.
      do ict = 1, par%nct
        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          st_cdens_ct(ict) = this%st_nn(idx_dir)%get_ptr()
        else
          st_cdens_ct(ict) = this%st_em%get_ptr()
        end if
      end do
      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = [this%st_nn(idx_dir)%get_ptr(), (st_cdens_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! per-Schottky-contact Jacobian on schottky_bc(ict): diagonal +A_ct only at matching contact
    allocate(st_sbc(1 + par%nct))
    do ict = 1, par%nct
      if (par%contacts(ict)%type /= CT_SCHOTTKY) cycle

      idep_sbc = this%depend(this%schottky_bc_sel(ict))

      st_sbc(1) = this%st_em%get_ptr()  ! interior — no contribution
      do ict2 = 1, par%nct
        if (ict2 == ict) then
          st_sbc(1 + ict2) = this%st_dir%get_ptr()
        else
          st_sbc(1 + ict2) = this%st_em%get_ptr()
        end if
      end do
      this%jaco_schottky_bc(ict)%p => this%init_jaco_f(idep_sbc, &
        & st = st_sbc, const = .true., dtime = .false.)
    end do

    ! cdens jacobian entries: surface contributions at interior + Schottky vertices
    do idx_dir = 1, par%g%idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        surf = par%tr_surf(idx_dir)%get(idx)

        ict = par%ict%get(idx1)
        if (ict == 0 .or. par%contacts(max(ict,1))%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%add(idx1, idx,  surf)
        end if
        ict = par%ict%get(idx2)
        if (ict == 0 .or. par%contacts(max(ict,1))%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%add(idx2, idx, -surf)
        end if
      end do
    end do

    ! density time-derivative volume factor at every transport vertex
    do i = 1, par%transport_vct(0)%n
      idx1 = par%transport_vct(0)%get_idx(i)
      call this%jaco_dens_t%add(idx1, idx1, par%tr_vol%get(idx1))
    end do
    do ict = 1, par%nct
      if (par%contacts(ict)%type /= CT_SCHOTTKY) cycle
      do i = 1, par%transport_vct(ict)%n
        idx1 = par%transport_vct(ict)%get_idx(i)
        call this%jaco_dens_t%add(idx1, idx1, par%tr_vol%get(idx1))
      end do
    end do

    ! generation-recombination volume factor
    if (par%smc%incomp_ion) then
      do i = 1, par%ionvert(ci)%n
        idx1 = par%ionvert(ci)%get_idx(i)
        call this%jaco_genrec%add(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! schottky_bc Jacobian: +A_ct at every Schottky vertex
    do ict = 1, par%nct
      if (par%contacts(ict)%type /= CT_SCHOTTKY) cycle
      do i = 1, par%transport_vct(ict)%n
        idx1 = par%transport_vct(ict)%get_idx(i)
        call this%jaco_schottky_bc(ict)%p%add(idx1, idx1, par%get_ct_surf(ict, idx1))
      end do
    end do

    ! boundary conditions: Dirichlet for ohmic/gate; b=0 (default) for Schottky
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      do i = 1, par%transport_vct(ict)%n
        j = j + 1
        idx1 = par%transport_vct(ict)%get_idx(i)

        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          ! Schottky: residual = jaco_dens_t*dens + jaco_cdens*cdens + jaco_schottky_bc*sbc;  b=0
          cycle
        end if

        ! Ohmic/Gate: F = 1.0 * dens(idx) - b(j)
        call this%jaco_dens%add(idx1, idx1, 1.0)
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
    !! evaluate continuity equation: F = sum_k jaco_k * dep_k - b
    class(continuity), intent(inout) :: this

    integer           :: idx_dir, ict
    real, allocatable :: tmp(:)

    allocate (tmp(this%f%n))
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)

    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do

    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    do ict = 1, this%par%nct
      if (associated(this%jaco_schottky_bc(ict)%p)) &
        call this%jaco_schottky_bc(ict)%p%matr%mul_vec(this%schottky_bc_sel(ict)%get(), tmp, fact_y = 1.0)
    end do

    call this%f%set(tmp - this%b)
  end subroutine

end module
