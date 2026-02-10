module continuity_m

  use contact_m,         only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use electric_field_m,  only: electric_field
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use res_equation_m,    only: res_equation
  use schottky_m,        only: get_normal_dir, schottky_velocity, schottky_n0b, schottky_tunneling
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil, stencil_ptr
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

    ! Schottky contact support
    type(vselector), allocatable :: efield(:)
      !! electric field components (only allocated when Schottky contacts present)
    integer, allocatable :: schottky_normal(:)
      !! normal direction for each contact (0 if not Schottky)
    real, allocatable :: v_surf(:)
      !! thermionic velocity per contact (for Robin BC)

    real, allocatable :: b(:)
      !! right-hand side

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, dens, cdens, genrec, efield)
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
    type(electric_field), optional, intent(in)  :: efield(:)
      !! electric field components (required when Schottky contacts are present)

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j, dir
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status, has_schottky
    real                 :: surf, F, dF
    type(stencil_ptr), allocatable :: st_dens_ct(:), st_dens_t_ct(:), st_cdens_ct(:)

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    call this%equation_init(CR_NAME(dens%ci)//"continuity")
    this%par => par
    this%ci  = ci

    idx_dim = par%g%idx_dim
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), icdens(idx_dim))
    allocate (this%cdens(idx_dim), this%st_nn(idx_dim), this%jaco_cdens(idx_dim))

    ! check for Schottky contacts
    has_schottky = any(par%contacts(1:par%nct)%type == CT_SCHOTTKY)

    ! initialize Schottky support if needed
    if (has_schottky) then
      if (.not. present(efield)) call program_error("E-field required for Schottky contacts")

      allocate(this%schottky_normal(par%nct))
      allocate(this%v_surf(par%nct), source=0.0)
      do ict = 1, par%nct
        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          this%schottky_normal(ict) = get_normal_dir(par, ict)
          this%v_surf(ict) = schottky_velocity(par, ci, ict)
        else
          this%schottky_normal(ict) = 0
        end if
      end do

      allocate(this%efield(par%g%dim))
      do dir = 1, par%g%dim
        call this%efield(dir)%init(efield(dir), [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
      end do
    end if

    ! init variable selectors
    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    do idx_dir = 1, idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
    end do
    if (par%smc%incomp_ion) call this%genrec%init(genrec, par%ionvert(ci))

    ! init residuals
    call this%init_f(this%dens)

    ! init stencils
    call this%st_dir%init(par%g)
    do idx_dir = 1, idx_dim
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
    end do
    call this%st_em%init()

    ! setup contact-specific stencils
    allocate(st_dens_ct(par%nct), st_dens_t_ct(par%nct), st_cdens_ct(par%nct))
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        st_dens_ct(ict)   = this%st_dir%get_ptr()   ! Schottky: updated per iteration
        st_dens_t_ct(ict) = this%st_dir%get_ptr()   ! Schottky: has time derivative
      else
        st_dens_ct(ict)   = this%st_dir%get_ptr()   ! Ohmic/Gate: Dirichlet
        st_dens_t_ct(ict) = this%st_em%get_ptr()    ! Ohmic/Gate: no time term
      end if
    end do

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians
    ! jaco_dens: non-constant when Schottky contacts present (updated each iteration)
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (st_dens_ct(ict), ict = 1, par%nct)], &
      & const = .not. has_schottky, dtime = .false.)
    this%jaco_dens_t => this%init_jaco_f(idens, &
      & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
      & const = .true., dtime = .true. )
    do idx_dir = 1, idx_dim
      ! Schottky contacts need edge contributions (current divergence at boundary)
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

    ! set current density jacobian entries
    do idx_dir = 1, idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        surf = par%tr_surf(idx_dir)%get(idx)

        ! set values for interior vertices and Schottky contacts
        ! (Ohmic/Gate use Dirichlet BC, so no current divergence term)
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

    ! density time derivative factor
    do i = 1, par%transport_vct(0)%n
      idx1 = par%transport_vct(0)%get_idx(i)
      call this%jaco_dens_t%add(idx1, idx1, par%tr_vol%get(idx1))
    end do

    ! Schottky contacts also have time derivative terms
    if (has_schottky) then
      do ict = 1, par%nct
        if (par%contacts(ict)%type /= CT_SCHOTTKY) cycle
        do i = 1, par%transport_vct(ict)%n
          idx1 = par%transport_vct(ict)%get_idx(i)
          call this%jaco_dens_t%add(idx1, idx1, par%tr_vol%get(idx1))
        end do
      end do
    end if

    ! generation-recombination
    if (par%smc%incomp_ion) then
      do i = 1, par%ionvert(ci)%n
        idx1 = par%ionvert(ci)%get_idx(i)
        call this%jaco_genrec%add(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! boundary conditions
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      do i = 1, par%transport_vct(ict)%n
        j = j + 1
        idx1 = par%transport_vct(ict)%get_idx(i)

        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          ! Schottky: placeholder jaco_dens entry (updated per iteration in eval)
          call this%jaco_dens%add(idx1, idx1, 1.0)
          ! b(j) = 0 for Schottky (residual computed in eval)
        else
          ! Ohmic/Gate: Dirichlet BC
          call this%jaco_dens%add(idx1, idx1, 1.0)
          if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
            this%b(j) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
          else
            call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
            this%b(j) = par%smc%edos(ci) * F
          end if
        end if
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer              :: idx_dir, ict, i, j, normal_dir
    integer, allocatable :: idx(:)
    real, allocatable    :: tmp(:), dens_arr(:), efield_arr(:)
    real                 :: E_normal, dens, A_ct
    real                 :: n0B, J_tn, dJ_tn_ddens, v_surf_ict

    allocate (tmp(this%f%n))
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)

    ! zero out jaco_dens contribution at Schottky contacts (before adding cdens)
    if (allocated(this%efield)) then
      j = this%par%transport_vct(0)%n
      do ict = 1, this%par%nct
        do i = 1, this%par%transport_vct(ict)%n
          j = j + 1
          if (this%par%contacts(ict)%type == CT_SCHOTTKY) tmp(j) = 0.0
        end do
      end do
    end if

    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    ! Robin BC for Schottky contacts
    if (allocated(this%efield)) then
      allocate(idx(this%par%g%idx_dim))
      dens_arr = this%dens%get()

      ! reset non-constant Jacobian entries before re-adding
      call this%jaco_dens%reset(const = .false., nonconst = .true.)

      j = this%par%transport_vct(0)%n
      do ict = 1, this%par%nct
        do i = 1, this%par%transport_vct(ict)%n
          j = j + 1
          idx = this%par%transport_vct(ict)%get_idx(i)

          if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
            normal_dir = this%schottky_normal(ict)
            if (normal_dir > 0) then
              efield_arr = this%efield(normal_dir)%get()
              E_normal = efield_arr(j)
              dens = dens_arr(j)

              A_ct = this%par%get_ct_surf(ict, idx)
              v_surf_ict = this%v_surf(ict)

              ! equilibrium injection density (with optional IFBL)
              call schottky_n0b(this%par, this%ci, ict, E_normal, n0B)

              ! tunneling current and its derivative
              call schottky_tunneling(this%par, this%ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)

              ! Jacobian: dF/d(dens) = v_surf*A_ct - A_ct*dJ_tn/ddens
              call this%jaco_dens%add(idx, idx, v_surf_ict * A_ct - A_ct * dJ_tn_ddens)

              ! residual: F = div(J) - A_ct * [v_surf*(n - n0B) + J_tn]
              tmp(j) = tmp(j) + v_surf_ict * A_ct * dens - v_surf_ict * A_ct * n0B - A_ct * J_tn
            end if
          else
            ! Ohmic/Gate: Dirichlet BC (re-add constant 1.0 after reset)
            call this%jaco_dens%add(idx, idx, 1.0)
          end if
        end do
      end do

      ! materialize non-constant Jacobian entries
      call this%jaco_dens%set_matr(const = .false., nonconst = .true.)
    end if

    call this%f%set(tmp - this%b)
  end subroutine

end module
