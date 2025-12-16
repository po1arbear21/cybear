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
  use res_equation_m,    only: res_equation
  use schottky_m,        only: get_normal_dir, schottky_velocity, schottky_n0b, schottky_tunneling
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL, CR_ELEC, CR_HOLE
  use contact_m,         only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use normalization_m,   only: norm, denorm
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
    type(vselector)              :: iref
      !! main variable: quasi-fermi potential
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities in 2 directions
    type(vselector)              :: genrec
      !! generation - recombination

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
    type(jacobian),     pointer     :: jaco_iref     => null()
      !! Direct Jacobian for iref (Schottky contacts only, bypasses chain rule)
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, stat, dens, iref, cdens, genrec, efield)
    !! initialize continuity equation
    class(continuity),              intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    logical,                        intent(in)  :: stat
      !! stationary? if true, use imref as main variable
    type(density),                  intent(in)  :: dens
      !! electron/hole density
    type(imref),                    intent(in)  :: iref
      !! electron/hole quasi-fermi potential
    type(current_density),          intent(in)  :: cdens(:)
      !! electron/hole current density
    type(generation_recombination), intent(in)  :: genrec
      !! generation-recombination rate
    type(electric_field), optional, intent(in)  :: efield(:)
      !! electric field components (for Schottky, optional)

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j, dir
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status, has_schottky
    real                 :: surf, F, dF
    type(stencil_ptr), allocatable :: st_dens_ct(:), st_dens_t_ct(:), st_cdens_ct(:)

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    if (stat) then
      call this%equation_init(CR_NAME(dens%ci)//"continuity_stat")
    else
      call this%equation_init(CR_NAME(dens%ci)//"continuity")
    end if
    this%par => par
    this%ci  = ci

    idx_dim = par%g%idx_dim
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), icdens(idx_dim))
    allocate (this%cdens(idx_dim), this%st_nn(idx_dim), this%jaco_cdens(idx_dim))

    ! Check for Schottky contacts early (needed for iref init)
    has_schottky = any(par%contacts(1:par%nct)%type == CT_SCHOTTKY)

    ! Initialize Schottky support if needed
    if (has_schottky) then
      if (.not. present(efield)) then
        call program_error("E-field required for Schottky contacts")
      end if

      ! Store normal directions and thermionic velocity for each contact
      allocate(this%schottky_normal(par%nct))
      allocate(this%accum_TE(par%nct), source=0.0)
      allocate(this%accum_TN(par%nct), source=0.0)
      allocate(this%v_surf(par%nct), source=0.0)
      do ict = 1, par%nct
        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          this%schottky_normal(ict) = get_normal_dir(par, ict)
          this%v_surf(ict) = schottky_velocity(par, ci, ict)
        else
          this%schottky_normal(ict) = 0
        end if
      end do

      ! Store E-field selectors
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

    ! init residuals using this%dens or this%iref as main variable
    if (stat) then
      call this%iref%init(iref, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
      call this%init_f(this%iref)
    else
      call this%init_f(this%dens)
    end if

    ! init stencils
    call this%st_dir%init(par%g)
    do idx_dir = 1, idx_dim
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
    end do
    call this%st_em%init()

    ! setup contact-specific stencils
    allocate(st_cdens_ct(par%nct), st_dens_ct(par%nct), st_dens_t_ct(par%nct))

    ! Setup contact-specific stencils by contact type
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        ! Schottky: use jaco_dens directly (updated each iteration)
        st_dens_ct(ict) = this%st_dir%get_ptr()
        st_dens_t_ct(ict) = this%st_dir%get_ptr()  ! Time derivative term for Schottky
      else
        st_dens_ct(ict) = this%st_dir%get_ptr()    ! Dirichlet BC for Ohmic/Gate
        st_dens_t_ct(ict) = this%st_em%get_ptr()   ! No time term for Ohmic/Gate
      end if
    end do

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians
    ! jaco_dens: diagonal stencil, non-constant for Schottky (updated each iteration)
      this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (st_dens_ct(ict), ict = 1, par%nct)], &
      & const = .false., dtime = .false.)

    if (.not. stat) then
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
    end if
    do idx_dir = 1, idx_dim
      ! setup contact-specific stencils for this direction
      do ict = 1, par%nct
        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          st_cdens_ct(ict) = this%st_nn(idx_dir)%get_ptr()  ! Edge contributions for Schottky
        else
          st_cdens_ct(ict) = this%st_em%get_ptr()           ! No edge contributions for Ohmic/Gate
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
        ! get indices of edge, first and second vertex
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get edge surface
        surf = par%tr_surf(idx_dir)%get(idx)

        ! set values for interior vertices and Schottky contacts
        ! (Ohmic/Gate use Dirichlet BC, so no current divergence term)
        ict = par%ict%get(idx1)
        if (ict == 0) then
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx, surf)
        else if (par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx, surf)
        end if
        ict = par%ict%get(idx2)
        if (ict == 0) then
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        else if (par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        end if
      end do
    end do

    ! density time derivative factor
    if (.not. stat) then
      do i = 1, par%transport_vct(0)%n
        idx1 = par%transport_vct(0)%get_idx(i)
        call this%jaco_dens_t%set(idx1, idx1, par%tr_vol%get(idx1))
      end do
    end if

    ! generation-recombination
    if (par%smc%incomp_ion) then
      do i = 1, par%ionvert(ci)%n
        idx1 = par%ionvert(ci)%get_idx(i)
        call this%jaco_genrec%set(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! boundary conditions: Dirichlet for Ohmic/Gate, handled in eval for Schottky
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      do i = 1, par%transport_vct(ict)%n
        j = j + 1
        idx1 = par%transport_vct(ict)%get_idx(i)

        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%jaco_dens%set(idx1, idx1, 1.0)
          ! Schottky: FVM flux BC with Jacobian via compensation
          ! jaco_dens entry set in eval(), residual compensated to undo mul_vec effect
          ! Equation: Σ(J_ij·A_ij) - J_sch·A_ct = 0, b=0

        else if (par%contacts(ict)%type == CT_OHMIC .or. par%contacts(ict)%type == CT_GATE) then
          ! Dirichlet BC for Ohmic/Gate contacts
          call this%jaco_dens%set(idx1, idx1, 1.0)
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
    real                 :: n0B, J_tn, dJ_tn_ddens, v_surf_ict, J_te

    allocate (tmp(this%f%n))
    allocate (idx(this%par%g%idx_dim))

    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)

    ! print "(A)", "DEBUG_CONTINUITY: Node structure"
    ! print "(A,I6)", "  n_interior = ", this%par%transport_vct(0)%n
    ! do ict = 1, this%par%nct
    !   print "(A,I3,A,I6,A,I3)", "  contact ", ict, ": n_nodes=", this%par%transport_vct(ict)%n, &
    !                            ", type=", this%par%contacts(ict)%type
    ! end do
    ! dens_arr = this%dens%get()
    ! print "(A)", "DEBUG_CONTINUITY: after jaco_dens*dens (tmp = J_dens * n)"
    ! print "(A)", "  node | dens | tmp | J_diag (=tmp/dens)"
    ! do i = 1, size(tmp)
    !   print "(A,I6,A,3ES14.6)", "  ", i, " |", denorm(dens_arr(i), "cm^-3"), tmp(i), tmp(i)/dens_arr(i)
    ! end do

    ! Zero out jaco_dens contribution for Schottky contacts (before adding cdens)
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

    ! Add Schottky current contribution and set Jacobian
    if (allocated(this%efield)) then
      dens_arr = this%dens%get()

      j = this%par%transport_vct(0)%n
      do ict = 1, this%par%nct
        do i = 1, this%par%transport_vct(ict)%n
          j = j + 1
          idx = this%par%transport_vct(ict)%get_idx(i)

          if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
            normal_dir = this%schottky_normal(ict)
            if (normal_dir > 0) then
              ! Get E-field and density
              efield_arr = this%efield(normal_dir)%get()
              E_normal = efield_arr(j)
              dens = dens_arr(j)

              ! Get contact surface area and thermionic velocity
              A_ct = this%par%get_ct_surf(ict, idx)
              v_surf_ict = this%v_surf(ict)

              ! Robin BC + Tunneling approach:
              ! J_schottky = v_surf*(n - n0B) + J_tn
              ! F = div(J) - A_ct * J_schottky = 0

              ! Compute equilibrium injection density n0B (with optional IFBL)
              call schottky_n0b(this%par, this%ci, ict, E_normal, n0B)

              ! Compute tunneling current and its derivative
              call schottky_tunneling(this%par, this%ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)

              ! ! Numerical verification of dJ_tn_ddens
              ! block
              !   real :: J_tn_plus, J_tn_minus, dJ_numerical, eps_fd, dummy
              !   real :: jaco_full, edge_flux_before
              !   eps_fd = 1e-6 * max(abs(dens), 1e-12)
              !   call schottky_tunneling(this%par, this%ci, ict, E_normal, dens + eps_fd, J_tn_plus, dummy)
              !   call schottky_tunneling(this%par, this%ci, ict, E_normal, dens - eps_fd, J_tn_minus, dummy)
              !   dJ_numerical = (J_tn_plus - J_tn_minus) / (2.0 * eps_fd)

              !   edge_flux_before = tmp(j)  ! Edge flux contribution BEFORE Schottky term
              !   jaco_full = v_surf_ict * A_ct - A_ct * dJ_tn_ddens

              !   print "(A,I5,A,I2)", "FD_CHECK j=", j, " ci=", this%ci
              !   print "(A,ES12.4,A,ES12.4,A,F8.4)", &
              !     "  dJ_tn/dp: analytical=", dJ_tn_ddens, " numerical=", dJ_numerical, &
              !     " ratio=", dJ_tn_ddens / (dJ_numerical + 1e-30)
              !   print "(A,ES12.4,A,ES12.4,A,ES12.4)", &
              !     "  v_surf*A_ct=", v_surf_ict * A_ct, " -A_ct*dJ_tn/dp=", -A_ct * dJ_tn_ddens, &
              !     " jaco_set=", jaco_full
              !   print "(A,ES12.4)", "  edge_flux_before_schottky=", edge_flux_before
              ! end block

              ! TE component for debugging (J_tn already computed above)
              J_te = -v_surf_ict * (dens - n0B)

              ! ! Debug output
              ! print "(A,I3,A,I3)", "DEBUG_SCHOTTKY_ROBIN: contact=", ict, " vertex=", i
              ! print "(A,2ES14.6)", "  dens, n0B = ", denorm(dens, "cm^-3"), denorm(n0B, "cm^-3")
              ! print "(A,ES14.6,A)", "  v_surf = ", denorm(v_surf_ict, "cm/s"), " cm/s"
              ! print "(A,ES14.6,A)", "  J_TE = -v*(n-n0B) = ", denorm(J_te, "A/cm^2"), " A/cm²"
              ! print "(A,ES14.6,A)", "  J_TN = ", denorm(J_tn, "A/cm^2"), " A/cm²"
              ! print "(A,ES14.6,A)", "  J_total = ", denorm(J_te + J_tn, "A/cm^2"), " A/cm²"
              ! print "(A,ES14.6,A)", "  dJ/ddens = ", v_surf_ict * A_ct - A_ct * dJ_tn_ddens

              ! Set Jacobian entry: dF/d(dens) = +v_surf*A_ct - A_ct*dJ_tn/ddens
              call this%jaco_dens%set(idx, idx, v_surf_ict * A_ct - A_ct * dJ_tn_ddens)

              ! ! Debug: print what we're setting vs what might already be there
              ! if (j == 3365) then
              !   block
              !     real :: jaco_schottky, residual_F, expected_dx
              !     jaco_schottky = v_surf_ict * A_ct - A_ct * dJ_tn_ddens
              !     residual_F = tmp(j) + v_surf_ict * A_ct * dens - v_surf_ict * A_ct * n0B - A_ct * J_tn
              !     expected_dx = residual_F / jaco_schottky  ! If diagonal-only

              !     print "(A)", "=== JACOBIAN DEBUG j=3365 ==="
              !     print "(A,ES12.4)", "  jaco_schottky (what we set) = ", jaco_schottky
              !     print "(A,ES12.4)", "  Residual F = ", residual_F
              !     print "(A,ES12.4)", "  Expected dx = F/J = ", expected_dx
              !     print "(A,ES12.4)", "  (Jacobian test showed actual dx = -4.10E-14)"
              !     print "(A)", "  If expected ≈ actual but still overshoots → coupling issue"
              !   end block
              ! end if

              ! Add residual: F = div(J) - A_ct * [v_surf*(n - n0B) + J_tn]
              ! The div(J) is already in tmp from jaco_cdens
              ! jaco_dens contribution was zeroed earlier, so we add the full BC term:
              tmp(j) = tmp(j) + v_surf_ict * A_ct * dens - v_surf_ict * A_ct * n0B - A_ct * J_tn

            end if
          else
            ! Ohmic/Gate: Dirichlet BC, jaco_dens = 1.0
            call this%jaco_dens%set(idx, idx, 1.0)
          end if
        end do
      end do

      ! Materialize the non-constant Jacobian entries
      call this%jaco_dens%set_matr(const = .false., nonconst = .true.)
    end if

    call this%f%set(tmp - this%b)
    ! print *, "Max residual at vertex:", maxloc(abs(tmp - this%b)), maxval(abs(tmp - this%b))
  end subroutine

end module
