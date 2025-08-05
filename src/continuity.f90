module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL, CR_ELEC, CR_HOLE
  use contact_m,         only: CT_SCHOTTKY, CT_OHMIC, CT_GATE
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

    real, allocatable :: b(:)
      !! right-hand side
    real, allocatable :: b_edge(:,:)
      !! edge contributions to RHS for Robin BC (edge_index, direction)

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

  subroutine continuity_init(this, par, stat, dens, iref, cdens, genrec)
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

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    integer              :: ict1, ict2
    logical              :: status
    real                 :: surf, F, dF, dx, S, D, n0
    type(stencil_ptr), allocatable :: st_cdens(:)

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
    if (.not. stat) then
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
    end if
    ! Build stencil arrays based on contact types
    do idx_dir = 1, idx_dim
      ! Create stencil array: index 0 for uncontacted, 1:nct for contacts
      allocate(st_cdens(0:par%nct))
      
      ! Uncontacted vertices always use near_neighb stencil
      st_cdens(0) = this%st_nn(idx_dir)%get_ptr()
      
      ! For each contact, choose stencil based on contact type
      do ict = 1, par%nct
        select case (par%contacts(ict)%type)
        case (CT_SCHOTTKY)
          ! Schottky contacts need edge contributions for Robin BC
          st_cdens(ict) = this%st_nn(idx_dir)%get_ptr()
        case (CT_OHMIC, CT_GATE)
          ! Ohmic/Gate contacts use empty stencil (no edge contributions)
          st_cdens(ict) = this%st_em%get_ptr()
        case default
          ! Default to empty for unknown contact types
          st_cdens(ict) = this%st_em%get_ptr()
        end select
      end do
      
      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = st_cdens, const = .true., dtime = .false.)
        
      deallocate(st_cdens)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! allocate edge RHS storage for Robin BC
    allocate (this%b_edge(maxval([(par%transport(IDX_EDGE,idx_dir)%n, idx_dir = 1, idx_dim)]), idx_dim), &
              source = 0.0)

    ! set current density jacobian entries
    do idx_dir = 1, idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        ! get indices of edge, first and second vertex
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get edge surface and length
        surf = par%tr_surf(idx_dir)%get(idx)
        dx = par%g%get_len(idx, idx_dir)

        ! check contact types
        ict1 = par%ict%get(idx1)
        ict2 = par%ict%get(idx2)

        ! handle different vertex combinations
        if (ict1 == 0 .and. ict2 == 0) then
          ! both uncontacted: standard flux continuity
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        
        else if (ict1 > 0 .and. ict2 == 0) then
          if (par%contacts(ict1)%type == CT_SCHOTTKY) then
            ! Schottky contact (vertex 1) to interior (vertex 2): Robin BC
            n0 = calculate_equilibrium_density(par, ci, ict1)
            S = calculate_thermionic_velocity(par, ci, ict1)
            D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! D = mobility in normalized units
            
            ! Robin BC: -D(n2 - n1)/dx = S(n1 - n0)
            ! Rearrange: D/dx * n2 - (D/dx + S) * n1 = -S * n0
            call this%jaco_cdens(idx_dir)%p%set(idx1, idx, -(D/dx + S) * surf)  ! contact vertex
            call this%jaco_cdens(idx_dir)%p%set(idx2, idx,  (D/dx) * surf)      ! interior vertex
            
            ! store RHS contribution
            this%b_edge(i, idx_dir) = -S * n0 * surf
          end if
          
        else if (ict2 > 0 .and. ict1 == 0) then
          if (par%contacts(ict2)%type == CT_SCHOTTKY) then
            ! Interior (vertex 1) to Schottky contact (vertex 2): Robin BC
            n0 = calculate_equilibrium_density(par, ci, ict2)
            S = calculate_thermionic_velocity(par, ci, ict2)
            D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! D = mobility in normalized units
            
            ! Robin BC with opposite orientation
            call this%jaco_cdens(idx_dir)%p%set(idx1, idx, (D/dx) * surf)       ! interior vertex
            call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -(D/dx + S) * surf)  ! contact vertex
            
            ! store RHS contribution (opposite sign due to orientation)
            this%b_edge(i, idx_dir) = S * n0 * surf
          end if
          
        else
          ! other cases: keep standard behavior for uncontacted vertices
          if (ict1 == 0) call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
          if (ict2 == 0) call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
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

    ! dirichlet conditions
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      select case (par%contacts(ict)%type)
      case (CT_SCHOTTKY)
        ! Skip Dirichlet BC for Schottky contacts - handled by Robin BC in edge assembly
        j = j + par%transport_vct(ict)%n
      case default  ! CT_OHMIC, CT_GATE
        ! Apply Dirichlet BC for Ohmic/Gate contacts
        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)
          call this%jaco_dens%set(idx1, idx1, 1.0)
          ! Ohmic contact: infinite injection (equilibrium density)
          if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
            this%b(j) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
          else
            call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
            this%b(j) = par%smc%edos(ci) * F
          end if
        end do
      end select
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer           :: idx_dir, i, j, ict1, ict2, j1, j2
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    logical           :: status
    real, allocatable :: tmp(:)

    allocate (tmp(this%f%n))
    allocate (idx(this%par%g%idx_dim), idx1(this%par%g%idx_dim), idx2(this%par%g%idx_dim))
    
    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)
    
    ! Add edge contributions from Robin BC for Schottky contacts
    do idx_dir = 1, this%par%g%idx_dim
      do i = 1, this%par%transport(IDX_EDGE,idx_dir)%n
        if (abs(this%b_edge(i, idx_dir)) > 0.0) then  ! Only process if there's a contribution
          idx = this%par%transport(IDX_EDGE,idx_dir)%get_idx(i)
          call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
          call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)
          
          ict1 = this%par%ict%get(idx1)
          ict2 = this%par%ict%get(idx2)
          
          ! Get indices in solution vector using vselector
          if (ict1 == 0) then
            j1 = this%dens%itab%get(idx1)
            if (j1 > 0 .and. j1 <= this%par%transport_vct(0)%n) then
              tmp(j1) = tmp(j1) - this%b_edge(i, idx_dir)
            end if
          end if
          
          if (ict2 == 0) then
            j2 = this%dens%itab%get(idx2)
            if (j2 > 0 .and. j2 <= this%par%transport_vct(0)%n) then
              tmp(j2) = tmp(j2) - this%b_edge(i, idx_dir)
            end if
          end if
        end if
      end do
    end do
    
    call this%f%set(tmp - this%b)
  end subroutine

  function calculate_equilibrium_density(par, ci, ict) result(n0)
    !! Calculate equilibrium carrier density at Schottky contact
    !! n₀ = N_c * exp(-q*Φ_B/kT)
    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in) :: ict  ! contact index
    real                            :: n0
    
    real :: phi_B
    
    phi_B = par%contacts(ict)%barrier_height
    
    ! Use same pattern as existing ohmic contact calculation
    if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
      n0 = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(-CR_CHARGE(ci) * phi_B)
    else
      call par%smc%get_dist(-CR_CHARGE(ci) * phi_B, 0, n0, phi_B)  ! use phi_B as dummy for dF
      n0 = par%smc%edos(ci) * n0
    end if
  end function

  function calculate_thermionic_velocity(par, ci, ict) result(S)
    !! Calculate thermionic emission velocity
    !! S = (A*T²)/(q*N_c) * exp(-q*Φ_B/kT)
    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index
    integer,             intent(in) :: ict  ! contact index  
    real                            :: S
    
    real :: A_star, phi_B, T, N_c, n0
    
    A_star = par%contacts(ict)%richardson_const  ! A/cm²/K²
    phi_B = par%contacts(ict)%barrier_height     ! eV
    T = par%T                                    ! K
    N_c = par%smc%edos(ci)                      ! cm⁻³
    
    ! Get equilibrium density (includes exp(-q*Φ_B/kT) factor)  
    n0 = calculate_equilibrium_density(par, ci, ict)
    
    ! S = A*T²/(q*N_c) * exp(-q*Φ_B/kT) = A*T²*n0/(q*N_c²)
    ! But need to be careful about charge sign and units
    S = A_star * T**2 * n0 / (abs(CR_CHARGE(ci)) * N_c)
    
    ! Note: This gives S in units of cm/s if A* is in A/cm²/K² and other units are consistent
  end function

end module
