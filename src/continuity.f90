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
    type(vselector)              :: iref
      !! main variable: quasi-fermi potential
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities in 2 directions
    type(vselector)              :: genrec
      !! generation - recombination

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
    logical              :: status
    real                 :: surf, F, dF

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
    print "(A)", "DEBUG: Initializing variable selectors for continuity equation"
    print "(A,I0)", "DEBUG: Carrier index ci = ", ci
    print "(A,I0)", "DEBUG: Number of contacts nct = ", par%nct
    print "(A,I0)", "DEBUG: Grid dimensions idx_dim = ", idx_dim

    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])

    do idx_dir = 1, idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
      print "(A,I0,A)", "DEBUG: Initialized current density selector for direction ", idx_dir, " with edge transport region"
    end do

    if (par%smc%incomp_ion) then
      call this%genrec%init(genrec, par%ionvert(ci))
      print "(A,I0)", "DEBUG: Incomplete ionization enabled - initialized generation/recombination selector for carrier ", ci
    else
      print "(A)", "DEBUG: Incomplete ionization disabled - skipping genrec initialization"
    end if


    print "(A)", "************** Check from Here **************"

    ! init residuals using this%dens or this%iref as main variable
    print "(A)", "DEBUG: Initializing residual calculation method"
    if (stat) then
      print "(A)", "DEBUG: Stationary mode - using quasi-Fermi potential (iref) as main variable"
      call this%iref%init(iref, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
      call this%init_f(this%iref)
      print "(A)", "DEBUG: Initialized residual function with iref selector"
    else
      print "(A)", "DEBUG: Transient mode - using carrier density (dens) as main variable"
      call this%init_f(this%dens)
      print "(A)", "DEBUG: Initialized residual function with density selector"
    end if

    ! init stencils
    print "(A)", "DEBUG: Initializing finite difference stencils"
    print "(A,A)", "DEBUG: Grid name = ", par%g%name
    print "(A,I0)", "DEBUG: Grid spatial dimensions (dim) = ", par%g%dim
    print "(A,I0)", "DEBUG: Grid index dimensions (idx_dim) = ", par%g%idx_dim
    print "(A,I0)", "DEBUG: Cell vertices per cell = ", par%g%cell_nvert
    print "(A,I0)", "DEBUG: Maximum edges per cell = ", par%g%max_cell_nedge

    print "(A)", "DEBUG: About to initialize st_dir (dirichlet_stencil)"
    print "(A)", "DEBUG: st_dir will create identity mapping for diagonal Jacobian terms"
    call this%st_dir%init(par%g)
    print "(A)", "DEBUG: st_dir initialized with following properties:"
    print "(A,I0)", "DEBUG: st_dir%nmax (max dependencies per point) = ", this%st_dir%nmax
    if (allocated(this%st_dir%perm)) then
      print "(A,*(I0,:,','))", "DEBUG: st_dir%perm (permutation array) = ", this%st_dir%perm
    else
      print "(A)", "DEBUG: st_dir%perm not allocated"
    end if
    if (allocated(this%st_dir%off1)) then
      print "(A,*(I0,:,','))", "DEBUG: st_dir%off1 (start offsets) = ", this%st_dir%off1
      print "(A,*(I0,:,','))", "DEBUG: st_dir%off2 (end offsets) = ", this%st_dir%off2
    else
      print "(A)", "DEBUG: st_dir offset arrays not allocated"
    end if
    print "(A)", "DEBUG: st_dir purpose: handles identity mapping (idx -> idx) for time derivatives"

    do idx_dir = 1, idx_dim
      print "(A,I0)", "DEBUG: About to initialize st_nn for spatial direction ", idx_dir
      print "(A)", "DEBUG: st_nn maps: IDX_VERTEX (grid points) -> IDX_EDGE (flux locations)"
      print "(A,I0)", "DEBUG: Source: IDX_VERTEX = ", IDX_VERTEX
      print "(A,I0)", "DEBUG: Target: IDX_EDGE = ", IDX_EDGE
      print "(A,I0)", "DEBUG: Direction = ", idx_dir
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
      print "(A,I0,A,I0)", "DEBUG: st_nn(", idx_dir, ")%nmax = ", this%st_nn(idx_dir)%nmax
      print "(A,I0,A,I0)", "DEBUG: st_nn(", idx_dir, ")%idx1_type = ", this%st_nn(idx_dir)%idx1_type
      print "(A,I0,A,I0)", "DEBUG: st_nn(", idx_dir, ")%idx1_dir = ", this%st_nn(idx_dir)%idx1_dir
      print "(A,I0,A,I0)", "DEBUG: st_nn(", idx_dir, ")%idx2_type = ", this%st_nn(idx_dir)%idx2_type
      print "(A,I0,A,I0)", "DEBUG: st_nn(", idx_dir, ")%idx2_dir = ", this%st_nn(idx_dir)%idx2_dir
      print "(A,I0,A)", "DEBUG: st_nn(", idx_dir, ") purpose: finds neighbor vertices for each edge flux"
    end do

    call this%st_em%init()
    print "(A)", "DEBUG: Initialized empty stencil st_em"

    ! dependencies
    print "(A)", "DEBUG: Setting up equation dependencies"
    idens = this%depend(this%dens)
    print "(A,I0)", "DEBUG: Density dependency index = ", idens
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
      print "(A,I0,A,I0)", "DEBUG: Current density dependency index for direction ", idx_dir, " = ", icdens(idx_dir)
    end do
    if (par%smc%incomp_ion) then
      igenrec = this%depend(this%genrec)
      print "(A,I0)", "DEBUG: Generation/recombination dependency index = ", igenrec
    else
      print "(A)", "DEBUG: No generation/recombination dependency (complete ionization)"
    end if

    ! init jacobians
    print "(A)", "DEBUG: ==================== JACOBIAN INITIALIZATION ===================="
    print "(A)", "DEBUG: Creating Jacobian matrices for continuity equation"
    print "(A)", "DEBUG: Jacobian represents ∂(continuity_residual)/∂(variables)"

    ! Density Jacobian - ∂(continuity)/∂(density)
    print "(A,I0,A,I0)", "DEBUG: init_jaco_f idens=", idens, " with stencil array size=", 1+par%nct
    print "(A,I0,A,I0)", "DEBUG: v1 (result): ntab=", this%f%ntab, " nval=", this%f%nval
    print "(A,I0,A,I0)", "DEBUG: v2 (depend): ntab=", this%vdep%d(idens)%p%ntab, " nval=", this%vdep%d(idens)%p%nval
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (this%st_dir%get_ptr(), ict = 1, par%nct)], &
      & const = .true., dtime = .false.)
    print "(A,I0,A,I0)", "DEBUG: jaco_dens allocated: matr blocks=", size(this%jaco_dens%matr%const,1), " x ", size(this%jaco_dens%matr%const,2)
    print "(A,L1,A,L1,A,L1)", "DEBUG: block(1,1) const=", this%jaco_dens%matr%const(1,1), " zero=", this%jaco_dens%matr%zero(1,1), " dense=", this%jaco_dens%matr%dense(1,1)
    if (allocated(this%jaco_dens%matr%s)) then
      print "(A,I0,A,I0)", "DEBUG: sparse blocks allocated: ", size(this%jaco_dens%matr%s,1), " x ", size(this%jaco_dens%matr%s,2)
    end if
    if (allocated(this%jaco_dens%st)) then
      print "(A,I0)", "DEBUG: stencils stored: ", size(this%jaco_dens%st)
    end if

    if (.not. stat) then
      print "(A)", "DEBUG: Creating transient jaco_dens_t with dtime=.true."
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
      print "(A,I0,A,I0)", "DEBUG: jaco_dens_t allocated: matr blocks=", size(this%jaco_dens_t%matr%const,1), " x ", size(this%jaco_dens_t%matr%const,2)
      print "(A,L1)", "DEBUG: jaco_dens_t%const(1,1) = ", this%jaco_dens_t%matr%const(1,1)
    else
      print "(A)", "DEBUG: stat=.true. - skipping transient jacobian"
    end if

    ! Current density Jacobians - ∂(continuity)/∂(J_x,J_y,J_z)
    do idx_dir = 1, idx_dim
      print "(A,I0,A,I0,A,I0)", "DEBUG: jaco_cdens(", idx_dir, ") icdens=", icdens(idx_dir), " st_nn%nmax=", this%st_nn(idx_dir)%nmax
      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = [this%st_nn(idx_dir)%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
      print "(A,I0,A,I0,A,I0)", "DEBUG: jaco_cdens(", idx_dir, ") allocated: blocks=", &
        & size(this%jaco_cdens(idx_dir)%p%matr%const,1), " x ", size(this%jaco_cdens(idx_dir)%p%matr%const,2)
      if (allocated(this%jaco_cdens(idx_dir)%p%sd)) then
        print "(A,I0,A,I0)", "DEBUG: jaco_cdens(", idx_dir, ") sparse_data allocated: ", size(this%jaco_cdens(idx_dir)%p%sd)
      end if
      print "(A,I0,A,L1)", "DEBUG: jaco_cdens(", idx_dir, ")%dense(1,1)=", this%jaco_cdens(idx_dir)%p%matr%dense(1,1)
    end do

    if (par%smc%incomp_ion) then
      print "(A,I0)", "DEBUG: Creating jaco_genrec with igenrec=", igenrec
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
      print "(A,I0,A,I0)", "DEBUG: jaco_genrec allocated: blocks=", size(this%jaco_genrec%matr%const,1), " x ", size(this%jaco_genrec%matr%const,2)
      print "(A,L1)", "DEBUG: jaco_genrec%const(1,1)=", this%jaco_genrec%matr%const(1,1)
    else
      print "(A)", "DEBUG: par%smc%incomp_ion=.false. - skipping genrec jacobian"
    end if

    print "(A,I0,A,I0,A,I0,A)", "DEBUG: Total Jacobians created: jaco_dens + ", &
      & merge(1,0,.not.stat), " jaco_dens_t + ", idx_dim, " jaco_cdens + ", merge(1,0,par%smc%incomp_ion), " jaco_genrec"

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
    call this%f%set(tmp - this%b)
  end subroutine

end module
