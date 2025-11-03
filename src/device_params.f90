m4_include(util/macro.f90.inc)

module device_params_m

  use bin_search_m,     only: bin_search
  use contact_m,        only: CT_OHMIC, CT_GATE, CT_SCHOTTKY, contact
  use container_m,      only: container, STORAGE_WRITE
  use error_m,          only: assert_failed, program_error
  use galene_m,         only: gal_file, gal_block, GALDATA_CHAR, GALDATA_INT, GALDATA_REAL
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid, grid_ptr
  use grid_data_m,      only: allocate_grid_data0_real, allocate_grid_data1_real, allocate_grid_data2_real, &
  &                         allocate_grid_data3_real, allocate_grid_data0_int, allocate_grid_data1_int, &
  &                         grid_data_int, grid_data_real
  use grid_generator_m, only: DIR_NAME, generate_1D_grid m4_ifdef({m4_triangle},{, generate_triangle_grid})
  use grid_table_m,     only: grid_table
  use grid1D_m,         only: grid1D
  use input_m,          only: input_file
  use map_m,            only: map_string_int, mapnode_string_int
  use normalization_m,  only: norm, denorm
  use region_m,         only: region_ptr, region_poisson, region_transport, region_doping, region_mobility, region_contact
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON, CR_NAME, DOP_NAME, DOS_PARABOLIC, DOS_PARABOLIC_TAIL, DOS_GAUSS, DIST_MAXWELL, DIST_FERMI, DIST_FERMI_REG, semiconductor
  use string_m,         only: string, new_string
  use tensor_grid_m,    only: tensor_grid
  use triang_grid_m,    only: triang_grid
  m4_ifdef({m4_triangle},{  use triangle_m,       only: triangulation})
  use util_m,           only: int2str, split_string

  implicit none

  private
  public device_params

  type device_params
    !! device geometry and material parameters

    real :: T
      !! temperature in K

    integer             :: ci0, ci1
      !! enabled carrier index range (maximal: CR_ELEC..CR_HOLE)
    type(semiconductor) :: smc
      !! semiconductor charge carrier parameters

    real :: curr_fact
      !! current factor for converting A/m² or A/m to A

    class(grid_data_real), allocatable :: eps(:,:)
      !! electrostatic permittivity (idx_type, idx_dir)
    class(grid_data_real), allocatable :: surf(:)
      !! adjoint volume surfaces
    class(grid_data_real), allocatable :: dop(:,:,:)
      !! donor/acceptor concentration (idx_type, idx_dir, dcon/acon)
    class(grid_data_real), allocatable :: ii_E_dop(:)
      !! incomplete ionization donor/acceptor dopant energy level (dcon/acon)
    class(grid_data_real), allocatable :: mob0(:,:,:)
      !! zero-field mobility (idx_type, idx_dir, carrier index)
    class(grid_data_real), allocatable :: tr_surf(:)
      !! adjoint volume surfaces in transport region
    class(grid_data_real), allocatable :: tr_vol
      !! adjoint volumes for transport region
    class(grid_data_real), allocatable :: vol
      !! adjoint volumes

    type(grid_table),      allocatable :: poisson(:,:)
      !! poisson grid tables (idx_type, idx_dir)
    type(grid_table),      allocatable :: oxide(:,:)
      !! oxide grid tables (idx_type, idx_dir)
    type(grid_table),      allocatable :: transport(:,:)
      !! transport grid tables (idx_type, idx_dir)
    type(grid_table)                   :: dopvert(2)
      !! uncontacted transport vertices with doping > 0 (donor, acceptor)
    type(grid_table)                   :: ionvert(2)
      !! uncontacted transport vertices with partial ionization enabled
    type(contact),         allocatable :: contacts(:)
      !! device contacts
    type(map_string_int)               :: contact_map
      !! get contact index by name
    class(grid_data_int),  allocatable :: ict
      !! get contact index by grid indices
    integer                            :: nct
      !! number of contacts
    type(grid_table),      allocatable :: contacted(:)
      !! contacted vertices (1:nct)
    type(grid_table),      allocatable :: poisson_vct(:)
      !! poisson vertices grouped by contacts (0:nct)
    type(grid_table),      allocatable :: oxide_vct(:)
      !! oxide vertices grouped by contacts (0:nct)
    type(grid_table),      allocatable :: transport_vct(:)
      !! transport vertices grouped by contacts (0:nct)

    type(string)              :: gtype
      !! grid type ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    type(grid1D)              :: g1D(3)
      !! x, y, z grids
    type(tensor_grid)         :: tg
      !! tensor grid
    m4_ifdef({m4_triangle},{type(triangulation)       :: tr})
      !! triangles generated from regions
    type(triang_grid)         :: gtr
      !! triangle grid
    class(grid),  pointer     :: g => null()
      !! pointer to grid that is actually used
    character(:), allocatable :: idx_dir_name(:)
      !! name of index direction (e.g. ["x", "y"] for gtype = "xy"; or ["xy", "z"] for gtype = "tr_xyz")

    type(region_ptr),       allocatable :: reg(:)
      !! pointer to all regions
    type(region_poisson),   allocatable :: reg_poiss(:)
      !! poisson regions
    type(region_transport), allocatable :: reg_trans(:)
      !! transport regions
    type(region_doping),    allocatable :: reg_dop(:)
      !! doping regions
    type(region_mobility),  allocatable :: reg_mob(:)
      !! mobility regions
    type(region_contact),   allocatable :: reg_ct(:)
      !! contact regions

    logical              :: gal
      !! use GALENE III file
    type(gal_file)       :: gal_fl
      !! GALENE III file
    type(map_string_int) :: gal_mat
      !! GALENE III materials/contacts
    integer              :: gal_ox, gal_si
      !! GALENE III silicon and oxide material indices
    real                 :: gal_epsox, gal_epssi
      !! GALENE III permittivities
    real                 :: gal_phims
      !! GALENE III phims for gate contacts
  contains
    procedure :: init     => device_params_init
    procedure :: destruct => device_params_destruct
    procedure :: get_ct_surf

    procedure, private :: init_transport_params  => device_params_init_transport_params
    procedure, private :: init_regions           => device_params_init_regions
    procedure, private :: init_grid              => device_params_init_grid
    procedure, private :: init_poisson           => device_params_init_poisson
    procedure, private :: init_transport         => device_params_init_transport
    procedure, private :: init_doping            => device_params_init_doping
    procedure, private :: init_mobility          => device_params_init_mobility
    procedure, private :: init_contacts          => device_params_init_contacts
  end type

  real, parameter :: SURF_TOL = 1e-32
    !! use if adjoint surface is zero

contains

  subroutine device_params_init(this, file, T)
    !! initialize device parameter object
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file
      !! device file
    real,                         intent(in)  :: T
      !! temperature in K

    integer                   :: gal_sid, i, status, ci, idx_dir
    type(string)              :: gal_filename
    type(string), allocatable :: s(:)
    type(gal_block), pointer  :: gblock
    type(container)           :: ctnr

    this%T = T

    call this%init_transport_params(file)

    gal_sid = -1
    call file%get_section("galene", gal_sid, status = status)
    m4_assert(status <= 0) ! at most one galene sections allowed
    this%gal = (status == 0)
    if (this%gal) then
      ! only transport params and galene section allowed
      m4_assert(file%sections%n == 2)

      ! load galene file
      call file%get(gal_sid, "file", gal_filename)
      call this%gal_fl%init(gal_filename%s)

      ! use xy tensor grid
      this%gtype%s = "xy"

      ! get materials
      gblock => this%gal_fl%get_block("name(material)")
      s = split_string(gblock%cdata, [" "])
      call this%gal_mat%init()
      do i = 1, size(s)
        call this%gal_mat%insert(s(i), i)
      end do
      this%gal_ox = this%gal_mat%get(new_string("OX"))
      this%gal_si = this%gal_mat%get(new_string("SIL"))
      call file%get(gal_sid, "epsox", this%gal_epsox)
      call file%get(gal_sid, "epssi", this%gal_epssi)
      call file%get(gal_sid, "phims", this%gal_phims)
    else
      ! get grid type ("x", "xy", "xyz", "tr_xy", "tr_xyz")
      call file%get("grid", "gtype", this%gtype)
      select case (this%gtype%s)
      case ("x", "xy", "xyz", "tr_xy", "tr_xyz")
      case default
        call program_error("Invalid grid type: "//this%gtype%s)
      end select

      ! get regions
      call this%init_regions(file)
    end if

    ! init rest
    call this%init_grid(file)
    call this%init_poisson()
    call this%init_transport()
    call this%init_doping(file)
    call this%init_mobility()
    call this%init_contacts()

    ! output
    call ctnr%open("device.fbs", flag = STORAGE_WRITE)
    ! grid
    call ctnr%save(this%g)
    ! permittivity on cells
    call ctnr%save("perm", this%eps(IDX_CELL, 0), unit = "eps0")
    do ci = this%ci0, this%ci1
      ! doping on cells
      call ctnr%save(DOP_NAME(ci)//"con_cell", this%dop(IDX_CELL, 0, ci), unit = "cm^-3")
      ! doping on vertices
      call ctnr%save(DOP_NAME(ci)//"con_vert", this%dop(IDX_VERTEX, 0, ci), unit = "cm^-3")
      ! dopant ionization energy on vertices
      if (this%smc%incomp_ion) call ctnr%save(DOP_NAME(ci)//"conE_vert", this%ii_E_dop(ci), unit = "eV")
      do idx_dir = 1, this%g%idx_dim
        ! mobility on edges
        call ctnr%save(CR_NAME(ci)//"mob"//DIR_NAME(idx_dir), this%mob0(IDX_EDGE, idx_dir, ci), unit = "cm^2/V/s")
      end do
    end do
    call ctnr%close()
  end subroutine

  subroutine device_params_destruct(this)
    !! destruct device parameter object
    class(device_params), intent(inout) :: this

    call this%contact_map%destruct()
  end subroutine

  subroutine device_params_init_transport_params(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)    :: file

    integer      :: sid
    logical      :: elec, hole, status
    type(string) :: dos, dist

    ! find transport parameters section id
    call file%get_section("transport parameters", sid)

    ! general parameters
    call file%get(sid, "electrons", elec)
    call file%get(sid, "holes",     hole)
    this%ci0 = CR_ELEC
    this%ci1 = CR_HOLE
    if (.not. elec) this%ci0 = CR_HOLE
    if (.not. hole) this%ci1 = CR_ELEC
    call file%get(sid, "N_c0",       this%smc%edos(CR_ELEC))
    call file%get(sid, "N_v0",       this%smc%edos(CR_HOLE))
    this%smc%edos(:) = this%smc%edos(:) * this%T ** 1.5
    call file%get(sid, "E_gap",      this%smc%band_gap)
    this%smc%band_edge(CR_ELEC) =   0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))
    this%smc%band_edge(CR_HOLE) = - 0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))

    ! Debug: Show band edge positions (Ec0 and Ev0) relative to intrinsic level
    print "(A)", "DEBUG: Band edge positions relative to Ei = 0:"
    print "(A,ES12.5,A,ES12.5,A)", "  Ec0 = ", denorm(this%smc%band_edge(CR_ELEC), "eV"), &
                                    " eV (normalized: ", this%smc%band_edge(CR_ELEC), ")"
    print "(A,ES12.5,A,ES12.5,A)", "  Ev0 = ", denorm(this%smc%band_edge(CR_HOLE), "eV"), &
                                    " eV (normalized: ", this%smc%band_edge(CR_HOLE), ")"

    ! density of states
    call file%get(sid, "dos", dos, status)
    if (.not. status) then
      dos%s = "parabolic"
      print "(A)", "Assume parabolic density of states"
    end if
    select case (dos%s)
    case ("parabolic")
      this%smc%dos = DOS_PARABOLIC
    case ("parabolic_tail")
      this%smc%dos = DOS_PARABOLIC_TAIL
    case ("gauss")
      this%smc%dos = DOS_GAUSS
    case default
      call program_error("unknown dos '" // dos%s // "'")
    end select
    if (this%smc%dos == DOS_PARABOLIC_TAIL) then
      call file%get(sid, "dos_t0", this%smc%dos_params(1))
    elseif (this%smc%dos == DOS_GAUSS) then
      call file%get(sid, "dos_sigma", this%smc%dos_params(1))
    end if

    ! distribution
    call file%get(sid, "dist", dist, status)
    if (.not. status) then
      dist%s = "maxwell"
      print "(A)", "Assume Maxwell-Boltzmann distribution density"
    end if
    select case (dist%s)
    case ("maxwell")
      this%smc%dist = DIST_MAXWELL
    case ("fermi")
      this%smc%dist = DIST_FERMI
    case ("fermi_reg")
      this%smc%dist = DIST_FERMI_REG
    case default
      call program_error("unknown distribution density '" // dist%s // "'")
    end select
    if (this%smc%dist == DIST_FERMI_REG) then
      call file%get(sid, "reg_A", this%smc%dist_params(1))
      call file%get(sid, "reg_B", this%smc%dist_params(2))
    end if

    ! mobility
    call file%get(sid, "mob",        this%smc%mob)
    call file%get(sid, "alpha",      this%smc%alpha)
    call file%get(sid, "beta",       this%smc%beta)
    call file%get(sid, "mob_min",    this%smc%mob_min)
    call file%get(sid, "mob_max",    this%smc%mob_max)
    call file%get(sid, "N_ref",      this%smc%N_ref)
    call file%get(sid, "v_sat",      this%smc%v_sat)

    call file%get(sid, "curr_fact",  this%curr_fact)

    ! incomplete ionization
    call file%get(sid, "incomp_ion", this%smc%incomp_ion)
    call file%get(sid, "ii_tau",     this%smc%ii_tau)
    call file%get(sid, "ii_E_dop0",  this%smc%ii_E_dop0)
    call file%get(sid, "ii_g",       this%smc%ii_g)
    call file%get(sid, "ii_N_crit",  this%smc%ii_N_crit)
    call file%get(sid, "ii_dop_th",  this%smc%ii_dop_th)

    ! make sure parameters are valid
    m4_assert(this%ci0 <= this%ci1)
    m4_assert(size(this%smc%alpha)     == 2)
    m4_assert(size(this%smc%beta)      == 2)
    m4_assert(size(this%smc%mob_min)   == 2)
    m4_assert(size(this%smc%mob_max)   == 2)
    m4_assert(size(this%smc%N_ref)     == 2)
    m4_assert(size(this%smc%v_sat)     == 2)
    m4_assert(size(this%smc%ii_tau)    == 2)
    m4_assert(size(this%smc%ii_E_dop0) == 2)
    m4_assert(size(this%smc%ii_g)      == 2)
  end subroutine

  subroutine device_params_init_regions(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)    :: file

    integer, allocatable :: sids(:)
    integer              :: si, i, j

    print "(A)", "init_regions"

    ! initialize poisson regions
    call file%get_sections("poisson", sids)
    allocate (this%reg_poiss(size(sids)))
    do si = 1, size(sids)
      call this%reg_poiss(si)%init(file, sids(si), this%gtype)
    end do

    ! initialize transport regions
    call file%get_sections("transport", sids)
    allocate (this%reg_trans(size(sids)))
    do si = 1, size(sids)
      call this%reg_trans(si)%init(file, sids(si), this%gtype)
    end do

    ! initialize doping regions
    call file%get_sections("doping", sids)
    allocate (this%reg_dop(size(sids)))
    do si = 1, size(sids)
      call this%reg_dop(si)%init(file, sids(si), this%gtype)
    end do

    ! initialize mobility regions
    call file%get_sections("mobility", sids)
    allocate (this%reg_mob(size(sids)))
    do si = 1, size(sids)
      call this%reg_mob(si)%init(file, sids(si), this%gtype)
    end do

    ! initialize contact regions
    call file%get_sections("contact", sids)
    allocate (this%reg_ct(size(sids)))
    do si = 1, size(sids)
      call this%reg_ct(si)%init(file, sids(si), this%gtype)
    end do

    ! create region pointer array
    allocate (this%reg(size(this%reg_poiss) + size(this%reg_trans) + size(this%reg_dop) + size(this%reg_ct)))
    i = 0
    do j = 1, size(this%reg_poiss)
      i = i + 1
      this%reg(i)%p => this%reg_poiss(j)
    end do
    do j = 1, size(this%reg_trans)
      i = i + 1
      this%reg(i)%p => this%reg_trans(j)
    end do
    do j = 1, size(this%reg_dop)
      i = i + 1
      this%reg(i)%p => this%reg_dop(j)
    end do
    do j = 1, size(this%reg_ct)
      i = i + 1
      this%reg(i)%p => this%reg_ct(j)
    end do
  end subroutine

  subroutine device_params_init_grid(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)    :: file

    integer        :: dim, dir, ngptr
    logical        :: load, status
    type(grid_ptr) :: gptr(3)

    print "(A)", "init_grid"

    ! grid dimension
    select case(this%gtype%s)
    case("x")
      dim = 1
    case("xy")
      dim = 2
    case("xyz")
      dim = 3
    case("tr_xy")
      dim = 2
    case("tr_xyz")
      dim = 3
    end select

    ! load or generate ?
    if (this%gal) then
      load = .true.
    else
      call file%get("grid", "load", load, status)
      if (.not. status) load = .false.
    end if
    ngptr = 0

    ! generate or load grid
    select case(this%gtype%s)
    case("x", "xy", "xyz")
      do dir = 1, dim
        call generate_1D_grid(file, load, this%gal, this%gal_fl, dir, this%reg, this%g1D(dir), gptr, ngptr)

        print "(A,I0)", "#("//DIR_NAME(dir)//"): ", this%g1D(dir)%n
      end do
      allocate (character(1) :: this%idx_dir_name(dim))
      do dir = 1, dim
        this%idx_dir_name(dir) = DIR_NAME(dir)
      end do

    case("tr_xy", "tr_xyz")
      m4_ifdef({m4_triangle},{      call generate_triangle_grid(file, load, this%reg, this%tr, this%gtr, gptr, ngptr)})
      if (dim == 3) then
        call generate_1D_grid(file, load, .false., this%gal_fl, 3, this%reg, this%g1D(3), gptr, ngptr)
      end if
      allocate (character(2) :: this%idx_dir_name(dim-1))
      this%idx_dir_name(1) = "xy"
      if (dim == 3) this%idx_dir_name(2) = "z"

    end select

    ! create tensor grid if necessary and set this%g pointer
    if (ngptr == 1) then
      this%g => gptr(1)%p
    else
      call this%tg%init("grid", gptr(1:ngptr))
      this%g => this%tg
    end if
  end subroutine

  subroutine device_params_init_poisson(this)
    class(device_params), intent(inout) :: this

    integer                  :: dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable     :: idx(:), idx2(:)
    logical                  :: status
    real                     :: mid(2), p(2,3), eps
    real,    allocatable     :: surf(:,:), vol(:)
    type(string)             :: gtype_tmp
    type(gal_block), pointer :: gblock

    print "(A)", "init_poisson"

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! allocate temp arrays
    allocate (idx(idx_dim), idx2(idx_dim))
    allocate (surf(this%g%max_cell_nedge,idx_dim), vol(this%g%cell_nvert))

    ! allocate/initialize grid data
    call allocate_grid_data2_real(this%eps,  idx_dim, [1, 0], [4, idx_dim])
    call allocate_grid_data1_real(this%surf, idx_dim, 1, idx_dim)
    call allocate_grid_data0_real(this%vol,  idx_dim)
    call this%eps(IDX_CELL,0)%init(this%g, IDX_CELL, 0)
    do idx_dir = 1, idx_dim
      call this%eps(IDX_EDGE,idx_dir)%init(this%g, IDX_EDGE, idx_dir)
      call this%surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
    end do
    call this%vol%init(this%g, IDX_VERTEX, 0)

    ! allocate/initialize poisson grid tables
    allocate (this%poisson(4,0:idx_dim))
    call this%poisson(IDX_VERTEX,0)%init("poisson_v", this%g, IDX_VERTEX, 0)
    call this%poisson(IDX_CELL,  0)%init("poisson_c", this%g, IDX_CELL,   0)
    do idx_dir = 1, idx_dim
      call this%poisson(IDX_EDGE,idx_dir)%init("poisson_e"//trim(this%idx_dir_name(idx_dir)), this%g, IDX_EDGE, idx_dir)
    end do

    ! set poisson cells
    if (this%gal) then
      gblock => this%gal_fl%get_block("desc(volume)")
      if (size(gblock%idata) /= (this%g1D(1)%n-1) * (this%g1D(2)%n-1)) then
        call program_error("Number of quads in GALENE file does not match grid")
      end if
      k = 0
      do j = 1, this%g1D(2)%n-1; do i = 1, this%g1D(1)%n - 1
        k = k + 1
        if (gblock%idata(k) == 0) cycle
        idx = [i, j]
        call this%poisson(IDX_CELL,0)%flags%set(idx, .true.)
        if (gblock%idata(k) == this%gal_si) then
          call this%eps(IDX_CELL,0)%set(idx, this%gal_epssi)
        elseif (gblock%idata(k) == this%gal_ox) then ! OX
          call this%eps(IDX_CELL,0)%set(idx, this%gal_epsox)
        else
          call program_error("(" // int2str(i) // "," // int2str(j) // "): unexpected material " // int2str(gblock%idata(k)))
        end if
      end do; end do
    else
      i0 = 1
      i1 = 1
      do ri = 1, size(this%reg_poiss)
        select case (this%gtype%s)
        case ("x", "xy", "xyz")
          ! get bounds
          do idx_dir = 1, idx_dim
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(ri)%xyz(idx_dir,1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(ri)%xyz(idx_dir,2)) - 1
          end do

          ! loop over cells in line-segment (x), rectangle (xy) or cuboid (xyz)
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            idx = ijk(1:idx_dim)

            ! enable poisson for cell and set permittivity
            call this%poisson(IDX_CELL,0)%flags%set(idx, .true.)
            call this%eps(IDX_CELL,0)%set(idx, this%reg_poiss(ri)%eps)
          end do; end do; end do

        case ("tr_xy", "tr_xyz")
          ! get z bounds
          if (dim == 3) then
            i0(3) = bin_search(this%g1D(3)%x, this%reg_poiss(ri)%xyz(3,1))
            i1(3) = bin_search(this%g1D(3)%x, this%reg_poiss(ri)%xyz(3,2)) - 1
          end if

          ! enable poisson in region
          gtype_tmp = new_string("tr_xy")
          do i = 1, this%gtr%ncell
            ! check if triangle is in region using centroid
            call this%gtr%get_cell([i], p)
            mid(1) = (p(1,1) + p(1,2) + p(1,3)) / 3
            mid(2) = (p(2,1) + p(2,2) + p(2,3)) / 3
            if (.not. this%reg_poiss(ri)%point_test(gtype_tmp, mid)) cycle

            ! loop over z direction
            do k = i0(3), i1(3)
              ijk(1:2) = [i, k]
              idx      = ijk(1:idx_dim)

              ! enable poisson for cell and set permittivity
              call this%poisson(IDX_CELL,0)%flags%set(idx, .true.)
              call this%eps(IDX_CELL,0)%set(idx, this%reg_poiss(ri)%eps)
            end do
          end do
        end select
      end do
    end if
    call this%poisson(IDX_CELL,0)%init_final()

    ! process vertices, edges
    do i = 1, this%poisson(IDX_CELL,0)%n
      idx = this%poisson(IDX_CELL,0)%get_idx(i)

      ! get adjoint cell parts
      call this%g%get_adjoint(idx, surf = surf, vol = vol)
      where (surf == 0) surf = SURF_TOL

      ! enable poisson for vertices and update adjoint volumes
      do j = 1, this%g%cell_nvert
        call this%g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx, j, idx2, status)
        call this%poisson(IDX_VERTEX,0)%flags%set(idx2, .true.)
        call this%vol%update(idx2, vol(j))
      end do

      ! enable poisson for edges and update adjoint surfaces + permittivity
      eps = this%eps(IDX_CELL,0)%get(idx)
      do idx_dir = 1, idx_dim
        do j = 1, this%g%cell_nedge(idx_dir)
          call this%g%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx2, .true.)
          call this%surf(idx_dir)%update(idx2, surf(j,idx_dir))
          call this%eps(IDX_EDGE,idx_dir)%update(idx2, surf(j,idx_dir) * eps)
        end do
      end do
    end do

    ! finish vertex + edge tables
    call this%poisson(IDX_VERTEX,0)%init_final()
    do idx_dir = 1, idx_dim
      call this%poisson(IDX_EDGE,idx_dir)%init_final()
    end do

    ! weighted average of permittivity on adjoint surfaces
    do idx_dir = 1, idx_dim
      do i = 1, this%poisson(IDX_EDGE,idx_dir)%n
        idx = this%poisson(IDX_EDGE,idx_dir)%get_idx(i)
        call this%eps(IDX_EDGE,idx_dir)%set(idx, this%eps(IDX_EDGE,idx_dir)%get(idx) / this%surf(idx_dir)%get(idx))
      end do
    end do

    print "(A,I0)", "#(Poisson vertices): ", this%poisson(IDX_VERTEX,0)%n
  end subroutine

  subroutine device_params_init_transport(this)
    class(device_params), intent(inout) :: this

    integer                  :: dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable     :: idx(:), idx2(:)
    logical                  :: status
    real                     :: mid(2), p(2,3)
    real,    allocatable     :: surf(:,:), vol(:)
    type(string)             :: gtype_tmp
    type(gal_block), pointer :: gblock

    print "(A)", "init_transport"

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! allocate temp arrays
    allocate (idx(idx_dim), idx2(idx_dim))
    allocate (surf(this%g%max_cell_nedge,idx_dim), vol(this%g%cell_nvert))

    ! allocate/initialize grid data
    call allocate_grid_data1_real(this%tr_surf, idx_dim, 1, idx_dim)
    call allocate_grid_data0_real(this%tr_vol, idx_dim)
    do idx_dir = 1, idx_dim
      call this%tr_surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
    end do
    call this%tr_vol%init(this%g, IDX_VERTEX, 0)

    ! allocate/initialize oxide and transport grid tables
    allocate (this%oxide(4,0:this%g%idx_dim))
    allocate (this%transport(4,0:this%g%idx_dim))
    call this%oxide(IDX_VERTEX,0)%init("oxide_v", this%g, IDX_VERTEX, 0)
    call this%oxide(IDX_CELL,  0)%init("oxide_c", this%g, IDX_CELL,   0)
    this%oxide(IDX_VERTEX,0)%flags = this%poisson(IDX_VERTEX,0)%flags
    this%oxide(IDX_CELL,  0)%flags = this%poisson(IDX_CELL,  0)%flags
    call this%transport(IDX_VERTEX,0)%init("transport_v", this%g, IDX_VERTEX, 0)
    call this%transport(IDX_CELL,  0)%init("transport_c", this%g, IDX_CELL,   0)
    do idx_dir = 1, idx_dim
      call this%oxide(    IDX_EDGE,idx_dir)%init(    "oxide_e"//trim(this%idx_dir_name(idx_dir)), this%g, IDX_EDGE, idx_dir)
      this%oxide(IDX_EDGE,idx_dir)%flags = this%poisson(IDX_EDGE,idx_dir)%flags
      call this%transport(IDX_EDGE,idx_dir)%init("transport_e"//trim(this%idx_dir_name(idx_dir)), this%g, IDX_EDGE, idx_dir)
    end do

    ! set transport cells
    if (this%gal) then
      gblock => this%gal_fl%get_block("desc(volume)")
      if (size(gblock%idata) /= (this%g1D(1)%n-1) * (this%g1D(2)%n-1)) then
        call program_error("Number of quads in GALENE file does not match grid")
      end if
      k = 0
      do j = 1, this%g1D(2)%n-1; do i = 1, this%g1D(1)%n - 1
        k = k + 1
        if (gblock%idata(k) /= this%gal_si) cycle
        idx = [i, j]
        call this%oxide(    IDX_CELL,0)%flags%set(idx, .false.)
        call this%transport(IDX_CELL,0)%flags%set(idx, .true.)
      end do; end do
    else
      i0 = 1
      i1 = 1
      do ri = 1, size(this%reg_trans)
        select case (this%gtype%s)
        case ("x", "xy", "xyz")
          ! get bounds
          do idx_dir = 1, idx_dim
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(ri)%xyz(idx_dir,1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(ri)%xyz(idx_dir,2)) - 1
          end do

          ! loop over cells in line-segment (x), rectangle (xy) or cuboid (xyz)
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            idx = ijk(1:idx_dim)

            ! enable transport for cell
            call this%oxide(    IDX_CELL,0)%flags%set(idx, .false.)
            call this%transport(IDX_CELL,0)%flags%set(idx, .true.)
          end do; end do; end do

        case ("tr_xy", "tr_xyz")
          ! get z bounds
          if (dim == 3) then
            i0(3) = bin_search(this%g1D(3)%x, this%reg_trans(ri)%xyz(3,1))
            i1(3) = bin_search(this%g1D(3)%x, this%reg_trans(ri)%xyz(3,2)) - 1
          end if

          ! enable transport in region
          gtype_tmp = new_string("tr_xy")
          do i = 1, this%gtr%ncell
            ! check if triangle is in region using centroid
            call this%gtr%get_cell([i], p)
            mid(1) = (p(1,1) + p(1,2) + p(1,3)) / 3
            mid(2) = (p(2,1) + p(2,2) + p(2,3)) / 3
            if (.not. this%reg_trans(ri)%point_test(gtype_tmp, mid)) cycle

            ! loop over z direction
            do k = i0(3), i1(3)
              ijk(1:2) = [i, k]
              idx      = ijk(1:idx_dim)

              ! enable transport for cell
              call this%oxide(    IDX_CELL,0)%flags%set(idx, .false.)
              call this%transport(IDX_CELL,0)%flags%set(idx, .true.)
            end do
          end do
        end select
      end do
    end if
    call this%oxide(    IDX_CELL,0)%init_final()
    call this%transport(IDX_CELL,0)%init_final()

    ! process vertices, edges
    do i = 1, this%transport(IDX_CELL,0)%n
      idx = this%transport(IDX_CELL,0)%get_idx(i)

      ! get adjoint cell parts
      call this%g%get_adjoint(idx, surf = surf, vol = vol)
      where (surf == 0) surf = SURF_TOL

      ! enable transport for vertices and update adjoint volumes
      do j = 1, this%g%cell_nvert
        call this%g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx, j, idx2, status)
        call this%oxide(    IDX_VERTEX,0)%flags%set(idx2, .false.)
        call this%transport(IDX_VERTEX,0)%flags%set(idx2, .true.)
        call this%tr_vol%update(idx2, vol(j))
      end do

      ! enable transport for edges and update adjoint surfaces
      do idx_dir = 1, idx_dim
        do j = 1, this%g%cell_nedge(idx_dir)
          call this%g%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          call this%oxide(    IDX_EDGE,idx_dir)%flags%set(idx2, .false.)
          call this%transport(IDX_EDGE,idx_dir)%flags%set(idx2, .true.)
          call this%tr_surf(idx_dir)%update(idx2, surf(j,idx_dir))
        end do
      end do
    end do

    ! finish vertex + edge tables
    call this%oxide(    IDX_VERTEX,0)%init_final()
    call this%transport(IDX_VERTEX,0)%init_final()
    do idx_dir = 1, idx_dim
      call this%oxide(    IDX_EDGE,idx_dir)%init_final()
      call this%transport(IDX_EDGE,idx_dir)%init_final()
    end do

    print "(A,I0)", "#(Transport vertices): ", this%transport(IDX_VERTEX,0)%n
  end subroutine

  subroutine device_params_init_doping(this, file)
    class(device_params), intent(inout) :: this
    type(input_file),     intent(in)    :: file

    integer                  :: ci, dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable     :: idx(:), idx2(:)
    logical                  :: status, load
    real                     :: cdop, dop(2), p(2,3), mid(2), ii_E_dop, ni
    real,    allocatable     :: surf(:,:), vol(:), tmp(:), ddop(:), tdop(:)
    type(string)             :: gtype_tmp
    type(gal_block), pointer :: gblock

    print "(A)", "init_doping"

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! allocate temp arrays
    allocate (idx(idx_dim), idx2(idx_dim))
    allocate (surf(this%g%max_cell_nedge,idx_dim), vol(this%g%cell_nvert))

    ! allocate/initialize grid data
    call allocate_grid_data3_real(this%dop,  this%g%idx_dim, [1, 0, DOP_DCON], [4, this%g%idx_dim, DOP_ACON])
    call allocate_grid_data1_real(this%ii_E_dop, idx_dim, DOP_DCON, DOP_ACON)
    do ci = DOP_DCON, DOP_ACON
      call this%dop(IDX_VERTEX,0,ci)%init(this%g, IDX_VERTEX, 0)
      call this%dop(IDX_CELL,  0,ci)%init(this%g, IDX_CELL,   0)
      do idx_dir = 1, idx_dim
        call this%dop(IDX_EDGE,idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)
      end do
    end do

    load = this%gal
    if (load) then ! load from GALENE III file
      gblock => this%gal_fl%get_block("norm fac")
      ni = norm(gblock%rdata(1), "1/cm^3")

      gblock => this%gal_fl%get_block("blank.dd")
      ddop = gblock%rdata
      gblock => this%gal_fl%get_block("blank.td")
      tdop = gblock%rdata

      call this%dop(IDX_VERTEX,0,DOP_DCON)%set(max(0.0, ni * 0.5 * (tdop + ddop)))
      call this%dop(IDX_VERTEX,0,DOP_ACON)%set(max(0.0, ni * 0.5 * (tdop - ddop)))
    else
      call file%get_section("load doping", ri, status = i)
      if (i > 0) call program_error("multiple 'load doping' sections found")
      if (i == 0) then ! load from input file
        if (size(this%reg_dop) > 0) call program_error("found 'load doping' AND 'doping' sections")
        load = .true.

        ! load doping on vertices
        call file%get(ri, "dcon", tmp, status = status)
        if (status) call this%dop(IDX_VERTEX,0,DOP_DCON)%set(tmp)
        call file%get(ri, "acon", tmp, status = status)
        if (status) call this%dop(IDX_VERTEX,0,DOP_ACON)%set(tmp)
      end if
    end if

    if (load) then
      ! set doping in cells
      do i = 1, this%transport(IDX_CELL,0)%n
        idx = this%transport(IDX_CELL,0)%get_idx(i)

        ! get adjoint cell volume parts, normalize to sum (= cell volume)
        call this%g%get_adjoint(idx, vol = vol)
        vol = vol / sum(vol)

        ! get doping for cell (mean)
        do j = 1, this%g%cell_nvert
          call this%g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx, j, idx2, status)
          do ci = DOP_DCON, DOP_ACON
            call this%dop(IDX_CELL,0,ci)%update(idx, vol(j) * this%dop(IDX_VERTEX,0,ci)%get(idx2))
          end do
        end do
      end do
    else ! regions
      i0   = 1
      i1   = 1
      do ri = 1, size(this%reg_dop)
        dop = this%reg_dop(ri)%dop

        ! set doping in cells
        select case (this%gtype%s)
        case ("x", "xy", "xyz")
          ! get bounds
          do idx_dir = 1, idx_dim
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(ri)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(ri)%xyz(idx_dir, 2)) - 1
          end do

          ! loop over cells in line-segment (x), rectangle (xy) or cuboid (xyz)
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            idx = ijk(1:idx_dim)

            ! set doping in cells
            do ci = DOP_DCON, DOP_ACON
              if (this%reg_dop(ri)%set_dop(ci)) call this%dop(IDX_CELL,0,ci)%set(idx, dop(ci))
            end do
          end do; end do; end do

        case ("tr_xy", "tr_xyz")
          ! get z bounds
          if (dim == 3) then
            i0(3) = bin_search(this%g1D(3)%x, this%reg_dop(ri)%xyz(3,1))
            i1(3) = bin_search(this%g1D(3)%x, this%reg_dop(ri)%xyz(3,2)) - 1
          end if

          ! loop over triangle cells
          gtype_tmp = new_string("tr_xy")
          do i = 1, this%gtr%ncell
            ! check if triangle is in region using centroid
            call this%gtr%get_cell([i], p)
            mid(1) = (p(1,1) + p(1,2) + p(1,3)) / 3
            mid(2) = (p(2,1) + p(2,2) + p(2,3)) / 3
            if (.not. this%reg_dop(ri)%point_test(gtype_tmp, mid)) cycle

            ! loop over z direction
            do k = i0(3), i1(3)
              ijk(1:2) = [i, k]
              idx      = ijk(1:idx_dim)

              ! set doping in cells
              do ci = DOP_DCON, DOP_ACON
                if (this%reg_dop(ri)%set_dop(ci)) call this%dop(IDX_CELL,0,ci)%set(idx, dop(ci))
              end do
            end do
          end do

        end select
      end do

      ! set doping on vertices
      do i = 1, this%transport(IDX_CELL,0)%n
        idx = this%transport(IDX_CELL,0)%get_idx(i)

        ! get adjoint cell volume parts
        call this%g%get_adjoint(idx, vol = vol)

        ! update doping for all vertices belonging to cell
        do j = 1, this%g%cell_nvert
          call this%g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx, j, idx2, status)
          do ci = DOP_DCON, DOP_ACON
            dop(ci) = this%dop(IDX_CELL,0,ci)%get(idx)
            call this%dop(IDX_VERTEX,0,ci)%update(idx2, dop(ci) * vol(j) / this%tr_vol%get(idx2))
          end do
        end do
      end do
    end if

    ! set doping on edges
    do i = 1, this%transport(IDX_CELL,0)%n
      idx = this%transport(IDX_CELL,0)%get_idx(i)

      ! get adjoint cell surfaces
      call this%g%get_adjoint(idx, surf = surf)

      ! loop over edges in this cell, average doping over adjoint surface
      do idx_dir = 1, idx_dim
        do j = 1, this%g%cell_nedge(idx_dir)
          call this%g%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          do ci = DOP_DCON, DOP_ACON
            dop(ci) = this%dop(IDX_CELL,0,ci)%get(idx)
            call this%dop(IDX_EDGE,idx_dir,ci)%update(idx2, dop(ci) * surf(j,idx_dir) / this%tr_surf(idx_dir)%get(idx2))
          end do
        end do
      end do
    end do

    ! Pearson-Bardeen ionization model init
    do ci = DOP_DCON, DOP_ACON
      call this%ii_E_dop(ci)%init(this%g, IDX_VERTEX, 0)
      cdop = this%smc%ii_E_dop0(ci) / (this%smc%ii_N_crit(ci)**(1.0/3.0))
      do i = 1, this%transport(IDX_VERTEX,0)%n
        idx     = this%transport(IDX_VERTEX,0)%get_idx(i)
        dop(ci) = this%dop(IDX_VERTEX,0,ci)%get(idx)

        ii_E_dop = this%smc%ii_E_dop0(ci) - cdop * dop(ci)**(1.0/3.0) !max(this%smc%ii_E_dop0(ci) - cdop * dop(ci)**(1.0/3.0), 0.0)
        call this%ii_E_dop(ci)%set(idx, ii_E_dop)
      end do
    end do
  end subroutine

  subroutine device_params_init_mobility(this)
    class(device_params), intent(inout) :: this

    integer              :: ci, dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable :: idx(:), idx2(:)
    logical              :: status
    real                 :: alpha, dop(2), mid(2), mob0(2), mob_min, mob_max, N_ref, p(2,3)
    real,    allocatable :: surf(:,:)
    type(string)         :: gtype_tmp

    print "(A)", "init_mobility"

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! allocate temp arrays
    allocate (idx(idx_dim), idx2(idx_dim))
    allocate (surf(this%g%max_cell_nedge,idx_dim))

    call allocate_grid_data3_real(this%mob0, this%g%idx_dim, [1, 0, this%ci0], [4, this%g%idx_dim, this%ci1])
    do ci = this%ci0, this%ci1
      call this%mob0(IDX_CELL,0,ci)%init(this%g, IDX_CELL, 0)
      do idx_dir = 1, idx_dim
        call this%mob0(IDX_EDGE,idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)
      end do
    end do

    ! init zero-field mobility based on doping
    do ci = this%ci0, this%ci1
      mob_min = this%smc%mob_min(ci)
      mob_max = this%smc%mob_max(ci)
      N_ref   = this%smc%N_ref(ci)
      alpha   = this%smc%alpha(ci)

      ! mobility in cells
      do i = 1, this%transport(IDX_CELL,0)%n
        idx           = this%transport(IDX_CELL,0)%get_idx(i)
        dop(DOP_DCON) = this%dop(IDX_CELL,0,DOP_DCON)%get(idx)
        dop(DOP_ACON) = this%dop(IDX_CELL,0,DOP_ACON)%get(idx)
        mob0(ci)      = mob_min + (mob_max - mob_min)/(1 + ((dop(DOP_DCON) + dop(DOP_ACON))/N_ref)**alpha)
        call this%mob0(IDX_CELL,0,ci)%set(idx, mob0(ci))
      end do
    end do

    ! additionally process mobility regions
    i0   = 1
    i1   = 1
    do ri = 1, size(this%reg_mob)
      mob0 = this%reg_mob(ri)%mob0

      ! set mobility in cells
      select case (this%gtype%s)
      case ("x", "xy", "xyz")
        ! get bounds
        do idx_dir = 1, idx_dim
          i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_mob(ri)%xyz(idx_dir, 1))
          i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_mob(ri)%xyz(idx_dir, 2)) - 1
        end do

        ! loop over cells in line-segment (x), rectangle (xy) or cuboid (xyz)
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          idx = ijk(1:idx_dim)

          ! overwrite mobility in cells
          do ci = this%ci0, this%ci1
            if (this%reg_mob(ri)%set_mob(ci)) call this%mob0(IDX_CELL,0,ci)%set(idx, mob0(ci))
          end do
        end do; end do; end do

      case ("tr_xy", "tr_xyz")
        ! get z bounds
        if (dim == 3) then
          i0(3) = bin_search(this%g1D(3)%x, this%reg_mob(ri)%xyz(3,1))
          i1(3) = bin_search(this%g1D(3)%x, this%reg_mob(ri)%xyz(3,2)) - 1
        end if

        ! loop over triangle cells
        gtype_tmp = new_string("tr_xy")
        do i = 1, this%gtr%ncell
          ! check if triangle is in region using centroid
          call this%gtr%get_cell([i], p)
          mid(1) = (p(1,1) + p(1,2) + p(1,3)) / 3
          mid(2) = (p(2,1) + p(2,2) + p(2,3)) / 3
          if (.not. this%reg_mob(ri)%point_test(gtype_tmp, mid)) cycle

          ! loop over z direction
          do k = i0(3), i1(3)
            ijk(1:2) = [i, k]
            idx      = ijk(1:idx_dim)

            ! overwrite mobility in cells
            do ci = this%ci0, this%ci1
              if (this%reg_mob(ri)%set_mob(ci)) call this%mob0(IDX_CELL,0,ci)%set(idx, mob0(ci))
            end do
          end do
        end do

      end select
    end do

    ! set mobility on edges
    do i = 1, this%transport(IDX_CELL,0)%n
      idx = this%transport(IDX_CELL,0)%get_idx(i)

      ! get adjoint cell surfaces
      call this%g%get_adjoint(idx, surf = surf)

      ! loop over edges in this cell, average mobility over adjoint surface
      do idx_dir = 1, idx_dim
        do j = 1, this%g%cell_nedge(idx_dir)
          call this%g%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          do ci = this%ci0, this%ci1
            mob0(ci) = this%mob0(IDX_CELL,0,ci)%get(idx)
            call this%mob0(IDX_EDGE,idx_dir,ci)%update(idx2, mob0(ci) * surf(j,idx_dir) / this%tr_surf(idx_dir)%get(idx2))
          end do
        end do
      end do
    end do
  end subroutine

  subroutine device_params_init_contacts(this)
    class(device_params), intent(inout) :: this

    integer                           :: ci, dim, i, i0(3), i1(3), ict, ict0, idx_dim, idx_dir, ijk(3), j, k, k0, k1, kk, nct, ri, nv_poiss, nv_conti
    integer, allocatable              :: idx(:), gimat(:)
    real                              :: dop(2), ii_E_dop(2), p(2)
    type(string)                      :: name
    type(string), allocatable         :: gmat(:)
    type(mapnode_string_int), pointer :: node
    type(gal_block), pointer          :: gblock1, gblock2, gblock3

    print "(A)", "init_contacts"

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! get all contact names and number of contacts (nct)
    call this%contact_map%init()
    nct = 0
    if (this%gal) then
      allocate (gmat(this%gal_mat%n), gimat(this%gal_mat%n))
      call this%gal_mat%to_array(keys = gmat, values = gimat)
      do i = 1, size(gmat)
        if ((gmat(i)%s == "SIL") .or. (gmat(i)%s == "OX")) cycle

        ! new contact
        nct = nct + 1
        call this%contact_map%insert(gmat(i), nct)
      end do
    else
      do ri = 1, size(this%reg_ct)
        ! search for contact name in map
        node => this%contact_map%find(this%reg_ct(ri)%name)

        ! do nothing if name already exists
        if (associated(node)) cycle

        ! new contact
        nct = nct + 1
        call this%contact_map%insert(this%reg_ct(ri)%name, nct)
      end do
    end if

    this%nct = nct
    allocate (this%contacts(nct))

    ! allocate/initialize grid data
    allocate (idx(idx_dim))
    call allocate_grid_data0_int(this%ict, idx_dim)
    call this%ict%init(this%g, IDX_VERTEX, 0)

    ! allocate/initialize grid tables
    allocate (this%contacted(nct), this%poisson_vct(0:nct), this%oxide_vct(0:nct), this%transport_vct(0:nct))
    call this%poisson_vct(  0)%init("poisson_VCT0",   this%g, IDX_VERTEX, 0)
    call this%oxide_vct(    0)%init("oxide_VCT0",     this%g, IDX_VERTEX, 0)
    call this%transport_vct(0)%init("transport_VCT0", this%g, IDX_VERTEX, 0)
    this%poisson_vct(  0)%flags = this%poisson(  IDX_VERTEX,0)%flags
    this%oxide_vct(    0)%flags = this%oxide(    IDX_VERTEX,0)%flags
    this%transport_vct(0)%flags = this%transport(IDX_VERTEX,0)%flags

    ! init contacts
    if (this%gal) then
      gblock1 => this%gal_fl%get_block("#vertex(contact,poisson)")
      gblock2 => this%gal_fl%get_block("#vertex(contact,continuity)")
      gblock3 => this%gal_fl%get_block("vertex(contact,poisson)")
      k1 = 0
      do i = 1, size(gmat)
        do j = 1, size(gimat)
          if (gimat(j) == i) exit
        end do
        if ((gmat(j)%s == "SIL") .or. (gmat(j)%s == "OX")) cycle
        ict = this%contact_map%get(gmat(j))

        !i: gal contact order
        !j: gimat contact order
        !ict: this contact order

        ! number of vertices
        nv_poiss = gblock1%idata(i)
        nv_conti = gblock2%idata(i)

        this%contacts(ict)%name = gmat(j)%s
        if (nv_conti == 0) then
          this%contacts(ict)%type = CT_GATE
        else
          this%contacts(ict)%type = CT_OHMIC
        end if
        this%contacts(ict)%phims = this%gal_phims

        call this%contacted(    ict)%init("contacted_"//name%s, this%g, IDX_VERTEX, 0)
        call this%poisson_vct(  ict)%init("poisson_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        call this%oxide_vct(    ict)%init("oxide_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        call this%transport_vct(ict)%init("transport_VCT_"//name%s, this%g, IDX_VERTEX, 0)

        k0 = k1 + 1
        k1 = k1 + nv_poiss
        do k = k0, k1
          kk = gblock3%idata(k)
          idx(1) = mod(kk-1, this%g1D(1)%n) + 1
          idx(2) = (kk - idx(1)) / this%g1D(1)%n + 1
          call this%ict%set(idx, ict)
          call this%contacted(ict)%flags%set(idx, .true.)
        end do
      end do
    else
      i0 = 1
      i1 = 1
      do ri = 1, size(this%reg_ct)
        ! get contact name, index and type
        name = this%reg_ct(ri)%name
        ict  = this%contact_map%get(name)

        ! create new contact if name is encountered for the first time
        if (.not. allocated(this%contacts(ict)%name)) then
          this%contacts(ict)%name  = name%s
          this%contacts(ict)%type  = this%reg_ct(ri)%type
          this%contacts(ict)%phims = this%reg_ct(ri)%phims

          ! Transfer Schottky-specific parameters if applicable
          if (this%reg_ct(ri)%type == CT_SCHOTTKY) then
            this%contacts(ict)%phi_b = this%reg_ct(ri)%phi_b
            this%contacts(ict)%A_richardson = this%reg_ct(ri)%A_richardson
            this%contacts(ict)%ifbl = this%reg_ct(ri)%ifbl
            print "(A,L1)", "DEBUG: Contact " // name%s // " ifbl = ", this%contacts(ict)%ifbl
            ! Transfer tunneling parameters
            this%contacts(ict)%tunneling = this%reg_ct(ri)%tunneling
            this%contacts(ict)%m_tunnel = this%reg_ct(ri)%m_tunnel
            print "(A,L1)", "DEBUG: Contact " // name%s // " tunneling = ", this%contacts(ict)%tunneling
            print "(A,F6.3)", "DEBUG: Contact " // name%s // " m_tunnel = ", this%contacts(ict)%m_tunnel
          end if

          call this%contacted(    ict)%init("contacted_"//name%s, this%g, IDX_VERTEX, 0)
          call this%poisson_vct(  ict)%init("poisson_VCT_"//name%s, this%g, IDX_VERTEX, 0)
          call this%oxide_vct(    ict)%init("oxide_VCT_"//name%s, this%g, IDX_VERTEX, 0)
          call this%transport_vct(ict)%init("transport_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        elseif (this%contacts(ict)%type /= this%reg_ct(ri)%type) then
          call program_error("contact "//name%s//" given multiple times with different types")
        end if

        select case (this%gtype%s)
        case("x", "xy", "xyz")
          ! get bounds
          do idx_dir = 1, idx_dim
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(ri)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(ri)%xyz(idx_dir, 2))
          end do

          ! update vertex tables
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk  = [i, j, k]
            idx  = ijk(1:idx_dim)
            ict0 = this%ict%get(idx)
            if ((ict0 /= 0) .and. (ict0 /= ict)) then
              call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
            end if
            call this%ict%set(idx, ict)
            call this%contacted(ict)%flags%set(idx, .true.)
          end do; end do; end do

        case ("tr_xy", "tr_xyz")
          ! get z bounds
          if (dim == 3) then
            i0(3) = bin_search(this%g1D(3)%x, this%reg_ct(ri)%xyz(3,1))
            i1(3) = bin_search(this%g1D(3)%x, this%reg_ct(ri)%xyz(3,2))
          end if

          ! loop over triangle grid vertices
          do i = 1, this%gtr%nvert
            ! test if vertex belongs to contact surface
            call this%gtr%get_vertex([i], p)
            if (.not. this%reg_ct(ri)%point_test(new_string("tr_xy"), p)) cycle

            ! loop over z direction
            do k = i0(3), i1(3)
              ijk(1:2) = [i, k]
              idx      = ijk(1:idx_dim)
              ict0     = this%ict%get(idx)
              if ((ict0 /= 0) .and. (ict0 /= ict)) then
                call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
              end if
              call this%ict%set(idx, ict)
              call this%contacted(ict)%flags%set(idx, .true.)
            end do
          end do
        end select
      end do
    end if

    ! finish initialization of contacts
    do ict = 1, nct
      call this%contacted(ict)%init_final()

      ! set phims for ohmic contacts (not for Schottky)
      if (this%contacts(ict)%type == CT_OHMIC) then
        ! find vertex indices in transport region
        do i = 1, this%contacted(ict)%n
          idx = this%contacted(ict)%get_idx(i)
          if (this%transport(IDX_VERTEX,0)%flags%get(idx)) exit
        end do
        if (i > this%contacted(ict)%n) call program_error("Ohmic contact "//this%contacts(ict)%name//" not at transport region")

        ! calculate phims using charge neutrality
        dop(DOP_DCON) = this%dop(IDX_VERTEX,0,DOP_DCON)%get(idx)
        dop(DOP_ACON) = this%dop(IDX_VERTEX,0,DOP_ACON)%get(idx)
        ii_E_dop(DOP_DCON) = this%ii_E_dop(DOP_DCON)%get(idx)
        ii_E_dop(DOP_ACON) = this%ii_E_dop(DOP_ACON)%get(idx)
        call this%contacts(ict)%set_phims_ohmic(this%ci0, this%ci1, dop, ii_E_dop, this%smc)
      elseif (this%contacts(ict)%type == CT_SCHOTTKY) then
        ! set phims for Schottky contacts using Sentaurus Device approach
        call this%contacts(ict)%set_phims_schottky(this%smc)
      end if

      ! update contact vertex tables
      do i = 1, this%contacted(ict)%n
        idx = this%contacted(ict)%get_idx(i)
        if (this%poisson(IDX_VERTEX,0)%flags%get(idx)) then
          call this%poisson_vct(ict)%flags%set(idx, .true. )
          call this%poisson_vct(  0)%flags%set(idx, .false.)
        end if
        if (this%oxide(IDX_VERTEX,0)%flags%get(idx)) then
          call this%oxide_vct(ict)%flags%set(idx, .true. )
          call this%oxide_vct(  0)%flags%set(idx, .false.)
        end if
        if (this%transport(IDX_VERTEX,0)%flags%get(idx)) then
          call this%transport_vct(ict)%flags%set(idx, .true. )
          call this%transport_vct(  0)%flags%set(idx, .false.)
        end if
      end do
      call this%poisson_vct(  ict)%init_final()
      call this%oxide_vct(    ict)%init_final()
      call this%transport_vct(ict)%init_final()
    end do
    call this%poisson_vct(  0)%init_final()
    call this%oxide_vct(    0)%init_final()
    call this%transport_vct(0)%init_final()

    ! doping vertices
    do ci = DOP_DCON, DOP_ACON
      call this%dopvert(ci)%init("dop"//DOP_NAME(ci)//"_v", this%g, IDX_VERTEX, 0)
      call this%ionvert(ci)%init("ion"//DOP_NAME(ci)//"_v", this%g, IDX_VERTEX, 0)
    end do
    do i = 1, this%transport_vct(0)%n
      idx = this%transport_vct(0)%get_idx(i)

      do ci = DOP_DCON, DOP_ACON
        dop(ci) = this%dop(IDX_VERTEX,0,ci)%get(idx)
        if (dop(ci) > 0) then
          call this%dopvert(ci)%flags%set(idx, .true.)
          if (dop(ci) < this%smc%ii_dop_th(ci)) then
            call this%ionvert(ci)%flags%set(idx, .true.)
          end if
        end if
      end do
    end do
    do ci = DOP_DCON, DOP_ACON
      call this%dopvert(ci)%init_final()
      call this%ionvert(ci)%init_final()
    end do
  end subroutine


  function get_ct_surf(this, ict, idx_vert) result(ct_surf)
    !! Calculate boundary surface area for a contact vertex
    !! This is the area of the metal-semiconductor interface at this vertex

    class(device_params), intent(in) :: this
    integer,              intent(in) :: ict
    integer,              intent(in) :: idx_vert(:)
    real                             :: ct_surf

    integer :: l, m, idx_dir, n_verts
    integer :: idx_neighb(size(idx_vert)), idx_surf(size(idx_vert)), idx_surf_vert(size(idx_vert))
    logical :: status, surf_contacted
    real    :: p1(3), p2(3), p3(3), p4(3), edge_len
    real    :: surf_verts(3,4), vec1(3), vec2(3), surf_area

    ! Initialize contact surface
    ct_surf = 0.0

    select case (this%g%dim)
    case (1)
      ! For 1D, use normalized cross-sectional area
      ct_surf = 1.0

    case (2)
      ! For 2D: sum half of edge lengths to neighboring contact vertices
      do l = 1, this%g%get_max_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0)
        call this%g%get_neighb(IDX_VERTEX, 0, IDX_VERTEX, 0, idx_vert, l, idx_neighb, status)
        if (.not. status) exit

        ! Check if neighbor is in the same contact
        if (this%ict%get(idx_neighb) == ict) then
          ! Get vertex positions
          call this%g%get_vertex(idx_vert, p1(1:2))
          call this%g%get_vertex(idx_neighb, p2(1:2))

          ! Calculate edge length and add half to this vertex's surface
          edge_len = sqrt((p2(1) - p1(1))**2 + (p2(2) - p1(2))**2)
          ct_surf = ct_surf + 0.5 * edge_len
        end if
      end do

    case (3)
      ! For 3D: sum portions of adjacent surface areas that are fully contacted
      do idx_dir = 1, 3
        do l = 1, this%g%get_max_neighb(IDX_VERTEX, 0, IDX_FACE, idx_dir)
          call this%g%get_neighb(IDX_VERTEX, 0, IDX_FACE, idx_dir, idx_vert, l, idx_surf, status)
          if (.not. status) exit

          ! Check if all vertices of the surface belong to the same contact
          surf_contacted = .true.
          n_verts = 0

          do m = 1, this%g%get_max_neighb(IDX_FACE, idx_dir, IDX_VERTEX, 0)
            call this%g%get_neighb(IDX_FACE, idx_dir, IDX_VERTEX, 0, idx_surf, m, idx_surf_vert, status)
            if (.not. status) exit

            n_verts = n_verts + 1
            if (n_verts > 4) exit  ! Should not happen for tensor grids

            ! Check if this surface vertex is in the same contact
            if (this%ict%get(idx_surf_vert) /= ict) then
              surf_contacted = .false.
              exit
            end if

            ! Store vertex coordinates for area calculation
            call this%g%get_vertex(idx_surf_vert, surf_verts(:, n_verts))
          end do

          ! If surface is fully contacted, add its area contribution
          if (surf_contacted .and. n_verts == 4) then
            ! Calculate area of rectangular surface using cross product
            ! Vectors from first vertex to second and third
            vec1 = surf_verts(:, 2) - surf_verts(:, 1)
            vec2 = surf_verts(:, 3) - surf_verts(:, 1)

            ! Area = |vec1 x vec2| (parallelogram area)
            surf_area = sqrt((vec1(2)*vec2(3) - vec1(3)*vec2(2))**2 + &
                           (vec1(3)*vec2(1) - vec1(1)*vec2(3))**2 + &
                           (vec1(1)*vec2(2) - vec1(2)*vec2(1))**2)

            ! Each vertex gets 1/4 of the surface area
            ct_surf = ct_surf + 0.25 * surf_area
          end if
        end do
      end do

    case default
      call program_error("Invalid dimension for contact surface calculation")
    end select
  end function

end module
