m4_include(util/macro.f90.inc)

module device_params_m

  use bin_search_m,     only: bin_search
  use contact_m,        only: CT_OHMIC, CT_GATE, contact
  use error_m,          only: assert_failed, program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid, grid_ptr
  use grid_data_m,      only: allocate_grid_data0_real, allocate_grid_data1_real, allocate_grid_data2_real, &
    &                         allocate_grid_data3_real, allocate_grid_data0_int, allocate_grid_data1_int, &
    &                         grid_data_int, grid_data_real
  use grid_generator_m, only: DIR_NAME, generate_1D_grid, generate_triangle_grid
  use grid_table_m,     only: grid_table
  use grid1D_m,         only: grid1D
  use input_m,          only: input_file
  use map_m,            only: map_string_int, mapnode_string_int
  use region_m,         only: region_ptr, region_poisson, region_transport, region_doping, region_contact
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON, semiconductor
  use string_m,         only: string, new_string
  use tensor_grid_m,    only: tensor_grid
  use triang_grid_m,    only: triang_grid
  use triangle_m,       only: triangulation

  implicit none

  private
  public device_params

  type device_params
    !! device geometry and material parameters

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
      !! donator/acceptor concentration (idx_type, idx_dir, dcon/acon)
    class(grid_data_real), allocatable :: asb(:)
      !! Altermatt-Schenk factor of bound states (dcon/acon)
    class(grid_data_real), allocatable :: edop(:)
      !! donator/acceptor dopant energy level (dcon/acon)
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
    type(triangulation)       :: tr
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
    type(region_contact),   allocatable :: reg_ct(:)
      !! contact regions

  contains
    procedure :: init     => device_params_init
    procedure :: destruct => device_params_destruct

    procedure, private :: init_transport_params  => device_params_init_transport_params
    procedure, private :: init_regions           => device_params_init_regions
    procedure, private :: init_grid              => device_params_init_grid
    procedure, private :: init_poisson           => device_params_init_poisson
    procedure, private :: init_transport         => device_params_init_transport
    procedure, private :: init_doping            => device_params_init_doping
    procedure, private :: init_contacts          => device_params_init_contacts
  end type

  real, parameter :: SURF_TOL = 1e-32
    !! use if adjoint surface is zero

contains

  subroutine device_params_init(this, file)
    !! initialize device parameter object
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file
      !! device file

    ! get grid type ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    call file%get("grid", "gtype", this%gtype)
    select case (this%gtype%s)
    case ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    case default
      call program_error("Invalid grid type: "//this%gtype%s)
    end select

    ! process file sections
    call this%init_transport_params(file)
    call this%init_regions(file)
    call this%init_grid(file)

    ! process regions
    call this%init_poisson()
    call this%init_transport()
    call this%init_doping(file)
    call this%init_contacts()
  end subroutine

  subroutine device_params_destruct(this)
    !! destruct device parameter object
    class(device_params), intent(inout) :: this

    call this%contact_map%destruct()
  end subroutine

  subroutine device_params_init_transport_params(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)    :: file

    integer :: sid
    logical :: elec, hole

    ! find transport parameters section id
    call file%get_section("transport parameters", sid)

    ! load parameters
    call file%get(sid, "electrons", elec)
    call file%get(sid, "holes",     hole)
    this%ci0 = CR_ELEC
    this%ci1 = CR_HOLE
    if (.not. elec) this%ci0 = CR_HOLE
    if (.not. hole) this%ci1 = CR_ELEC
    call file%get(sid, "N_c",        this%smc%edos(CR_ELEC))
    call file%get(sid, "N_v",        this%smc%edos(CR_HOLE))
    call file%get(sid, "E_gap",      this%smc%band_gap)
    this%smc%n_intrin = sqrt(this%smc%edos(CR_ELEC) * this%smc%edos(CR_HOLE) * exp(-this%smc%band_gap))
    this%smc%band_edge(CR_ELEC) =   0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))
    this%smc%band_edge(CR_HOLE) = - 0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))
    call file%get(sid, "mass",       this%smc%mass)
    call file%get(sid, "degen",      this%smc%degen)
    call file%get(sid, "alpha",      this%smc%alpha)
    call file%get(sid, "beta",       this%smc%beta)
    call file%get(sid, "mob_min",    this%smc%mob_min)
    call file%get(sid, "mob_max",    this%smc%mob_max)
    call file%get(sid, "N_ref",      this%smc%N_ref)
    call file%get(sid, "v_sat",      this%smc%v_sat)
    call file%get(sid, "curr_fact",  this%curr_fact)
    call file%get(sid, "incomp_ion", this%smc%incomp_ion)
    call file%get(sid, "edop0",      this%smc%edop)
    call file%get(sid, "asc",        this%smc%asc)
    call file%get(sid, "asd",        this%smc%asd)
    call file%get(sid, "N_asb",      this%smc%N_asb)
    call file%get(sid, "N_asr",      this%smc%N_asr)
    call file%get(sid, "dop_degen",  this%smc%g_dop)

    ! make sure parameters are valid
    m4_assert(this%ci0 <= this%ci1)
    m4_assert(size(this%smc%mass) == 2)
    m4_assert(size(this%smc%alpha) == 2)
    m4_assert(size(this%smc%beta) == 2)
    m4_assert(size(this%smc%mob_min) == 2)
    m4_assert(size(this%smc%mob_max) == 2)
    m4_assert(size(this%smc%N_ref) == 2)
    m4_assert(size(this%smc%v_sat) == 2)
    m4_assert(size(this%smc%edop) == 2)
    m4_assert(size(this%smc%asc) == 2)
    m4_assert(size(this%smc%asd) == 2)
    m4_assert(size(this%smc%N_asb) == 2)
    m4_assert(size(this%smc%N_asr) == 2)
    m4_assert(size(this%smc%g_dop) == 2)
  end subroutine

  subroutine device_params_init_regions(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)    :: file

    integer, allocatable :: sids(:)
    integer              :: si, i, j

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
    call file%get("grid", "load", load, status)
    if (.not. status) load = .false.
    ngptr = 0

    ! generate or load grid
    select case(this%gtype%s)
    case("x", "xy", "xyz")
      do dir = 1, dim
        call generate_1D_grid(file, load, dir, this%reg, this%g1D(dir), gptr, ngptr)
      end do
      allocate (character(1) :: this%idx_dir_name(dim))
      do dir = 1, dim
        this%idx_dir_name(dir) = DIR_NAME(dir)
      end do

    case("tr_xy", "tr_xyz")
      call generate_triangle_grid(file, load, this%reg, this%tr, this%gtr, gptr, ngptr)
      if (dim == 3) then
        call generate_1D_grid(file, load, 3, this%reg, this%g1D(3), gptr, ngptr)
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

    integer              :: dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable :: idx(:), idx2(:)
    logical              :: status
    real                 :: mid(2), p(2,3)
    real,    allocatable :: surf(:,:), vol(:)
    type(string)         :: gtype_tmp

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

    ! process regions, focus on cells
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
      do idx_dir = 1, idx_dim
        do j = 1, this%g%cell_nedge(idx_dir)
          call this%g%get_neighb(IDX_CELL, 0, IDX_EDGE, idx_dir, idx, j, idx2, status)
          call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx2, .true.)
          call this%surf(idx_dir)%update(idx2, surf(j,idx_dir))
          call this%eps(IDX_EDGE,idx_dir)%update(idx2, surf(j,idx_dir) * this%eps(IDX_CELL,0)%get(idx))
        end do
      end do
    end do

    ! weighted average of permittivity on adjoint surfaces
    do idx_dir = 1, idx_dim
      call this%eps(IDX_EDGE,idx_dir)%set(this%eps(IDX_EDGE,idx_dir)%get() / this%surf(idx_dir)%get())
    end do

    ! finish vertex + edge tables
    call this%poisson(IDX_VERTEX,0)%init_final()
    do idx_dir = 1, idx_dim
      call this%poisson(IDX_EDGE,idx_dir)%init_final()
    end do
  end subroutine

  subroutine device_params_init_transport(this)
    class(device_params), intent(inout) :: this

    integer              :: dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable :: idx(:), idx2(:)
    logical              :: status
    real                 :: mid(2), p(2,3)
    real,    allocatable :: surf(:,:), vol(:)
    type(string)         :: gtype_tmp

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

    ! process regions, focus on cells
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
  end subroutine

  subroutine device_params_init_doping(this, file)
    class(device_params), intent(inout) :: this
    type(input_file),     intent(in)    :: file

    integer              :: ci, dim, i, i0(3), i1(3), idx_dim, idx_dir, ijk(3), j, k, ri
    integer, allocatable :: idx(:), idx2(:)
    logical              :: status
    real                 :: dop(2), p(2,3), mid(2), asb, edop, mob0, mob_min, mob_max, N_ref, alpha
    real,    allocatable :: surf(:,:), vol(:), tmp(:)
    type(string)         :: gtype_tmp

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! allocate temp arrays
    allocate (idx(idx_dim), idx2(idx_dim))
    allocate (surf(this%g%max_cell_nedge,idx_dim), vol(this%g%cell_nvert))

    ! allocate/initialize grid data
    call allocate_grid_data3_real(this%dop,  this%g%idx_dim, [1, 0, DOP_DCON], [4, this%g%idx_dim, DOP_ACON])
    call allocate_grid_data3_real(this%mob0, this%g%idx_dim, [1, 0, this%ci0], [4, this%g%idx_dim, this%ci1])
    call allocate_grid_data1_real(this%asb,  idx_dim, DOP_DCON, DOP_ACON)
    call allocate_grid_data1_real(this%edop, idx_dim, DOP_DCON, DOP_ACON)
    do ci = DOP_DCON, DOP_ACON
      call this%dop(IDX_VERTEX,0,ci)%init(this%g, IDX_VERTEX, 0)
      call this%dop(IDX_CELL,  0,ci)%init(this%g, IDX_CELL,   0)
      do idx_dir = 1, idx_dim
        call this%dop(IDX_EDGE, idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)
      end do
    end do

    ! load doping?
    call file%get_section("load doping", ri, status = i)
    if (i > 0) call program_error("multiple 'load doping' sections found")
    if (i == 0) then ! load
      if (size(this%reg_dop) > 0) call program_error("found 'load doping' AND 'doping' sections")

      ! load doping on vertices
      call file%get(ri, "dcon", tmp, status = status)
      if (status) call this%dop(IDX_VERTEX,0,DOP_DCON)%set(tmp)
      call file%get(ri, "acon", tmp, status = status)
      if (status) call this%dop(IDX_VERTEX,0,DOP_ACON)%set(tmp)

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
              call this%dop(IDX_CELL,0,ci)%set(idx, dop(ci))
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
                call this%dop(IDX_CELL,0,ci)%set(idx, dop(ci))
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

    ! Altermatt-Schenk ionization model init
    do ci = DOP_DCON, DOP_ACON
      call this%edop(ci)%init(this%g, IDX_VERTEX, 0)
      call this%asb( ci)%init(this%g, IDX_VERTEX, 0)
      do i = 1, this%transport(IDX_VERTEX,0)%n
        idx     = this%transport(IDX_VERTEX,0)%get_idx(i)
        dop(ci) = this%dop(IDX_VERTEX,0,ci)%get(idx)

        edop = this%smc%edop(ci) / (1 + (dop(ci) / this%smc%N_asr(ci))**this%smc%asc(ci))
        asb  = 1 / (1 + (dop(ci) / this%smc%N_asb(ci))**this%smc%asd(ci))

        call this%asb(ci)%set( idx, asb)
        call this%edop(ci)%set(idx, edop)
      end do
    end do

    ! init zero-field mobility
    do ci = this%ci0, this%ci1
      mob_min = this%smc%mob_min(ci)
      mob_max = this%smc%mob_max(ci)
      N_ref   = this%smc%N_ref(ci)
      alpha   = this%smc%alpha(ci)

      ! mobility in cells
      call this%mob0(IDX_CELL,0,ci)%init(this%g, IDX_CELL, 0)
      do i = 1, this%transport(IDX_CELL,0)%n
        idx           = this%transport(IDX_CELL,0)%get_idx(i)
        dop(DOP_DCON) = this%dop(IDX_CELL,0,DOP_DCON)%get(idx)
        dop(DOP_ACON) = this%dop(IDX_CELL,0,DOP_ACON)%get(idx)
        mob0          = mob_min + (mob_max - mob_min)/(1 + ((dop(DOP_DCON) + dop(DOP_ACON))/N_ref)**alpha)
        call this%mob0(IDX_CELL,0,ci)%set(idx, mob0)
      end do

      ! mobility on edges
      do idx_dir = 1, this%g%idx_dim
        call this%mob0(IDX_EDGE,idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)
        do i = 1, this%transport(IDX_EDGE,idx_dir)%n
          idx           = this%transport(IDX_EDGE,idx_dir)%get_idx(i)
          dop(DOP_DCON) = this%dop(IDX_EDGE,idx_dir,DOP_DCON)%get(idx)
          dop(DOP_ACON) = this%dop(IDX_EDGE,idx_dir,DOP_ACON)%get(idx)
          mob0          = mob_min + (mob_max - mob_min)/(1 + ((dop(DOP_DCON) + dop(DOP_ACON))/N_ref)**alpha)
          call this%mob0(IDX_EDGE,idx_dir,ci)%set(idx, mob0)
        end do
      end do
    end do
  end subroutine

  subroutine device_params_init_contacts(this)
    class(device_params), intent(inout) :: this

    integer                           :: dim, i, i0(3), i1(3), ict, ict0, idx_dim, idx_dir, ijk(3), j, k, nct, ri
    integer, allocatable              :: idx(:)
    real                              :: dop(2), p(2)
    type(string)                      :: name
    type(mapnode_string_int), pointer :: node

    ! abbreviations
    dim     = this%g%dim
    idx_dim = this%g%idx_dim

    ! get all contact names and number of contacts (nct)
    call this%contact_map%init()
    nct = 0
    do ri = 1, size(this%reg_ct)
      ! search for contact name in map
      node => this%contact_map%find(this%reg_ct(ri)%name)

      ! do nothing if name already exists
      if (associated(node)) cycle

      ! new contact
      nct = nct + 1
      call this%contact_map%insert(this%reg_ct(ri)%name, nct)
    end do
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
          do k = i0(1), i1(1)
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

    ! finish initialization of contacts
    do ict = 1, nct
      call this%contacted(ict)%init_final()

      ! set phims for ohmic contacts
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
        call this%contacts(ict)%set_phims_ohmic(this%ci0, this%ci1, dop, this%smc)
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
  end subroutine

end module
