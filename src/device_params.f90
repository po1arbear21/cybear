m4_include(util/macro.f90.inc)

module device_params_m

  use bin_search_m,     only: bin_search
  use contact_m,        only: CT_OHMIC, CT_GATE, contact
  use error_m,          only: assert_failed, program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid, grid_ptr
  use grid_data_m,      only: allocate_grid_data0_real, allocate_grid_data1_real, allocate_grid_data2_real, &
    &                         allocate_grid_data3_real, allocate_grid_data0_int, allocate_grid_data1_int, &
    &                         grid_data_int, grid_data_real, grid_data2_int, grid_data2_real
  use grid_generator_m, only: DIR_NAME, generate_cartesian_grid, generate_triangle_grid
  use grid_table_m,     only: grid_table
  use grid1D_m,         only: grid1D
  use input_m,          only: input_file
  use map_m,            only: map_string_int, mapnode_string_int
  use math_m,           only: linspace
  use qsort_m,          only: qsort
  use region_m,         only: region, region_ptr, region_poisson, region_transport, region_doping, region_contact
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, CR_NAME, CR_CHARGE, semiconductor
  use string_m,         only: string, new_string
  use tensor_grid_m,    only: tensor_grid
  use triang_grid_m,    only: triang_grid
  use triangle_m,       only: triangulation
  use vector_m,         only: vector_real

  implicit none

  type device_params
    !! device geometry and material parameters

    integer           :: ci0, ci1
      !! enabled carrier index range (maximal: CR_ELEC..CR_HOLE)

    type(semiconductor) :: smc
    real              :: curr_fact
    integer           :: dim
      !! grid dimension

    class(grid_data_real), allocatable  :: eps(:,:)
      !! electrostatic permittivity (idx_type, idx_dir)
    class(grid_data_real), allocatable :: surf(:)
      !! adjoint volume surfaces
    class(grid_data_real), allocatable :: dop(:,:,:)
      !! donator/acceptor concentration (idx_type, idx_dir, dcon/acon)
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
    type(grid_table),      allocatable :: poisson_vct(:)
      !! poisson vertices grouped by contacts (0:size(contacts))
    type(grid_table),      allocatable :: oxide_vct(:)
      !! oxide vertices grouped by contacts (0:size(contacts))
    type(grid_table),      allocatable :: transport_vct(:)
      !! transport vertices grouped by contacts (0:size(contacts))

    type(string) :: gtype
      !! gitter name("x", "xy", "xyz", "tr_xy", "tr_xyz")
    type(region_ptr),       allocatable :: reg(:)
      !! pointer of regions
    type(region_poisson),   allocatable :: reg_poiss(:)
      !! poisson regions
    type(region_transport), allocatable :: reg_trans(:)
      !! transport regions
    type(region_doping),    allocatable :: reg_dop(:)
      !! doping regions
    type(region_contact),   allocatable :: reg_ct(:)
      !! contact regions
    type(grid1D)                        :: g1D(3)
      !! x, y, z grids
    type(tensor_grid)                   :: tg
      !! tensor grid
    type(triangulation)                 :: tr
      !! add polygons
    type(triang_grid)                   :: gtr
      !! triangle grid
    class(grid),            pointer     :: g => null()
    character(32),          allocatable :: idx_dir_name(:)
      !! "xyz"

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

contains

  subroutine device_params_init(this, file)
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file

    call this%init_transport_params(file)

    ! get grid type ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    call file%get("grid", "gtype", this%gtype)
    select case(this%gtype%s)
    case("x")
      this%dim = 1
    case("xy")
      this%dim = 2
    case("xyz")
      this%dim = 3
    case("tr_xy")
      this%dim = 2
    case("tr_xyz")
      this%dim = 3
    end select

    call this%init_regions(file)
    call this%init_grid(file)

    call this%init_poisson()
    call this%init_transport()
    call this%init_doping()
    call this%init_contacts()
  end subroutine

  subroutine device_params_destruct(this)
    !! destruct device parameter object
    class(device_params), intent(inout) :: this

    call this%contact_map%destruct()
  end subroutine

  subroutine device_params_init_transport_params(this, file)
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file

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
    call file%get(sid, "N_c",       this%smc%edos(CR_ELEC))
    call file%get(sid, "N_v",       this%smc%edos(CR_HOLE))
    call file%get(sid, "E_gap",     this%smc%band_gap)
    this%smc%n_intrin = sqrt(this%smc%edos(CR_ELEC) * this%smc%edos(CR_HOLE) * exp(-this%smc%band_gap))
    this%smc%band_edge(CR_ELEC) =   0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))
    this%smc%band_edge(CR_HOLE) = - 0.5 * this%smc%band_gap + 0.5 * log(this%smc%edos(CR_ELEC) / this%smc%edos(CR_HOLE))
    call file%get(sid, "mass",      this%smc%mass)
    call file%get(sid, "degen",     this%smc%degen)
    call file%get(sid, "alpha",     this%smc%alpha)
    call file%get(sid, "beta",      this%smc%beta)
    call file%get(sid, "mob_min",   this%smc%mob_min)
    call file%get(sid, "mob_max",   this%smc%mob_max)
    call file%get(sid, "N_ref",     this%smc%N_ref)
    call file%get(sid, "v_sat",     this%smc%v_sat)
    call file%get(sid, "curr_fact", this%curr_fact)

    ! make sure parameters are valid
    m4_assert(this%ci0 <= this%ci1)
    m4_assert(size(this%smc%mass) == 2)
    m4_assert(size(this%smc%alpha) == 2)
    m4_assert(size(this%smc%beta) == 2)
    m4_assert(size(this%smc%mob_min) == 2)
    m4_assert(size(this%smc%mob_max) == 2)
    m4_assert(size(this%smc%N_ref) == 2)
    m4_assert(size(this%smc%v_sat) == 2)
  end subroutine

  subroutine device_params_init_regions(this, file)
    class(device_params), target, intent(inout) :: this
    type(input_file),             intent(in)  :: file

    integer, allocatable :: sids(:)
    integer              :: si, i, j

    ! get poisson sections
    call file%get_sections("poisson", sids)
    ! initialize poisson regions
    allocate (this%reg_poiss(size(sids)))
    do si = 1, size(sids)
      call this%reg_poiss(si)%init(file, sids(si), this%gtype)
    end do

    ! get transport sections
    call file%get_sections("transport", sids)
    ! initialize transport regions
    allocate (this%reg_trans(size(sids)))
    do si = 1, size(sids)
      call this%reg_trans(si)%init(file, sids(si), this%gtype)
    end do

    ! get doping sections
    call file%get_sections("doping", sids)

    ! initialize doping regions
    allocate (this%reg_dop(size(sids)))
    do si = 1, size(sids)
      call this%reg_dop(si)%init(file, sids(si), this%gtype)
    end do

    ! get contact sections
    call file%get_sections("contact", sids)

    ! initialize contact regions
    allocate (this%reg_ct(size(sids)))
    do si = 1, size(sids)
      call this%reg_ct(si)%init(file, sids(si), this%gtype)
    end do

    allocate (this%reg(size(this%reg_poiss) + size(this%reg_trans) &
        &           + size(this%reg_dop) + size(this%reg_ct)))
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

    real, allocatable :: max_dxyz(:), max_areadz(:)
    real              :: max_dz, max_area
    integer           :: i

    select case(this%gtype%s)
    case("x", "xy", "xyz")
      allocate (max_dxyz(this%dim))
      do i = 1, this%dim
        call file%get("grid", "max_d"//DIR_NAME(i), max_dxyz(i))
      end do
      call generate_cartesian_grid(this%dim, this%reg, max_dxyz, this%g1D, this%g, this%tg)
      allocate (this%idx_dir_name(this%dim))
      do i = 1, this%dim
        this%idx_dir_name(i) = DIR_NAME(i)
      end do

    case("tr_xy", "tr_xyz")
      allocate (max_areadz(this%dim - 1))
      call file%get("grid", "max_area", max_area)
      max_areadz(1) = max_area
      if (this%dim == 3) then
        call file%get("grid", "max_dz", max_dz)
        max_areadz(2) = max_dz
      end if
      call generate_triangle_grid(this%dim, this%reg, max_areadz, this%tr, this%g, this%g1D, this%gtr, this%tg)
      allocate (this%idx_dir_name(this%dim-1))
      this%idx_dir_name(1) = "xy"
      if (this%dim == 3) this%idx_dir_name(2) = "z"
    end select
  end subroutine

  subroutine device_params_init_poisson(this)
    class(device_params), target, intent(inout) :: this

    integer,      allocatable  :: idx_c(:), idx_f(:), idx_v(:), idx_e(:), idx_v1(:), idx_v2(:)
    character(:), allocatable  :: table_name0, table_name
    integer  :: idx_dir, idx_dir2, idx_type, idx_dir0(4), idx_dir1(4), i0(3), i1(3), ii, jj, kk, i, j, k, si, ijk(3)
    real     :: vol, len, surf, edge(3), p(2,3), trsurf(3), mid(2)

    idx_dir0 = [ 0, 1, 1, 0]
    idx_dir1 = [ 0, this%g%idx_dim, this%g%idx_dim, 0]

    ! allocate grid data
    call allocate_grid_data2_real(this%eps, this%g%idx_dim, [1, 0], [4, this%g%idx_dim])
    call this%eps(IDX_CELL,0)%init(this%g, IDX_CELL, 0)
    call allocate_grid_data1_real(this%surf, this%g%idx_dim, 1, this%g%idx_dim)
    do idx_dir = 1, this%g%idx_dim
      call this%eps(IDX_EDGE,idx_dir)%init(this%g, IDX_EDGE, idx_dir)
      call this%surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
    end do
    call allocate_grid_data0_real(this%vol, this%g%idx_dim)
    call this%vol%init(this%g, IDX_VERTEX, 0)

    ! initialize poisson grid tables
    allocate (this%poisson(4,0:this%g%idx_dim))
    allocate (idx_v(this%g%idx_dim), idx_v1(this%g%idx_dim), idx_v2(this%g%idx_dim), idx_e(this%g%idx_dim), idx_f(this%g%idx_dim), idx_c(this%g%idx_dim))
    do idx_type = 1, 4
      table_name0 = "poisson_"//IDX_NAME(idx_type)(1:1)
      do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
        if (idx_dir > 0) then
          table_name = table_name0//this%idx_dir_name(idx_dir)
        else
          table_name = table_name0
        end if
        call this%poisson(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
      end do
    end do

    do si = 1, size(this%reg_poiss)
      select case (this%gtype%s)
      case ("x", "xy", "xyz")
        ! get bounds
        do idx_dir = 1, this%dim
          i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,1))
          i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,2)) - 1
        end do

        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%dim
            idx_c(idx_dir) = ijk(idx_dir)
          end do

          ! set permittivity for cell
          call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

          ! enable poisson for vertices and update adjoint volumes
          vol = this%g%get_vol(idx_c)
          do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
            ijk = [i+ii, j+jj, k+kk]
            do idx_dir = 1, this%dim
              idx_v(idx_dir) = ijk(idx_dir)
            end do
            call this%poisson(IDX_VERTEX,0)%flags%set(idx_v, .true.)
            call this%vol%update(idx_v, 0.125*vol)
          end do; end do; end do

          do idx_dir = 1, this%dim
            surf = this%g%get_surf(idx_c, idx_dir)

            ! edges
            do jj = 0, 1; do ii = 0, 1
              if (idx_dir == 1) then
                ijk = [i, j+ii, k+jj]
              elseif (idx_dir == 2) then
                ijk = [i+ii, j, k+jj]
              else
                ijk = [i+ii, j+jj, k]
              end if
              do idx_dir2 = 1, this%dim
                idx_e(idx_dir2) = ijk(idx_dir2)
              end do
              call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx_e, .true.)
              call this%surf(idx_dir)%update(idx_e, 0.25*surf)
              call this%eps(IDX_EDGE,idx_dir)%update(idx_e, 0.25*surf*this%reg_poiss(si)%eps)
            end do; end do

            ! faces
            do ii = 0, 1
              ijk = [i, j, k]
              ijk(idx_dir) = ijk(idx_dir) + ii
              do idx_dir2 = 1, this%dim
                idx_f(idx_dir2) = ijk(idx_dir2)
              end do
              call this%poisson(IDX_FACE,idx_dir)%flags%set(idx_f, .true.)
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do; end do; end do
        end do

      case("tr_xy", "tr_xyz")
        ! get z bounds
        if (this%dim == 3) then
          i0(1) = bin_search(this%g1D(3)%x, this%reg_poiss(si)%xyz(3,1))
          i1(1) = bin_search(this%g1D(3)%x, this%reg_poiss(si)%xyz(3,2)) - 1
        else
          i0(1) = 1
          i1(1) = 1
        end if

        ! if (this%dim == 3) mid(3) = (i1(1) + i0(1)) / 2
        ! enable poisson in region
        do i = 1, this%gtr%ncell
          call this%gtr%get_cell([i], p)
          mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
          mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
          if (.not. this%reg_poiss(si)%point_test(new_string("tr_xy"), mid)) cycle
          call this%gtr%adjoint(i, edge, trsurf)

          do k = i0(1), i1(1)
            idx_c(1) = i
            if (this%dim == 3) then
              idx_c(2) = k

              ! z edge length
              len = this%g1D(3)%get_len([k], 1)
            else
              ! factor of 2 cancels 0.5 in volume/surface calculation
              len = 2
            end if

            ! set permittivity for cell
            call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

            do ii = 1, 3
              idx_e(1)  = this%gtr%cell2edge(ii,i)
              idx_v1(1) = this%gtr%edge2vert(1,idx_e(1))
              idx_v2(1) = this%gtr%edge2vert(2,idx_e(1))
              vol = 0.125 * len * edge(ii) * trsurf(ii)
              do kk = 0, this%dim-2
                if (this%dim == 3) then
                  idx_v1(2) = k + kk
                  idx_v2(2) = k + kk
                  idx_e( 2) = k + kk
                end if
                call this%poisson(IDX_VERTEX,0)%flags%set(idx_v1, .true.)
                call this%poisson(IDX_VERTEX,0)%flags%set(idx_v2, .true.)
                call this%vol%update(idx_v1, vol)
                call this%vol%update(idx_v2, vol)
                call this%poisson(IDX_EDGE,1)%flags%set(idx_e, .true.)
                if (trsurf(ii)==0) then
                  call this%surf(1)%update(idx_e, 1e-16)
                else
                  call this%surf(1)%update(idx_e, trsurf(ii))
                end if
                call this%eps(IDX_EDGE,1)%update(idx_e, 0.5 * len * trsurf(ii) * this%reg_poiss(si)%eps)
              end do

              ! faces perpendicular to triangle
              idx_f(1) = idx_e(1)
              if (this%dim == 3) then
                idx_f(2) = k
              end if
              call this%poisson(IDX_FACE,1)%flags%set(idx_f, .true.)

              if (this%dim == 3) then
                ! edge perpendicular to triangle through first vertex
                idx_e = [this%gtr%edge2vert(1,this%gtr%cell2edge(ii,i)), k]
                call this%poisson(IDX_EDGE,2)%flags%set(idx_e, .true.)
                if (trsurf(ii) == 0) then
                  call this%surf(2)%update(idx_e, 0.25 * edge(ii) * 1e-16)
                else
                  call this%surf(2)%update(idx_e, 0.25 * edge(ii) * trsurf(ii))
                end if
                call this%eps(IDX_EDGE,2)%update(idx_e, 0.25 * edge(ii) * trsurf(ii) * this%reg_poiss(si)%eps)

                ! edge perpendicular to triangle through second vertex
                idx_e = [this%gtr%edge2vert(2,this%gtr%cell2edge(ii,i)), k]
                call this%poisson(IDX_EDGE,2)%flags%set(idx_e, .true.)
                if (trsurf(ii) == 0) then
                  call this%surf(2)%update(idx_e, 0.25 * edge(ii) * 1e-16)
                else
                  call this%surf(2)%update(idx_e, 0.25 * edge(ii) * trsurf(ii))
                end if
                call this%eps(IDX_EDGE,2)%update(idx_e, 0.25 * edge(ii) * trsurf(ii) * this%reg_poiss(si)%eps)

                ! faces in triangle
                do kk = 0, 1
                  idx_f = [i, k + kk]
                  call this%poisson(IDX_FACE,2)%flags%set(idx_f, .true.)
                end do
              end if
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do
        end do
      end select
    end do

    ! weighted average of permittivity on adjoint surfaces
    do idx_dir = 1, this%g%idx_dim
      call this%eps(IDX_EDGE,idx_dir)%set(this%eps(IDX_EDGE,idx_dir)%get() / this%surf(idx_dir)%get())
    end do

    ! finish initialization of poisson grid tables
    do idx_type = 1, 4
      do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
        call this%poisson(idx_type, idx_dir)%init_final()
      end do
    end do
  end subroutine

  subroutine device_params_init_transport(this)
    class(device_params), target, intent(inout) :: this

    integer,       allocatable :: idx_e(:), idx_c(:), idx_v(:), idx_f(:)
    integer,       allocatable :: idx_v1(:), idx_v2(:)
    character(:),  allocatable :: table_name0, table_name
    integer                    :: i0(3), i1(3), idx_type, i, j, k, ii, jj, kk, si, idx_dir, idx_dir2, idx_dir0(4), idx_dir1(4), ijk(3)
    real                       :: len, surf, vol, p(2,3), mid(2), edge(3), trsurf(3)

    idx_dir0 = [ 0, 1, 1, 0]
    idx_dir1 = [ 0, this%g%idx_dim, this%g%idx_dim, 0]
    ! allocate grid data
    call allocate_grid_data1_real(this%tr_surf, this%g%idx_dim, 1, this%g%idx_dim)
    call allocate_grid_data0_real(this%tr_vol, this%g%idx_dim)
    do idx_dir = 1, this%g%idx_dim
      call this%tr_surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
    end do
    call this%tr_vol%init(this%g, IDX_VERTEX, 0)

    ! initialize oxide and transport grid tables
    allocate (this%oxide(4,0:this%g%idx_dim))
    allocate (this%transport(4,0:this%g%idx_dim))
    allocate (idx_v(this%g%idx_dim), idx_v1(this%g%idx_dim), idx_v2(this%g%idx_dim), idx_e(this%g%idx_dim), idx_f(this%g%idx_dim), idx_c(this%g%idx_dim))
    do idx_type = 1, 4
      table_name0 = "oxide_"//IDX_NAME(idx_type)(1:1)
      do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
        if (idx_dir > 0) then
          table_name = table_name0//this%idx_dir_name(idx_dir)
        else
          table_name = table_name0
        end if
        call this%oxide(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        this%oxide(idx_type, idx_dir)%flags = this%poisson(idx_type, idx_dir)%flags
      end do
      table_name0 = "transport_"//IDX_NAME(idx_type)(1:1)
      do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
        if (idx_dir > 0) then
          ! dir = this%ig1D(idx_dir)
          table_name = table_name0//DIR_NAME(idx_dir)
        else
          table_name = table_name0
        end if
        call this%transport(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
      end do
    end do

    do si = 1, size(this%reg_trans)
      select case (this%gtype%s)
      case("x", "xy", "xyz")
        ! get bounds
        do idx_dir = 1, this%dim
          i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 1))
          i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 2)) - 1
        end do

        ! enable transport in region
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%dim
            idx_c(idx_dir) = ijk(idx_dir)
          end do

          ! enable transport for vertices and update adjoint volumes
          vol = this%g%get_vol(idx_c)
          do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
            ijk = [i+ii, j+jj, k+kk]
            do idx_dir = 1, this%g%idx_dim
              idx_v(idx_dir) = ijk(idx_dir)
            end do
            call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v, .false.)
            call this%transport(IDX_VERTEX,0)%flags%set(idx_v, .true.)
            call this%tr_vol%update(idx_v, 0.125*vol)
          end do; end do; end do

          do idx_dir = 1, this%dim
            surf = this%g%get_surf(idx_c, idx_dir)
            ! edges
            do jj = 0, 1; do ii = 0, 1
              if (idx_dir == 1) then
                ijk = [i, j+ii, k+jj]
              elseif (idx_dir == 2) then
                ijk = [i+ii, j, k+jj]
              else
                ijk = [i+ii, j+jj, k]
              end if
              do idx_dir2 = 1, this%dim
                idx_e(idx_dir2) = ijk(idx_dir2)
              end do
              ! enable transport for edges
              call this%oxide(    IDX_EDGE,idx_dir)%flags%set(idx_e, .false.)
              call this%transport(IDX_EDGE,idx_dir)%flags%set(idx_e, .true. )
              call this%tr_surf(idx_dir)%update(idx_e, 0.25*surf)
            end do; end do

            ! faces
            do ii = 0, 1
              ijk = [i, j, k]
              ijk(idx_dir) = ijk(idx_dir) + ii
              do idx_dir2 = 1, this%dim
                idx_f(idx_dir2) = ijk(idx_dir2)
              end do
              call this%oxide(    IDX_FACE,idx_dir)%flags%set(idx_f, .false.)
              call this%transport(IDX_FACE,idx_dir)%flags%set(idx_f, .true. )
            end do
          end do

          ! enable transport for cell
          call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
          call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
        end do; end do; end do

      case ("tr_xy", "tr_xyz")
        ! get z bounds
        if (this%dim == 3) then
          i0(1) = bin_search(this%g1D(3)%x, this%reg_trans(si)%xyz(3,1))
          i1(1) = bin_search(this%g1D(3)%x, this%reg_trans(si)%xyz(3,2)) - 1
        else
          i0(1) = 1
          i1(1) = 1
        end if

        ! enable transport in region
        do i = 1, this%gtr%ncell
          call this%gtr%get_cell([i], p)
          mid(1) = (p(1, 1) + p(1, 2) + p(1, 3))/3
          mid(2) = (p(2, 1) + p(2, 2) + p(2, 3))/3
          if (.not. this%reg_trans(si)%point_test(new_string("tr_xy"), mid)) cycle
          call this%gtr%adjoint(i, edge, trsurf)

          do k = i0(1), i1(1)
            idx_c(1) = i
            if (this%dim == 3) then
              idx_c(2) = k
              ! z edge length
              len = this%g1D(3)%get_len([k], 1)
            else
              ! factor of 2 cancels 0.5 in volume/surface calculation
              len = 2
            end if

            do ii = 1, 3
              ! edges in triangle
              idx_e(1) = this%gtr%cell2edge(ii,i)
              idx_v1(1) = this%gtr%edge2vert(1,idx_e(1))
              idx_v2(1) = this%gtr%edge2vert(2,idx_e(1))
              vol = 0.125 * len * edge(ii) *trsurf(ii)
              do kk = 0, this%dim-2
                if (this%dim == 3) then
                  idx_v1(2) = k + kk
                  idx_v2(2) = k + kk
                  idx_e(2) = k + kk
                end if
                ! enable transport for vertices and update adjoint volumes
                call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v1, .false.)
                call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v2, .false.)
                call this%transport(IDX_VERTEX,0)%flags%set(idx_v1, .true.)
                call this%transport(IDX_VERTEX,0)%flags%set(idx_v2, .true.)
                call this%tr_vol%update(idx_v1, vol)
                call this%tr_vol%update(idx_v2, vol)

                ! enable transport for edges in triangle
                call this%oxide(    IDX_EDGE,1)%flags%set(idx_e, .false.)
                call this%transport(IDX_EDGE,1)%flags%set(idx_e, .true. )
                call this%tr_surf(1)%update(idx_e, 0.5 * len * trsurf(ii))
              end do

              ! faces perpendicular to triangle
              idx_f(1) = idx_e(1)
              if (this%dim == 3) then
                idx_f(2) = k
              end if
              call this%oxide(    IDX_FACE,1)%flags%set(idx_f, .false.)
              call this%transport(IDX_FACE,1)%flags%set(idx_f, .true. )

              if (this%dim == 3) then
                ! edge perpendicular to triangle through first vertex
                idx_e = [this%gtr%edge2vert(1,this%gtr%cell2edge(ii,i)), k]
                call this%oxide(    IDX_EDGE,2)%flags%set(idx_e, .false.)
                call this%transport(IDX_EDGE,2)%flags%set(idx_e, .true. )
                call this%tr_surf(2)%update(idx_e, 0.25 * edge(ii) *trsurf(ii))

                ! edge perpendicular to triangle through second vertex
                idx_e = [this%gtr%edge2vert(2,this%gtr%cell2edge(ii,i)), k]
                call this%oxide(    IDX_EDGE,2)%flags%set(idx_e, .false.)
                call this%transport(IDX_EDGE,2)%flags%set(idx_e, .true. )
                call this%tr_surf(2)%update(idx_e, 0.25 * edge(ii) *trsurf(ii))

                ! faces in triangle
                do kk = 0, 1
                  idx_f = [i, k+kk]
                  call this%oxide(    IDX_FACE,2)%flags%set(idx_f, .false.)
                  call this%transport(IDX_FACE,2)%flags%set(idx_f, .true. )
                end do
              end if
            end do
          end do

          ! enable transport for cell
          call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
          call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
        end do
      end select
    end do

    ! finish initialization of oxide, transport grid tables
    do idx_type = 1, 4
      do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
        call this%oxide(    idx_type, idx_dir)%init_final()
        call this%transport(idx_type, idx_dir)%init_final()
      end do
    end do
  end subroutine

  subroutine device_params_init_doping(this)
    class(device_params), target, intent(inout) :: this

    integer, allocatable :: idx_e(:), idx_c(:), idx_v(:), idx_v1(:), idx_v2(:), idx(:)
    integer              :: i0(3), i1(3), i, j, k, si, ii, jj, kk, ci, idx_dir, idx_dir2, idx_dir0(4), idx_dir1(4), idx_type, ijk(3)
    real                 :: mob_min, mob_max, N_ref, alpha, tr_surf, dop(2), acon, dcon, tr_vol(2), vol, len, surf, edge(3), trsurf(3), p(2,3), mid(2)

    idx_dir0 = [ 0, 1, 1, 0]
    idx_dir1 = [ 0, this%g%idx_dim, this%g%idx_dim, 0]
    ! allocate grid data
    call allocate_grid_data3_real(this%dop, this%g%idx_dim, [1, 0, CR_ELEC], [4, this%g%idx_dim, CR_HOLE])
    call allocate_grid_data3_real(this%mob0, this%g%idx_dim, [1, 0, this%ci0], [4, this%g%idx_dim, this%ci1])
    allocate (idx_v(this%g%idx_dim), idx_v1(this%g%idx_dim), idx_v2(this%g%idx_dim), idx_e(this%g%idx_dim), idx_c(this%g%idx_dim))
    do ci = CR_ELEC, CR_HOLE
      do idx_type = 1, 4
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          call this%dop(idx_type,idx_dir,ci)%init(this%g, idx_type, idx_dir)
        end do
      end do
    end do

    do si = 1, size(this%reg_dop)
      dop(CR_ELEC) = this%reg_dop(si)%dcon
      dop(CR_HOLE) = this%reg_dop(si)%acon

      select case(this%gtype%s)
      case("x", "xy", "xyz")
        do idx_dir = 1, this%dim
          i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 1))
          i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 2)) - 1
        end do

        ! set doping in region
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%dim
            idx_c(idx_dir) = ijk(idx_dir)
          end do
          vol = this%g%get_vol(idx_c)

          ! set doping on vertices
          do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
            ijk = [i+ii, j+jj, k+kk]
            do idx_dir = 1, this%dim
              idx_v(idx_dir) = ijk(idx_dir)
            end do
            tr_vol(1) = this%tr_vol%get(idx_v)
            do ci = CR_ELEC, CR_HOLE
              call this%dop(IDX_VERTEX,0,ci)%update(idx_v, 0.125*vol/tr_vol(1)*dop(ci))
            end do
          end do; end do; end do

          ! set doping for edges
          do idx_dir = 1, this%dim
            surf    = this%g%get_surf(idx_c, idx_dir)
            tr_surf = this%tr_surf(idx_dir)%get(idx_c)
            do jj = 0, 1; do ii = 0, 1
              if (idx_dir == 1) then
                ijk = [i, j+ii, k+jj]
              elseif (idx_dir == 2) then
                ijk = [i+ii, j, k+jj]
              else
                ijk = [i+ii, j+jj, k]
              end if
              do idx_dir2 = 1, this%dim
                idx_e(idx_dir2) = ijk(idx_dir2)
              end do
              do ci = CR_ELEC, CR_HOLE
                call this%dop(IDX_EDGE,idx_dir,ci)%update(idx_e, 0.25*surf/tr_surf*dop(ci))
              end do
            end do; end do
          end do

          ! set doping in cells
          do ci = CR_ELEC, CR_HOLE
            call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
          end do
        end do; end do; end do

      case("tr_xy", "tr_xyz")
        ! get z bounds
        if (this%dim == 3) then
          i0(1) = bin_search(this%g1D(3)%x, this%reg_dop(si)%xyz(3,1))
          i1(1) = bin_search(this%g1D(3)%x, this%reg_dop(si)%xyz(3,2)) - 1
        else
          i0(1) = 1
          i1(1) = 1
        end if

        ! set doping in region
        do i = 1, this%gtr%ncell
          call this%gtr%get_cell([i], p)
          mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
          mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
          if (.not. this%reg_dop(si)%point_test(new_string("tr_xy"), mid)) cycle
          call this%gtr%adjoint(i, edge, trsurf)

          do k = i0(1), i1(1)
            idx_c(1) = i
            if (this%dim == 3) then
              idx_c(2) = k
              ! z edge length
              len = this%g1D(3)%get_len([k], 1)
            else
              ! factor of 2 cancels 0.5 in volume/surface calculation
              len = 2
            end if

            do ii = 1, 3
              ! set doping on vertices
              idx_e(1)  = this%gtr%cell2edge(ii,i)
              idx_v1(1) = this%gtr%edge2vert(1,idx_e(1))
              idx_v2(1) = this%gtr%edge2vert(2,idx_e(1))
              vol = 0.125 * len * edge(ii) * trsurf(ii)
              do kk = 0, this%dim-2
                if (this%dim == 3) then
                  idx_v1(2) = k + kk
                  idx_v2(2) = k + kk
                  idx_e(2) = k + kk
                end if
                tr_vol(1) = this%tr_vol%get(idx_v1)
                tr_vol(2) = this%tr_vol%get(idx_v2)
                do ci = CR_ELEC, CR_HOLE
                  call this%dop(IDX_VERTEX,0,ci)%update(idx_v1, vol/tr_vol(1)*dop(ci))
                  call this%dop(IDX_VERTEX,0,ci)%update(idx_v2, vol/tr_vol(2)*dop(ci))
                end do

                ! set doping for edges in triangles
                tr_surf = this%tr_surf(1)%get(idx_e)
                do ci = CR_ELEC, CR_HOLE
                  if (tr_surf == 0) then
                    call this%dop(IDX_EDGE,1,ci)%update(idx_e, 0.5 * len * trsurf(ii)/1e-16*dop(ci))
                  else
                    call this%dop(IDX_EDGE,1,ci)%update(idx_e, 0.5 * len * trsurf(ii)/tr_surf*dop(ci))
                  end if
                end do
              end do

              if (this%dim == 3) then
                ! edges perpendicular to triangle through first vertex
                idx_e = [this%gtr%edge2vert(1,this%gtr%cell2edge(ii,i)), k]
                tr_surf = this%tr_surf(2)%get(idx_e)
                do ci = CR_ELEC, CR_HOLE
                  if (tr_surf == 0) then
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/1e-16*dop(ci))
                  else
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/tr_surf*dop(ci))
                  end if
                end do

                ! edges perpendicular to triangle through second vertex
                idx_e = [this%gtr%edge2vert(2,this%gtr%cell2edge(ii,i)), k]
                tr_surf = this%tr_surf(2)%get(idx_e)
                do ci = CR_ELEC, CR_HOLE
                  if (tr_surf == 0) then
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/1e-16*dop(ci))
                  else
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/tr_surf*dop(ci))
                  end if
                end do
              end if
            end do
          end do

          ! set doping in cells
          do ci = CR_ELEC, CR_HOLE
            call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
          end do
        end do
      end select
    end do

    ! init zero-field mobility
    acon = 0
    dcon = 0
    do ci = this%ci0, this%ci1
      mob_min = this%smc%mob_min(ci)
      mob_max = this%smc%mob_max(ci)
      N_ref   = this%smc%N_ref(ci)
      alpha   = this%smc%alpha(ci)

      ! mobility in cells
      call this%mob0(IDX_CELL,0,ci)%init(this%g, IDX_CELL, 0)
      do i = 1, this%transport(IDX_CELL,0)%n
        idx = this%transport(IDX_CELL,0)%get_idx(i)
        if (this%ci0 == CR_ELEC) dcon = this%dop(IDX_CELL,0,CR_ELEC)%get(idx)
        if (this%ci1 == CR_HOLE) acon = this%dop(IDX_CELL,0,CR_HOLE)%get(idx)
        call this%mob0(IDX_CELL,0,ci)%set(idx, mob_min + (mob_max - mob_min)/(1 + ((acon + dcon)/N_ref)**alpha))
      end do

      ! mobility on edges
      do idx_dir = 1, this%g%idx_dim
        call this%mob0(IDX_EDGE,idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)
        do i = 1, this%transport(IDX_EDGE,idx_dir)%n
          idx  = this%transport(IDX_EDGE,idx_dir)%get_idx(i)
          if (this%ci0 == CR_ELEC) dcon = this%dop(IDX_EDGE,idx_dir,CR_ELEC)%get(idx)
          if (this%ci1 == CR_HOLE) acon = this%dop(IDX_EDGE,idx_dir,CR_HOLE)%get(idx)
          call this%mob0(IDX_EDGE,idx_dir,ci)%set(idx, mob_min + (mob_max - mob_min)/(1 + ((acon + dcon)/N_ref)**alpha))
        end do
      end do
    end do
  end subroutine

  subroutine device_params_init_contacts(this)
    class(device_params), target, intent(inout) :: this

    type(mapnode_string_int), pointer     :: node
    integer,                  allocatable :: idx_v(:)
    integer       :: ict, ict0, nct, ct_type
    type(string)  :: name
    integer       :: i0(3), i1(3), i, j, k, si, idx_dir, ijk(3)
    real          :: p(2,3), dcon, acon

    ! allocate grid data
    call allocate_grid_data0_int(this%ict, this%g%idx_dim)
    allocate (idx_v(this%g%idx_dim))
    ! get all contact names
    call this%contact_map%init()
    nct = 0
    do si = 1, size(this%reg_ct)
      node => this%contact_map%find(this%reg_ct(si)%name)
      if (.not. associated(node)) then
        ! new contact
        nct = nct + 1
        call this%contact_map%insert(this%reg_ct(si)%name, nct)
      end if
    end do

    ! init grid tables
    allocate (this%poisson_vct(  0:nct))
    allocate (this%oxide_vct(    0:nct))
    allocate (this%transport_vct(0:nct))
    call this%poisson_vct(  0)%init("poisson_VCT0",   this%g, IDX_VERTEX, 0)
    call this%oxide_vct(    0)%init("oxide_VCT0",     this%g, IDX_VERTEX, 0)
    call this%transport_vct(0)%init("transport_VCT0", this%g, IDX_VERTEX, 0)
    this%poisson_vct(  0)%flags = this%poisson(  IDX_VERTEX,0)%flags
    this%oxide_vct(    0)%flags = this%oxide(    IDX_VERTEX,0)%flags
    this%transport_vct(0)%flags = this%transport(IDX_VERTEX,0)%flags
    call this%ict%init(this%g, IDX_VERTEX, 0)

    ! init contacts
    allocate (this%contacts(nct))
    do si = 1, size(this%reg_ct)
      ! get name
      name = this%reg_ct(si)%name
      ict = this%contact_map%get(name)

      ! get type
      if (this%reg_ct(si)%type%s == "ohmic") then
        ct_type = CT_OHMIC
      elseif (this%reg_ct(si)%type%s == "gate") then
        ct_type = CT_GATE
      else
        call program_error("unknown contact type "//this%reg_ct(si)%type%s)
      end if

      if (.not. allocated(this%contacts(ict)%name)) then
        ! new contact
        this%contacts(ict)%name  = name%s
        this%contacts(ict)%type  = ct_type
        call this%poisson_vct(  ict)%init("poisson_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        call this%oxide_vct(    ict)%init("oxide_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        call this%transport_vct(ict)%init("transport_VCT_"//name%s, this%g, IDX_VERTEX, 0)
      elseif (this%contacts(ict)%type /= ct_type) then
        call program_error("contact "//name%s//" given multiple times with different types")
      end if

      select case (this%gtype%s)
      case("x", "xy", "xyz")
        ! bounds
        i0 = 1
        i1 = 1
        do idx_dir = 1, this%dim
          i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 1))
          i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 2))
        end do

        ! phims
        if (ct_type == CT_OHMIC) then
          outer: do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            do idx_dir = 1, this%dim
              idx_v(idx_dir) = ijk(idx_dir)
            end do
            if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit outer
          end do; end do; end do outer
          if ((i > i1(1)) .or. (j > i1(2)).or. (k > i1(3))) call program_error("Ohmic contact "//name%s//" not in transport region")
          dcon = this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v)
          acon = this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)
          call this%contacts(ict)%set_phims_ohmic(CR_ELEC, CR_HOLE, dcon, acon, this%smc)
        else
          this%contacts(ict)%phims = this%reg_ct(si)%phims
        end if

        ! update vertex tables
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%dim
            idx_v(idx_dir) = ijk(idx_dir)
          end do
          ict0 = this%ict%get(idx_v)
          if ((ict0 /= 0) .and. (ict0 /= ict)) then
            call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
          end if
          call this%ict%set(idx_v, ict)
          if (this%poisson(IDX_VERTEX,0)%flags%get(idx_v)) then
            call this%poisson_vct(ict)%flags%set(idx_v, .true. )
            call this%poisson_vct(  0)%flags%set(idx_v, .false.)
          end if
          if (this%oxide(IDX_VERTEX,0)%flags%get(idx_v)) then
            call this%oxide_vct(ict)%flags%set(idx_v, .true. )
            call this%oxide_vct(  0)%flags%set(idx_v, .false.)
          end if
          if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) then
            call this%transport_vct(ict)%flags%set(idx_v, .true. )
            call this%transport_vct(  0)%flags%set(idx_v, .false.)
          end if
        end do; end do; end do

      case("tr_xy", "tr_xyz")
        ! get z bounds
        if (this%dim == 3) then
          i0(1) = bin_search(this%g1D(3)%x, this%reg_ct(si)%xyz(3,1))
          i1(1) = bin_search(this%g1D(3)%x, this%reg_ct(si)%xyz(3,2)) - 1
        else
          i0(1) = 1
          i1(1) = 1
        end if

        do i = 1, this%gtr%nvert
          call this%gtr%get_vertex([i], p(:,1))
          if (.not. this%reg_ct(si)%point_test(new_string("tr_xy"), p(:,1))) cycle

          do k = i0(1), i1(1)
            idx_v(1) = i
            if (this%dim == 3) then
              idx_v(2) = k
            end if

            ! phims
            if (ct_type == CT_OHMIC) then
              dcon = this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v)
              acon = this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)
              call this%contacts(ict)%set_phims_ohmic(CR_ELEC, CR_HOLE, dcon, acon, this%smc)
            else
              this%contacts(ict)%phims = this%reg_ct(si)%phims
            end if

            ! update vertex tables
            ict0 = this%ict%get(idx_v)
            if ((ict0 /= 0) .and. (ict0 /= ict)) then
              call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
            end if
            call this%ict%set(idx_v, ict)
            if (this%poisson(IDX_VERTEX,0)%flags%get(idx_v)) then
              call this%poisson_vct(ict)%flags%set(idx_v, .true. )
              call this%poisson_vct(  0)%flags%set(idx_v, .false.)
            end if
            if (this%oxide(IDX_VERTEX,0)%flags%get(idx_v)) then
              call this%oxide_vct(ict)%flags%set(idx_v, .true. )
              call this%oxide_vct(  0)%flags%set(idx_v, .false.)
            end if
            if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) then
              call this%transport_vct(ict)%flags%set(idx_v, .true. )
              call this%transport_vct(  0)%flags%set(idx_v, .false.)
            end if
          end do
        end do
      end select
    end do

    ! finish initialization of grid tables
    call this%poisson_vct(  0)%init_final()
    call this%oxide_vct(    0)%init_final()
    call this%transport_vct(0)%init_final()
    do ict = 1, nct
      call this%poisson_vct(  ict)%init_final()
      call this%oxide_vct(    ict)%init_final()
      call this%transport_vct(ict)%init_final()
    end do
  end subroutine

end module
