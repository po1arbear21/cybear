m4_include(util/macro.f90.inc)

module device_params_m

  use bin_search_m,  only: bin_search
  use contact_m,     only: CT_OHMIC, CT_GATE, contact
  use error_m,       only: assert_failed, program_error
  use grid_m,        only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid, grid_ptr
  use grid_data_m,   only: allocate_grid_data0_real, allocate_grid_data1_real, allocate_grid_data2_real, allocate_grid_data3_real, &
  &                        allocate_grid_data0_int, allocate_grid_data1_int, grid_data_int, grid_data_real, grid_data2_int, grid_data2_real
  use grid_table_m,  only: grid_table
  use grid1D_m,      only: grid1D
  use input_m,       only: input_file
  use map_m,         only: map_string_int, mapnode_string_int
  use math_m,        only: linspace
  use qsort_m,       only: qsort
  use region_m,      only: region, region_poisson, region_transport, region_doping, region_contact
  use string_m,      only: string, new_string
  use tensor_grid_m, only: tensor_grid
  use triang_grid_m, only: triang_grid
  use triangle_m,    only: triangulation
  use vector_m,      only: vector_real

  implicit none

  ! parameters
  character(*), parameter :: DIR_NAME(3)  = [ "x", "y", "z"]
  integer,      parameter :: CR_ELEC      = 1
  integer,      parameter :: CR_HOLE      = 2
  character(*), parameter :: CR_NAME(2)   = [ "n", "p" ]
  real,         parameter :: CR_CHARGE(2) = [ -1.0, +1.0 ]

  type device_params
    !! device geometry and material parameters

    real              :: n_intrin
      !! intrinsic carrier density
    integer           :: ci0, ci1
      !! enabled carrier index range (maximal: CR_ELEC..CR_HOLE)
    real, allocatable :: mass(:)
      !! effective mass
    real, allocatable :: alpha(:)
      !! Caughey-Thomas alpha parameter
    real, allocatable :: beta(:)
      !! Caughey-Thomas beta parameter
    real, allocatable :: mob_min(:)
      !! Caughey-Thomas minimal mobility
    real, allocatable :: mob_max(:)
      !! Caughey-Thomas maximal mobility
    real, allocatable :: N_ref(:)
      !! Caughey-Thomas reference density
    real, allocatable :: v_sat(:)
      !! Caughey-Thomas saturation velocity
    real              :: curr_fact

    class(grid_data_real), allocatable  :: eps(:,:)
      !! electrostatic permittivity (idx_type, idx_dir)
    class(grid_data_real), allocatable :: surf(:)
      !! adjoint volume surfaces
    class(grid_data_real), allocatable :: dop(:,:,:)
      !! acceptor/donator concentration (idx_type, idx_dir, carrier index)
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
  contains
    procedure :: init     => device_params_init
    procedure :: destruct => device_params_destruct
  end type

contains

  subroutine device_params_init(this, file)
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file

    integer,      allocatable :: sids(:), idx_v(:), idx_v1(:), idx_v2(:), idx_e(:), idx_f(:), idx_c(:), idx(:)
    integer                   :: si, i, j, k, ii, jj, kk, i0(3), i1(3), idx_dir, idx_type, idx_dir0(4), idx_dir1(4), ci
    real                      :: vol, mid(2), surf, trsurf(3), edge(3), p(2,3), len, R
    character(:), allocatable :: table_name0, table_name
    character(32), allocatable :: idx_dir_name(:)

    call init_transport_params()

    ! get region name ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    call file%get("grid", "gtype", this%gtype)
print *, "A"
    call init_regions()
print *, "B"
    call init_grid()
print *, "C"

    idx_dir0 = [ 0, 1, 1, 0]
    idx_dir1 = [ 0, this%g%idx_dim, this%g%idx_dim, 0]
    allocate (idx_v(this%g%idx_dim), idx_v1(this%g%idx_dim), idx_v2(this%g%idx_dim), idx_e(this%g%idx_dim), idx_f(this%g%idx_dim), idx_c(this%g%idx_dim))

    call init_poisson()
    print *, "D"
    call init_transport()
    print *, "E"
    call init_doping()
    print *, "F"
    call init_contacts()
    print *, "G"

  contains

    subroutine init_transport_params()
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
      call file%get(sid, "n_intrin",  this%n_intrin)
      call file%get(sid, "mass",      this%mass)
      call file%get(sid, "alpha",     this%alpha)
      call file%get(sid, "beta",      this%beta)
      call file%get(sid, "mob_min",   this%mob_min)
      call file%get(sid, "mob_max",   this%mob_max)
      call file%get(sid, "N_ref",     this%N_ref)
      call file%get(sid, "v_sat",     this%v_sat)
      call file%get(sid, "curr_fact", this%curr_fact)

      ! make sure parameters are valid
      m4_assert(this%ci0 <= this%ci1)
      m4_assert(size(this%mass) == 2)
      m4_assert(size(this%alpha) == 2)
      m4_assert(size(this%beta) == 2)
      m4_assert(size(this%mob_min) == 2)
      m4_assert(size(this%mob_max) == 2)
      m4_assert(size(this%N_ref) == 2)
      m4_assert(size(this%v_sat) == 2)
    end subroutine

    subroutine init_regions()

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
    end subroutine

    subroutine init_grid()
      real              :: max_dx, max_dy, max_dz, max_area
      real, allocatable :: xyz(:), x0(:), y0(:), z0(:)
      type(grid_ptr)    :: gptr(3)

      select case(this%gtype%s)
      case("x")
        ! get max_dx
        call file%get("grid", "max_dx", max_dx)

        ! get x bounds
        call get_bounds(1, x0)
        call generate_axis(x0, max_dx, xyz)

        ! create grids
        call this%g1D(1)%init(DIR_NAME(1), xyz)
        this%g => this%g1D(1)

        allocate (idx_dir_name(1))
        idx_dir_name(1) = "x"

      case("xy")
        call file%get("grid", "max_dx", max_dx)
        call get_bounds(1, x0)
        call generate_axis(x0, max_dx, xyz)
        call this%g1D(1)%init(DIR_NAME(1), xyz)
        gptr(1) = this%g1D(1)%get_ptr()

        call file%get("grid", "max_dy", max_dy)
        call get_bounds(2, y0)
        if (allocated(xyz)) deallocate(xyz)
        call generate_axis(y0, max_dy, xyz)
        call this%g1D(2)%init(DIR_NAME(2), xyz)
        gptr(2) = this%g1D(2)%get_ptr()

        ! tensor grid
        call this%tg%init("grid", gptr(1: 2))
        this%g => this%tg

        allocate (idx_dir_name(2))
        idx_dir_name(1) = "x"
        idx_dir_name(2) = "y"

      case("xyz")
        call file%get("grid", "max_dx", max_dx)
        call get_bounds(1, x0)
        call generate_axis(x0, max_dx, xyz)
        call this%g1D(1)%init(DIR_NAME(1), xyz)
        gptr(1) = this%g1D(1)%get_ptr()

        call file%get("grid", "max_dy", max_dy)
        call get_bounds(2, y0)
        if (allocated(xyz)) deallocate(xyz)
        call generate_axis(y0, max_dy, xyz)
        call this%g1D(2)%init(DIR_NAME(2), xyz)
        gptr(2) = this%g1D(2)%get_ptr()

        call file%get("grid", "max_dz", max_dz)
        call get_bounds(3, z0)
        if (allocated(xyz)) deallocate(xyz)
        call generate_axis(z0, max_dz, xyz)
        call this%g1D(3)%init(DIR_NAME(3), xyz)
        gptr(3) = this%g1D(3)%get_ptr()

        ! tensor grid
        call this%tg%init("grid", gptr(1: 3))
        this%g => this%tg

        allocate (idx_dir_name(3))
        idx_dir_name(1) = "x"
        idx_dir_name(2) = "y"
        idx_dir_name(3) = "z"


      case("tr_xy")
        call file%get("grid", "max_area", max_area)

        call this%tr%init()
        do i = 1, size(this%reg_poiss)
          call this%tr%add_polygon(this%reg_poiss(i)%xyz(1,:), this%reg_poiss(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_trans)
          call this%tr%add_polygon(this%reg_trans(i)%xyz(1,:), this%reg_trans(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_dop)
          call this%tr%add_polygon(this%reg_dop(i)%xyz(1,:), this%reg_dop(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_ct)
          call this%tr%add_polygon(this%reg_ct(i)%xyz(1,:), this%reg_ct(i)%xyz(2,:), closed = .false.)
        end do
        call this%tr%triangulate(max_area = max_area)
        call this%tr%get_grid("tr", this%gtr)
        call this%gtr%output_plotmtv("grid.csv", "nm")

        this%g => this%gtr

        allocate (idx_dir_name(1))
        idx_dir_name(1) = "xy"

      case("tr_xyz")
        call file%get("grid", "max_area", max_area)

        call this%tr%init()
        do i = 1, size(this%reg_poiss)
          call this%tr%add_polygon(this%reg_poiss(i)%xyz(1,:), this%reg_poiss(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_trans)
          call this%tr%add_polygon(this%reg_trans(i)%xyz(1,:), this%reg_trans(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_dop)
          call this%tr%add_polygon(this%reg_dop(i)%xyz(1,:), this%reg_dop(i)%xyz(2,:), closed = .true.)
        end do
        do i = 1, size(this%reg_ct)
          call this%tr%add_polygon(this%reg_ct(i)%xyz(1,:), this%reg_ct(i)%xyz(2,:), closed = .false.)
        end do
        call this%tr%triangulate(max_area = max_area)
        call this%tr%get_grid("tr", this%gtr)

        gptr(1) = this%gtr%get_ptr()
        call file%get("grid", "max_dz", max_dz)
        call get_bounds(3, z0)
        if (allocated(xyz)) deallocate(xyz)
        call generate_axis(z0, max_dz, xyz)
        call this%g1D(3)%init(DIR_NAME(3), xyz)
        gptr(2) = this%g1D(3)%get_ptr()

        ! tensor grid
        call this%tg%init("grid", gptr(1: 2))
        this%g => this%tg

        allocate (idx_dir_name(2))
        idx_dir_name(1) = "xy"
        idx_dir_name(2) = "z"
      end select

    end subroutine

    subroutine get_bounds(dim, x0)

      integer,           intent(in)  :: dim
      real, allocatable, intent(out) :: x0(:)

      allocate (x0(2*(size(this%reg_poiss) + size(this%reg_trans) &
        &           + size(this%reg_dop) + size(this%reg_ct))))
      j = 0
      do i = 1, size(this%reg_poiss)
        x0(j+1:j+2) = this%reg_poiss(i)%xyz(dim,1:2)
        j = j + 2
      end do
      do i = 1, size(this%reg_trans)
        x0(j+1:j+2) = this%reg_trans(i)%xyz(dim,1:2)
        j = j + 2
      end do
      do i = 1, size(this%reg_dop)
        x0(j+1:j+2) = this%reg_dop(i)%xyz(dim,1:2)
        j = j + 2
      end do
      do i = 1, size(this%reg_ct)
        x0(j+1:j+2) = this%reg_ct(i)%xyz(dim,1:2)
        j = j + 2
      end do
    end subroutine

    subroutine generate_axis(x0, dx, x)
      real,              intent(in)  :: x0(:)
      real,              intent(in)  :: dx
      real, allocatable, intent(out) :: x(:)

      real, allocatable :: xtmp(:), xx(:)
      integer           :: n, m
      type(vector_real) :: xv

      real, parameter :: TOL = 1e-10

      allocate (xtmp(size(x0)))
      xtmp = x0
      call qsort(xtmp)
      n = size(xtmp)

      ! remove duplicates
      j = 1
      do i = 2 , n
        if (abs(xtmp(i)-xtmp(j)) > TOL) then
          j = j + 1
          xtmp(j) = xtmp(i)
        end if
      end do
      n = j

      call xv%init(0, c = 2 * ceiling((xtmp(n) - xtmp(1)) / dx))
      call xv%push(xtmp(1))
      do i = 1, n-1
        m  = ceiling((xtmp(i+1)-xtmp(i))/dx) + 1
        if (allocated(xx)) deallocate (xx)
        xx = linspace(xtmp(i), xtmp(i+1), m)
        call xv%push(xx(2:m))
      end do
      x= xv%to_array()
    end subroutine

    subroutine init_poisson()
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
      do idx_type = 1, 4
        table_name0 = "poisson_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) then
            table_name = table_name0//idx_dir_name(idx_dir)
          else
            table_name = table_name0
          end if
          call this%poisson(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        end do
      end do

      do si = 1, size(this%reg_poiss)
        select case (this%gtype%s)
        case ("x")
          ! get bounds
          i0(1) = bin_search(this%g1D(1)%x, this%reg_poiss(si)%xyz(1,1))
          i1(1) = bin_search(this%g1D(1)%x, this%reg_poiss(si)%xyz(1,2)) - 1

          ! enable poisson in region
          do i = i0(1), i1(1)
            idx_c(1) = i
            ! set permittivity for cell
            call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

            ! enable poisson for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do ii = 0, 1
              idx_v = i+ii
              call this%poisson(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%vol%update(idx_v, 0.5*vol)
            end do

            ! edges
            surf = this%g%get_surf(idx_c, 1)
            call this%poisson(IDX_EDGE,1)%flags%set(idx_c, .true.)
            call this%surf(1)%update(idx_c, surf)
            call this%eps(IDX_EDGE,1)%update(idx_c, surf*this%reg_poiss(si)%eps)

            ! faces
            do ii = 0, 1
              idx_f(1) = i+ii
              call this%poisson(IDX_FACE,1)%flags%set(idx_f, .true.)
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do

        case ("xy")
          ! get bounds
          do idx_dir = 1, 2
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,2)) - 1
          end do

          ! enable poisson in region
          do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j]

            ! set permittivity for cell
            call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

            ! enable poisson for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj]
              call this%poisson(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%vol%update(idx_v, 0.25*vol)
            end do; end do

            do idx_dir = 1, 2
              surf = this%g%get_surf(idx_c, idx_dir)

              ! edges
              do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii]
                elseif (idx_dir == 2) then
                  idx_e = [i+ii, j]
                end if

                call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx_e, .true.)
                call this%surf(idx_dir)%update(idx_e, 0.5*surf)
                call this%eps(IDX_EDGE,idx_dir)%update(idx_e, 0.5*surf*this%reg_poiss(si)%eps)
              end do

              ! faces
              do ii = 0, 1
                idx_f = [i,j]
                idx_f(idx_dir) = idx_f(idx_dir) + ii

                call this%poisson(IDX_FACE,idx_dir)%flags%set(idx_f, .true.)
              end do
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do; end do

        case ("xyz")
          ! get bounds
          do idx_dir = 1, 3
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_poiss(si)%xyz(idx_dir,2)) - 1
          end do

          ! enable poisson in region
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j, k]

            ! set permittivity for cell
            call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

            ! enable poisson for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do kk = 0, 1;do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj, k+kk]
              call this%poisson(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%vol%update(idx_v, 0.125*vol)
            end do; end do; end do

            do idx_dir = 1, 3
              surf = this%g%get_surf(idx_c, idx_dir)

              ! edges
              do jj = 0, 1; do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii, k+jj]
                elseif (idx_dir == 2) then
                  idx_e = [i+ii, j, k+jj]
                else
                  idx_e = [i+ii, j+jj, k]
                end if

                call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx_e, .true.)
                call this%surf(idx_dir)%update(idx_e, 0.25*surf)
                call this%eps(IDX_EDGE,idx_dir)%update(idx_e, 0.25*surf*this%reg_poiss(si)%eps)
              end do; end do

              ! faces
              do ii = 0, 1
                idx_f = [i, j, k]
                idx_f(idx_dir) = idx_f(idx_dir) + ii

                call this%poisson(IDX_FACE,idx_dir)%flags%set(idx_f, .true.)
              end do
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do; end do; end do

        case("tr_xy")
          do i = 1, this%gtr%ncell
            idx_c = i
            call this%gtr%get_cell(idx_c, p)

            mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
            mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
            if (.not. this%reg_poiss(si)%point_test(this%gtype, mid)) cycle

            call adjoint_triangle(idx_c, edge, trsurf)

            ! set permittivity for cell
            call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

            do ii = 1, 3
              idx_e  = this%gtr%cell2edge(ii,i)
              idx_v1 = this%gtr%edge2vert(1,idx_e(1))
              idx_v2 = this%gtr%edge2vert(2,idx_e(1))

              vol = 0.25 * edge(ii) * trsurf(ii)

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
              call this%eps(IDX_EDGE,1)%update(idx_e, trsurf(ii)*this%reg_poiss(si)%eps)

              idx_f = idx_e
              call this%poisson(IDX_FACE,1)%flags%set(idx_f, .true.)
            end do

            ! enable poisson for cell
            call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do

        case("tr_xyz")
          ! get bounds
            i0 = bin_search(this%g1D(3)%x, this%reg_poiss(si)%xyz(3,1))
            i1 = bin_search(this%g1D(3)%x, this%reg_poiss(si)%xyz(3,2)) - 1

            ! enable poisson in region
            do i = 1, this%gtr%ncell
              call this%gtr%get_cell([i], p)

              mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
              mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
              if (.not. this%reg_poiss(si)%point_test(new_string("tr_xy"), mid)) cycle

              call adjoint_triangle([i], edge, trsurf)
              ! ! get edge lengths
              ! do ii = 1, 3
              !   idx_e(1)  = this%gtr%cell2edge(ii,i)
              !   idx_v1(1) = this%gtr%edge2vert(1,idx_e(1))
              !   idx_v2(1) = this%gtr%edge2vert(2,idx_e(1))

              !   call this%gtr%get_vertex(idx_v1(1:1), p(:,1))
              !   call this%gtr%get_vertex(idx_v2(1:1), p(:,2))

              !   edge(ii) = sqrt((p(1,2) - p(1,1))**2 + (p(2,2) - p(2,1))**2)
              ! end do

              ! ! get radius of circumscribed circle
              ! R = (edge(2) * edge(1)*edge(3)) / sqrt(  (edge(1)+edge(2)+edge(3))*(edge(1)+edge(2)-edge(3))*(edge(1)-edge(2)+edge(3))*(-edge(1)+edge(2)+edge(3) ) )

              ! ! adjoint surface parts
              ! do ii = 1, 3
              !   trsurf(ii) = sqrt(R**2 - (edge(ii)/2)**2)
              ! end do

              do k = i0(1), i1(1)
                idx_c = [i, k]

                ! z edge length
                len = this%g1D(3)%get_len([k], 1)

                ! set permittivity for cell
                call this%eps(IDX_CELL,0)%set(idx_c, this%reg_poiss(si)%eps)

                ! enable poisson for vertices and update adjoint volumes
                do ii = 1, 3
                  ! edges in triangle
                  idx_e(1) = this%gtr%cell2edge(ii,i)
                  vol = 0.125 * edge(ii) * trsurf(ii) * len
                  do kk = 0, 1
                    idx_v1 = [this%gtr%edge2vert(1,idx_e(1)), k+kk]
                    idx_v2 = [this%gtr%edge2vert(2,idx_e(1)), k+kk]
                    call this%poisson(IDX_VERTEX,0)%flags%set(idx_v1, .true.)
                    call this%poisson(IDX_VERTEX,0)%flags%set(idx_v2, .true.)
                    call this%vol%update(idx_v1, vol)
                    call this%vol%update(idx_v2, vol)

                    idx_e(2) = k + kk
                    call this%poisson(IDX_EDGE,1)%flags%set(idx_e, .true.)
                    if (trsurf(ii) == 0) then
                      call this%surf(1)%update(idx_e, 0.5 * len * 1e-16)
                    else
                      call this%surf(1)%update(idx_e, 0.5 * len * trsurf(ii))
                    end if
                    call this%eps(IDX_EDGE,1)%update(idx_e, 0.5 * len * trsurf(ii) * this%reg_poiss(si)%eps)
                  end do

                  ! faces perpendicular to triangle
                  idx_f = [idx_e(1), k]
                  call this%poisson(IDX_FACE,1)%flags%set(idx_f, .true.)

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
                end do
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

    subroutine init_transport()
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
      do idx_type = 1, 4
        table_name0 = "oxide_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) then
            table_name = table_name0//idx_dir_name(idx_dir)
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
        case ("x")
          ! get bounds
          i0(1) = bin_search(this%g1D(1)%x, this%reg_trans(si)%xyz(1,1))
          i1(1) = bin_search(this%g1D(1)%x, this%reg_trans(si)%xyz(1,2)) - 1

          ! enable transport in region
          do i = i0(1), i1(1)
            idx_c(1) = i

            ! enable transport for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do ii = 0, 1
              idx_v(1) = i+ii
              call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v, .false.)
              call this%transport(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%tr_vol%update(idx_v, 0.5*vol)
            end do

            surf = this%g%get_surf(idx_c, 1)


            ! edges
            idx_e = idx_c
            ! enable transport for edges
            call this%oxide(    IDX_EDGE,1)%flags%set(idx_e, .false.)
            call this%transport(IDX_EDGE,1)%flags%set(idx_e, .true. )
            call this%tr_surf(1)%update(idx_e, surf)

            ! faces
            do ii = 0, 1
              idx_f(1) = i+ii

              call this%oxide(    IDX_FACE,1)%flags%set(idx_f, .false.)
              call this%transport(IDX_FACE,1)%flags%set(idx_f, .true. )
            end do

            ! enable transport for cell
            call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
            call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do

        case ("xy")
          do idx_dir = 1, 2
            ! get bounds
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 2)) - 1
          end do

          ! enable transport in region
          do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j]

            ! enable transport for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj]

              call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v, .false.)
              call this%transport(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%tr_vol%update(idx_v, 0.25*vol)
            end do; end do

            do idx_dir = 1, 2
              surf = this%g%get_surf(idx_c, idx_dir)
              ! edges
              do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii]
                elseif (idx_dir == 2) then
                  idx_e = [i+ii, j]
                end if

                ! enable transport for edges
                call this%oxide(    IDX_EDGE,idx_dir)%flags%set(idx_e, .false.)
                call this%transport(IDX_EDGE,idx_dir)%flags%set(idx_e, .true. )
                call this%tr_surf(idx_dir)%update(idx_e, 0.5*surf)
              end do

              ! faces
              do ii = 0, 1
                idx_f = [i, j]
                idx_f(idx_dir) = idx_f(idx_dir) + ii

                call this%oxide(    IDX_FACE,idx_dir)%flags%set(idx_f, .false.)
                call this%transport(IDX_FACE,idx_dir)%flags%set(idx_f, .true. )
              end do
            end do

            ! enable transport for cell
            call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
            call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do; end do

        case ("xyz")
          do idx_dir = 1, 3
            ! get bounds
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_trans(si)%xyz(idx_dir, 2)) - 1
          end do

          ! enable transport in region
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j, k]

            ! enable transport for vertices and update adjoint volumes
            vol = this%g%get_vol(idx_c)
            do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj, k+kk]
              call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v, .false.)
              call this%transport(IDX_VERTEX,0)%flags%set(idx_v, .true.)
              call this%tr_vol%update(idx_v, 0.125*vol)
            end do; end do; end do

            do idx_dir = 1, 3
              surf = this%g%get_surf(idx_c, idx_dir)
              ! edges
              do jj = 0, 1; do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii, k+jj]
                elseif (idx_dir == 2) then
                  idx_e = [i+ii, j, k+jj]
                else
                  idx_e = [i+ii, j+jj, k]
                end if

                ! enable transport for edges
                call this%oxide(    IDX_EDGE,idx_dir)%flags%set(idx_e, .false.)
                call this%transport(IDX_EDGE,idx_dir)%flags%set(idx_e, .true. )
                call this%tr_surf(idx_dir)%update(idx_e, 0.25*surf)
              end do; end do

              ! faces
              do ii = 0, 1
                idx_f = [i, j, k]
                idx_f(idx_dir) = idx_f(idx_dir) + ii

                call this%oxide(    IDX_FACE,idx_dir)%flags%set(idx_f, .false.)
                call this%transport(IDX_FACE,idx_dir)%flags%set(idx_f, .true. )
              end do
            end do

            ! enable transport for cell
            call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
            call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
          end do; end do; end do

        case ("tr_xy")
          ! enable transport in region
          do i = 1, this%gtr%ncell
            idx_c = i
            call this%gtr%get_cell(idx_c, p)

            mid(1) = (p(1, 1) + p(1, 2) + p(1, 3))/3
            mid(2) = (p(2, 1) + p(2, 2) + p(2, 3))/3
            if (.not. this%reg_trans(si)%point_test(this%gtype, mid)) cycle

            call adjoint_triangle(idx_c, edge, trsurf)

            ! enable transport for vertices and update adjoint volumes
            do ii = 1, 3
              idx_e = this%gtr%cell2edge(ii,i)
              idx_v1 = this%gtr%edge2vert(1,idx_e)
              idx_v2 = this%gtr%edge2vert(2,idx_e)

              vol = 0.25 * edge(ii) *trsurf(ii)

              call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v1, .false.)
              call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v2, .false.)
              call this%transport(IDX_VERTEX,0)%flags%set(idx_v1, .true.)
              call this%transport(IDX_VERTEX,0)%flags%set(idx_v2, .true.)
              call this%tr_vol%update(idx_v1, vol)
              call this%tr_vol%update(idx_v2, vol)

              ! enable transport for edges
              call this%oxide(    IDX_EDGE,1)%flags%set(idx_e, .false.)
              call this%transport(IDX_EDGE,1)%flags%set(idx_e, .true. )
              call this%tr_surf(1)%update(idx_e, trsurf(ii))

              ! faces
              idx_f = idx_e
              call this%oxide(    IDX_FACE,1)%flags%set(idx_f, .false.)
              call this%transport(IDX_FACE,1)%flags%set(idx_f, .true. )
            end do

            ! enable transport for cell
            call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
            call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)

          end do
        case ("tr_xyz")
          ! get bounds
          i0 = bin_search(this%g1D(3)%x, this%reg_trans(si)%xyz(3,1))
          i1 = bin_search(this%g1D(3)%x, this%reg_trans(si)%xyz(3,2)) - 1

          ! enable transport in region
          do i = 1, this%gtr%ncell
            call this%gtr%get_cell([i], p)

            mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
            mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
            if (.not. this%reg_trans(si)%point_test(new_string("tr_xy"), mid)) cycle

            call adjoint_triangle([i], edge, trsurf)

            do k = i0(1), i1(1)
              idx_c = [i, k]

              ! z edge length
              len = this%g1D(3)%get_len([k], 1)

              ! enable transport for vertices and update adjoint volumes
              do ii = 1, 3
                ! edges in triangle
                idx_e(1) = this%gtr%cell2edge(ii,i)
                vol = 0.125 * edge(ii) * trsurf(ii) * len
                do kk = 0, 1
                  idx_v1 = [this%gtr%edge2vert(1,idx_e(1)), k+kk]
                  idx_v2 = [this%gtr%edge2vert(2,idx_e(1)), k+kk]
                  call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v1, .false.)
                  call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v2, .false.)
                  call this%transport(IDX_VERTEX,0)%flags%set(idx_v1, .true.)
                  call this%transport(IDX_VERTEX,0)%flags%set(idx_v2, .true.)
                  call this%tr_vol%update(idx_v1, vol)
                  call this%tr_vol%update(idx_v2, vol)

                  ! enable transport for edges in triangle
                  idx_e(2) = k + kk
                  call this%oxide(    IDX_EDGE,1)%flags%set(idx_e, .false.)
                  call this%transport(IDX_EDGE,1)%flags%set(idx_e, .true. )
                  call this%tr_surf(1)%update(idx_e, 0.5 * len *trsurf(ii))
                end do

                ! faces perpendicular to triangle
                idx_f = [idx_e(1), k]
                call this%oxide(    IDX_FACE,1)%flags%set(idx_f, .false.)
                call this%transport(IDX_FACE,1)%flags%set(idx_f, .true. )

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
              end do
            end do
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

    subroutine init_doping()
      real :: mob_min, mob_max, N_ref, alpha, tr_surf, dop(2), acon, dcon, tr_vol(2)

      ! allocate grid data
      call allocate_grid_data3_real(this%dop, this%g%idx_dim, [1, 0, this%ci0], [4, this%g%idx_dim, this%ci1])
      call allocate_grid_data3_real(this%mob0, this%g%idx_dim, [1, 0, this%ci0], [4, this%g%idx_dim, this%ci1])
      do ci = this%ci0, this%ci1
        do idx_type = 1, 4
          do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
            call this%dop(idx_type,idx_dir,ci)%init(this%g, idx_type, idx_dir)
          end do
        end do
      end do

      do si = 1, size(this%reg_dop)
        if (this%ci0 == CR_ELEC) dop(CR_ELEC) = this%reg_dop(si)%dcon
        if (this%ci1 == CR_HOLE) dop(CR_HOLE) = this%reg_dop(si)%acon

        select case(this%gtype%s)
        case("x")
          i0(1) = bin_search(this%g1D(1)%x, this%reg_dop(si)%xyz(1,1))
          i1(1) = bin_search(this%g1D(1)%x, this%reg_dop(si)%xyz(1,2)) - 1

          ! set doping in region
          do i = i0(1), i1(1)
            idx_c(1) = i
            vol = this%g%get_vol(idx_c)

            ! set doping on vertices
            do ii = 0, 1
              idx_v  = i+ii
              tr_vol(1) = this%tr_vol%get(idx_v)
              do ci = this%ci0, this%ci1
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v, 0.5*vol/tr_vol(1)*dop(ci))
              end do
            end do

            ! set doping for edges
            surf    = this%g%get_surf(idx_c, 1)

            idx_e = idx_c
            tr_surf = this%tr_surf(1)%get(idx_e)
            do ci = this%ci0, this%ci1
              call this%dop(IDX_EDGE,1,ci)%update(idx_e, surf/tr_surf*dop(ci))
            end do

            ! set doping in cells
            do ci = this%ci0, this%ci1
              call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
            end do
          end do

        case("xy")
          do idx_dir = 1, 2
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 2)) - 1
          end do

          ! set doping in region
          do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j]
            vol = this%g%get_vol(idx_c)

            ! set doping on vertices
            do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj]
              tr_vol(1) = this%tr_vol%get(idx_v)
              do ci = this%ci0, this%ci1
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v, 0.25*vol/tr_vol(1)*dop(ci))
              end do
            end do; end do

            ! set doping for edges
            do idx_dir = 1, 2
              surf    = this%g%get_surf(idx_c, idx_dir)

              do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii]
                elseif (idx_dir == 2) then
                  idx_e = [i+ii, j]
                end if
                tr_surf = this%tr_surf(idx_dir)%get(idx_e)
                do ci = this%ci0, this%ci1
                  call this%dop(IDX_EDGE,idx_dir,ci)%update(idx_e, 0.5*surf/tr_surf*dop(ci))
                end do
              end do
            end do

            ! set doping in cells
            do ci = this%ci0, this%ci1
              call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
            end do
          end do; end do

        case("xyz")
          do idx_dir = 1, 3
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_dop(si)%xyz(idx_dir, 2)) - 1
          end do

          ! set doping in region
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_c = [i, j, k]
            vol = this%g%get_vol(idx_c)

            ! set doping on vertices
            do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
              idx_v = [i+ii, j+jj, k+kk]
              tr_vol(1) = this%tr_vol%get(idx_v)
              do ci = this%ci0, this%ci1
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v, 0.125*vol/tr_vol(1)*dop(ci))
              end do
            end do; end do; end do

            ! set doping for edges
            do idx_dir = 1, 3
              surf    = this%g%get_surf(idx_c, idx_dir)

              do jj = 0, 1; do ii = 0, 1
                if (idx_dir == 1) then
                  idx_e = [i, j+ii, k+jj]
                  elseif (idx_dir == 2) then
                    idx_e = [i+ii, j, k+jj]
                  else
                    idx_e = [i+ii, j+jj, k]
                  end if
                  tr_surf = this%tr_surf(idx_dir)%get(idx_c)
                do ci = this%ci0, this%ci1
                  call this%dop(IDX_EDGE,idx_dir,ci)%update(idx_e, 0.25*surf/tr_surf*dop(ci))
                end do
              end do; end do
            end do

            ! set doping in cells
            do ci = this%ci0, this%ci1
              call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
            end do
          end do; end do; end do

        case("tr_xy")
          do i = 1, this%gtr%ncell
            idx_c = i
            call this%gtr%get_cell(idx_c, p)

            mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
            mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
            if (.not. this%reg_dop(si)%point_test(this%gtype, mid)) cycle

            call adjoint_triangle(idx_c, edge, trsurf)

            ! set doping in region
            do ii = 1, 3
              idx_e  = this%gtr%cell2edge(ii,i)
              idx_v1 = this%gtr%edge2vert(1,idx_e(1))
              idx_v2 = this%gtr%edge2vert(2,idx_e(1))
              ! idx_v = this%gtr%cell2vert(ii,i)
              ! set doping on vertices
              vol = 0.25 * edge(ii) * trsurf(ii)
              tr_vol(1) = this%tr_vol%get(idx_v1)
              tr_vol(2) = this%tr_vol%get(idx_v2)

              do ci = this%ci0, this%ci1
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v1, vol/tr_vol(1)*dop(ci))
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v2, vol/tr_vol(2)*dop(ci))
              end do

              ! set doping for edges
              tr_surf = this%tr_surf(1)%get(idx_e)

              do ci = this%ci0, this%ci1
                if (tr_surf == 0) then
                  call this%dop(IDX_EDGE,1,ci)%update(idx_e, trsurf(ii)/1e-16*dop(ci))
                else
                  call this%dop(IDX_EDGE,1,ci)%update(idx_e, trsurf(ii)/tr_surf*dop(ci))
                end if
              end do

            end do

            ! set doping in cells
            do ci = this%ci0, this%ci1
              call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
            end do
          end do

        case("tr_xyz")
          ! get bounds
          i0 = bin_search(this%g1D(3)%x, this%reg_dop(si)%xyz(3,1))
          i1 = bin_search(this%g1D(3)%x, this%reg_dop(si)%xyz(3,2)) - 1

          ! set doping in region
          do i = 1, this%gtr%ncell
            call this%gtr%get_cell([i], p)
            mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
            mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
            if (.not. this%reg_dop(si)%point_test(new_string("tr_xy"), mid)) cycle

            call adjoint_triangle([i], edge, trsurf)

            do k = i0(1), i1(1)
              idx_c = [i, k]

              ! z edge length
              len = this%g1D(3)%get_len([k], 1)

              do ii = 1, 3
                idx_e(1) = this%gtr%cell2edge(ii,i)
                vol = 0.125 * edge(ii) * trsurf(ii) * len
                ! set doping on vertices
                do kk = 0, 1
                  idx_v1 = [this%gtr%edge2vert(1,idx_e(1)), k+kk]
                  idx_v2 = [this%gtr%edge2vert(2,idx_e(1)), k+kk]
                  tr_vol(1) = this%tr_vol%get(idx_v1)
                  tr_vol(2) = this%tr_vol%get(idx_v2)
                  do ci = this%ci0, this%ci1
                    call this%dop(IDX_VERTEX,0,ci)%update(idx_v1, vol/tr_vol(1)*dop(ci))
                    call this%dop(IDX_VERTEX,0,ci)%update(idx_v2, vol/tr_vol(2)*dop(ci))
                  end do

                  ! set doping for edges in triangles
                  idx_e(2) = k + kk
                  tr_surf = this%tr_surf(1)%get(idx_e)
                  do ci = this%ci0, this%ci1
                    if (tr_surf == 0) then
                      call this%dop(IDX_EDGE,1,ci)%update(idx_e, 0.5 * len * trsurf(ii)/1e-16*dop(ci))
                    else
                      call this%dop(IDX_EDGE,1,ci)%update(idx_e, 0.5 * len * trsurf(ii)/tr_surf*dop(ci))
                    end if
                  end do
                end do

                ! edges perpendicular to triangle through first vertex
                idx_e = [this%gtr%edge2vert(1,this%gtr%cell2edge(ii,i)), k]
                tr_surf = this%tr_surf(2)%get(idx_e)
                do ci = this%ci0, this%ci1
                  if (tr_surf == 0) then
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/1e-16*dop(ci))
                  else
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/tr_surf*dop(ci))
                  end if
                end do

                ! edges perpendicular to triangle through second vertex
                idx_e = [this%gtr%edge2vert(2,this%gtr%cell2edge(ii,i)), k]
                tr_surf = this%tr_surf(2)%get(idx_e)
                do ci = this%ci0, this%ci1
                  if (tr_surf == 0) then
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/1e-16*dop(ci))
                  else
                    call this%dop(IDX_EDGE,2,ci)%update(idx_e, 0.25 * edge(ii) * trsurf(ii)/tr_surf*dop(ci))
                  end if
                end do
              end do
            end do

            ! set doping in cells
            do ci = this%ci0, this%ci1
              call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))
            end do
          end do
        end select
      end do

      ! init zero-field mobility
      acon = 0
      dcon = 0
      do ci = this%ci0, this%ci1
        mob_min = this%mob_min(ci)
        mob_max = this%mob_max(ci)
        N_ref   = this%N_ref(ci)
        alpha   = this%alpha(ci)

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

    subroutine init_contacts()
      integer                           :: ict, ict0, nct, ct_type
      type(mapnode_string_int), pointer :: node
      type(string)                      :: name

      ! allocate grid data
      call allocate_grid_data0_int(this%ict, this%g%idx_dim)
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
        case("x")
          ! bounds
          i0(1) = bin_search(this%g1D(1)%x, this%reg_ct(si)%xyz(1,1))
          i1(1) = bin_search(this%g1D(1)%x, this%reg_ct(si)%xyz(1,2))

          ! phims
          if (ct_type == CT_OHMIC) then
            outer: do i = i0(1), i1(1)
              idx_v = i
              if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit outer
            end do outer
            if (i > i1(1)) call program_error("Ohmic contact "//name%s//" not in transport region")
            if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
              this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
            elseif (this%ci0 == CR_ELEC) then
              this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
            elseif (this%ci0 == CR_HOLE) then
              this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
            end if
          else
            this%contacts(ict)%phims = this%reg_ct(si)%phims
          end if

          ! update vertex tables
          do i = i0(1), i1(1)
            idx_v = i
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

        case("xy")
          ! bounds
          i0 = 1
          i1 = 1
          do idx_dir = 1, 2
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 2))
          end do

          ! phims
          if (ct_type == CT_OHMIC) then
            oute: do j = i0(2), i1(2); do i = i0(1), i1(1)
              idx_v = [i, j]
              if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit oute
            end do; end do oute
            if ((i > i1(1)) .or. (j > i1(2))) call program_error("Ohmic contact "//name%s//" not in transport region")
            if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
              this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
            elseif (this%ci0 == CR_ELEC) then
              this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
            elseif (this%ci0 == CR_HOLE) then
              this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
            end if
          else
            this%contacts(ict)%phims = this%reg_ct(si)%phims
          end if

          ! update vertex tables
          do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_v = [i, j]
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
          end do; end do

        case("xyz")
          ! bounds
          i0 = 1
          i1 = 1
          do idx_dir = 1, 3
            i0(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 1))
            i1(idx_dir) = bin_search(this%g1D(idx_dir)%x, this%reg_ct(si)%xyz(idx_dir, 2))
          end do

          ! phims
          if (ct_type == CT_OHMIC) then
            out: do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
              idx_v = [i, j, k]
              if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit out
            end do; end do; end do out
            if ((i > i1(1)) .or. (j > i1(2)).or. (k > i1(3))) call program_error("Ohmic contact "//name%s//" not in transport region")
            if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
              this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
            elseif (this%ci0 == CR_ELEC) then
              this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
            elseif (this%ci0 == CR_HOLE) then
              this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
            end if
          else
            this%contacts(ict)%phims = this%reg_ct(si)%phims
          end if

          ! update vertex tables
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            idx_v = [i, j, k]
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

        case("tr_xy")
          do i = 1, this%gtr%nvert
            idx_v = i
            call this%gtr%get_vertex(idx_v, p(:,1))

            if (.not. this%reg_ct(si)%point_test(this%gtype, p(:,1))) cycle

            ! phims
            if (ct_type == CT_OHMIC) then
             ! if ((i > i1(1)) .or. (j > i1(2))) call program_error("Ohmic contact "//name%s//" not in transport region")
              if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
                this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
              elseif (this%ci0 == CR_ELEC) then
                this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
              elseif (this%ci0 == CR_HOLE) then
                this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
              end if
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

         case("tr_xyz")
          ! get bounds
          i0(3) = bin_search(this%g1D(3)%x, this%reg_ct(si)%xyz(3,1))
          i1(3) = bin_search(this%g1D(3)%x, this%reg_ct(si)%xyz(3,2)) - 1

          do i = 1, this%gtr%nvert
            call this%gtr%get_vertex([i], p(:,1))
            if (.not. this%reg_ct(si)%point_test(new_string("tr_xy"), p(:,1))) cycle

            do k = i0(3), i1(3)
              idx_v = [i, k]

              ! phims
              if (ct_type == CT_OHMIC) then
                ! if ((i > i1(1)) .or. (j > i1(2))) call program_error("Ohmic contact "//name%s//" not in transport region")
                if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
                  this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
                elseif (this%ci0 == CR_ELEC) then
                  this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
                elseif (this%ci0 == CR_HOLE) then
                  this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
                end if
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

            ! call this%gtr%get_cell([i], p)

            ! mid(1) = (p(1,1) + p(1,2) + p(1,3))/3
            ! mid(2) = (p(2,1) + p(2,2) + p(2,3))/3
            ! if (.not. this%reg_ct(si)%point_test(new_string("tr_xy"), mid)) cycle

            ! ! phims
            ! if (ct_type == CT_OHMIC) then
            !   outer_trxyz: do ii = 1, 3; do k = i0(1), i1(1)
            !     idx_v = [this%gtr%cell2vert(ii,i), k]
            !     if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit outer_trxyz
            !   end do; end do outer_trxyz
            !   ! if ((i > i1(1)) .or. (j > i1(2))) call program_error("Ohmic contact "//name%s//" not in transport region")
            !   if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
            !     this%contacts(ict)%phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
            !   elseif (this%ci0 == CR_ELEC) then
            !     this%contacts(ict)%phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
            !   elseif (this%ci0 == CR_HOLE) then
            !     this%contacts(ict)%phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
            !   end if
            ! else
            !   this%contacts(ict)%phims = this%reg_ct(si)%phims
            ! end if

            ! ! update vertex tables
            ! do ii = 1, 3; do k = i0(1), i1(1)
            !   idx_v = [this%gtr%cell2vert(ii,i), k]
            !   ict0 = this%ict%get(idx_v)
            !   if ((ict0 /= 0) .and. (ict0 /= ict)) then
            !     call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
            !   end if
            !   call this%ict%set(idx_v, ict)
            !   if (this%poisson(IDX_VERTEX,0)%flags%get(idx_v)) then
            !     call this%poisson_vct(ict)%flags%set(idx_v, .true. )
            !     call this%poisson_vct(  0)%flags%set(idx_v, .false.)
            !   end if
            !   if (this%oxide(IDX_VERTEX,0)%flags%get(idx_v)) then
            !     call this%oxide_vct(ict)%flags%set(idx_v, .true. )
            !     call this%oxide_vct(  0)%flags%set(idx_v, .false.)
            !   end if
            !   if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) then
            !     call this%transport_vct(ict)%flags%set(idx_v, .true. )
            !     call this%transport_vct(  0)%flags%set(idx_v, .false.)
            !   end if
            ! end do; end do

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

    subroutine adjoint_triangle(indx_c, len, surf)
      integer, intent(in)  :: indx_c(:)
        !! triangle cell indices
      real,    intent(out) :: len(:)
        !! output: 3 triangle egdes lengths
      real,    intent(out) :: surf(:)
        !! output: adjoint surfaces per edge (3)

      real :: R
        !! circumscribed radius
      real :: p(2,3)
        !! vertex coordinates, size = (dim=2, 3)

      call this%gtr%get_cell(indx_c, p)

      ! get edge lengths
      do ii = 1, 3
        idx_e  = this%gtr%cell2edge(ii,i)
        idx_v1 = this%gtr%edge2vert(1,idx_e(1))
        idx_v2 = this%gtr%edge2vert(2,idx_e(1))

        call this%gtr%get_vertex(idx_v1(1:1), p(:,1))
        call this%gtr%get_vertex(idx_v2(1:1), p(:,2))

        len(ii) = sqrt((p(1,2) - p(1,1))**2 + (p(2,2) - p(2,1))**2)
      end do

      ! circumradius
      R = (len(2) * len(1)*len(3)) / sqrt(  (len(1)+len(2)+len(3))*(len(1)+len(2)-len(3))*(len(1)-len(2)+len(3))*(-len(1)+len(2)+len(3) ) )

      ! adjoint surface parts
      do ii = 1, 3
        if (R <= len(ii)*0.5) then
          surf(ii) = 0
        else
          surf(ii) = sqrt(R**2 - (len(ii)*0.5)**2)
        end if
      end do

    end subroutine

  end subroutine

  subroutine device_params_destruct(this)
    !! destruct device parameter object
    class(device_params), intent(inout) :: this

    call this%contact_map%destruct()
  end subroutine
end module
