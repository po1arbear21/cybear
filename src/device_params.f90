m4_include(util/macro.f90.inc)

module device_params_m

  use bin_search_m,  only: bin_search
  use contact_m,     only: CT_OHMIC, CT_GATE, contact
  use grid_m,        only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid, grid_ptr
  use grid_data_m,   only: allocate_grid_data0_real, allocate_grid_data1_real, allocate_grid_data2_real, allocate_grid_data3_real, &
    &                      allocate_grid_data0_int, allocate_grid_data1_int, grid_data_int, grid_data_real, grid_data2_int, grid_data2_real
  use grid_table_m,  only: grid_table
  use grid1D_m,      only: grid1D
  use error_m,       only: assert_failed, program_error
  use input_m,       only: input_file
  use math_m,        only: linspace
  use map_m,         only: map_string_int, mapnode_string_int
  use string_m,      only: string
  use tensor_grid_m, only: tensor_grid

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

    type(grid1D)      :: g1D(3)
      !! x, y, z grids
    integer           :: ng1D
      !! number of enabled 1D grids
    integer           :: ig1D(3)
      !! grid index translation table
    type(tensor_grid) :: g
      !! tensor grid

    class(grid_data_real), allocatable :: surf(:)
      !! adjoint volume surfaces
    class(grid_data_real), allocatable :: tr_surf(:)
      !! adjoint volume surfaces in transport region
    class(grid_data_real), allocatable :: vol
      !! adjoint volumes
    class(grid_data_real), allocatable :: tr_vol
      !! adjoint volumes for transport region

    class(grid_data_real), allocatable  :: eps(:,:)
      !! electrostatic permittivity (idx_type, idx_dir)

    type(grid_table), allocatable :: poisson(:,:)
      !! poisson grid tables (idx_type, idx_dir)
    type(grid_table), allocatable :: oxide(:,:)
      !! oxide grid tables (idx_type, idx_dir)
    type(grid_table), allocatable :: transport(:,:)
      !! transport grid tables (idx_type, idx_dir)

    class(grid_data_real), allocatable :: dop(:,:,:)
      !! acceptor/donator concentration (idx_type, idx_dir, carrier index)

    class(grid_data_real), allocatable :: mob0(:,:,:)
      !! zero-field mobility (idx_type, idx_dir, carrier index)

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
  contains
    procedure :: init     => device_params_init
    procedure :: destruct => device_params_destruct
  end type

contains

  subroutine device_params_init(this, file)
    !! initialize device parameters using device input file
    class(device_params), target, intent(out) :: this
    type(input_file),             intent(in)  :: file
      !! device input file

    integer                   :: i, i0(3), i1(3), ii, ijk(3), j, jj, k, kk, si, ci, n(3), dir, dir2
    integer                   :: idx_type, idx_dir, idx_dir0(4), idx_dir1(4), idx_dir2
    integer,      allocatable :: sids(:), idx_v(:), idx_e(:), idx_f(:), idx_c(:), idx(:)
    real                      :: vol, surf, xyz_bounds(2,3)
    real,         allocatable :: bounds(:)
    character(:), allocatable :: table_name0, table_name

    call init_transport_params()
    call init_grid()

    idx_dir0 = [ 0, 1, 1, 0]
    idx_dir1 = [ 0, this%g%idx_dim, this%g%idx_dim, 0]
    allocate (idx_v(this%g%idx_dim), idx_e(this%g%idx_dim), idx_f(this%g%idx_dim), idx_c(this%g%idx_dim))

    call init_poisson()
    call init_transport()
    call init_doping()
    call init_contacts()

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

    subroutine init_grid()
      integer               :: si_gen, si_load, status_gen, status_load, dir
      real, allocatable     :: xyz(:)
      type(grid_ptr)        :: gptr(3)
      logical               :: status

      ! find generate/load grid section ids
      call file%get_section("generate grid", si_gen, status = status_gen)
      call file%get_section("load grid", si_load, status = status_load)

      ! error checking
      if ((status_gen >= 0) .and. (status_load >= 0)) call program_error("found both generate grid and load grid sections")
      if ((status_gen <  0) .and. (status_load <  0)) call program_error("found neither generate nor load grid section")
      if (status_gen  > 0) call program_error("found multiple generate grid sections")
      if (status_load > 0) call program_error("found multiple load grid sections")

      this%ng1D = 0
      do dir = 1, 3
        if (status_gen == 0) then
          ! load grid parameters
          call file%get(si_gen, "n"//DIR_NAME(dir), n(dir), status = status)
          if (.not. status) cycle
          call file%get(si_gen, DIR_NAME(dir)//"bounds", bounds)

          ! generate x, y values
          if (allocated(xyz)) deallocate(xyz)
          xyz = linspace(bounds(1), bounds(2), n(dir))
        elseif (status_load == 0) then

          ! load grid values directly
          call file%get(si_load, DIR_NAME(dir), xyz, status = status)
          if (.not. status) cycle
        end if

        this%ng1D = this%ng1D + 1
        this%ig1D(this%ng1D) = dir

        ! create grid
        call this%g1D(dir)%init(DIR_NAME(dir), xyz)
        gptr(this%ng1D) = this%g1D(dir)%get_ptr()
      end do

      call this%g%init("grid", gptr(1:this%ng1D))
    end subroutine

    subroutine init_poisson()
      real :: eps

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
            dir = this%ig1D(idx_dir)
            table_name = table_name0//DIR_NAME(dir)
          else
            table_name = table_name0
          end if
          call this%poisson(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        end do
      end do

      ! get poisson sections
      call file%get_sections("poisson", sids)

      ! process all poisson sections
      do si = 1, size(sids)
        ! get bounds
        i0 = 1
        i1 = 1
        do idx_dir = 1, this%g%idx_dim
          dir = this%ig1D(idx_dir)
          call file%get(sids(si), DIR_NAME(dir)//"bounds", bounds)
          xyz_bounds(:,dir) = bounds
          i0(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(1,dir))
          i1(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(2,dir)) - 1
        end do

        ! get permittivity
        call file%get(sids(si), "eps", eps)

        ! enable poisson in region
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            idx_c(idx_dir) = ijk(dir)
          end do

          ! set permittivity for cell
          call this%eps(IDX_CELL,0)%set(idx_c, eps)

          ! enable poisson for vertices and update adjoint volumes
          vol = this%g%get_vol(idx_c)
          do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
            ijk = [i+ii, j+jj, k+kk]
            do idx_dir = 1, this%g%idx_dim
              dir = this%ig1D(idx_dir)
              idx_v(idx_dir) = ijk(dir)
            end do
            call this%poisson(IDX_VERTEX,0)%flags%set(idx_v, .true.)
            call this%vol%update(idx_v, 0.125*vol)
          end do; end do; end do

          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            surf = this%g%get_surf(idx_c, idx_dir)

            ! edges
            do jj = 0, 1; do ii = 0, 1
              if (dir == 1) then
                ijk = [i, j+ii, k+jj]
              elseif (dir == 2) then
                ijk = [i+ii, j, k+jj]
              else
                ijk = [i+ii, j+jj, k]
              end if
              do idx_dir2 = 1, this%g%idx_dim
                dir2 = this%ig1D(idx_dir2)
                idx_e(idx_dir2) = ijk(dir2)
              end do

              call this%poisson(IDX_EDGE,idx_dir)%flags%set(idx_e, .true.)
              call this%surf(idx_dir)%update(idx_e, 0.25*surf)
              call this%eps(IDX_EDGE,idx_dir)%update(idx_e, 0.25*surf*eps)
            end do; end do

            ! faces
            do ii = 0, 1
              ijk = [i, j, k]
              ijk(dir) = ijk(dir) + ii
              do idx_dir2 = 1, this%g%idx_dim
                dir2 = this%ig1D(idx_dir2)
                idx_f(idx_dir2) = ijk(dir2)
              end do

              call this%poisson(IDX_FACE,idx_dir)%flags%set(idx_f, .true.)
            end do
          end do

          ! enable poisson for cell
          call this%poisson(IDX_CELL,0)%flags%set(idx_c, .true.)
        end do; end do; end do
      end do

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
            dir = this%ig1D(idx_dir)
            table_name = table_name0//DIR_NAME(dir)
          else
            table_name = table_name0
          end if
          call this%oxide(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
          this%oxide(idx_type, idx_dir)%flags = this%poisson(idx_type, idx_dir)%flags
        end do
        table_name0 = "transport_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) then
            dir = this%ig1D(idx_dir)
            table_name = table_name0//DIR_NAME(dir)
          else
            table_name = table_name0
          end if
          call this%transport(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        end do
      end do

      ! get transport sections
      call file%get_sections("transport", sids)

      ! process all sections
      do si = 1, size(sids)
        do idx_dir = 1, this%g%idx_dim
          dir = this%ig1D(idx_dir)
          call file%get(sids(si), DIR_NAME(dir)//"bounds", bounds)
          xyz_bounds(:,dir) = bounds
          ! search bounds
          i0(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(1,dir))
          i1(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(2,dir)) - 1
        end do

        ! enable transport in region
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            idx_c(idx_dir) = ijk(dir)
          end do

          ! enable transport for vertices and update adjoint volumes
          vol = this%g%get_vol(idx_c)
          do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
            ijk = [i+ii, j+jj, k+kk]
            do idx_dir = 1, this%g%idx_dim
              dir = this%ig1D(idx_dir)
              idx_v(idx_dir) = ijk(dir)
            end do
            call this%oxide(    IDX_VERTEX,0)%flags%set(idx_v, .false.)
            call this%transport(IDX_VERTEX,0)%flags%set(idx_v, .true.)
            call this%tr_vol%update(idx_v, 0.125*vol)
          end do; end do; end do

          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            surf = this%g%get_surf(idx_c, idx_dir)

            ! edges
            do jj = 0, 1; do ii = 0, 1
              if (dir == 1) then
                ijk = [i, j+ii, k+jj]
              elseif (dir == 2) then
                ijk = [i+ii, j, k+jj]
              else
                ijk = [i+ii, j+jj, k]
              end if
              do idx_dir2 = 1, this%g%idx_dim
                dir2 = this%ig1D(idx_dir2)
                idx_e(idx_dir2) = ijk(dir2)
              end do

              ! enable transport for edges
              call this%oxide(    IDX_EDGE,idx_dir)%flags%set(idx_e, .false.)
              call this%transport(IDX_EDGE,idx_dir)%flags%set(idx_e, .true. )
              call this%tr_surf(idx_dir)%update(idx_e, 0.25*surf)
            end do; end do

            ! faces
            do ii = 0, 1
              ijk = [i, j, k]
              ijk(dir) = ijk(dir) + ii
              do idx_dir2 = 1, this%g%idx_dim
                dir2 = this%ig1D(idx_dir2)
                idx_f(idx_dir2) = ijk(dir2)
              end do

              call this%oxide(    IDX_FACE,idx_dir)%flags%set(idx_f, .false.)
              call this%transport(IDX_FACE,idx_dir)%flags%set(idx_f, .true. )
            end do
          end do

          ! enable transport for cell
          call this%oxide(    IDX_CELL,0)%flags%set(idx_c, .false.)
          call this%transport(IDX_CELL,0)%flags%set(idx_c, .true.)
        end do; end do; end do
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
      integer           :: si_load, status
      real              :: dcon, acon, mob_min, mob_max, N_ref, alpha, tr_vol, tr_surf, dop(2)
      real, allocatable :: dop_ar(:)

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

      ! get load doping and doping region sections
      call file%get_section("load doping", si_load, status = status)
      call file%get_sections("doping", sids)

      ! error checking
      if (status > 0) call program_error("found multiple load doping sections")
      if ((status == 0) .and. (size(sids) > 0)) call program_error("found both load and regular doping sections")

      if (status == 0) then
        ! load doping (on vertices)
        if (this%ci0 == CR_ELEC) then
          call file%get(si_load, "dcon", dop_ar)
          call this%dop(IDX_VERTEX,0,CR_ELEC)%set(dop_ar)
        end if
        if (this%ci1 == CR_HOLE) then
          call file%get(si_load, "acon", dop_ar)
          call this%dop(IDX_VERTEX,0,CR_HOLE)%set(dop_ar)
        end if

        ! set doping for edges and cells
        do ci = this%ci0, this%ci1
          do i = 1, this%transport(IDX_CELL,0)%n
            idx_c = this%transport(IDX_CELL,0)%get_idx(i)

            ! cell doping
            dop(ci) = 0
            do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
              ijk = [ii, jj, kk]
              do idx_dir = 1, this%g%idx_dim
                dir = this%ig1D(idx_dir)
                idx_v(idx_dir) = idx_c(idx_dir) + ijk(dir)
              end do
              dop(ci) = dop(ci) + 0.125 * this%dop(IDX_VERTEX,0,ci)%get(idx_v)
            end do; end do; end do
            call this%dop(IDX_CELL,0,ci)%set(idx_c, dop(ci))

            ! edges
            do idx_dir = 1, this%g%idx_dim
              dir     = this%ig1D(idx_dir)
              surf    = this%g%get_surf(idx_c, idx_dir)
              tr_surf = this%tr_surf(idx_dir)%get(idx_c)
              do jj = 0, 1; do ii = 0, 1
                if (dir == 1) then
                  ijk = [0, ii, jj]
                elseif (dir == 2) then
                  ijk = [ii, 0, jj]
                else
                  ijk = [ii, jj, 0]
                end if
                do idx_dir2 = 1, this%g%idx_dim
                  dir2 = this%ig1D(idx_dir2)
                  idx_e(idx_dir2) = idx_c(idx_dir2) + ijk(dir2)
                end do
                call this%dop(IDX_EDGE,idx_dir,ci)%update(idx_e, 0.25*surf/tr_surf*dop(ci))
              end do; end do
            end do
          end do
        end do
      else

        ! process all doping sections
        do si = 1, size(sids)
          i0 = 1
          i1 = 1
          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            call file%get(sids(si), DIR_NAME(dir)//"bounds", bounds)
            xyz_bounds(:,dir) = bounds
            i0(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(1,dir))
            i1(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(2,dir)) - 1
          end do
          if (this%ci0 == CR_ELEC) call file%get(sids(si), "dcon", dop(CR_ELEC))
          if (this%ci1 == CR_HOLE) call file%get(sids(si), "acon", dop(CR_HOLE))

          ! set doping in region
          do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            do idx_dir = 1, this%g%idx_dim
              dir = this%ig1D(idx_dir)
              idx_c(idx_dir) = ijk(dir)
            end do
            vol = this%g%get_vol(idx_c)

            ! set doping on vertices
            do kk = 0, 1; do jj = 0, 1; do ii = 0, 1
              ijk = [i+ii, j+jj, k+kk]
              do idx_dir = 1, this%g%idx_dim
                dir = this%ig1D(idx_dir)
                idx_v(idx_dir) = ijk(dir)
              end do
              tr_vol = this%tr_vol%get(idx_v)
              do ci = this%ci0, this%ci1
                call this%dop(IDX_VERTEX,0,ci)%update(idx_v, 0.125*vol/tr_vol*dop(ci))
              end do
            end do; end do; end do

            ! set doping for edges
            do idx_dir = 1, this%g%idx_dim
              dir = this%ig1D(idx_dir)
              surf    = this%g%get_surf(idx_c, idx_dir)
              tr_surf = this%tr_surf(idx_dir)%get(idx_c)

              do jj = 0, 1; do ii = 0, 1
                if (dir == 1) then
                  ijk = [i, j+ii, k+jj]
                elseif (dir == 2) then
                  ijk = [i+ii, j, k+jj]
                else
                  ijk = [i+ii, j+jj, k]
                end if
                do idx_dir2 = 1, this%g%idx_dim
                  dir2 = this%ig1D(idx_dir2)
                  idx_e(idx_dir2) = ijk(dir2)
                end do
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
        end do
      end if

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
      real                              :: phims
      type(string)                      :: name, typename
      type(mapnode_string_int), pointer :: node

      ! allocate grid data
      call allocate_grid_data0_int(this%ict, this%g%idx_dim)

      ! get contact sections
      call file%get_sections("contact", sids)

      ! get all contact names
      call this%contact_map%init()
      nct = 0
      do si = 1, size(sids)
        call file%get(sids(si), "name", name)
        node => this%contact_map%find(name)
        if (.not. associated(node)) then
          ! new contact
          nct = nct + 1
          call this%contact_map%insert(name, nct)
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
      do si = 1, size(sids)

        ! lookup name
        call file%get(sids(si), "name", name)
        ict = this%contact_map%get(name)

        ! get type
        call file%get(sids(si), "type", typename)
        if (typename%s == "ohmic") then
          ct_type = CT_OHMIC
        elseif (typename%s == "gate") then
          ct_type = CT_GATE
        else
          call program_error("unknown contact type "//typename%s)
        end if

        ! bounds
        i0 = 1
        i1 = 1
        do idx_dir = 1, this%g%idx_dim
          dir = this%ig1D(idx_dir)
          call file%get(sids(si), DIR_NAME(dir)//"bounds", bounds)
          xyz_bounds(:,dir) = bounds
          i0(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(1,dir))
          i1(dir) = bin_search(this%g1D(dir)%x, xyz_bounds(2,dir))
        end do

        ! phims
        if (ct_type == CT_OHMIC) then
          outer: do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
            ijk = [i, j, k]
            do idx_dir = 1, this%g%idx_dim
              dir = this%ig1D(idx_dir)
              idx_v(idx_dir) = ijk(dir)
            end do
            if (this%transport(IDX_VERTEX,0)%flags%get(idx_v)) exit outer
          end do; end do; end do outer
          if ((i > i1(1)) .or. (j > i1(2)).or. (k > i1(3))) call program_error("Ohmic contact "//name%s//" not in transport region")
          if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
            phims = asinh(0.5 * (this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) - this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v)) / this%n_intrin)
          elseif (this%ci0 == CR_ELEC) then
            phims = log(this%dop(IDX_VERTEX,0,CR_ELEC)%get(idx_v) / this%n_intrin)
          elseif (this%ci0 == CR_HOLE) then
            phims = -log(this%dop(IDX_VERTEX,0,CR_HOLE)%get(idx_v) / this%n_intrin)
          end if
        else
          call file%get(sids(si), "phims", phims)
        end if

        if (.not. allocated(this%contacts(ict)%name)) then
          ! new contact
          this%contacts(ict)%name  = name%s
          this%contacts(ict)%type  = ct_type
          this%contacts(ict)%phims = phims
          call this%poisson_vct(  ict)%init("poisson_VCT_"//name%s, this%g, IDX_VERTEX, 0)
          call this%oxide_vct(    ict)%init("oxide_VCT_"//name%s, this%g, IDX_VERTEX, 0)
          call this%transport_vct(ict)%init("transport_VCT_"//name%s, this%g, IDX_VERTEX, 0)
        elseif (this%contacts(ict)%type /= ct_type) then
          call program_error("contact "//name%s//" given multiple times with different types")
        end if

        ! update vertex tables
        do k = i0(3), i1(3); do j = i0(2), i1(2); do i = i0(1), i1(1)
          ijk = [i, j, k]
          do idx_dir = 1, this%g%idx_dim
            dir = this%ig1D(idx_dir)
            idx_v(idx_dir) = ijk(dir)
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

  end subroutine

  subroutine device_params_destruct(this)
    !! destruct device parameter object
    class(device_params), intent(inout) :: this

    call this%contact_map%destruct()
  end subroutine

end module
