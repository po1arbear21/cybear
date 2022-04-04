#include "util/macro.f90.inc"

module device_params_m

  use bin_search_m,  only: bin_search
  use contact_m,     only: CT_OHMIC, CT_GATE, contact
  use grid_m,        only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME, grid_data2_real, grid_data2_int, grid_table
  use grid1D_m,      only: grid1D
  use error_m,       only: assert_failed, program_error
  use input_m,       only: input_file
  use math_m,        only: linspace
  use map_m,         only: map_string_int, mapnode_string_int
  use string_m,      only: string
  use tensor_grid_m, only: tensor_grid

  implicit none

  ! parameters
  character(*), parameter :: DIR_NAME(2)  = [ "x", "y" ]
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

    type(grid1D)      :: gx, gy
      !! x and y grids
    type(tensor_grid) :: g
      !! 2D tensor grid based on gx and gy

    type(grid_data2_real) :: surf(2)
      !! adjoint volume surfaces
    type(grid_data2_real) :: tr_surf(2)
      !! adjoint volume surfaces in transport region
    type(grid_data2_real) :: vol
      !! adjoint volumes
    type(grid_data2_real) :: tr_vol
      !! adjoint volumes for transport region

    type(grid_data2_real) :: eps
      !! electrostatic permittivity

    type(grid_table) :: poisson(4,0:2)
      !! poisson grid tables (idx_type, idx_dir)
    type(grid_table) :: oxide(4,0:2)
      !! oxide grid tables (idx_type, idx_dir)
    type(grid_table) :: transport(4,0:2)
      !! transport grid tables (idx_type, idx_dir)

    type(grid_data2_real) :: acon(4,0:2), dcon(4,0:2)
      !! acceptor/donator concentration (idx_type, idx_dir)

    type(grid_data2_real) :: mob0(2,2)
      !! zero-field mobility on edges (direction, carrier index)

    type(contact),    allocatable :: contacts(:)
      !! device contacts
    type(map_string_int)          :: contact_map
      !! get contact index by name
    type(grid_data2_int)          :: ict
      !! get contact index by grid indices
    type(grid_table), allocatable :: poisson_vct(:)
      !! poisson vertices grouped by contacts (0:size(contacts))
    type(grid_table), allocatable :: oxide_vct(:)
      !! oxide vertices grouped by contacts (0:size(contacts))
    type(grid_table), allocatable :: transport_vct(:)
      !! transport vertices grouped by contacts (0:size(contacts))
  contains
    procedure :: init     => device_params_init
    procedure :: destruct => device_params_destruct
  end type

contains

  subroutine device_params_init(this, file)
    !! initialize device parameters using device input file
    class(device_params), intent(out) :: this
    type(input_file),     intent(in)  :: file
      !! device input file

    integer, parameter :: idx_dir0(4) = [ 0, 1, 1, 0]
    integer, parameter :: idx_dir1(4) = [ 0, 2, 2, 0]
    integer, parameter :: eye(2,2) = reshape([1, 0, 0, 1], [2,2])

    integer                   :: i, i0, i1, ii, j, j0, j1, jj, idx(2), si, ci, nx, ny, idx_type, idx_dir
    integer,      allocatable :: sids(:)
    real                      :: vol, surf
    real,         allocatable :: xbounds(:), ybounds(:)
    character(:), allocatable :: table_name

    call init_transport_params()
    call init_grid()
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

      ! make sure parameters are valid
      ASSERT(this%ci0 <= this%ci1)
      ASSERT(size(this%mass) == 2)
      ASSERT(size(this%alpha) == 2)
      ASSERT(size(this%beta) == 2)
      ASSERT(size(this%mob_min) == 2)
      ASSERT(size(this%mob_max) == 2)
      ASSERT(size(this%N_ref) == 2)
      ASSERT(size(this%v_sat) == 2)
    end subroutine

    subroutine init_grid()
      ! find grid section id
      call file%get_section("grid", si)

      ! load grid parameters
      call file%get(si, "nx", nx)
      call file%get(si, "xbounds", xbounds)
      call file%get(si, "ny", ny)
      call file%get(si, "ybounds", ybounds)

      ! init x, y, xy grids
      call this%gx%init(linspace(xbounds(1), xbounds(2), nx))
      call this%gy%init(linspace(ybounds(1), ybounds(2), ny))
      call this%g%init([this%gx%get_ptr(), this%gy%get_ptr()])
    end subroutine

    subroutine init_poisson()
      real :: eps

      ! allocate grid data
      call this%eps%init(this%g, IDX_CELL, 0)
      do idx_dir = 1, 2
        call this%surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
      end do
      call this%vol%init(this%g, IDX_VERTEX, 0)

      ! initialize poisson grid tables
      do idx_type = 1, 4
        table_name = "poisson_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) table_name = table_name//DIR_NAME(idx_dir)
          call this%poisson(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        end do
      end do

      ! get poisson sections
      call file%get_sections("poisson", sids)

      ! process all poisson sections
      do si = 1, size(sids)
        call file%get(sids(si), "xbounds", xbounds)
        call file%get(sids(si), "ybounds", ybounds)
        call file%get(sids(si), "eps",     eps    )

        ! search bounds
        i0 = bin_search(this%gx%x, xbounds(1))
        i1 = bin_search(this%gx%x, xbounds(2)) - 1
        j0 = bin_search(this%gy%x, ybounds(1))
        j1 = bin_search(this%gy%x, ybounds(2)) - 1

        ! enable poisson in region
        do j = j0, j1; do i = i0, i1
          ! set permittivity for cell
          call this%eps%set([i,j], eps)

          ! enable poisson for vertices and update adjoint volumes
          vol = this%g%get_vol([i,j])
          do jj = j, j+1; do ii = i, i+1
            call this%poisson(IDX_VERTEX,0)%flags%set([ii,jj], .true.)
            call this%vol%update([ii,jj], 0.25*vol)
          end do; end do

          do idx_dir = 1, 2
            ! enable poisson for edges
            call this%poisson(IDX_EDGE,idx_dir)%flags%set([i,j],                  .true.)
            call this%poisson(IDX_EDGE,idx_dir)%flags%set([i,j]+eye(:,3-idx_dir), .true.)

            ! enable poisson for faces
            call this%poisson(IDX_FACE,idx_dir)%flags%set([i,j],                  .true.)
            call this%poisson(IDX_FACE,idx_dir)%flags%set([i,j]+eye(:,  idx_dir), .true.)

            ! update adjoint surfaces
            surf = this%g%get_surf([i,j], idx_dir)
            call this%surf(idx_dir)%update([i,j],                  0.5*surf)
            call this%surf(idx_dir)%update([i,j]+eye(:,3-idx_dir), 0.5*surf)
          end do

          ! enable poisson for cell
          call this%poisson(IDX_CELL,0)%flags%set([i, j], .true.)
        end do; end do
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
      do idx_dir = 1, 2
        call this%tr_surf(idx_dir)%init(this%g, IDX_EDGE, idx_dir)
      end do
      call this%tr_vol%init(this%g, IDX_VERTEX, 0)

      ! initialize oxide and transport grid tables
      do idx_type = 1, 4
        table_name = "oxide_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) table_name = table_name//DIR_NAME(idx_dir)
          call this%oxide(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
          this%oxide(idx_type, idx_dir)%flags = this%poisson(idx_type, idx_dir)%flags
        end do
        table_name = "transport_"//IDX_NAME(idx_type)(1:1)
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          if (idx_dir > 0) table_name = table_name//DIR_NAME(idx_dir)
          call this%transport(idx_type, idx_dir)%init(table_name, this%g, idx_type, idx_dir)
        end do
      end do

      ! get transport sections
      call file%get_sections("transport", sids)

      ! process all sections
      do si = 1, size(sids)
        call file%get(sids(si), "xbounds", xbounds)
        call file%get(sids(si), "ybounds", ybounds)

        ! search bounds
        i0 = bin_search(this%gx%x, xbounds(1))
        i1 = bin_search(this%gx%x, xbounds(2)) - 1
        j0 = bin_search(this%gy%x, ybounds(1))
        j1 = bin_search(this%gy%x, ybounds(2)) - 1

        ! enable transport in region
        do j = j0, j1; do i = i0, i1
          ! enable transport for vertices and update adjoint volumes
          vol = this%g%get_vol([i,j])
          do jj = j, j+1; do ii = i, i+1
            call this%oxide(    IDX_VERTEX,0)%flags%set([ii,jj], .false.)
            call this%transport(IDX_VERTEX,0)%flags%set([ii,jj], .true.)
            call this%tr_vol%update([ii,jj], 0.25*vol)
          end do; end do

          do idx_dir = 1, 2
            ! enable transport for edges
            call this%oxide(    IDX_EDGE,idx_dir)%flags%set([i,j],                  .false.)
            call this%transport(IDX_EDGE,idx_dir)%flags%set([i,j],                  .true. )
            call this%oxide(    IDX_EDGE,idx_dir)%flags%set([i,j]+eye(:,3-idx_dir), .false.)
            call this%transport(IDX_EDGE,idx_dir)%flags%set([i,j]+eye(:,3-idx_dir), .true. )

            ! enable transport for faces
            call this%oxide(    IDX_FACE,idx_dir)%flags%set([i,j],                  .false.)
            call this%transport(IDX_FACE,idx_dir)%flags%set([i,j],                  .true. )
            call this%oxide(    IDX_FACE,idx_dir)%flags%set([i,j]+eye(:,  idx_dir), .false.)
            call this%transport(IDX_FACE,idx_dir)%flags%set([i,j]+eye(:,  idx_dir), .true. )

            ! update adjoint surfaces
            surf = this%g%get_surf([i,j], idx_dir)
            call this%tr_surf(idx_dir)%update([i,j],                  0.5*surf)
            call this%tr_surf(idx_dir)%update([i,j]+eye(:,3-idx_dir), 0.5*surf)
          end do

          ! enable transport for cell
          call this%oxide(    IDX_CELL,0)%flags%set([i, j], .false.)
          call this%transport(IDX_CELL,0)%flags%set([i, j], .true.)
        end do; end do
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
      real :: acon, dcon, tr_vol, tr_surf, mob_min, mob_max, N_ref, alpha

      ! allocate grid data
      do idx_type = 1, 4
        do idx_dir = idx_dir0(idx_type), idx_dir1(idx_type)
          call this%acon(idx_type,idx_dir)%init(this%g, idx_type, idx_dir)
          call this%dcon(idx_type,idx_dir)%init(this%g, idx_type, idx_dir)
        end do
      end do

      ! get doping sections
      call file%get_sections("doping", sids)

      ! process all sections
      do si = 1, size(sids)
        call file%get(sids(si), "xbounds", xbounds)
        call file%get(sids(si), "ybounds", ybounds)
        if (this%ci1 == CR_HOLE) call file%get(sids(si), "acon", acon)
        if (this%ci0 == CR_ELEC) call file%get(sids(si), "dcon", dcon)

        ! search bounds
        i0 = bin_search(this%gx%x, xbounds(1))
        i1 = bin_search(this%gx%x, xbounds(2)) - 1
        j0 = bin_search(this%gy%x, ybounds(1))
        j1 = bin_search(this%gy%x, ybounds(2)) - 1

        ! set doping in region
        do j = j0, j1; do i = i0, i1
          vol = this%g%get_vol([i,j])

          ! set doping on vertices
          do jj = j, j+1; do ii = i, i+1
            tr_vol = this%tr_vol%get([ii,jj])
            if (this%ci1 == CR_HOLE) call this%acon(IDX_VERTEX,0)%update([ii,jj], 0.25*vol/tr_vol*acon)
            if (this%ci1 == CR_ELEC) call this%dcon(IDX_VERTEX,0)%update([ii,jj], 0.25*vol/tr_vol*dcon)
          end do; end do

          ! set doping for edges
          do idx_dir = 1, 2
            surf    = this%g%get_surf([i,j], idx_dir)
            tr_surf = this%tr_surf(idx_dir)%get([i,j])
            if (this%ci1 == CR_HOLE) then
              call this%acon(IDX_EDGE,idx_dir)%update([i,j],                  0.5*surf/tr_surf*acon)
              call this%acon(IDX_EDGE,idx_dir)%update([i,j]+eye(:,3-idx_dir), 0.5*surf/tr_surf*acon)
            end if
            if (this%ci0 == CR_ELEC) then
              call this%dcon(IDX_EDGE,idx_dir)%update([i,j],                  0.5*surf/tr_surf*dcon)
              call this%dcon(IDX_EDGE,idx_dir)%update([i,j]+eye(:,3-idx_dir), 0.5*surf/tr_surf*dcon)
            end if
          end do

          ! set doping in cells
          if (this%ci1 == CR_HOLE) call this%acon(IDX_CELL,0)%set([i,j], acon)
          if (this%ci0 == CR_ELEC) call this%dcon(IDX_CELL,0)%set([i,j], dcon)
        end do; end do
      end do

      ! init zero-field mobility
      acon = 0
      dcon = 0
      do ci = this%ci0, this%ci1
        mob_min = this%mob_min(ci)
        mob_max = this%mob_max(ci)
        N_ref   = this%N_ref(ci)
        alpha   = this%alpha(ci)

        do idx_dir = 1, 2
          call this%mob0(idx_dir,ci)%init(this%g, IDX_EDGE, idx_dir)

          do i = 1, this%transport(IDX_EDGE,idx_dir)%n
            idx  = this%transport(IDX_EDGE,idx_dir)%get_idx(i)
            if (this%ci1 == CR_HOLE) acon = this%acon(IDX_EDGE,idx_dir)%get(idx)
            if (this%ci0 == CR_ELEC) dcon = this%dcon(IDX_EDGE,idx_dir)%get(idx)
            call this%mob0(idx_dir,ci)%set(idx, mob_min + (mob_max - mob_min)/(1 + ((acon + dcon)/N_ref)**alpha))
          end do
        end do
      end do
    end subroutine

    subroutine init_contacts()
      integer                           :: ict, ict0, nct, ct_type
      character(:), allocatable         :: typename
      real                              :: phims
      type(string)                      :: name
      type(mapnode_string_int), pointer :: node

      ! get contact sections
      call file%get_sections("contact", sids)

      ! get all contact names
      call this%contact_map%init()
      nct = 0
      do si = 1, size(sids)
        call file%get(sids(si), "name", name%s)
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
        call file%get(sids(si), "name", name%s)
        ict = this%contact_map%get(name)

        ! get type
        call file%get(sids(si), "type", typename)
        if (typename == "ohmic") then
          ct_type = CT_OHMIC
        elseif (typename == "gate") then
          ct_type = CT_GATE
        else
          call program_error("unknown contact type "//typename)
        end if

        ! bounds
        call file%get(sids(si), "xbounds", xbounds)
        call file%get(sids(si), "ybounds", ybounds)
        i0 = bin_search(this%gx%x, xbounds(1))
        i1 = bin_search(this%gx%x, xbounds(2))
        j0 = bin_search(this%gy%x, ybounds(1))
        j1 = bin_search(this%gy%x, ybounds(2))

        ! phims
        if (ct_type == CT_OHMIC) then
          outer: do j = j0, j1; do i = i0, i1
              if (this%transport(IDX_VERTEX,0)%flags%get([i,j])) exit outer
          end do; end do outer
          if ((i > i1) .or. (j > j1)) call program_error("Ohmic contact "//name%s//" not in transport region")
          if ((this%ci0 == CR_ELEC) .and. (this%ci1 == CR_HOLE)) then
            phims = asinh(0.5 * (this%dcon(IDX_VERTEX,0)%get([i,j]) - this%acon(IDX_VERTEX,0)%get([i,j])) / this%n_intrin)
          elseif (this%ci0 == CR_ELEC) then
            phims = log(this%dcon(IDX_VERTEX,0)%get([i,j]) / this%n_intrin)
          elseif (this%ci0 == CR_HOLE) then
            phims = -log(this%acon(IDX_VERTEX,0)%get([i,j]) / this%n_intrin)
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
        do j = j0, j1; do i = i0, i1
          idx = [i, j]
          ict0 = this%ict%get(idx)
          if ((ict0 /= 0) .and. (ict0 /= ict)) then
            call program_error("contacts "//this%contacts(ict)%name//" and "//this%contacts(ict0)%name//" are overlapping")
          end if
          call this%ict%set(idx, ict)
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
        end do; end do
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
