m4_include(util/macro.f90.inc)

module region_m

  use contact_m,       only: CT_OHMIC, CT_GATE, CT_SCHOTTKY, CT_REALOHMIC
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: grid, IDX_CELL, IDX_VERTEX
  use grid_table_m,    only: grid_table
  use input_m,         only: input_file
  use semiconductor_m, only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON
  use string_m,        only: string

  implicit none

  type, abstract :: region
    type(string) :: gtype
      !! grid type ("x", "xy", "xyz", "tr_xy", "tr_xyz")
    real, allocatable :: xyz(:,:)
      !! bounds in x, y, z xyz:(dir, 2) tr:(dir, >2)
    type(grid_table) :: cell_table
      !! grid table of cells that belong to region
    type(grid_table) :: vert_table
      !! grid table of vertices that belong to region, including vertices on the boundary
  contains
    procedure :: init             => region_init
    procedure :: init_grid_tables => region_init_grid_tables
    procedure :: point_test       => region_point_test
  end type

  type, extends(region) :: region_poisson
    real :: eps
      !! permittivity
  end type

  type, extends(region) :: region_transport
    type(string) :: material
      !! material name (matches a [semiconductor].name); empty when set_material is .false.
    logical      :: set_material = .false.
      !! true if material= key was present in the [transport] section
    integer      :: material_id
      !! resolved index into device_params%mat(:); set by init_region_materials post-pass
  end type

  type, extends(region) :: region_doping
    real :: dop(2)
      !! donor/acceptor concentration
    logical :: set_dop(2)
      !! enable/disable set doping for each carrier
  end type

  type, extends(region) :: region_mobility
    real    :: mob0(2)
      !! zero-field mobility
    logical :: set_mob(2)
      !! enable/disable set mobility for each carrier
  end type

  type, extends(region) :: region_contact
    type(string) :: name
      !! contact name
    integer      :: type
      !! contact type (CT_OHMIC, CT_GATE, CT_SCHOTTKY, or CT_REALOHMIC)
    real         :: phims
      !! metal-semiconductor workfunction difference (one value per contact)
    real         :: phi_b
      !! Schottky barrier height (normalized to kT)
    real         :: A_richardson_n
      !! Richardson constant for electrons (A/cm^2/K^2)
    real         :: A_richardson_p
      !! Richardson constant for holes (A/cm^2/K^2)
    logical      :: ifbl
      !! image force barrier lowering flag
    logical      :: tunneling
      !! enable Tsu-Esaki tunneling model
    real         :: m_tunnel_n
      !! tunneling effective mass ratio for electrons (m*/m0)
    real         :: m_tunnel_p
      !! tunneling effective mass ratio for holes (m*/m0)
    real         :: vrec
      !! metal-semiconductor recombination velocity (real ohmic contact)
  ! contains
  !   procedure :: point_test => region_contact_point_test
  end type

  type region_ptr
    class(region), pointer :: p => null()
  end type

contains

  subroutine region_init(this, file, sid, gtype)
    class(region),    intent(out) :: this
    type(input_file), intent(in)  :: file
    integer,          intent(in)  :: sid
    type(string),     intent(in)  :: gtype

    real, allocatable :: tmpx(:), tmpy(:), tmpz(:)
    logical           :: st
    type(string)      :: type

    this%gtype = gtype

    ! get bounds
    select case(gtype%s)
    case ("x")
      call file%get(sid, "x", tmpx)
      if (size(tmpx) /= 2) call program_error("region x bounds must be an array of size 2")

      allocate (this%xyz(1,2))
      this%xyz(1,:) = tmpx

    case ("xy")
      call file%get(sid, "x", tmpx)
      if (size(tmpx) /= 2) call program_error("region x bounds must be an array of size 2")
      call file%get(sid, "y", tmpy)
      if (size(tmpy) /= 2) call program_error("region y bounds must be an array of size 2")

      allocate (this%xyz(2,2))
      this%xyz(1,:) = tmpx
      this%xyz(2,:) = tmpy

    case("xyz")
      call file%get(sid, "x", tmpx)
      if (size(tmpx) /= 2) call program_error("region x bounds must be an array of size 2")
      call file%get(sid, "y", tmpy)
      if (size(tmpy) /= 2) call program_error("region y bounds must be an array of size 2")
      call file%get(sid, "z", tmpz)
      if (size(tmpz) /= 2) call program_error("region z bounds must be an array of size 2")

      allocate (this%xyz(3,2))
      this%xyz(1,:) = tmpx
      this%xyz(2,:) = tmpy
      this%xyz(3,:) = tmpz

    case("tr_xy")
      call file%get(sid, "x", tmpx)
      if (size(tmpx) < 2) call program_error("region x bounds must be an array of at least size 2")
      call file%get(sid, "y", tmpy)
      if (size(tmpx) /= size(tmpy)) call program_error("region x and y bounds must be of same size")

      allocate (this%xyz(2,size(tmpx)))
      this%xyz(1,:) = tmpx
      this%xyz(2,:) = tmpy

    case("tr_xyz")
      call file%get(sid, "x", tmpx)
      if (size(tmpx) < 2) call program_error("region x bounds must be an array of at least size 2")
      call file%get(sid, "y", tmpy)
      if (size(tmpx) /= size(tmpy)) call program_error("region x and y bounds must be of same size")
      call file%get(sid, "z", tmpz)
      if (size(tmpz) /= 2) call program_error("region z bounds must be an array of size 2")

      allocate (this%xyz(3,size(tmpx)))
      this%xyz(1,:) = tmpx
      this%xyz(2,:) = tmpy
      this%xyz(3,1:2) = tmpz

    end select

    select type(this)
    class is (region_poisson)
      ! get permittivity
      call file%get(sid, "eps", this%eps)

    class is (region_transport)
      ! optional material reference; resolution to material_id happens in
      ! device_params%init_region_materials once the catalog is built.
      call file%get(sid, "material", this%material, status = this%set_material)

    class is (region_doping)
      ! get ND, NA
      call file%get(sid, "dcon", this%dop(DOP_DCON), status = this%set_dop(DOP_DCON))
      call file%get(sid, "acon", this%dop(DOP_ACON), status = this%set_dop(DOP_ACON))

    class is (region_mobility)
      ! get ND, NA
      call file%get(sid, "nmob0", this%mob0(CR_ELEC), status = this%set_mob(CR_ELEC))
      call file%get(sid, "pmob0", this%mob0(CR_HOLE), status = this%set_mob(CR_HOLE))

    class is (region_contact)
      ! get name, type and phims
      call file%get(sid, "name",  this%name)
      call file%get(sid, "type",  type)
      if (type%s == "ohmic") then
        this%type = CT_OHMIC
      elseif (type%s == "gate") then
        this%type = CT_GATE
      elseif (type%s == "schottky") then
        this%type = CT_SCHOTTKY
        call file%get(sid, "phi_b", this%phi_b, status = st)
        if (.not. st) call program_error("Schottky contact '"//this%name%s//"': phi_b is required")
        call file%get(sid, "A_richardson_n", this%A_richardson_n, status = st)
        if (.not. st) call program_error("Schottky contact '"//this%name%s//"': A_richardson_n is required")
        call file%get(sid, "A_richardson_p", this%A_richardson_p, status = st)
        if (.not. st) call program_error("Schottky contact '"//this%name%s//"': A_richardson_p is required")
        call file%get(sid, "ifbl", this%ifbl, status = st)
        if (.not. st) this%ifbl = .false.
        call file%get(sid, "tunneling", this%tunneling, status = st)
        if (.not. st) this%tunneling = .false.
        call file%get(sid, "m_tunnel_n", this%m_tunnel_n, status = st)
        if (.not. st) this%m_tunnel_n = 1.0
        call file%get(sid, "m_tunnel_p", this%m_tunnel_p, status = st)
        if (.not. st) this%m_tunnel_p = 1.0
      elseif (type%s == "realohmic") then
        this%type = CT_REALOHMIC
        call file%get(sid, "vrec", this%vrec)
      else
        call program_error("unknown contact type "//type%s)
      end if
      call file%get(sid, "phims", this%phims, status = st)
      if (.not. st) this%phims = 0

    end select

  end subroutine

  subroutine region_init_grid_tables(this, g)
    !! initialize cell and vertex grid table
    class(region), intent(inout) :: this
    class(grid),   intent(in)    :: g
      !! full grid of the device

    integer              :: i, j
    integer, allocatable :: idx_bnd(:,:), idx(:), idx1(:)
    logical              :: stat
    real,    allocatable :: p(:,:), mid(:)
    real,    parameter   :: TOL = 1e-10 ! tolerance when checking if vertices are inside region

    allocate (idx(g%idx_dim), idx1(g%idx_dim), idx_bnd(2, g%idx_dim), p(g%dim, g%cell_nvert), mid(g%dim))

    ! determine which cells belong to region
    call this%cell_table%init("region_cells", g, IDX_CELL, 0)
    ! iterate over all cell indices and check if they belong to region
    call g%get_idx_bnd(IDX_CELL, 0, idx_bnd)
    idx = idx_bnd(1, :)
    do while (idx(g%idx_dim) <= idx_bnd(2, g%idx_dim))
      ! cell center = average of vertex coordinates, this works because cells are convex
      call g%get_cell(idx, p)
      mid = sum(p, dim = 2) / g%cell_nvert
      ! add cell to grid table if it is in region
      ! sufficient to check cell center because region boundary should be aligned with grid/cell boundaries
      if (this%point_test(mid, 0.0)) call this%cell_table%flags%set(idx, .true.)

      ! go to next index
      idx(1) = idx(1) + 1
      do j = 1, g%idx_dim-1
        if (idx(j) <= idx_bnd(2,j)) exit
        idx(j  ) = idx_bnd(1,j)
        idx(j+1) = idx(j+1) + 1
      end do
    end do
    call this%cell_table%init_final()

    ! determine which vertices belong to region
    call this%vert_table%init("region_vertices", g, IDX_VERTEX, 0)
    if (this%cell_table%n > 0) then
      ! region with finite extension: loop over cells in region and add their vertices
      do i = 1, this%cell_table%n
        idx1 = this%cell_table%get_idx(i)
        do j = 1, g%cell_nvert
          call g%get_neighb(IDX_CELL, 0, IDX_VERTEX, 0, idx1, j, idx, stat)
          m4_assert(stat)
          call this%vert_table%flags%set(idx, .true.)
        end do
      end do
    else
      ! region with zero extension: iterate over all vertex indices and check if they belong to region
      call g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
      idx = idx_bnd(1, :)
      do while (idx(g%idx_dim) <= idx_bnd(2, g%idx_dim))
        call g%get_vertex(idx, mid)
        ! add vertex to grid table if it is in region
        if (this%point_test(mid, TOL)) call this%vert_table%flags%set(idx, .true.)

        ! go to next index
        idx(1) = idx(1) + 1
        do j = 1, g%idx_dim-1
          if (idx(j) <= idx_bnd(2,j)) exit
          idx(j  ) = idx_bnd(1,j)
          idx(j+1) = idx(j+1) + 1
        end do
      end do
    end if
    call this%vert_table%init_final()
  end subroutine

  function region_point_test(this, point, tol) result(t)
    !! check if point is inside region (triangular grid: count ray intersections with boundary)
    class(region), intent(in) :: this
    real,          intent(in) :: point(:)
      !! point to check
    real,          intent(in) :: tol
      !! tolerance for point to be considered inside region

    logical :: t

    select case(this%gtype%s)
    case ("x")
      t = (this%xyz(1,1) - tol <= point(1)) .and. (this%xyz(1,2) + tol >= point(1))

    case ("xy")
      t = (this%xyz(1,1) - tol <= point(1)) .and. (this%xyz(1,2) + tol >= point(1)) .and. &
        & (this%xyz(2,1) - tol <= point(2)) .and. (this%xyz(2,2) + tol >= point(2))

    case ("xyz")
      t = (this%xyz(1,1) - tol <= point(1)) .and. (this%xyz(1,2) + tol >= point(1)) .and. &
        & (this%xyz(2,1) - tol <= point(2)) .and. (this%xyz(2,2) + tol >= point(2)) .and. &
        & (this%xyz(3,1) - tol <= point(3)) .and. (this%xyz(3,2) + tol >= point(3))

    case("tr_xy", "tr_xyz")
      ! TODO: implement point test for triangular regions (old and possibly erroneous implementation commented out below)
      call program_error("point test for triangular regions not implemented/sufficiently tested yet")

    end select
  end function

  ! function region_point_test(this, gtype, point) result(t)
  !   !! check if point is inside region (triangular grid: count ray intersections with boundary)
  !   class(region), intent(in) :: this
  !   type(string),  intent(in) :: gtype
  !     !! grid type
  !   real,          intent(in) :: point(:)

  !   real, parameter :: e2 = 2.7182818284590452353602874713/2, d(2) = [cos(e2), sin(e2)]
  !     !! avoid d || region boundary by using strange angle (FIXME: zero determinant still possible in edge case)

  !   logical :: t
  !   real    :: k, s, p0x, p0y, p1x, p1y, p2x, p2y, dx, dy, ex, ey, fx, fy, det
  !   integer :: n_intersection, i, j

  !   select case(gtype%s)
  !   case ("x")
  !     t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1))

  !   case ("xy")
  !     t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1)) .and. &
  !       & (this%xyz(2,1) <= point(2)) .and. (this%xyz(2,2) >= point(2))

  !   case ("xyz")
  !     t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1)) .and. &
  !       & (this%xyz(2,1) <= point(2)) .and. (this%xyz(2,2) >= point(2)) .and. &
  !       & (this%xyz(3,1) <= point(3)) .and. (this%xyz(3,2) >= point(3))

  !   case("tr_xy")
  !     n_intersection = 0
  !     p0x = point(1)
  !     p0y = point(2)
  !     dx = d(1)
  !     dy = d(2)
  !     do i = 1, size(this%xyz, dim = 2)
  !       j = mod(i, size(this%xyz, dim = 2)) + 1
  !       p1x = this%xyz(1, i)
  !       p1y = this%xyz(2, i)
  !       p2x = this%xyz(1, j)
  !       p2y = this%xyz(2, j)
  !       ex = p1x - p2x
  !       ey = p1y - p2y
  !       fx = p1x - p0x
  !       fy = p1y - p0y
  !       det = dx*ey - dy*ex
  !       if (det == 0) call program_error("zero determinant, d must be changed")

  !       ! check whether intersection point lies within segment
  !       k = -(ex*fy - ey*fx)/det
  !       s = (dx*fy - dy*fx)/det
  !       if ((s >= 0) .and. (s <= 1) .and. (k >= 0)) n_intersection = n_intersection +1
  !     end do

  !     ! point inside region if number of intersections uneven
  !     t = (mod(n_intersection,2) == 1)

  !   case("tr_xyz")
  !     n_intersection = 0
  !     p0x = point(1)
  !     p0y = point(2)
  !     dx = d(1)
  !     dy = d(2)
  !     do i = 1, size(this%xyz, dim = 2)
  !       j = mod(i, size(this%xyz, dim = 2)) + 1
  !       p1x = this%xyz(1, i)
  !       p1y = this%xyz(2, i)
  !       p2x = this%xyz(1, j)
  !       p2y = this%xyz(2, j)
  !       ex = p1x - p2x
  !       ey = p1y - p2y
  !       fx = p1x - p0x
  !       fy = p1y - p0y
  !       det = dx*ey - dy*ex
  !       if (det == 0) call program_error("zero determinant, d must be changed")

  !       ! check whether intersection point lies within segment
  !       k = -(ex*fy - ey*fx)/det
  !       s = (dx*fy - dy*fx)/det
  !       if ((s>=0) .and. (s<=1) .and. (k>=0)) n_intersection = n_intersection +1
  !     end do

  !     ! point inside region if number of intersections uneven and within z boundaries
  !     t = (mod(n_intersection,2) == 1) .and. (this%xyz(3,1) <= point(3)) .and. (this%xyz(3,2) >= point(3))

  !   end select
  ! end function

  ! function region_contact_point_test(this, gtype, point) result(t)
  !   !! check if point is inside contact region (triangular grids: check "tube" around polygonal chain)
  !   class(region_contact), intent(in) :: this
  !   type(string),          intent(in) :: gtype
  !     !! grid type
  !   real,                  intent(in) :: point(:)
  !    !! to check

  !   logical         :: t, t1
  !   real, parameter :: TOL = 1e-10
  !   real            :: x, y, dx, dy, p0x, p0y, p1x, p1y, p2x, p2y, ex, ey, fx, fy, det, k, s, L, d, e
  !   integer         :: i, j

  !   select case(gtype%s)
  !   case ("x")
  !     t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1))

  !   case ("xy")
  !     t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1)) .and. &
  !       & (this%xyz(2,1) - TOL <= point(2)) .and. (this%xyz(2,2) + TOL >= point(2))

  !   case ("xyz")
  !     t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1)) .and. &
  !       & (this%xyz(2,1) - TOL <= point(2)) .and. (this%xyz(2,2) + TOL >= point(2)) .and. &
  !       & (this%xyz(3,1) - TOL <= point(3)) .and. (this%xyz(3,2) + TOL >= point(3))

  !   case("tr_xy")
  !     p0x = point(1)
  !     p0y = point(2)
  !     x = 1
  !     t = .true.
  !     do i = 1, size(this%xyz, dim = 2)-1
  !       j = i + 1
  !       p1x = this%xyz(1, i)
  !       p1y = this%xyz(2, i)
  !       p2x = this%xyz(1, j)
  !       p2y = this%xyz(2, j)
  !       ex = p1x - p2x
  !       ey = p1y - p2y
  !       fx = p1x - p0x
  !       fy = p1y - p0y
  !       L = sqrt(ex**2 + ey**2)
  !       if (ey==0) then
  !         dx = 0
  !         dy = 1
  !       elseif (ex == 0) then
  !         dx = 1
  !         dy = 0
  !       else
  !         y = - x*ex / ey
  !         dx = x/sqrt(x**2 + y**2)
  !         dy = y/sqrt(x**2 + y**2)
  !       end if
  !       det = dx*ey - dy*ex
  !       k = -(ex*fy - ey*fx)/det
  !       s = (dx*fy - dy*fx)/det
  !       d = abs(k)
  !       e = s*L
  !       if ((d <= TOL) .and. (e >= -TOL) .and. (e <= L+TOL)) return
  !     end do
  !     t = .false.

  !   case("tr_xyz")
  !     p0x = point(1)
  !     p0y = point(2)
  !     x = 1
  !     do i = 1, size(this%xyz, dim = 2)
  !       j = mod(i, size(this%xyz, dim = 2)) + 1
  !       p1x = this%xyz(1, i)
  !       p1y = this%xyz(2, i)
  !       p2x = this%xyz(1, j)
  !       p2y = this%xyz(2, j)
  !       ex = p1x - p2x
  !       ey = p1y - p2y
  !       fx = p1x - p0x
  !       fy = p1y - p0y
  !       L = sqrt(ex**2 + ey**2)
  !       y = - x*ex / ey
  !       dx = x/sqrt(x**2 + y**2)
  !       dy = y/sqrt(x**2 + y**2)
  !       det = dx*ey - dy*ex
  !       k = -(ex*fy - ey*fx)/det
  !       s = (dx*fy - dy*fx)/det
  !       d = abs(k)
  !       e = sqrt(fx**2 + fy**2 - k**2)
  !       t1 = ((d <= TOL) .and. (e >= -TOL) .and. (e <= L+TOL) )
  !       if (t1) t = .true.
  !     end do
  !     t = (t  .and. (this%xyz(3,1) - TOL <= point(3)) .and. (this%xyz(3,2) + TOL >= point(3)) )

  !   end select
  ! end function

end module
