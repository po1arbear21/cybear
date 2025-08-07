module region_m

  use contact_m,       only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use error_m,         only: program_error
  use input_m,         only: input_file
  use math_m,          only: PI
  use semiconductor_m, only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON
  use string_m,        only: string

  implicit none

  type, abstract :: region
    real, allocatable :: xyz(:,:)
      !! bounds in x, y, z xyz:(dir, 2) tr:(dir, >2)
  contains
    procedure :: init       => region_init
    procedure :: point_test => region_point_test
  end type

  type, extends(region) :: region_poisson
    real :: eps
      !! permittivity
  end type

  type, extends(region) :: region_transport
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
      !! contact type (CT_OHMIC, CT_GATE, or CT_SCHOTTKY)
    real         :: phims
      !! metal-semiconductor workfunction difference (one value per contact)
    real         :: barrier_height
      !! Schottky barrier height Φ_B (eV)
    real         :: richardson_const
      !! Richardson constant A* (A/cm²/K²)
    real         :: surf_recomb_vel(2)
      !! surface recombination velocities S_n, S_p (cm/s)
    logical      :: tunneling_enabled
      !! enable field emission tunneling
  contains
    procedure :: point_test => region_contact_point_test
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
      else
        call program_error("unknown contact type "//type%s)
      end if
      call file%get(sid, "phims", this%phims, status = st)
      if (.not. st) this%phims = 0
      
      ! get Schottky-specific parameters
      if (this%type == CT_SCHOTTKY) then
        call file%get(sid, "barrier_height", this%barrier_height, status = st)
        if (.not. st) this%barrier_height = 0.7  ! default 0.7 eV
        call file%get(sid, "richardson_const", this%richardson_const, status = st)
        if (.not. st) this%richardson_const = 112.0  ! default for Si
        print *, "DEBUG region_contact_init: Contact", trim(this%name%s), "type=SCHOTTKY"
        print *, "  barrier_height =", this%barrier_height, "eV"
        print *, "  richardson_const =", this%richardson_const, "A/cm²/K²"
        call file%get(sid, "surf_recomb_vel_n", this%surf_recomb_vel(1), status = st)
        if (.not. st) this%surf_recomb_vel(1) = 1.0e7  ! default 10^7 cm/s
        call file%get(sid, "surf_recomb_vel_p", this%surf_recomb_vel(2), status = st)
        if (.not. st) this%surf_recomb_vel(2) = 1.0e7  ! default 10^7 cm/s
        call file%get(sid, "tunneling", this%tunneling_enabled, status = st)
        if (.not. st) this%tunneling_enabled = .false.  ! default off
      end if

    end select

  end subroutine

  function region_point_test(this, gtype, point) result(t)
    !! check if point is inside region (triangular grid: count ray intersections with boundary)
    class(region), intent(in) :: this
    type(string),  intent(in) :: gtype
      !! grid type
    real,          intent(in) :: point(:)

    real, parameter :: e2 = 2.7182818284590452353602874713/2, d(2) = [cos(e2), sin(e2)]
      !! avoid d || region boundary by using strange angle (FIXME: zero determinant still possible in edge case)

    logical :: t
    real    :: k, s, p0x, p0y, p1x, p1y, p2x, p2y, dx, dy, ex, ey, fx, fy, det
    integer :: n_intersection, i, j

    select case(gtype%s)
    case ("x")
      t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1))

    case ("xy")
      t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1)) .and. &
        & (this%xyz(2,1) <= point(2)) .and. (this%xyz(2,2) >= point(2))

    case ("xyz")
      t = (this%xyz(1,1) <= point(1)) .and. (this%xyz(1,2) >= point(1)) .and. &
        & (this%xyz(2,1) <= point(2)) .and. (this%xyz(2,2) >= point(2)) .and. &
        & (this%xyz(3,1) <= point(3)) .and. (this%xyz(3,2) >= point(3))

    case("tr_xy")
      n_intersection = 0
      p0x = point(1)
      p0y = point(2)
      dx = d(1)
      dy = d(2)
      do i = 1, size(this%xyz, dim = 2)
        j = mod(i, size(this%xyz, dim = 2)) + 1
        p1x = this%xyz(1, i)
        p1y = this%xyz(2, i)
        p2x = this%xyz(1, j)
        p2y = this%xyz(2, j)
        ex = p1x - p2x
        ey = p1y - p2y
        fx = p1x - p0x
        fy = p1y - p0y
        det = dx*ey - dy*ex
        if (det == 0) call program_error("zero determinant, d must be changed")

        ! check whether intersection point lies within segment
        k = -(ex*fy - ey*fx)/det
        s = (dx*fy - dy*fx)/det
        if ((s >= 0) .and. (s <= 1) .and. (k >= 0)) n_intersection = n_intersection +1
      end do

      ! point inside region if number of intersections uneven
      t = (mod(n_intersection,2) == 1)

    case("tr_xyz")
      n_intersection = 0
      p0x = point(1)
      p0y = point(2)
      dx = d(1)
      dy = d(2)
      do i = 1, size(this%xyz, dim = 2)
        j = mod(i, size(this%xyz, dim = 2)) + 1
        p1x = this%xyz(1, i)
        p1y = this%xyz(2, i)
        p2x = this%xyz(1, j)
        p2y = this%xyz(2, j)
        ex = p1x - p2x
        ey = p1y - p2y
        fx = p1x - p0x
        fy = p1y - p0y
        det = dx*ey - dy*ex
        if (det == 0) call program_error("zero determinant, d must be changed")

        ! check whether intersection point lies within segment
        k = -(ex*fy - ey*fx)/det
        s = (dx*fy - dy*fx)/det
        if ((s>=0) .and. (s<=1) .and. (k>=0)) n_intersection = n_intersection +1
      end do

      ! point inside region if number of intersections uneven and within z boundaries
      t = (mod(n_intersection,2) == 1) .and. (this%xyz(3,1) <= point(3)) .and. (this%xyz(3,2) >= point(3))

    end select
  end function

  function region_contact_point_test(this, gtype, point) result(t)
    !! check if point is inside contact region (triangular grids: check "tube" around polygonal chain)
    class(region_contact), intent(in) :: this
    type(string),          intent(in) :: gtype
      !! grid type
    real,                  intent(in) :: point(:)
     !! to check

    logical         :: t, t1
    real, parameter :: TOL = 1e-10
    real            :: x, y, dx, dy, p0x, p0y, p1x, p1y, p2x, p2y, ex, ey, fx, fy, det, k, s, L, d, e
    integer         :: i, j

    select case(gtype%s)
    case ("x")
      t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1))

    case ("xy")
      t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1)) .and. &
        & (this%xyz(2,1) - TOL <= point(2)) .and. (this%xyz(2,2) + TOL >= point(2))

    case ("xyz")
      t = (this%xyz(1,1) - TOL <= point(1)) .and. (this%xyz(1,2) + TOL >= point(1)) .and. &
        & (this%xyz(2,1) - TOL <= point(2)) .and. (this%xyz(2,2) + TOL >= point(2)) .and. &
        & (this%xyz(3,1) - TOL <= point(3)) .and. (this%xyz(3,2) + TOL >= point(3))

    case("tr_xy")
      p0x = point(1)
      p0y = point(2)
      x = 1
      t = .true.
      do i = 1, size(this%xyz, dim = 2)-1
        j = i + 1
        p1x = this%xyz(1, i)
        p1y = this%xyz(2, i)
        p2x = this%xyz(1, j)
        p2y = this%xyz(2, j)
        ex = p1x - p2x
        ey = p1y - p2y
        fx = p1x - p0x
        fy = p1y - p0y
        L = sqrt(ex**2 + ey**2)
        if (ey==0) then
          dx = 0
          dy = 1
        elseif (ex == 0) then
          dx = 1
          dy = 0
        else
          y = - x*ex / ey
          dx = x/sqrt(x**2 + y**2)
          dy = y/sqrt(x**2 + y**2)
        end if
        det = dx*ey - dy*ex
        k = -(ex*fy - ey*fx)/det
        s = (dx*fy - dy*fx)/det
        d = abs(k)
        e = s*L
        if ((d <= TOL) .and. (e >= -TOL) .and. (e <= L+TOL)) return
      end do
      t = .false.

    case("tr_xyz")
      p0x = point(1)
      p0y = point(2)
      x = 1
      do i = 1, size(this%xyz, dim = 2)
        j = mod(i, size(this%xyz, dim = 2)) + 1
        p1x = this%xyz(1, i)
        p1y = this%xyz(2, i)
        p2x = this%xyz(1, j)
        p2y = this%xyz(2, j)
        ex = p1x - p2x
        ey = p1y - p2y
        fx = p1x - p0x
        fy = p1y - p0y
        L = sqrt(ex**2 + ey**2)
        y = - x*ex / ey
        dx = x/sqrt(x**2 + y**2)
        dy = y/sqrt(x**2 + y**2)
        det = dx*ey - dy*ex
        k = -(ex*fy - ey*fx)/det
        s = (dx*fy - dy*fx)/det
        d = abs(k)
        e = sqrt(fx**2 + fy**2 - k**2)
        t1 = ((d <= TOL) .and. (e >= -TOL) .and. (e <= L+TOL) )
        if (t1) t = .true.
      end do
      t = (t  .and. (this%xyz(3,1) - TOL <= point(3)) .and. (this%xyz(3,2) + TOL >= point(3)) )

    end select
  end function

end module
