m4_include(util/macro.f90.inc)

module grid_generator_m

  use error_m,       only: assert_failed, program_error
  use galene_m,      only: gal_file, gal_block
  use grid_m,        only: grid, grid_ptr
  use grid1D_m,      only: grid1D
  use input_m,       only: input_file
  use tensor_grid_m, only: tensor_grid
  use triang_grid_m, only: triang_grid
  use triangle_m,    only: triangulation
  use region_m,      only: region, region_ptr, region_poisson, region_transport, region_doping, region_contact
  use vector_m,      only: vector_real
  use math_m,        only: linspace, expm1
  use qsort_m,       only: qsort
use normalization_m
  implicit none

  private
  public  :: DIR_NAME, generate_1D_grid, generate_triangle_grid

   ! parameters
  character(*), parameter :: DIR_NAME(3)  = [ "x", "y", "z"]

  type refinement
    integer :: dir
    real    :: x(2)
      !! interval
    real    :: dx(2)
      !! delta x at interval end points
    real    :: e
      !! exponent
  contains
    procedure :: init   => refinement_init
    procedure :: get_dx => refinement_get_dx
    procedure :: get_x  => refinement_get_x
  end type

contains

  subroutine generate_1D_grid(file, load, gal, gal_fl, dir, reg, g1D, gptr, ngptr)
    !! create 1D grid x, y, or z
    type(input_file),     intent(in)    :: file
      !! device file
    logical,              intent(in)    :: load
      !! load or generate?
    logical,              intent(in)    :: gal
      !! load from GALENE III file
    type(gal_file),       intent(in)    :: gal_fl
      !! GALENE III file
    integer,              intent(in)    :: dir
      !! axis direction
    type(region_ptr),     intent(in)    :: reg(:)
      !! pointers to all regions
    type(grid1D), target, intent(out)   :: g1D
      !! output newly generated 1D grid
    type(grid_ptr),       intent(inout) :: gptr(:)
      !! save pointer to g1D in gptr (used to create tensor grid afterwards)
    integer,              intent(inout) :: ngptr
      !! number of already used gptr entries (gets increased by 1)

    real                          :: max_dx
    real,             allocatable :: x0(:), x(:)
    type(gal_block),  pointer     :: gblock
    type(refinement), allocatable :: ref(:)

    ! load or generate axis
    if (load) then
      if (gal) then
        gblock => gal_fl%get_block(DIR_NAME(dir) // "-coord")
        x = gblock%rdata
      else
        call file%get("grid", DIR_NAME(dir), x)
      end if
    else
      m4_assert(.not. gal)
      call file%get("grid", "max_d"//DIR_NAME(dir), max_dx)

      ! generate grid points
      ref = get_refinements(file, dir)
      x0  = get_bounds(dir, reg, ref)
      x   = generate_axis(x0, max_dx, ref)
    end if

    ! initialize grid
    call g1D%init(DIR_NAME(dir), x)

    ! update pointer list
    ngptr         =  ngptr + 1
    gptr(ngptr)%p => g1D
  end subroutine

  subroutine generate_triangle_grid(file, load, reg, tr, gtr, gptr, ngptr)
    !! create triangle grid
    type(input_file),          intent(in)    :: file
      !! device file
    logical,                   intent(in)    :: load
      !! load or generate? (load not implemented yet)
    type(region_ptr),          intent(in)    :: reg(:)
      !! pointers to all regions
    type(triangulation),       intent(out)   :: tr
      !! output triangulation
    type(triang_grid), target, intent(out)   :: gtr
      !! output corresponding triangle grid
    type(grid_ptr),            intent(inout) :: gptr(:)
      !! save pointer to g1D in gptr (used to create tensor grid afterwards)
    integer,                   intent(inout) :: ngptr
      !! number of already used gptr entries (gets increased by 1)

    integer :: i
    real    :: max_area

    ! load or generate triangle grid
    if (load) then
      call program_error("loading triangle grid not implemented yet")
    else
      call tr%init()
      do i = 1, size(reg)
        call tr%add_polygon(reg(i)%p%xyz(1,:), reg(i)%p%xyz(2,:), closed = .true.)
      end do
      call file%get("grid", "max_area", max_area)
      call tr%triangulate(max_area = max_area)
      call tr%get_grid("tr", gtr)
    end if

    ! update pointer list
    ngptr = ngptr + 1
    gptr(ngptr)%p => gtr
  end subroutine

  function get_refinements(file, dim) result(ref)
    type(input_file), intent(in)  :: file
    integer,          intent(in)  :: dim
    type(refinement), allocatable :: ref(:)

    integer                       :: i, j
    integer,          allocatable :: sids(:), perm(:)
    real,             allocatable :: x1(:)
    type(refinement), allocatable :: tmp(:)

    call file%get_sections("refine grid", sids)
    allocate (tmp(size(sids)))
    j = 0
    do i = 1, size(sids)
      j = j + 1
      call tmp(j)%init(file, sids(i))
      if (tmp(j)%dir /= dim) j = j - 1
    end do

    x1 = tmp(1:j)%x(1)
    allocate (perm(j))
    call qsort(x1, perm=perm)
    ref = tmp(perm)

    do i = 2, j
      if (ref(i)%x(1) < ref(i-1)%x(2)) call program_error("refinements are overlapping")
    end do
  end function

  function get_bounds(dim, reg, ref) result(x0)
    integer,          intent(in) :: dim
    type(region_ptr), intent(in) :: reg(:)
    type(refinement), intent(in) :: ref(:)
    real, allocatable            :: x0(:)

    integer :: i, j

    ! collect bounds of all regions and refinements in array
    allocate (x0(2*(size(reg) + size(ref))))
    j = 0
    do i = 1, size(reg)
      x0(j+1:j+2) = reg(i)%p%xyz(dim,1:2)
      j = j + 2
    end do
    do i = 1, size(ref)
      x0(j+1:j+2) = ref(i)%x
      j = j + 2
    end do
  end function

  function generate_axis(x0, dx, ref) result(x)
    real,             intent(in) :: x0(:)
    real,             intent(in) :: dx
    type(refinement), intent(in) :: ref(:)
    real, allocatable            :: x(:)

    real, parameter :: TOL = 1e-10

    real              :: xm
    real, allocatable :: xtmp(:), xx(:)
    integer           :: n, m, i, j, iref
    logical           :: refine
    type(vector_real) :: xv

    ! sort region bounds
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
    iref = 1
    do i = 1, n-1
      refine = .false.
      if (iref <= size(ref)) then
        xm = 0.5*(xtmp(i) + xtmp(i+1))
        if ((ref(iref)%x(1) <= xm) .and. (ref(iref)%x(2) >= xm)) refine = .true.
      end if

      if (refine) then
        ! spacing from refinement
        if ((xtmp(i) > ref(iref)%x(1)) .or. (xtmp(i+1) < ref(iref)%x(2))) then
          call program_error("refinement and region bounds are overlapping")
        end if
        xx   = ref(iref)%get_x()
        iref = iref + 1
      else
        ! linear spacing
        m  = ceiling((xtmp(i+1)-xtmp(i))/dx) + 1
        xx = linspace(xtmp(i), xtmp(i+1), m)
      end if
      call xv%push(xx(2:size(xx)))
    end do
    x= xv%to_array()
  end function

  subroutine refinement_init(this, file, sid)
    class(refinement), intent(out) :: this
    type(input_file),  intent(in)  :: file
    integer,           intent(in)  :: sid

    integer           :: dir
    logical           :: status
    real, allocatable :: tmp(:)

    this%dir = 0
    do dir = 1, 3
      call file%get(sid, DIR_NAME(dir), tmp, status = status)
      if (.not. status) cycle
      if (this%dir > 0) call program_error("Found multiple directions in refine grid section")
      this%dir = dir

      if (size(tmp) /= 2) call program_error("coordinates must have size 2")
      this%x = tmp

      call file%get(sid, "d"//DIR_NAME(dir), tmp, status = status)
      if (.not. status) call program_error("expected delta array, e.g. dx")
      if (size(tmp) /= 2) call program_error("delta array must have size 2")
      this%dx = tmp

      call file%get(sid, "e", this%e)
    end do
  end subroutine

  function refinement_get_dx(this, x) result(dx)
    class(refinement), intent(in) :: this
    real,              intent(in) :: x
    real                          :: dx

    real :: t

    ! FIXME: check
    t = (x - this%x(1)) / (this%x(2) - this%x(1))
    if (this%e == 0) then
      dx = this%dx(1) + (this%dx(2) - this%dx(1)) * t
    else
      dx = this%dx(1) * exp(this%e * t) + (this%dx(2) - this%dx(1) * exp(this%e)) * expm1(this%e * t) / expm1(this%e)
    end if
  end function

  function refinement_get_x(this) result(x)
    class(refinement), intent(in) :: this
    real, allocatable             :: x(:)

    integer           :: i0, i1, di, n
    type(vector_real) :: xvec

    call xvec%init(0, c = 16)

    if (this%dx(1) <= this%dx(2)) then
      i0 = 1
      i1 = 2
      di = 1
    else
      i0 = 2
      i1 = 1
      di = -1
    end if

    call xvec%push(this%x(i0))
    do while (di * xvec%back() < di * this%x(i1))
      call xvec%push(xvec%back() + di * this%get_dx(xvec%back()))
    end do
    n = xvec%n
    x = xvec%to_array()

    ! scale
    x = this%x(i0) + (x - this%x(i0)) * (1 + (this%x(i1) - x(n)) * (x - this%x(i0)) / (x(n) - this%x(i0))**2)

    if (di < 0) x = x(n:1:-1)
  end function

end module
