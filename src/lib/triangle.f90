m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_triangle},0,-1))

module triangle_m

  use, intrinsic :: iso_c_binding

  use error_m,        only: program_error
  use grid_m,         only: grid, grid_ptr
  use hashmap_m,      only: hashmap_int
  use triang_grid_m,  only: triang_grid
  use util_m,         only: f2cstring, int2str
  use vector_m,       only: vector_int, vector_real

  implicit none

  private
  public triangulation

  type, bind(c) :: triangulateio
    type(c_ptr)    :: pointlist = c_null_ptr
    type(c_ptr)    :: pointattributelist = c_null_ptr
    type(c_ptr)    :: pointmarkerlist = c_null_ptr

    integer(c_int) :: numberofpoints
    integer(c_int) :: numberofpointattributes

    type(c_ptr)    :: trianglelist = c_null_ptr
    type(c_ptr)    :: triangleattributelist = c_null_ptr
    type(c_ptr)    :: trianglearealist = c_null_ptr
    type(c_ptr)    :: neighborlist = c_null_ptr
    integer(c_int) :: numberoftriangles
    integer(c_int) :: numberofcorners
    integer(c_int) :: numberoftriangleattributes

    type(c_ptr)    :: segmentlist = c_null_ptr
    type(c_ptr)    :: segmentmarkerlist = c_null_ptr
    integer(c_int) :: numberofsegments

    type(c_ptr)    :: holelist       = c_null_ptr
    integer(c_int) :: numberofholes

    type(c_ptr)    :: regionlist     = c_null_ptr
    integer(c_int) :: numberofregions

    type(c_ptr)    :: edgelist       = c_null_ptr
    type(c_ptr)    :: edgemarkerlist = c_null_ptr
    type(c_ptr)    :: normlist       = c_null_ptr
    integer(c_int) :: numberofedges
  end type triangulateio

  interface
    ! void triangulate(char *triswitches, struct triangulateio *in, struct triangulateio *out, struct triangulateio *vorout)
    subroutine triangulate(triswitches, in, out, vorout) bind(c)
      import :: triangulateio, c_char, c_ptr
      character(kind=c_char), intent(in) :: triswitches(*)
      type(c_ptr), value,     intent(in) :: in
      type(c_ptr), value,     intent(in) :: out
      type(c_ptr), value,     intent(in) :: vorout
    end subroutine

    subroutine trifree(memptr) bind(c)
      import
      type(c_ptr),value, intent(in) :: memptr
    end subroutine

  end interface

  type triangulation
    type(vector_real) :: xy
      !! points
    type(vector_int)  :: segments
      !! segments with point index
    type(vector_int)  :: triangles
      !! triangles (point indices)
    type(vector_real) :: holes
      !! hole points
    type(hashmap_int) :: hmap
      !! hash map for segments

  contains
    procedure :: init        => triangulation_init
    procedure :: add_point   => triangulation_add_point
    procedure :: add_segment => triangulation_add_segment
    procedure :: add_polygon => triangulation_add_polygon
    procedure :: triangulate => triangulation_triangulate
    procedure :: get_grid    => triangulation_get_grid


  end type

contains

  subroutine triangulation_init(this, c_xy, c_segments)
    !! initialize triangulation
    class(triangulation), intent(out) :: this
    integer, optional,    intent(in)  :: c_xy, c_segments
      !! initial capacity

    integer :: c_xy_, c_segments_

    ! deal with optional inputs
    c_xy_ = 32
    if (present(c_xy)) c_xy_ = c_xy
    c_segments_ = 32
    if (present(c_segments)) c_segments_ = c_segments

    ! init vectors
    call this%xy%init(0, c = c_xy_)
    call this%segments%init(0, c = c_segments_)
    call this%triangles%init(0, c = 32)
    call this%holes%init(0, c=32)
    call this%hmap%init()
  end subroutine

  function triangulation_add_point(this, x, y, fast) result(i)
    !! append new point(x,y) to end of vector xy
    class(triangulation), intent(inout) :: this
    real,                 intent(in)    :: x, y
      !! point coordinates
    integer                             :: i
      !! index of the new point
    logical, optional,    intent(in)    :: fast
      !! speed up inclusion of point by not checking whether a point is part of a segment
      !! default: false

    integer :: j, k, ix, iy
    logical :: status, fast_
    real    :: x1, y1, x2, y2

    real, parameter :: TOL = 1e-10

    fast_ = .false.
    if (present(fast)) fast_ = fast

    ! check if the point already exists in xy
    ix = nint(x*1e10)
    iy = nint(y*1e10)
    call this%hmap%get([ix, iy], j, status=status)
    if (status) then
      i = (j + 1) / 2
      return
    end if

    ! check if the point is part of an existing segment with tolerance
    if (.not. fast_) then
      do j = 1, this%segments%n, 2
        x1 = this%xy%d(2 * this%segments%d(j) - 1)
        y1 = this%xy%d(2 * this%segments%d(j))
        x2 = this%xy%d(2 * this%segments%d(j+1) - 1)
        y2 = this%xy%d(2 * this%segments%d(j+1))

        ! check if the point lies on the segment
        if (is_pnt_on_segment(x1, y1, x2, y2, x, y, TOL)) then
            ! add the new point (x, y)
            call this%xy%push(x)
            call this%xy%push(y)

            ! get the index of the new point
            i = this%xy%n / 2

            ! replace the current segment with two new segments
            k = this%segments%d(j+1)
            this%segments%d(j+1) = i  ! update the second endpoint of the original segment

            ! add a new segment with the new point and the second original endpoint
            call this%segments%push(i)
            call this%segments%push(k)

            return
        end if
      end do
  end if

    ! add the new point (x, y) if it doesn't exist and isn't on a segment
    call this%xy%push(x)
    call this%xy%push(y)
    call this%hmap%set([ix, iy], this%xy%n)

    ! return the index of the new point
    i = this%xy%n / 2
  end function

  logical function is_pnt_on_segment(x1, y1, x2, y2, x, y, tol) result(res)
    !! check if point (x, y) is on the segment (x1, y1) - (x2, y2) with tolerance tol
    real, intent(in) :: x1, y1, x2, y2
      !! segment: (x1, y1) -- (x2, y2)
    real, intent(in) :: x, y
      !! point: (x, y)
    real, intent(in) :: tol
      !! absolute tolerance

    real :: cross, cdot, len2

    ! check if point (x, y) is on the segment (x1, y1) - (x2, y2) with tolerance tol
    cross = (y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)
    if (abs(cross) > tol) then
        res = .false.
        return
    end if

    cdot = (x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)
    if (cdot < 0.0) then
        res = .false.
        return
    end if

    len2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
    if (cdot > len2) then
        res = .false.
        return
    end if

    res = .true.
  end function

  subroutine triangulation_add_segment(this, i1, i2)
    !! append new segment to end of vector segments
    class(triangulation), intent(inout) :: this
    integer,              intent(in)    :: i1, i2
      !! endpoints index of the new segment

    integer :: n
    logical, allocatable :: exist(:)

    ! check if the segment already exists
    n = this%segments%n
    exist = (i1 == this%segments%d(1:n-1:2) .and. i2 == this%segments%d(2:n:2))
    if (any(exist)) return
    exist = (i2 == this%segments%d(1:n-1:2) .and. i1 == this%segments%d(2:n:2))
    if (any(exist)) return

    ! push segment
    call this%segments%push(i1)
    call this%segments%push(i2)
  end subroutine

  subroutine triangulation_add_polygon(this, x, y, closed, hole, fast)
    !! append polygon to vector of points and segments
    class(triangulation), intent(inout) :: this
    real,                 intent(in)    :: x(:)
      !! x coordinate of points
    real,                 intent(in)    :: y(:)
      !! y coordinate of points
    logical, optional,    intent(in)    :: closed
      !! true, enclosure
    real,    optional,    intent(in)    :: hole(:)
      !! hole points
    logical, optional,    intent(in)    :: fast
      !! speed up inclusion of point by not checking whether a point is part of a segment
      !! default: false

    integer              :: i, n
    logical              :: closed_, hole_
    integer, allocatable :: p(:)

    ! number of points
    n = size(x)
    ! indices of points
    allocate(p(n))

    p(1) = this%add_point(x(1), y(1), fast=fast)
    do i = 2, n
      ! get index of the point
      p(i) = this%add_point(x(i), y(i), fast=fast)
      ! add segment
      call this%add_segment(p(i-1), p(i))
    end do

    closed_ = .false.
    if (present(closed)) closed_ = closed
    if (closed_) then
      call this%add_segment(p(n), p(1))
    end if

    hole_ = .false.
    if (present(hole)) hole_ = .true.
    if (hole_) then
      if (.not. closed_) call program_error("hole only possible for closed")
      call this%holes%push(hole(1))
      call this%holes%push(hole(2))
    end if
  end subroutine

  subroutine triangulation_triangulate(this, quiet, max_area, wpoly, max_steiner_pnts)
    class(triangulation), intent(inout) :: this
    logical, optional,    intent(in)    :: quiet
      !! if true -> no output. Default: true
    real,    optional,    intent(in)    :: max_area
      !! maximum area constraint
    logical, optional,    intent(in)    :: wpoly
      !! if true -> polygons are to be considered . Default: false
    integer, optional,    intent(in)    :: max_steiner_pnts
      !! maximum number of steiner points that may be added by Triangle. Default: Infinity

    character(32) :: str, buf
    logical       :: quiet_, wpoly_

    character(kind=c_char),   allocatable :: triswitches(:)
    real(     kind=c_double), pointer     :: pointlist(:)
    integer(  kind=c_int),    pointer     :: segmentlist(:), trianglelist(:)
    real(     kind=c_double), pointer     :: holelist(:)
    type(triangulateio),      target      :: c_in, c_out

    integer :: i
    real    :: xmin, xmax, ymin, ymax

    nullify (pointlist, segmentlist, trianglelist, holelist)

    xmin = minval(this%xy%d(1:this%xy%n-1:2))
    xmax = maxval(this%xy%d(1:this%xy%n-1:2))
    ymin = minval(this%xy%d(2:this%xy%n:2))
    ymax = maxval(this%xy%d(2:this%xy%n:2))

    ! set command line switches
    str = ""
    quiet_ = .true.
    wpoly_ = .true.
    if (present(quiet)) quiet_ = quiet
    if (present(wpoly)) wpoly_ = wpoly
    if (wpoly_) str = "pqD"
    if (quiet_) str = trim(str)//"Q"
    if (present(max_area)) then
      write (buf, "(A,F7.5)") "a", max_area / ((xmax - xmin) * (ymax - ymin))
      str = trim(str)//trim(buf)
    end if
    if (present(max_steiner_pnts)) str = trim(str)//"S"// int2str(max_steiner_pnts)
    triswitches = f2cstring(trim(str))

    ! set points
    allocate (pointlist(this%xy%n))
    do i = 1, this%xy%n/2
      pointlist(2*i-1) = (this%xy%d(2*i-1) - xmin) / (xmax - xmin)
      pointlist(2*i  ) = (this%xy%d(2*i  ) - ymin) / (ymax - ymin)
    end do
    c_in%pointlist               = c_loc(pointlist)
    c_in%numberofpoints          = int(this%xy%n / 2, kind = c_int)
    c_in%numberofpointattributes = 0

    ! set segmentlist
    if (this%segments%n > 0) then
      allocate (segmentlist(this%segments%n))
      do i=1, this%segments%n
        segmentlist(i) = int(this%segments%d(i), kind = c_int)
      end do
      c_in%segmentlist = c_loc(segmentlist)
    end if
    c_in%numberofsegments = int(this%segments%n / 2, kind = c_int)

    ! set holelist
    if (this%holes%n > 0) then
      allocate (holelist(this%holes%n))
      do i=1, this%holes%n/2
        holelist(2*i-1) = (this%holes%d(2*i-1) - xmin) / (xmax - xmin)
        holelist(2*i  ) = (this%holes%d(2*i  ) - ymin) / (ymax - ymin)
      end do
      c_in%holelist = c_loc(holelist)
    end if
    c_in%numberofholes = int(this%holes%n / 2, kind = c_int)

    ! regions
    c_in%numberofregions = 0
    call triangulate(triswitches, c_loc(c_in), c_loc(c_out), c_null_ptr)

    ! deallocate pointlist ....
    if (associated(pointlist)) deallocate(pointlist)
    if (associated(segmentlist)) deallocate(segmentlist)
    if (associated(holelist)) deallocate(holelist)

    ! extract points from c_out
    call this%xy%reset()
    call c_f_pointer(c_out%pointlist, pointlist, shape = [2*c_out%numberofpoints])
    do i = 1, c_out%numberofpoints
      call this%xy%push(pointlist(2*i-1)*(xmax - xmin) + xmin)
      call this%xy%push(pointlist(2*i  )*(ymax - ymin) + ymin)
    end do

    ! extract segements from c_out
    call this%segments%reset()
    call c_f_pointer(c_out%segmentlist, segmentlist, shape = [2*c_out%numberofsegments])
    if (associated(segmentlist)) then
      do i = 1, c_out%numberofsegments
        call this%segments%push(int(segmentlist(2*i-1)))
        call this%segments%push(int(segmentlist(2*i  )))
      end do
    end if

    ! extract triangles from c_out
    call this%triangles%reset()
    call c_f_pointer(c_out%trianglelist, trianglelist, shape = [3*c_out%numberoftriangles])
    do i = 1, c_out%numberoftriangles
      call this%triangles%push(int(trianglelist(3*i-2)))
      call this%triangles%push(int(trianglelist(3*i-1)))
      call this%triangles%push(int(trianglelist(3*i  )))
    end do

    ! trifree
    if (c_associated(c_out%pointlist)) call trifree(c_out%pointlist)
    if (c_associated(c_out%pointattributelist)) call trifree(c_out%pointattributelist)
    if (c_associated(c_out%pointmarkerlist)) call trifree(c_out%pointmarkerlist)
    if (c_associated(c_out%trianglelist)) call trifree(c_out%trianglelist)
    if (c_associated(c_out%triangleattributelist)) call trifree(c_out%triangleattributelist)
    if (c_associated(c_out%trianglearealist)) call trifree(c_out%trianglearealist)
    if (c_associated(c_out%neighborlist)) call trifree(c_out%neighborlist)
    if (c_associated(c_out%segmentlist)) call trifree(c_out%segmentlist)
    if (c_associated(c_out%segmentmarkerlist)) call trifree(c_out%segmentmarkerlist)
    if (c_associated(c_out%edgelist)) call trifree(c_out%edgelist)
    if (c_associated(c_out%edgemarkerlist)) call trifree(c_out%edgemarkerlist)
    if (c_associated(c_out%normlist)) call trifree(c_out%normlist)
  end subroutine

  subroutine triangulation_get_grid(this, name, g)
    class(triangulation), intent(inout) :: this
    character(*),         intent(in)    :: name
      !! grid name
    type(triang_grid),    intent(out)   :: g
      !! return triangle grid

    integer :: nvert, ncell

    ! get number of points
    nvert = this%xy%n / 2

    ! get number of triangles
    ncell = this%triangles%n / 3

    ! initialize triangle grid
    call g%init(name, reshape(this%xy%d(1:this%xy%n), [2, nvert]), reshape(this%triangles%d(1: this%triangles%n), [3, ncell]))
  end subroutine

end module

m4_divert(0)
