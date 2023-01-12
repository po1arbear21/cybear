m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_triangle},0,-1))

module triangle_m

  use, intrinsic :: iso_c_binding

  use error_m,        only: program_error
  use util_m,         only: f2cstring
  use vector_m,       only: vector_int, vector_real
  use grid_m,         only: grid, grid_ptr
  use triang_grid_m,  only: triang_grid

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

  contains
    procedure :: init        => triangulation_init
    procedure :: add_point   => triangulation_add_point
    procedure :: add_segment => triangulation_add_segment
    procedure :: add_polygon => triangulation_add_polygon
    procedure :: triangulate => triangulation_triangulate
    ! procedure :: refine      => triangulation_refine
    procedure :: get_grid    => triangulation_get_grid


  end type

contains

  subroutine triangulation_init(this)
    !! initialize triangulation
    class(triangulation), intent(out) :: this

    ! init vectors
    call this%xy%init(0, c = 32)
    call this%segments%init(0, c = 32)
    call this%triangles%init(0, c = 32)
    call this%holes%init(0, c=32)
  end subroutine

  function triangulation_add_point(this, x, y) result(i)
    !! append new point(x,y) to end of vector xy
    class(triangulation), intent(inout) :: this
    real,                 intent(in)    :: x, y
      !! new point(x,y)
    integer                             :: i
      !! return the index of the new point

    call this%xy%push(x)
    call this%xy%push(y)
    ! index
    i = this%xy%n / 2
  end function

  subroutine triangulation_add_segment(this, i1, i2)
    !! append new segment to end of vector segments
    class(triangulation), intent(inout) :: this
    integer,              intent(in)    :: i1, i2
     !! endpoints index of the new segment

    call this%segments%push(i1)
    call this%segments%push(i2)
  end subroutine

  subroutine triangulation_add_polygon(this, x, y, closed, hole)
    class(triangulation), intent(inout) :: this
    real,                 intent(in)    :: x(:)
      !! x coordinate of points
    real,                 intent(in)    :: y(:)
      !! y coordinate of points
    logical, optional,    intent(in)    :: closed
      !! true, enclosure
    real,    optional,    intent(in)    :: hole(:)
      !! hole points

    integer              :: i, n
    logical              :: closed_, hole_
    integer, allocatable :: p(:)

    ! number of points
    n = size(x)
    ! indices of points
    allocate(p(n))

    p(1) = this%add_point(x(1), y(1))
    do i = 2, n
      ! get index of the point
      p(i) = this%add_point(x(i), y(i))
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

  subroutine triangulation_triangulate(this, quiet, max_area)
    class(triangulation), intent(inout) :: this
    logical, optional,    intent(in)    :: quiet
    real,    optional,    intent(in)    :: max_area

    character(32) :: str, buf
    logical       :: quiet_

    character(kind=c_char),   allocatable :: triswitches(:)
    real(     kind=c_double), pointer     :: pointlist(:)
    integer(  kind=c_int),    pointer     :: segmentlist(:), trianglelist(:)
    real(     kind=c_double), pointer     :: holelist(:)
    type(triangulateio),      target      :: c_in, c_out

    integer :: i
    real    :: xmin, xmax, ymin, ymax

    xmin = minval(this%xy%d(1:this%xy%n-1:2))
    xmax = maxval(this%xy%d(1:this%xy%n-1:2))
    ymin = minval(this%xy%d(2:this%xy%n:2))
    ymax = maxval(this%xy%d(2:this%xy%n:2))

    str = "pqD"
    quiet_ = .true.
    if (present(quiet)) quiet_ = quiet
    if (quiet_) str = trim(str)//"Q"
    if (present(max_area)) then
      write (buf, "(A,F7.5)") "a", max_area / ((xmax - xmin) * (ymax - ymin))
      str = trim(str)//trim(buf)
    end if
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
    deallocate (pointlist, segmentlist, holelist)

    ! extract triangulation from c_out
    call this%xy%reset()
    call c_f_pointer(c_out%pointlist, pointlist, shape = [2*c_out%numberofpoints])
    do i = 1, c_out%numberofpoints
      call this%xy%push(pointlist(2*i-1)*(xmax - xmin) + xmin)
      call this%xy%push(pointlist(2*i  )*(ymax - ymin) + ymin)
    end do
    call this%segments%reset()
    call c_f_pointer(c_out%segmentlist, segmentlist, shape = [2*c_out%numberofsegments])
    do i = 1, c_out%numberofsegments
      call this%segments%push(int(segmentlist(2*i-1)))
      call this%segments%push(int(segmentlist(2*i  )))
    end do
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
