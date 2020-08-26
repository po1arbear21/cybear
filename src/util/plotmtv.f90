#include "macro.f90.inc"

module plotmtv_m

  use error_m

  implicit none

  private
  public write_plotmtv
  public plotmtv_opts

  type plotset_opts
    ! Plot Titles
    character(:), allocatable :: xlabel, ylabel, zlabel
    character(:), allocatable :: toplabel, subtitle, comment
    logical,      allocatable :: sidelabel
    integer,      allocatable :: sidelabellength, labeloffset

    ! Plot Aspect-ratio/Appearance
    logical, allocatable :: grid, equalscale, fitpage
    real,    allocatable :: xyratio

    ! Axis-Scales
    logical,      allocatable :: xautorange, yautorange, zautorange, autorange
    integer,      allocatable :: xticks, yticks, zticks
    logical,      allocatable :: xlog, ylog, zlog
    real,         allocatable :: xscale, yscale, zscale
    character(:), allocatable :: xticklabel, yticklabel, zticklabel
    logical,      allocatable :: innerticks

    ! Plot Boundaries
    real, allocatable :: xmin, ymin, zmin
    real, allocatable :: xmax, ymax, zmax
    real, allocatable :: vxmin, vymin, vzmin
    real, allocatable :: vxmax, vymax, vzmax

    ! Miscellaneous Options
    logical, allocatable :: xflip, yflip
    logical, allocatable :: xabs, yabs, zabs
    logical, allocatable :: overlay
  contains
    procedure :: write => plotset_opts_write
  end type

  type curve_opts
    ! Line Options
    character(:), allocatable :: linelabel
    integer,      allocatable :: linewidth, linetype, linecolor

    ! Marker Options
    integer, allocatable :: markertype, markercolor, markersize

    ! Fill Options
    integer, allocatable :: filltype, fillcolor
  contains
    procedure :: write => curve_opts_write
  end type

  type view3d_opts
    ! View Point
    real, allocatable :: eyepos_x, eyepos_y, eyepos_z
    real, allocatable :: viewcenter_x, viewcenter_y, viewcenter_z

    ! Axis Options
    real,    allocatable :: window_xmin, window_ymin
    real,    allocatable :: window_xmax, window_ymax
    logical, allocatable :: axislabel, axismove, axisscale
    real,    allocatable :: xaxisscale, yaxisscale, zaxisscale

    ! Miscellaneous Options
    logical, allocatable :: leftworld, hiddenline, paintcube, axisguides
  contains
    procedure :: write => view3d_opts_write
  end type

  type plotmtv_opts
    type(plotset_opts) :: ps
    type(curve_opts)   :: c
    type(view3d_opts)  :: v3

  contains
    procedure, private :: write => plotmtv_opts_write
  end type

  ! internal routines. can write any option, e.g. '% grid = True' or '% xlabel = "my x"'
  interface write_line
    module procedure :: write_line_char
    module procedure :: write_line_int
    module procedure :: write_line_log
    module procedure :: write_line_real
  end interface

contains

  subroutine write_plotmtv(fname, x, y, z, opts)
    !! writes data to file.
    !! options may be specified (otherwise default values from plotmtv are used).

    character(*),       intent(in)           :: fname
      !! file name, e.g. 'output/folder/my_data.asc'
    real,               intent(in)           :: x(:), y(:)
      !! data
    real,               intent(in), optional :: z(:)
      !! z-data is optional
    type(plotmtv_opts), intent(in), optional :: opts
      !! options for plotmtv

    integer :: iounit, ios, i

    ASSERT(size(x) == size(y))
    if (present(z)) then
      ASSERT(size(x) == size(z))
    end if

    call create_parent_dir(fname)

    open (newunit=iounit, file=fname, iostat=ios, action='WRITE')
    if (ios /= 0) call program_error("Error opening file")

    write (iounit, '(A)') '$ DATA=CURVE' // merge('3D', '2D', present(z))
    write (iounit, *)

    if (present(opts)) call opts%write(iounit)

    do i = 1, size(x)
      if (present(z)) then
        write (iounit, *) x(i), y(i), z(i)
      else
        write (iounit, *) x(i), y(i)
      end if
    end do

    close (unit=iounit, iostat=ios)
    if (ios /= 0) call program_error("Error closing file")
  end subroutine

  subroutine plotmtv_opts_write(this, iounit)
    !! write options to file.

    class(plotmtv_opts), intent(in) :: this
    integer,             intent(in) :: iounit
      !! file io unit

    call this%ps%write(iounit)
    call this%c%write( iounit)
    call this%v3%write(iounit)
  end subroutine

  subroutine curve_opts_write(this, iounit)
    !! write curve options to file.

    class(curve_opts), intent(in) :: this
    integer,           intent(in) :: iounit
      !! file io unit

    call write_line(iounit, "linelabel", this%linelabel)
    call write_line(iounit, "linewidth", this%linewidth)
    call write_line(iounit, "linetype",  this%linetype )
    call write_line(iounit, "linecolor", this%linecolor)

    call write_line(iounit, "markertype",  this%markertype )
    call write_line(iounit, "markercolor", this%markercolor)
    call write_line(iounit, "markersize",  this%markersize )

    call write_line(iounit, "filltype",  this%filltype )
    call write_line(iounit, "fillcolor", this%fillcolor)
  end subroutine

  subroutine view3d_opts_write(this, iounit)
    !! write view3d options to file.

    class(view3d_opts), intent(in) :: this
    integer,            intent(in) :: iounit
      !! file io unit

    call write_line(iounit, "eyepos.x",     this%eyepos_x    )
    call write_line(iounit, "eyepos.y",     this%eyepos_y    )
    call write_line(iounit, "eyepos.z",     this%eyepos_z    )
    call write_line(iounit, "viewcenter.x", this%viewcenter_x)
    call write_line(iounit, "viewcenter.y", this%viewcenter_y)
    call write_line(iounit, "viewcenter.z", this%viewcenter_z)

    call write_line(iounit, "window.xmin", this%window_xmin)
    call write_line(iounit, "window.ymin", this%window_ymin)
    call write_line(iounit, "window.xmax", this%window_xmax)
    call write_line(iounit, "window.ymax", this%window_ymax)
    call write_line(iounit, "axislabel",   this%axislabel  )
    call write_line(iounit, "axismove",    this%axismove   )
    call write_line(iounit, "axisscale",   this%axisscale  )
    call write_line(iounit, "xaxisscale",  this%xaxisscale )
    call write_line(iounit, "yaxisscale",  this%yaxisscale )
    call write_line(iounit, "zaxisscale",  this%zaxisscale )

    call write_line(iounit, "leftworld",  this%leftworld )
    call write_line(iounit, "hiddenline", this%hiddenline)
    call write_line(iounit, "paintcube",  this%paintcube )
    call write_line(iounit, "axisguides", this%axisguides)
  end subroutine

  subroutine plotset_opts_write(this, iounit)
    !! write plotset options to file.

    class(plotset_opts), intent(in) :: this
    integer,             intent(in) :: iounit
      !! file io unit

    call write_line(iounit, "xlabel",          this%xlabel         )
    call write_line(iounit, "ylabel",          this%ylabel         )
    call write_line(iounit, "zlabel",          this%zlabel         )
    call write_line(iounit, "toplabel",        this%toplabel       )
    call write_line(iounit, "subtitle",        this%subtitle       )
    call write_line(iounit, "comment",         this%comment        )
    call write_line(iounit, "sidelabel",       this%sidelabel      )
    call write_line(iounit, "sidelabellength", this%sidelabellength)
    call write_line(iounit, "labeloffset",     this%labeloffset    )

    call write_line(iounit, "grid",       this%grid      )
    call write_line(iounit, "equalscale", this%equalscale)
    call write_line(iounit, "fitpage",    this%fitpage   )
    call write_line(iounit, "xyratio",    this%xyratio   )

    call write_line(iounit, "xautorange",  this%xautorange)
    call write_line(iounit, "yautorange",  this%yautorange)
    call write_line(iounit, "zautorange",  this%zautorange)
    call write_line(iounit, "autorange",   this%autorange )
    call write_line(iounit, "xticks",      this%xticks    )
    call write_line(iounit, "yticks",      this%yticks    )
    call write_line(iounit, "zticks",      this%zticks    )
    call write_line(iounit, "xlog",        this%xlog      )
    call write_line(iounit, "ylog",        this%ylog      )
    call write_line(iounit, "zlog",        this%zlog      )
    call write_line(iounit, "xscale",      this%xscale    )
    call write_line(iounit, "yscale",      this%yscale    )
    call write_line(iounit, "zscale",      this%zscale    )
    call write_line(iounit, "xticklabel",  this%xticklabel)
    call write_line(iounit, "yticklabel",  this%yticklabel)
    call write_line(iounit, "zticklabel",  this%zticklabel)
    call write_line(iounit, "innerticks",  this%innerticks)

    call write_line(iounit, "xmin",  this%xmin )
    call write_line(iounit, "ymin",  this%ymin )
    call write_line(iounit, "zmin",  this%zmin )
    call write_line(iounit, "xmax",  this%xmax )
    call write_line(iounit, "ymax",  this%ymax )
    call write_line(iounit, "zmax",  this%zmax )
    call write_line(iounit, "vxmin", this%vxmin)
    call write_line(iounit, "vymin", this%vymin)
    call write_line(iounit, "vzmin", this%vzmin)
    call write_line(iounit, "vxmax", this%vxmax)
    call write_line(iounit, "vymax", this%vymax)
    call write_line(iounit, "vzmax", this%vzmax)

    call write_line(iounit, "xflip",   this%xflip  )
    call write_line(iounit, "yflip",   this%yflip  )
    call write_line(iounit, "xabs",    this%xabs   )
    call write_line(iounit, "yabs",    this%yabs   )
    call write_line(iounit, "zabs",    this%zabs   )
    call write_line(iounit, "overlay", this%overlay)
  end subroutine

  subroutine write_line_char(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    character(:), intent(in), allocatable :: value

    if (allocated(value)) write (iounit, '(A)') '% ' // name // ' = "' // value // '"'
  end subroutine

  subroutine write_line_int(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    integer,      intent(in), allocatable :: value

    if (allocated(value)) write (iounit, '(A, I10)') '% ' // name // ' = ', value
  end subroutine

  subroutine write_line_log(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    logical,      intent(in), allocatable :: value

    if (allocated(value)) write (iounit, '(A)') '% ' // name // ' = ' // merge('True ', 'False', value)
  end subroutine

  subroutine write_line_real(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    real,         intent(in), allocatable :: value

    if (allocated(value)) write (iounit, '(A, E15.5)') '% ' // name // ' = ', value
  end subroutine

  subroutine create_parent_dir(fname)
    !! makes sure that the parent directory exists

    character(*), intent(in) :: fname
      !! file name, e.g. 'output/folder/my_data.asc'

    integer :: i

    i = scan(fname, '/', back=.true.)
    if (i /= 0) call execute_command_line("mkdir -p " // fname(:i))
  end subroutine

end module
