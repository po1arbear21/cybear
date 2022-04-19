module plotmtv_m

  use error_m, only: program_error

  implicit none

  private
  public plotmtv_write
  public plotmtv
  public plotset_options
  public curve_options
  public view3d_options

  type plotset_options
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
    procedure :: write => plotset_options_write
  end type

  type curve_options
    ! Line Options
    character(:), allocatable :: linelabel
    integer,      allocatable :: linewidth, linetype, linecolor

    ! Marker Options
    integer, allocatable :: markertype, markercolor, markersize

    ! Fill Options
    integer, allocatable :: filltype, fillcolor
  contains
    procedure :: write => curve_options_write
  end type

  type view3d_options
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
    procedure :: write => view3d_options_write
  end type

  type plotmtv
    private
    logical :: file_open
      !! is a file opened?
    integer :: iounit
      !! file unit
    logical :: three_dim
      !! is it a 3d (true) or 2d (false) plot?

  contains
    procedure :: init         => plotmtv_init
    procedure :: write_header => plotmtv_write_header
    procedure :: write_curve  => plotmtv_write_curve
    procedure :: close        => plotmtv_close
  end type

  ! internal routines. can write any option, e.g. '% grid = True' or '% xlabel = "my x"'
  interface write_line
    module procedure :: write_line_char
    module procedure :: write_line_int
    module procedure :: write_line_log
    module procedure :: write_line_real
  end interface

contains

  subroutine plotmtv_write(fname, x, y, plotset_opts, curve_opts)
    !! simple 2d write routine.
    !! wrapper around more involved plotmtv type.

    character(*),                    intent(in) :: fname
      !! file name, e.g. "output/folder/my_data.asc"
    real,                            intent(in) :: x(:), y(:)
      !! data arrays (of same length!)
    type(plotset_options), optional, intent(in) :: plotset_opts
      !! plotset options
    type(curve_options),   optional, intent(in) :: curve_opts
      !! curve options

    type(plotmtv) :: p

    call p%init(fname)
    call p%write_header(plotset_opts = plotset_opts)
    call p%write_curve(x, y, opts = curve_opts)
    call p%close()
  end subroutine

  subroutine plotmtv_init(this, fname)
    !! inits plotmtv. creates file handler

    class(plotmtv), intent(out) :: this
    character(*),   intent(in)  :: fname
      !! file name, e.g. "output/folder/my_data.asc"

    integer :: i, ios

    i = scan(fname, "/", back = .true.)
    if (i /= 0) call execute_command_line("mkdir -p " // fname(:i))

    open (newunit = this%iounit, file = fname, iostat = ios, action = "WRITE")
    if (ios /= 0) call program_error("Error opening file")
    this%file_open = .true.
  end subroutine

  subroutine plotmtv_close(this)
    !! closes file handler

    class(plotmtv), intent(inout) :: this

    integer :: ios

    if (this%file_open) then
      close (unit = this%iounit, iostat = ios)
      if (ios /= 0) call program_error("Error closing file")

      this%file_open = .false.
    else
      print *, "Cannot close a file as none is open!"
    end if
  end subroutine

  subroutine plotmtv_write_header(this, three_dim, plotset_opts, view3d_opts, gl_curve_opts)
    class(plotmtv),        intent(inout)          :: this
    logical,               intent(in),   optional :: three_dim
      !! is it 3d (true) or 2d (false) plot?
      !! default: 2d
    type(plotset_options), intent(in),   optional :: plotset_opts
      !! plotset options
    type(view3d_options),  intent(in),   optional :: view3d_opts
      !! 3D view options
    type(curve_options),   intent(in),   optional :: gl_curve_opts
      !! global curve options

    if (.not. this%file_open) call program_error("call init beforehand to open a file!")

    ! 3d plot?
    this%three_dim = .false.
    if (present(three_dim)) this%three_dim = three_dim

    ! view3d only if 3d plot
    if (present(view3d_opts) .and. (.not. this%three_dim)) call program_error("view3d options mustnt be supplied for 3d plot!")

    write (this%iounit, "(A)") "$ DATA=CURVE" // merge("3D", "2D", this%three_dim)
    write (this%iounit, *)

    if (present(plotset_opts) ) call plotset_opts%write( this%iounit        )
    if (present(view3d_opts)  ) call view3d_opts%write(  this%iounit        )
    if (present(gl_curve_opts)) call gl_curve_opts%write(this%iounit, .true.)
  end subroutine

  subroutine plotmtv_write_curve(this, x, y, z, opts)
    !! write data for a curve, consisting of 2 or three data points per line.

    class(plotmtv),      intent(in)           :: this
    real,                intent(in)           :: x(:), y(:)
    real,                intent(in), optional :: z(:)
    type(curve_options), intent(in), optional :: opts

    integer :: i

    if (.not. this%file_open) call program_error("call init beforehand to open a file!")

    ! check header/curve dimensions are same
    if (this%three_dim .neqv. present(z)) call program_error("header specified a 2d/3d plot but curve is of opposite dimension!")

    ! check data have same lengths
    if (size(x) /= size(y)) call program_error("x, y arrays are of different lengths!")
    if (present(z)) then
      if (size(x) /= size(z)) call program_error("x, z arrays are of different lengths!")
    end if

    if (present(opts)) call opts%write(this%iounit, .false.)

    do i = 1, size(x)
      if (this%three_dim) then
        write (this%iounit, "(3ES24.16)") x(i), y(i), z(i)
      else
        write (this%iounit, "(2ES24.16)") x(i), y(i)
      end if
    end do

    write (this%iounit, *)
  end subroutine

  subroutine curve_options_write(this, iounit, global)
    !! write curve options to file.

    class(curve_options), intent(in) :: this
    integer,              intent(in) :: iounit
      !! file io unit
    logical,              intent(in) :: global
      !! is this a global curve option?

    character(:), allocatable :: gl_cstr

    ! preprend a "d" for global options
    allocate (character(0) :: gl_cstr)      ! remove gfortran warning
    gl_cstr = ""
    if (global) gl_cstr = "d"

    ! linelabel only for non-global curve options
    if (global) then
      if (allocated(this%linelabel)) call program_error("global curve options mustnt have a linelabel!")
    else
      call write_line(iounit, gl_cstr//"linelabel", this%linelabel)
    end if

    call write_line(iounit, gl_cstr//"linewidth", this%linewidth)
    call write_line(iounit, gl_cstr//"linetype",  this%linetype )
    call write_line(iounit, gl_cstr//"linecolor", this%linecolor)

    call write_line(iounit, gl_cstr//"markertype",  this%markertype )
    call write_line(iounit, gl_cstr//"markercolor", this%markercolor)
    call write_line(iounit, gl_cstr//"markersize",  this%markersize )

    call write_line(iounit, gl_cstr//"filltype",  this%filltype )
    call write_line(iounit, gl_cstr//"fillcolor", this%fillcolor)

    if (global) write (iounit, *)
  end subroutine

  subroutine view3d_options_write(this, iounit)
    !! write view3d options to file.

    class(view3d_options), intent(in) :: this
    integer,               intent(in) :: iounit
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

    write (iounit, *)
  end subroutine

  subroutine plotset_options_write(this, iounit)
    !! write plotset options to file.

    class(plotset_options), intent(in) :: this
    integer,                intent(in) :: iounit
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

    write (iounit, *)
  end subroutine

  subroutine write_line_char(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    character(:), intent(in), allocatable :: value

    if (allocated(value)) write (iounit, "(A)") "% " // name // ' = "' // value // '"'
  end subroutine

  subroutine write_line_int(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    integer,      intent(in), allocatable :: value

    if (allocated(value)) write (iounit, "(A,I0)") "% " // name // " = ", value
  end subroutine

  subroutine write_line_log(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    logical,      intent(in), allocatable :: value

    if (allocated(value)) write (iounit, "(A)") "% " // name // " = " // merge("True ", "False", value)
  end subroutine

  subroutine write_line_real(iounit, name, value)
    integer,      intent(in)              :: iounit
    character(*), intent(in)              :: name
    real,         intent(in), allocatable :: value

    if (allocated(value)) write (iounit, "(A,ES24.16)") "% " // name // " = ", value
  end subroutine

end module
