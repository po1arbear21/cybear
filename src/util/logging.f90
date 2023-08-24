m4_include(macro.f90.inc)

module logging_m

  use color_m,         only: COL_DEFAULT, COL_YELLOW, COL_MAGENTA
  use error_m,         only: program_error
  use iso_fortran_env, only: error_unit, output_unit
  use util_m,          only: get_memory_usage, is_letter, split_folder_file

  implicit none

  private
  public init_logging, finish_logging, add_level, logging
  public LVL_DEBUG, LVL_INFO, LVL_WARN, LVL_ERROR

  type datetime
    integer :: year  = 0
    integer :: month = 0
    integer :: day   = 0
    integer :: zone  = 0
    integer :: hour  = 0
    integer :: min   = 0
    integer :: sec   = 0
    integer :: msec  = 0
  contains
    procedure :: now => datetime_now
  end type

  integer, parameter :: TK_LIT    = 0
  integer, parameter :: TK_LVL    = 1
  integer, parameter :: TK_YEAR4  = 2
  integer, parameter :: TK_YEAR2  = 3
  integer, parameter :: TK_MONTH  = 4
  integer, parameter :: TK_DAY    = 5
  integer, parameter :: TK_ZONE   = 6
  integer, parameter :: TK_HOUR   = 7
  integer, parameter :: TK_MIN    = 8
  integer, parameter :: TK_SEC    = 9
  integer, parameter :: TK_MSEC3  = 10
  integer, parameter :: TK_MSEC2  = 11
  integer, parameter :: TK_MSEC1  = 12
  integer, parameter :: TK_FOLDER = 13
  integer, parameter :: TK_FILE   = 14
  integer, parameter :: TK_LINE   = 15
  integer, parameter :: TK_MEM    = 16
  integer, parameter :: TK_MSG    = 17

  type token
    integer :: istart
      !! start index in fmt
    integer :: istop
      !! stop index in fmt
    integer :: type
  end type

  m4_define({T},{token})
  m4_include(vector_def.f90.inc)

  type level
    character(:), allocatable :: name
      !! level name (e.g. "DEBUG", "INFO", ...)
    character(:), allocatable :: color
      !! color used when printing name on console
    character(:), allocatable :: fmt
      !! format string
    logical                   :: record_date
      !! record date
    logical                   :: record_time
      !! record time
    logical                   :: record_mem
      !! record memory usage
    logical                   :: print
      !! print to console
    logical                   :: error
      !! write to error unit?

    type(vector_token) :: tokens
      !! list of tokens in format string
  end type

  m4_define({T},{level})
  m4_include(vector_def.f90.inc)

  type(vector_level) :: levels

  integer, protected :: LVL_DEBUG = 0
  integer, protected :: LVL_INFO  = 0
  integer, protected :: LVL_WARN  = 0
  integer, protected :: LVL_ERROR = 0

  type entry
    integer                   :: ilvl = 0
      !! logging level index
    type(datetime)            :: time
      !! date and time
    character(:), allocatable :: folder
      !! source code folder
    character(:), allocatable :: file
      !! source code filename
    integer                   :: line = 0
      !! source code line number
    real                      :: mem = 0.0
      !! memory usage
    character(:), allocatable :: msg
      !! message
  contains
    procedure :: write => entry_write
  end type

  integer :: funit = 0

contains

  subroutine init_logging(file, record_date, record_time, record_mem, verbosity, fmt)
    !! initialize logging, create default levels
    character(*), optional, intent(in) :: file
      !! optional: output to log file
    logical,      optional, intent(in) :: record_date
      !! optional: record date for default levels? (default: true)
    logical,      optional, intent(in) :: record_time
      !! optional: record time for default levels? (default: true)
    logical,      optional, intent(in) :: record_mem
      !! optional: record memory usage for default levels? (default: false)
    integer,      optional, intent(in) :: verbosity
      !! optional: verbosity of console printing (default: 1) (1: all, 2: info,warn,error, 3: warn,error, 4: error)
    character(*), optional, intent(in) :: fmt
      !! optional: custom format string used for default levels

    character(:), allocatable :: fmt_
    integer                   :: iostat, verb
    logical                   :: record_date_, record_time_, record_mem_

    if (present(file)) then
      open (newunit = funit, file = file, status = "replace", action = "write", iostat = iostat)
      if (iostat /= 0) call program_error("unable to open file "//file)
    else
      funit = 0
    end if

    record_date_ = .true.
    if (present(record_date)) record_date_ = record_date

    record_time_ = .true.
    if (present(record_time)) record_time_ = record_time

    record_mem_ = .false.
    if (present(record_mem)) record_mem_ = record_mem

    verb = 1
    if (present(verbosity)) verb = verbosity

    if (present(fmt)) then
      fmt_ = fmt
    else
      fmt_ = "[%LVL] "
      if (record_date_) fmt_ = fmt_ // "%YYYY-%mm-%DD "
      if (record_time_) fmt_ = fmt_ // "%HH:%MM:%SS "
      fmt_ = fmt_ // "| %FILE:%LINE | "
      if (record_mem_) fmt_ = fmt_ // "%MEM GiB | "
      fmt_ = fmt_ // "%MSG"
    end if

    call levels%init(0, c = 8)
    call entries%init(0, c = 32)

    LVL_DEBUG = add_level("DEBUG",   COL_DEFAULT, fmt_, record_date_, record_time_, record_mem_, (verb <= 1), .false.)
    LVL_INFO  = add_level("INFO",    COL_DEFAULT, fmt_, record_date_, record_time_, record_mem_, (verb <= 2), .false.)
    LVL_WARN  = add_level("WARNING", COL_YELLOW,  fmt_, record_date_, record_time_, record_mem_, (verb <= 3), .false.)
    LVL_ERROR = add_level("ERROR",   COL_MAGENTA, fmt_, record_date_, record_time_, record_mem_, (verb <= 4), .true.)
  end subroutine

  subroutine finish_logging()
    !! finish logging, close logfile

    call levels%resize(0)
    call entries%resize(0)

    close (funit)
    funit = 0
  end subroutine

  function add_level(name, color, fmt, record_date, record_time, record_mem, print, error) result(ilvl)
    !! add new logging level
    character(*), intent(in) :: name
      !! level name
    character(*), intent(in) :: color
      !! console color
    character(*), intent(in) :: fmt
      !! format string
    logical,      intent(in) :: record_date
      !! record date?
    logical,      intent(in) :: record_time
      !! record time?
    logical,      intent(in) :: record_mem
      !! record memory usage?
    logical,      intent(in) :: print
      !! print to console?
    logical,      intent(in) :: error
      !! write to error unit?
    integer                  :: ilvl
      !! return index of new level

    integer     :: i, len
    type(level) :: lvl
    type(token) :: tk

    ! set level members
    lvl%name        = name
    lvl%color       = color
    lvl%fmt         = fmt
    lvl%record_date = record_date
    lvl%record_time = record_time
    lvl%record_mem  = record_mem
    lvl%print       = print
    lvl%error       = error

    ! get list of tokens
    call lvl%tokens%init(0, c = 32)
    len = len_trim(fmt)
    tk%istop = 0
    do while (tk%istop < len)
      tk%istart = tk%istop + 1
      if (fmt(tk%istart:tk%istart) == "%") then
        ! search for next non-letter character => end of token
        do i = tk%istart+1, len
          if (.not. is_letter(fmt(i:i))) exit
        end do
        tk%istop = i - 1

        select case (fmt(tk%istart:tk%istop))
        case ("%LVL")
          tk%type = TK_LVL
        case ("%YYYY")
          tk%type = TK_YEAR4
        case ("%YY")
          tk%type = TK_YEAR2
        case ("%mm")
          tk%type = TK_MONTH
        case ("%DD")
          tk%type = TK_DAY
        case ("%ZZZZ")
          tk%type = TK_ZONE
        case ("%HH")
          tk%type = TK_HOUR
        case ("%MM")
          tk%type = TK_MIN
        case ("%SS")
          tk%type = TK_SEC
        case ("%sss")
          tk%type = TK_MSEC3
        case ("%ss")
          tk%type = TK_MSEC2
        case ("%s")
          tk%type = TK_MSEC1
        case ("%FOLDER")
          tk%type = TK_FOLDER
        case ("%FILE")
          tk%type = TK_FILE
        case ("%LINE")
          tk%type = TK_LINE
        case ("%MEM")
          tk%type = TK_MEM
        case ("%MSG")
          tk%type = TK_MSG
        case default
          tk%type = TK_LIT
        end select
      else
        ! string literal until start of next token
        do i = tk%istart, len
          if (fmt(i:i) == "%") exit
        end do
        tk%istop = i - 1
        tk%type = TK_LIT
      end if

      call lvl%tokens%push(tk)
    end do

    ! save level in global vector
    call levels%push(lvl)

    ! return index
    ilvl = levels%n
  end function

  subroutine logging(ilvl, msg, file, line)
    integer,                intent(in) :: ilvl
    character(*),           intent(in) :: msg
    character(*), optional, intent(in) :: file
    integer,      optional, intent(in) :: line

    integer     :: ilvl_
    type(entry) :: e

    ilvl_ = ilvl
    if (ilvl == 0) ilvl_ = LVL_INFO

    associate (lvl => levels%d(ilvl_))
      e%ilvl = ilvl_
      if (lvl%record_date .or. lvl%record_time) call e%time%now()
      if (lvl%record_mem) e%mem = get_memory_usage()
      if (present(file)) then
        call split_folder_file(file, e%folder, e%file)
      end if
      if (present(line)) e%line = line
      e%msg = msg
      if (lvl%print) then
        if (lvl%error) then
          call e%write(error_unit)
        else
          call e%write(output_unit)
        end if
      end if
      if (funit /= 0) then
        call e%write(funit)
        flush (funit)
      end if
    end associate
  end subroutine

  m4_define({T},{token})
  m4_include(vector_imp.f90.inc)

  m4_define({T},{level})
  m4_include(vector_imp.f90.inc)

  subroutine datetime_now(this)
    !! set to current date and time
    class(datetime), intent(out) :: this

    integer :: values(8)

    ! get current date and time
    call date_and_time(values = values)

    ! extract values
    this%year  = values(1)
    this%month = values(2)
    this%day   = values(3)
    this%zone  = values(4)
    this%hour  = values(5)
    this%min   = values(6)
    this%sec   = values(7)
    this%msec  = values(8)
  end subroutine

  subroutine entry_write(this, unit)
    class(entry), intent(in) :: this
    integer,      intent(in) :: unit

    integer :: i

    associate (lvl => levels%d(this%ilvl))
      do i = 1, lvl%tokens%n
        associate (tk => lvl%tokens%d(i))
          select case (tk%type)
          case (TK_LIT)
            write (unit, "(A)", advance = "no") lvl%fmt(tk%istart:tk%istop)
          case (TK_LVL)
            if ((unit == error_unit) .or. (unit == output_unit)) then
              write (unit, "(A)", advance = "no") lvl%color
            end if
            write (unit, "(A)", advance = "no") lvl%name
            if ((unit == error_unit) .or. (unit == output_unit)) then
              write (unit, "(A)", advance = "no") COL_DEFAULT
            end if
          case (TK_YEAR4)
            write (unit, "(I4)", advance = "no") this%time%year
          case (TK_YEAR2)
            write (unit, "(I2)", advance = "no") mod(this%time%year, 100)
          case (TK_MONTH)
            write (unit, "(I2.2)", advance = "no") this%time%month
          case (TK_DAY)
            write (unit, "(I2.2)", advance = "no") this%time%day
          case (TK_ZONE)
            write (unit, "(I4.4)", advance = "no") this%time%zone
          case (TK_HOUR)
            write (unit, "(I2.2)", advance = "no") this%time%hour
          case (TK_MIN)
            write (unit, "(I2.2)", advance = "no") this%time%min
          case (TK_SEC)
            write (unit, "(I2.2)", advance = "no") this%time%sec
          case (TK_MSEC3)
            write (unit, "(I3.3)", advance = "no") this%time%msec
          case (TK_MSEC2)
            write (unit, "(I2.2)", advance = "no") this%time%msec / 10
          case (TK_MSEC1)
            write (unit, "(I1.1)", advance = "no") this%time%msec / 100
          case (TK_FOLDER)
            write (unit, "(A)", advance = "no") this%folder
          case (TK_FILE)
            write (unit, "(A)", advance = "no") this%file
          case (TK_LINE)
            write (unit, "(I0)", advance = "no") this%line
          case (TK_MEM)
            write (unit, "(F3.1)", advance = "no") this%mem
          case (TK_MSG)
            write (unit, "(A)", advance = "no") this%msg
          end select
        end associate
      end do
      write (unit, *)
    end associate
  end subroutine

end module
