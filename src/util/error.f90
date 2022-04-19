m4_include(macro.f90.inc)

module error_m

  use color_m, only: COL_DEFAULT, COL_MAGENTA
  m4_ifdef({m4_intel},use ifcore)

  implicit none

  private
  public error_msg, set_error_mode, program_error, program_has_error, program_get_errors, assert_failed

  logical :: error_mode = .false.
    !! save error messages to vector

  type error_msg
    !! error message
    character(:), allocatable :: msg
      !! message
    integer                   :: code
      !! additional error code
    character(:), allocatable :: file
      !! source code file
    integer                   :: line
      !! line in source code file
  end type

  m4_define({T},{error_msg})
  m4_include(vector_def.f90.inc)

  type(vector_error_msg) :: emsg

contains

  m4_define({T},{error_msg})
  m4_include(vector_imp.f90.inc)

  subroutine set_error_mode(mode)
    !! enable/disable saving of error messages
    logical, intent(in) :: mode
      !! false: do not save error messages, stop program instead; true: save error messages, do not stop program

    error_mode = mode
    if (mode .and. .not. allocated(emsg%d)) call emsg%init(0, c = 8)
  end subroutine

  subroutine program_error(msg, code, file, line)
    !! Depending on error_mode: Prints an error message and performs traceback or saves message for later

    character(*),           intent(in) :: msg
      !! error message
    integer,      optional, intent(in) :: code
      !! additional error code
    character(*), optional, intent(in) :: file
      !! source code file
    integer,      optional, intent(in) :: line
      !! source code line

    type(error_msg) :: em

    if (error_mode) then
      em%msg  = msg
      em%code = 0
      em%file = ""
      em%line = 0
      if (present(code)) em%code = code
      if (present(file)) em%file = file
      if (present(line)) em%line = line
      call emsg%push(em)
    else
      print "(A)", COL_MAGENTA//msg//COL_DEFAULT

      if (present(code)) then
        print "(A,I0)", "error code: ", code
      end if

      ! print traceback and stop program
      m4_ifdef({m4_intel},{call tracebackqq()})
      m4_ifdef({m4_gnu},{call backtrace()})
      call exit(1)
    end if
  end subroutine

  function program_has_error() result(f)
    !! check whether one or multiple errors occured
    logical :: f
      !! return error flag

    f = .false.
    if (error_mode) f = (emsg%n > 0)
  end function

  subroutine program_get_errors(em)
    !! return array with all error messages, clear error message vector
    type(error_msg), allocatable, intent(out) :: em(:)

    if (error_mode) then
      allocate (em(emsg%n), source = emsg%d(1:emsg%n))
      call emsg%resize(0)
    else
      allocate (em(0))
    end if
  end subroutine

  subroutine assert_failed(expr, file, line)
    character(*), intent(in) :: expr
    character(*), intent(in) :: file
    integer,      intent(in) :: line

    print *
    print "(A)", expr
    print "(A)", COL_MAGENTA//"Assertion failed!"//COL_DEFAULT
    print "(3A,I0)", "file: ", file, "; line: ", line
    print *

    m4_ifdef({m4_intel},{call tracebackqq()})
    m4_ifdef({m4_gnu},{call backtrace()})
    call exit(1)
  end subroutine

end module
