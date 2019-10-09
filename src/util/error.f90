module error_m
  use ifcore
  implicit none

contains

  subroutine program_error(msg, code)
    !! Prints an error message and performs traceback

    character(len=*),  intent(in) :: msg
      !! error message
    integer, optional, intent(in) :: code
      !! additional error code

    print *, msg

    if (present(code)) then
      print *, "error code: ", code
    end if

    ! print traceback and stop program
    call tracebackqq()
  end subroutine

  subroutine assert_failed(expr, file, line)
    character(len=*), intent(in) :: expr
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line

    print *
    print *, expr
    print *, "Assertion failed!"
    print "(3A, 1I6)", "file: ", file, "; line: ", line
    print *

    call tracebackqq()
  end subroutine

end module