m4_include(../../util/macro.f90.inc)

module test_logging_m

  use logging_m,   only: init_logging, finish_logging, logging, LVL_DEBUG, LVL_INFO, LVL_WARN
  use string_m,    only: new_string, string
  use test_case_m, only: test_case

  implicit none

  private
  public test_logging

contains

  subroutine test_logging()
    type(string)    :: line
    integer         :: funit, iostat
    type(test_case) :: tc

    call tc%init("logging")

    call init_logging(file = "/tmp/test_logging.log", record_date = .false., record_time = .false., record_mem = .false., verbosity = 5, fmt = "%LVL: %MSG")

    ! long version
    call logging(LVL_INFO, "logging test 1", file = "m4___file__", line = m4___line__)
    call logging(LVL_DEBUG, "logging test 2", file = "m4___file__", line = m4___line__)
    call logging(LVL_WARN, "logging test 3", file = "m4___file__", line = m4___line__)

    ! short version using macro
    m4_info("logging test 4")
    m4_info("logging test 5")
    m4_info("logging test 6")
    m4_error("logging test 7")

    call finish_logging()

    allocate (character(80) :: line%s)

    open (newunit = funit, file = "/tmp/test_logging.log", status = "old", action = "read", iostat = iostat)
    call tc%assert_eq(0, iostat, "file exists")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 1 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("INFO: logging test 1"), line, "line 1")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 2 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("DEBUG: logging test 2"), line, "line 2")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 3 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("WARNING: logging test 3"), line, "line 3")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 4 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("INFO: logging test 4"), line, "line 4")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 5 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("INFO: logging test 5"), line, "line 5")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 6 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("INFO: logging test 6"), line, "line 6")
    read (funit, "(A)", iostat = iostat) line%s
    call tc%assert_eq(0, iostat, "line 7 exists")
    if (iostat /= 0) goto 10
    call tc%assert_eq(new_string("ERROR: logging test 7"), line, "line 7")
    10 close (funit)

    call tc%finish()
  end subroutine

end module
