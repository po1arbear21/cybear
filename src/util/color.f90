module color_m
  implicit none

  ! console colors
  character(4), parameter :: COL_DEFAULT = achar(27)//"[0m"
  character(7), parameter :: COL_GRAY    = achar(27)//"[1;30m"
  character(7), parameter :: COL_RED     = achar(27)//"[1;31m"
  character(7), parameter :: COL_GREEN   = achar(27)//"[1;32m"
  character(7), parameter :: COL_YELLOW  = achar(27)//"[1;33m"
  character(7), parameter :: COL_BLUE    = achar(27)//"[1;34m"
  character(7), parameter :: COL_MAGENTA = achar(27)//"[1;35m"
  character(7), parameter :: COL_CYAN    = achar(27)//"[1;36m"
  character(7), parameter :: COL_WHITE   = achar(27)//"[1;37m"
end module
