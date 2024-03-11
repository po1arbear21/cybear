program cl_test_simple

  use cl_options_m, only: get_cl_options_simple
  use string_m,     only: string

  implicit none

  type(string) :: names(2), values(2)

  names(1)%s = "file"
  names(2)%s = "output"

  call get_cl_options_simple(names, values)

  print "(3A)", names(1)%s, ": ", values(1)%s
  print "(3A)", names(2)%s, ": ", values(2)%s

end program
