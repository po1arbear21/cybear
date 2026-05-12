program test

  use util_m, only: get_hostname

  use test_poisson_m

  implicit none

  print "(A)", "Start test on " // get_hostname()

  call test_poisson()

end program
