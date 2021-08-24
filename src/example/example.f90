program example

  use example_steady_state_m

  implicit none

  call init_configuration("src/example/example.inp")
  call init_dd()
  call init_full()
  call init_nlpe()

  call solve_nlpe()

  call output()

end program
