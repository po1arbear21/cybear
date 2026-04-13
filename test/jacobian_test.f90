program jacobian_test
  !! Taylor remainder tests for Jacobian validation of all equation systems

  use util_m,            only: get_hostname
  use device_m,          only: device
  use normalization_m,   only: init_normconst, norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE
  use approx_m,          only: approx_imref, approx_potential
  use taylor_remainder_m, only: taylor_test, taylor_test_per_block

  implicit none

  type(device) :: dev
  integer :: ci, ict
  real :: beam_y

  print "(A)", "Start test on " // get_hostname()

  ! initialize device
  call init_normconst(300.0)
  call dev%init("beam_minimal.ini", 300.0)

  ! set voltages to 0V
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = 0.0
  end do

  ! approximate initial conditions
  do ci = dev%par%ci0, dev%par%ci1
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)

  ! set beam position
  beam_y = 2000.0
  dev%beam_pos%x = norm(beam_y, 'nm')
  do ci = dev%par%ci0, dev%par%ci1
    call dev%calc_bgen(ci)%eval()
  end do
  print "(A,F8.1,A)", "  Beam position: y = ", beam_y, " nm"

  ! --- NLPE ---
  print "(A)", ""
  print "(A)", "========================================"
  print "(A)", " Taylor Test: NLPE (non-linear Poisson)"
  print "(A)", "========================================"
  call taylor_test(dev%sys_nlpe)
  call taylor_test_per_block(dev%sys_nlpe)

  ! --- Full Newton ---
  print "(A)", "========================================"
  print "(A)", " Taylor Test: Full Newton"
  print "(A)", "========================================"
  call taylor_test(dev%sys_full)
  call taylor_test_per_block(dev%sys_full)

  ! --- DD systems ---
  do ci = dev%par%ci0, dev%par%ci1
    print "(A)", "========================================"
    print "(A,I0,A)", " Taylor Test: DD (carrier ", ci, ")"
    print "(A)", "========================================"
    call taylor_test(dev%sys_dd(ci))
    call taylor_test_per_block(dev%sys_dd(ci))
  end do

end program
