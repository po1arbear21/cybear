program material_test
  !! Phase A material-catalog test runner.
  !!
  !! Main mode: load two_material_parse.ini, assert the catalog is populated
  !! correctly (smc(:), smc_map, reg_trans(:)%material_id).
  !!
  !! Subprocess mode (--load-only <ini>): just call dev%init on the INI and exit.
  !! Used by check_error to spawn the test binary against error-path INIs and
  !! verify each aborts with the expected message.

  use device_m,        only: device
  use normalization_m, only: init_normconst, denorm
  use string_m,        only: string

  implicit none

  type(device)              :: dev
  character(:), allocatable :: arg1, arg2, exe_path
  integer                   :: argc, n_pass, n_fail, mat_id
  logical                   :: ok

  argc = command_argument_count()

  ! ---- subprocess mode -----------------------------------------------------
  if (argc >= 2) then
    call get_arg(1, arg1)
    if (arg1 == "--load-only") then
      call get_arg(2, arg2)
      call init_normconst(300.0)
      call dev%init(arg2, 300.0)
      stop 0
    end if
  end if

  ! ---- main mode -----------------------------------------------------------
  n_pass = 0
  n_fail = 0

  print "(A)", "===================================="
  print "(A)", " Phase A — material catalog"
  print "(A)", "===================================="

  ! Happy-path catalog assertions
  call init_normconst(300.0)
  call dev%init("two_material_parse.ini", 300.0)

  call expect_int(size(dev%par%smc), 2,              "size(smc) == 2",                n_pass, n_fail)
  call expect_str(dev%par%smc(1)%name%s, "Si",       "smc(1).name == Si",             n_pass, n_fail)
  call expect_str(dev%par%smc(2)%name%s, "SiGe",     "smc(2).name == SiGe",           n_pass, n_fail)
  call expect_real(denorm(dev%par%smc(1)%chi, "eV"), 4.05, 1e-4, "smc(1).chi == 4.05 eV",         n_pass, n_fail)
  call expect_real(denorm(dev%par%smc(2)%chi, "eV"), 4.10, 1e-4, "smc(2).chi == 4.10 eV",         n_pass, n_fail)
  call expect_int(dev%par%smc_default, 1,            "smc_default == 1",              n_pass, n_fail)

  ! smc_map resolves names to indices
  call dev%par%smc_map%get(string("Si"),   mat_id, status = ok)
  call expect_log(ok,                                "smc_map[Si] resolves",          n_pass, n_fail)
  call expect_int(mat_id, 1,                         "smc_map[Si] == 1",              n_pass, n_fail)
  call dev%par%smc_map%get(string("SiGe"), mat_id, status = ok)
  call expect_log(ok,                                "smc_map[SiGe] resolves",        n_pass, n_fail)
  call expect_int(mat_id, 2,                         "smc_map[SiGe] == 2",            n_pass, n_fail)

  ! Transport regions resolve their material_id from the INI material= key
  call expect_int(size(dev%par%reg_trans), 2,         "2 transport regions",           n_pass, n_fail)
  call expect_int(dev%par%reg_trans(1)%material_id, 1,"reg_trans(1).material_id == 1", n_pass, n_fail)
  call expect_int(dev%par%reg_trans(2)%material_id, 2,"reg_trans(2).material_id == 2", n_pass, n_fail)

  ! Error-path subprocesses
  call get_arg(0, exe_path)
  call check_error(exe_path, "two_material_missing_name.ini", &
    "every block to declare name", n_pass, n_fail)
  call check_error(exe_path, "two_material_uniformity.ini", &
    "dist differs from default", n_pass, n_fail)
  call check_error(exe_path, "two_material_unknown.ini", &
    "unknown material 'Unobtainium'", n_pass, n_fail)

  print "(A)", "===================================="
  print "(A,I0,A,I0,A)", " Summary: ", n_pass, " passed, ", n_fail, " failed"
  print "(A)", "===================================="

  if (n_fail > 0) error stop 1

contains

  subroutine get_arg(n, val)
    integer,                   intent(in)  :: n
    character(:), allocatable, intent(out) :: val
    integer :: l
    call get_command_argument(n, length = l)
    allocate (character(l) :: val)
    call get_command_argument(n, val)
  end subroutine

  subroutine expect_int(actual, expected, label, np, nf)
    integer,      intent(in)    :: actual, expected
    character(*), intent(in)    :: label
    integer,      intent(inout) :: np, nf
    if (actual == expected) then
      np = np + 1
      print "(A,A)", "  PASS: ", label
    else
      nf = nf + 1
      print "(A,A,A,I0,A,I0,A)", "  FAIL: ", label, " (got ", actual, ", want ", expected, ")"
    end if
  end subroutine

  subroutine expect_str(actual, expected, label, np, nf)
    character(*), intent(in)    :: actual, expected, label
    integer,      intent(inout) :: np, nf
    if (actual == expected) then
      np = np + 1
      print "(A,A)", "  PASS: ", label
    else
      nf = nf + 1
      print "(A,A,A,A,A,A,A)", "  FAIL: ", label, " (got '", actual, "', want '", expected, "')"
    end if
  end subroutine

  subroutine expect_real(actual, expected, tol, label, np, nf)
    real,         intent(in)    :: actual, expected, tol
    character(*), intent(in)    :: label
    integer,      intent(inout) :: np, nf
    if (abs(actual - expected) <= tol) then
      np = np + 1
      print "(A,A)", "  PASS: ", label
    else
      nf = nf + 1
      print "(A,A,A,ES12.5,A,ES12.5,A)", "  FAIL: ", label, " (got ", actual, ", want ", expected, ")"
    end if
  end subroutine

  subroutine expect_log(actual, label, np, nf)
    logical,      intent(in)    :: actual
    character(*), intent(in)    :: label
    integer,      intent(inout) :: np, nf
    if (actual) then
      np = np + 1
      print "(A,A)", "  PASS: ", label
    else
      nf = nf + 1
      print "(A,A)", "  FAIL: ", label
    end if
  end subroutine

  subroutine check_error(exe, ini, expected_msg, np, nf)
    !! Spawn `<exe> --load-only <ini> 2> err.tmp` and verify it exits nonzero
    !! and that stderr contains `expected_msg`.
    character(*), intent(in)    :: exe, ini, expected_msg
    integer,      intent(inout) :: np, nf

    character(:), allocatable :: errfile, cmd
    integer                   :: ec, ios, u
    character(1024)           :: line
    logical                   :: found

    errfile = "err_" // ini // ".tmp"
    ! program_error writes to stdout, not stderr; capture both.
    cmd = trim(exe) // " --load-only " // ini // " > " // errfile // " 2>&1"
    call execute_command_line(cmd, exitstat = ec, wait = .true.)

    if (ec == 0) then
      nf = nf + 1
      print "(A,A,A)", "  FAIL: ", ini, " expected nonzero exit but got 0"
      return
    end if

    open (newunit = u, file = errfile, status = "old", action = "read", iostat = ios)
    if (ios /= 0) then
      nf = nf + 1
      print "(A,A,A)", "  FAIL: ", ini, " could not read stderr capture"
      return
    end if

    found = .false.
    do
      read (u, "(A)", iostat = ios) line
      if (ios /= 0) exit
      if (index(line, expected_msg) > 0) then
        found = .true.
        exit
      end if
    end do
    close (u)

    if (found) then
      np = np + 1
      print "(A,A,A,A,A)", "  PASS: ", ini, " aborted with '...", expected_msg, "...'"
    else
      nf = nf + 1
      print "(A,A,A,A,A)", "  FAIL: ", ini, " aborted but expected '", expected_msg, "' not in stderr"
    end if
  end subroutine

end program
