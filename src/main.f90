program main

use grid_m

  use approx_m,        only: approx_imref, approx_potential
  use device_m,        only: dev
  use error_m,         only: program_error
  use input_m,         only: input_file
  use input_src_m,     only: const_src
  use newton_m,        only: newton_opt
  use normalization_m, only: init_normconst, norm, denorm
  use steady_state_m,  only: steady_state

  implicit none

  character(:), allocatable :: filename, dev_filename
  real                      :: T
  type(input_file)          :: file
  type(newton_opt)          :: opt_nlpe, opt_dd, opt_gum, opt_full

  ! get input filename from command line arguments
  call parse_command_line_args()

  ! load input file
  call file%init(filename)

  ! get temperature and initialize normalization
  call file%get("", "temperature", T, normalize = .false.)
  call init_normconst(T)

  ! get device filename and initialize device
  call file%get("", "device", dev_filename)
  call dev%init(dev_filename)

  ! get iteration options
  call load_iteration_params("nlpe params",        dev%sys_nlpe%n,  opt_nlpe)
  opt_nlpe%dx_lim = norm(0.1, "V")
  call load_iteration_params("dd params",          dev%sys_dd(1)%n, opt_dd  )
  call load_iteration_params("gummel params",      1,               opt_gum )
  call load_iteration_params("full newton params", dev%sys_full%n,  opt_full)

  ! steady-state
  call solve_steady_state()

contains

  subroutine parse_command_line_args()
    integer :: nargs, length, status

    ! get number of command line arguments
    nargs = command_argument_count()
    if (nargs < 1) call program_error("Command filename expected!")
    if (nargs > 1) print *, "Warning: Command-line arguments after first one are ignored"

    ! read device filename
    allocate (character(32) :: filename)
    do while (.true.)
      call get_command_argument(1, value = filename, length = length, status = status)
      if (status == -1) then
        deallocate (filename)
        allocate (character(2*len(filename)) :: filename)
      elseif (status == 0) then
        filename = filename(1:length)
        exit
      else
        call program_error("Error at reading command-line argument")
      end if
    end do
  end subroutine

  subroutine load_iteration_params(section_name, n, opt)
    character(*),     intent(in)  :: section_name
      !! section name in input file
    integer,          intent(in)  :: n
      !! system size
    type(newton_opt), intent(out) :: opt
      !! output newton options

    integer :: max_it
    logical :: log
    real    :: rtol, atol

    call file%get(section_name, "max_it", max_it)
    call file%get(section_name, "rtol",   rtol  )
    call file%get(section_name, "atol",   atol  )
    call file%get(section_name, "log",    log   )

    call opt%init(n, atol = atol, rtol = rtol, max_it = max_it, log = log)
  end subroutine

  subroutine solve_steady_state()
    integer              :: si, ict
    integer, allocatable :: sids(:)

    type(const_src)    :: input
    type(steady_state) :: ss

    call file%get_sections("steady state", sids)
    do si = 1, size(sids)
      ! set voltages
      do ict = 1, size(dev%par%contacts)
        call file%get(sids(si), "V_"//dev%par%contacts(ict)%name, dev%volt(ict)%x)
      end do
      call input%init([(dev%volt(ict)%x, ict = 1, size(dev%par%contacts))])

      ! solve steady-state
      call ss%run(dev%sys_full, nopt = opt_full, input = input, gum = gummel)
    end do
  end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: i, ci
    real               :: error, err_pot, err_iref(2)
    real, allocatable  :: pot0(:), iref0(:,:)
    type(steady_state) :: ss_nlpe, ss_dd(2)

    allocate (pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n,2))

    ! approximate imrefs
    do ci = dev%par%ci0, dev%par%ci1
      call approx_imref(dev%par, dev%iref(ci), dev%volt)
    end do

    ! approximate potential
    call approx_potential(dev%par, dev%pot, dev%iref)

    ! gummel iteration
    i = 0
    error = huge(1.0)
    do while ((error > opt_gum%atol(1)) .and. (i < opt_gum%max_it))
      i = i + 1

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call ss_nlpe%run(dev%sys_nlpe, nopt = opt_nlpe)
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()
        call ss_dd(ci)%run(dev%sys_dd(ci), nopt = opt_dd)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

      ! log
      if (opt_gum%log) then
        print "(A,I6,ES24.16)", "Gummel: ", i, denorm(error, "V")
      end if
    end do

    if ((i > opt_gum%max_it) .and. (error > opt_gum%atol(1))) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if
  end subroutine

end program
