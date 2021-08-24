module example_steady_state_m

  use error_m,                   only: program_error
  use esystem_m,                 only: esystem
  use example_charge_density_m,  only: calc_e_dens, e_dens
  use example_contact_m,         only: contacts, init_contacts
  use example_continuity_m,      only: cont
  use example_current_density_m, only: curr, c_curr
  use example_density_m,         only: dens
  use example_device_m,          only: grd, init_device
  use example_imref_m,           only: c_dens, iref
  use example_poisson_m,         only: pois
  use example_potential_m,       only: pot
  use input_m,                   only: input_file
  use newton_m,                  only: newton_opt
  use normalization_m,           only: denorm, norm

  implicit none

  private
  public init_configuration
  public init_dd, init_full, init_nlpe
  public solve_full_newton, solve_gummel, solve_nlpe
  public output
  ! public sys_dd, sys_full, sys_nlpe

  type(esystem)    :: sys_dd, sys_full, sys_nlpe
  type(input_file) :: f
  type(newton_opt) :: opt_full, opt_nlpe

contains

  subroutine init_configuration(file_src)
    character(*), intent(in) :: file_src

    ! init input file
    call f%init(file_src)

    ! init devise
    call init_device(f)

    ! init contacts
    call init_contacts(f)

    ! init variables
    call pot%init()
    call dens%init()
    call e_dens%init()
    call iref%init()

    ! init equations
    call calc_e_dens%init()
    call pois%init()
    call c_dens%init()
    call cont%init()
  end subroutine

  subroutine init_dd()
    call sys_dd%init("drift-diffusion")
    call sys_dd%add_equation(cont)
    call sys_dd%add_equation(c_curr)
    call sys_dd%provide(pot)
    call sys_dd%init_final()
    call sys_dd%g%output("dd")
  end subroutine

  subroutine init_nlpe()
    integer :: max_it
    real    :: atol(1), rtol(1)

    call f%get("nlpe_parameters", "max_it", max_it)
    call f%get("nlpe_parameters", "atol",   atol(1))
    call f%get("nlpe_parameters", "rtol",   rtol(1))

    ! init optimizer
    call opt_nlpe%init(sys_nlpe%n, log = .true., max_it = max_it, atol = atol, rtol = rtol)
    opt_nlpe%dx_lim = norm(0.2, "V")

    ! init equation system
    call sys_nlpe%init("non-linear poisson")
    ! add related equations to the system
    call sys_nlpe%add_equation(pois)
    call sys_nlpe%add_equation(calc_e_dens)
    call sys_nlpe%add_equation(c_dens)
    ! provide variables
    call sys_nlpe%provide(iref)
    call sys_nlpe%provide(pois%volt)
    ! finalize the equation system
    call sys_nlpe%init_final()

    ! output the related graph to the esystem
    call sys_nlpe%g%output("nlpe")
  end subroutine

  subroutine init_full()
    integer :: i
    integer :: max_it
    real    :: atol(1), rtol(1)

    call f%get("nlpe_parameters", "max_it", max_it)
    call f%get("nlpe_parameters", "atol",   atol(1))
    call f%get("nlpe_parameters", "rtol",   rtol(1))

    ! init optimizer
    call opt_full%init(sys_nlpe%n, log = .true., max_it = max_it, atol = atol, rtol = rtol)
    opt_full%dx_lim = norm(0.2, "V")


    ! init equation system
    call sys_full%init("full newton")
    ! add related equations to the system
    call sys_full%add_equation(pois)
    call sys_full%add_equation(cont)
    call sys_full%add_equation(c_curr)
    ! provide variables
    do i = 1, size(contacts)
      call sys_full%provide(contacts(i)%volt, input = .true.)
    end do
    ! finalize the equation system
    call sys_full%init_final()

    ! output the related graph to the esystem
    call sys_full%g%output("full")
  end subroutine

  subroutine solve_nlpe()
    call sys_nlpe%solve(nopt = opt_nlpe)
  end subroutine

  subroutine solve_gummel()
    integer           :: max_it, i
    real              :: error, atol
    real, allocatable :: iref0(:), pot0(:)

    allocate(iref0(size(iref%get())), pot0(size(pot%get())))

    call f%get("gummel_parameters", "max_it", max_it)
    call f%get("gummel_parameters", "atol",   atol)

    i = 1
    do while (error > atol .and. i <= max_it)
      pot0 = pot%get()
      call sys_nlpe%solve()

      iref0 = iref%get()
      call sys_dd%solve()
      call iref%calc()

      error = max(maxval(abs(pot%get()-pot0)), maxval(abs(iref%get()-iref0)))
      print *, i, error
      i = i + 1
    end do

    if (i > max_it .and. error > atol) call program_error("max number of iterations reached")
  end subroutine

  subroutine solve_full_newton()
    call sys_full%solve(nopt = opt_full)
  end subroutine

  subroutine output()
    call dens%output_data("dens.csv")
    call pot%output_data( "pot.csv")
    call curr%output_data("curr.csv")
  end subroutine

end module
