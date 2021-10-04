module example_steady_state_m

  use error_m,                   only: program_error
  use esystem_m,                 only: esystem
  use example_charge_density_m,  only: calc_charge_dens, charge_dens
  use example_contact_m,         only: contacts, init_contacts
  use example_continuity_m,      only: contin
  use example_current_density_m, only: calc_current_dens, current_dens
  use example_density_m,         only: dens
  use example_device_m,          only: grd, init_device
  use example_imref_m,           only: calc_dens, calc_iref, iref
  use example_mobility_m,        only: calc_mobil, mobil
  use example_poisson_m,         only: pois
  use example_potential_m,       only: pot
  use example_ramo_m,            only: ramo_init
  use input_m,                   only: input_file
  use newton_m,                  only: newton_opt
  use normalization_m,           only: denorm, norm

  implicit none

  private
  public init_configuration
  public init_dd, init_full, init_nlpe
  public solve_full_newton, solve_gummel, solve_nlpe
  public output

  type(esystem)      :: sys_dd, sys_full, sys_nlpe
  type(input_file)   :: f
  type(newton_opt)   :: opt_dd, opt_full, opt_nlpe
  type(steady_state) :: steady_dd, steady_full, steady_nlpe

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
    call current_dens%init()
    call pot%init()
    call dens%init()
    call charge_dens%init()
    call iref%init()
    call mobil%init()

    ! init equations
    call calc_current_dens%init()
    call calc_charge_dens%init()
    call pois%init()
    call calc_dens%init()
    call calc_iref%init()
    call contin%init()
    call calc_mobil%init()

    ! init fundamental solution
    call ramo_init()
  end subroutine

  subroutine init_dd()
    integer :: max_it
    real    :: atol, rtol

    ! reading parameters
    call f%get("dd_parameters", "max_it", max_it)
    call f%get("dd_parameters", "atol",   atol)
    call f%get("dd_parameters", "rtol",   rtol)

    ! init equation system
    call sys_dd%init("drift-diffusion")

    ! add related equations to the system
    call sys_dd%add_equation(contin)
    call sys_dd%add_equation(calc_current_dens)
    call sys_dd%add_equation(calc_mobil)

    ! provide variables
    call sys_dd%provide(pot)
    call sys_dd%provide(iref)

    ! finalize the equation system
    call sys_dd%init_final()

    ! output the related graph to the equation system
    call sys_dd%g%output("dd")

    ! init solver options
    call opt_dd%init(sys_dd%n, log = .false., max_it = max_it)
    opt_dd%atol = atol
    opt_dd%rtol = rtol
  end subroutine

  subroutine init_nlpe()
    integer :: max_it, i
    real    :: atol, rtol

    ! reading parameters
    call f%get("nlpe_parameters", "max_it", max_it)
    call f%get("nlpe_parameters", "atol",   atol)
    call f%get("nlpe_parameters", "rtol",   rtol)

    ! init equation system
    call sys_nlpe%init("non-linear poisson")

    ! add related equations to the system
    call sys_nlpe%add_equation(pois)
    call sys_nlpe%add_equation(calc_charge_dens)
    call sys_nlpe%add_equation(calc_dens)

    ! provide variables
    call sys_nlpe%provide(iref)
    do i = 1, size(contacts)
      call sys_nlpe%provide(contacts(i)%volt)
    end do

    ! finalize the equation system
    call sys_nlpe%init_final()

    ! output the related graph to the equation system
    call sys_nlpe%g%output("nlpe")

    ! init solver options
    call opt_nlpe%init(sys_nlpe%n, log = .false., max_it = max_it)
    opt_nlpe%dx_lim = norm(0.2, "V")
    opt_nlpe%atol = atol
    opt_nlpe%rtol = rtol
  end subroutine

  subroutine init_full()
    integer :: i
    integer :: max_it
    real    :: atol, rtol

    ! reading parameters
    call f%get("nlpe_parameters", "max_it", max_it)
    call f%get("nlpe_parameters", "atol",   atol)
    call f%get("nlpe_parameters", "rtol",   rtol)

    ! init equation system
    call sys_full%init("full newton")

    ! add related equations to the system
    call sys_full%add_equation(pois)
    call sys_full%add_equation(contin)
    call sys_full%add_equation(calc_charge_dens)
    call sys_full%add_equation(calc_current_dens)
    call sys_full%add_equation(calc_mobil)
    call sys_full%add_equation(calc_iref)

    ! provide variables
    do i = 1, size(contacts)
      call sys_full%provide(contacts(i)%volt, input = .true.)
    end do

    ! finalize the equation system
    call sys_full%init_final()

    ! output the related graph to the equation system
    call sys_full%g%output("full")

    ! init solver options
    call opt_full%init(sys_full%n, log = .true., max_it = max_it)
    opt_full%dx_lim = norm(0.2, "V")
    opt_full%atol = atol
    opt_full%rtol = rtol
  end subroutine

  subroutine solve_nlpe()
    !! solve the non-linear poisson system
    call steady_nlpe%run(sys_nlpe, nopt = opt_nlpe)
  end subroutine

  subroutine solve_gummel()
    !! solve the gummel system
    integer           :: max_it, i
    real              :: error, atol
    real, allocatable :: iref0(:), pot0(:)

    allocate(iref0(size(iref%get())), pot0(size(pot%get())))

    ! reading parameters
    call f%get("gummel_parameters", "max_it", max_it)
    call f%get("gummel_parameters", "atol",   atol)

    i = 1
    error = huge(1.0)
    do while (error > atol .and. i <= max_it)
      pot0 = pot%get()
      call steady_nlpe%run(sys_nlpe, nopt = opt_nlpe)
      iref0 = iref%get()
      call steady_dd%run(  sys_dd,   nopt = opt_dd)
      call calc_iref%eval()

      error = max(maxval(abs(pot%get()-pot0)), maxval(abs(iref%get()-iref0)))
      print *, i, denorm(error, "V")
      i = i + 1
    end do

    if (i > max_it .and. error > atol) call program_error("max number of iterations reached")
  end subroutine

  subroutine solve_full_newton()
    !! solve the full newton system
    integer           :: i
    real, allocatable :: input(:,:)

    allocate(input(size(contacts), 1))
    input(:,1) = [(contacts(i)%volt%x, i = 1, size(contacts))]
    call steady_full%run(sys_full, nopt = opt_full, input = input)
  end subroutine

  subroutine output()
    !! output density, potential and current_density
    call dens%output_data("dens.csv")
    call pot%output_data( "pot.csv")
    call current_dens%output_data("current_dens.csv")
  end subroutine

end module
