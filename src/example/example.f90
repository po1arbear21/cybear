program example

  use normalization_m,          only: norm, denorm
  use newton_m,                 only: newton_opt
  use esystem_m,                only: esystem
  use example_density_m,        only: dens
  use example_poisson_m,        only: pois
  use example_imref_m,          only: iref, c_dens
  use example_charge_density_m, only: e_dens, calc_e_dens
  use example_contact_m,        only: init_contacts
  use example_device_m,         only: init_device, grd
  use example_potential_m,      only: pot
  use input_m,                  only: input_file

  implicit none

  type(input_file) :: f
  type(esystem)    :: sys
  type(newton_opt) :: opt
  integer :: i

  ! init input file
  call f%init("src/example/example.inp")

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

  ! init esystem
  call sys%init("non-linear poisson")
  call sys%add_equation(pois)
  call sys%add_equation(calc_e_dens)
  call sys%add_equation(c_dens)
  call sys%provide(iref)
  call sys%init_final()
  call sys%g%output("graph")

  call opt%init(sys%n, log = .true.)
  opt%dx_lim = norm(0.2, "V")

  call sys%solve(nopt = opt)

  do i = 1, size(pot%x)
    print "(3ES24.16)", denorm(grd%x(i), "nm"), denorm(pot%x(i), "V"), denorm(dens%x(i), "1/cm^3")
  end do

end program
