program example

  use esystem_m,                only: esystem
  use example_density_m,        only: dens
  use example_poisson_m,        only: pois
  use example_imref_m,          only: iref, c_dens
  use example_charge_density_m, only: e_dens, calc_e_dens
  use example_contact_m,        only: init_contacts
  use example_device_m,         only: init_device
  use example_potential_m,      only: pot
  use input_m,                  only: input_file

  implicit none

  type(input_file) :: f
  type(esystem)    :: sys

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

  call sys%g%output("output/degraph")
  call sys%init_final()
  call sys%solve()
end program
