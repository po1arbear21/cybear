program test_contact_normal
  !! Test program to verify contact normal direction detection

  use device_m,        only: device, dev
  use schottky_m,      only: get_schottky_contact_normal_dir
  use contact_m,       only: CT_SCHOTTKY
  use normalization_m, only: init_normconst, denorm

  implicit none

  integer :: ict, normal_dir
  character(len=100) :: config_file
  real :: T

  ! Get config file from command line or use default
  if (command_argument_count() > 0) then
    call get_command_argument(1, config_file)
  else
    ! Default to 2D test
    config_file = "/test/test_4contacts_2d.ini"
    print *, "No config file specified, using default: ", trim(config_file)
    print *, "Usage: test_contact_normal [config_file]"
    print *, "  e.g.: test_contact_normal test_6contacts_3d.ini"
    print *, ""
  end if

  print *, "Testing contact normal direction detection"
  print *, "Using config file: ", trim(config_file)
  print *, ""

  ! Initialize normalization constants FIRST
  T = 300.0  ! Temperature in Kelvin
  call init_normconst(T)

  ! Now initialize device
  call dev%init(config_file, T)

  ! Show device info
  print *, "Device dimensionality: ", dev%par%g%dim, "D"
  print *, "Number of contacts: ", dev%par%nct
  print *, ""

  do ict = 1, dev%par%nct
    print *, "Contact ", ict, ": ", trim(dev%par%contacts(ict)%name)
    print *, "  Type: ", dev%par%contacts(ict)%type

    if (dev%par%contacts(ict)%type == CT_SCHOTTKY) then
      normal_dir = get_schottky_contact_normal_dir(dev%par, ict)

      select case(normal_dir)
      case(1)
        print *, "  Normal direction: X"
      case(2)
        print *, "  Normal direction: Y"
      case(3)
        print *, "  Normal direction: Z"
      case default
        print *, "  Normal direction: NOT AXIS-ALIGNED (", normal_dir, ")"
      end select

      ! Get first vertex position to show contact location
      if (dev%par%transport_vct(ict)%n > 0) then
        block
          integer :: idx(dev%par%g%idx_dim)
          real :: pos(3)
          idx = dev%par%transport_vct(ict)%get_idx(1)
          call dev%par%g%get_vertex(idx, pos(1:dev%par%g%dim))
          print *, "  Sample position (nm): ", denorm(pos(1:dev%par%g%dim), "nm")
        end block
      end if
    else
      print *, "  (Not a Schottky contact, skipping normal detection)"
    end if
    print *, ""
  end do

  print *, "Test complete!"

end program
