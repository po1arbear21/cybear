program cl_test

  use cl_options_m, only: cl_option_descriptor, cl_option, get_cl_options

  implicit none

  integer                      :: idesc, i, j, nopt
  type(cl_option), allocatable :: opt(:)
  integer,         allocatable :: iopt(:), jopt(:)
  type(cl_option_descriptor) :: desc(4) = [ &
    cl_option_descriptor('f',     "file", .true.,  .false.,  .true.,  .true.), &
    cl_option_descriptor('o',   "output", .true.,  .false.,  .true.,  .true.), &
    cl_option_descriptor('p', "optional", .false.,  .true.,  .true., .false.), &
    cl_option_descriptor('d',    "debug", .false., .false., .false., .false.)  &
  ]

  call get_cl_options(desc, opt, iopt, jopt)

  print "(A)", "CL options sorted by occurence: "
  do i = 1, size(opt)
    select case (opt(i)%short)
    case ('f')
      print "(A)", "file: " // opt(i)%arg

    case ('o')
      print "(A)", "output: " // opt(i)%arg

    case ('d')
      print "(A)", "debug"

    case ('p')
      if (allocated(opt(i)%arg)) then
        print "(A)", "optional: " // opt(i)%arg
      else
        print "(A)", "optional"
      end if

    end select
  end do

  print *
  print "(A)", "CL options sorted by descriptor: "
  do idesc = 1, size(desc)
    nopt = iopt(idesc + 1) - iopt(idesc)
    print "(2A,I0,A)", trim(desc(idesc)%long), " (", nopt, "): "

    do i = iopt(idesc), iopt(idesc + 1) - 1
      j = jopt(i)
      select case (opt(j)%short)
      case ('f')
        print "(A)", "   " // opt(j)%arg

      case ('o')
        print "(A)", "   " // opt(j)%arg

      case ('d')
        print *

      case ('p')
        if (allocated(opt(j)%arg)) then
          print "(A)", "   " // opt(j)%arg
        else
          print *
        end if

      end select
    end do
  end do

end program
