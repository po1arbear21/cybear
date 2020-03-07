module test_random_m
  use test_case_m
  use random_m
  implicit none

contains

  subroutine test_random()
    type(test_case)    :: tc
    ! type(random)       :: rng
    ! integer, parameter :: N = 100000
    ! integer            :: i, seed(2), funit
    ! real               :: x(N)

    print "(1A)", "test_random"
    call tc%init("random")

    ! seed = [512323423456, 178234119012]
    ! call rng%init(seed, [835253408953, 1])
    ! call rng%next_reals(x)

    ! open (newunit = funit, file = "rand_pcg.asc", status = "replace", action = "write")
    ! do i = 1, N
    !   write (funit, "(1E64.56)") x(i)
    ! end do
    ! close (funit)

    ! ! seed
    ! seed = [123456, 789012]
    ! call rng1%init(seed, [0, 1])
    ! call rng2%init(seed, [0, 2])

    ! xi(1) = rng1%next_int()
    ! xi(2) = rng1%next_int()

    ! print *, xi(1)
    ! print *, xi(2)

    ! bound = 15
    ! do i = 1, 32
    !   xi(i) = rng1%next_int(bound = bound)
    !   print *, xi(i)
    ! end do

    ! call rng1%next_ints(xi, bound = bound)
    ! call rng2%next_ints(xi2, bound = bound)
    ! do i = 1, 32
    !   print *, xi(i), xi2(i)
    ! end do

    ! call rng1%next_reals(xr)
    ! do i = 1, 32
    !   print *, xr(i)
    ! end do

    call tc%finish()
  end subroutine

end module
