program degen_table_test

  use degen_table_m

  implicit none

  type(degen_table)  :: tab
  real, parameter    :: dd = 1e-3
  real, allocatable  :: dpot(:), eta1(:), eta2(:)
  real               :: cj, dcjdeta(2), dcjddpot, dcj1, dcj2
  real               :: dpotp, dpotm, eta2p, eta2m, eta1p, eta1m, cjp, cjm, tmp, dtmp(2)
  integer, parameter :: NE1 = 401, NE2 = 201, NPOT = 201
  integer            :: i, j, k, funit

  ! block
  !   use math_m, only: ber, dberdx, linspace
  !   integer, parameter :: NN = 10001
  !   real, allocatable  :: x(:)

  !   x = linspace(-100.0, 1000.0, NN)

  !   open (newunit = funit, file = "t", status = "replace", action = "write")
  !   do i = 1, NN
  !     write (funit, "(3ES25.16E3)") x(i), ber(x(i)), dberdx(x(i))
  !   end do
  !   close (funit)
  !   stop
  ! end block
  ! block
  !   use distributions_m, only: fermi_dirac_integral_1h
  !   use fermi_m,         only: fermi_dirac_integral_1h_reg
  !   real :: F, dF, FF
  !   eta1 = linspace(-200.0, 100.0, 1001)

  !   open (newunit = funit, file = "t", status = "replace", action = "write")
  !   do i = 1, 1001
  !     call fermi_dirac_integral_1h(eta1(i), F, dF)
  !     call fermi_dirac_integral_1h_reg(eta1(i), FF, dF)
  !     write (funit, "(3ES25.16E3)") eta1(i), F, FF
  !   end do
  !   close (funit)
  !   stop
  ! end block

  call tab%init(100, 100, 5.0, 1e4, 1e-3, 2e2, 500, 1e3, "/tmp")

  ! dcjdeta(1):
  ! eta1 =  -2.9599999999999994E+001
  ! eta2 =  -4.9949999999999996E+001
  ! dpot =  -5.0000000000000000E+001
  !   1.3462609379799645E-033  1.5046327690525280E-033  1.5046327690525280E-033

  ! dcjddpot:
  ! eta1 =   9.1750000000000096E+000
  ! eta2 =  -4.9949999999999996E+001
  ! dpot =  -5.0000000000000000E+001
  !   9.1227566001656282E-017  7.8098249117738503E-017  7.8099668652987895E-017

  ! dpot = linspace(-50.0, 50.0, NPOT)
  ! eta1 = linspace(-100.0, 10.0, NE1)
  ! eta2 = linspace(-100.0, 10.0, NE2)

  ! do k = 1, NPOT
  !   dpotp = dpot(k) + dd
  !   dpotm = dpot(k) - dd
  !   do j = 1, NE2
  !     eta2p = eta2(j) + dd
  !     eta2m = eta2(j) - dd
  !     do i = 1, NE1
  !       eta1p = eta1(i) + dd
  !       eta1m = eta1(i) - dd
  !       call tab%get([eta1(i), eta2(j)], dpot(k), cj, dcjdeta, dcjddpot)

  !       call tab%get([eta1p, eta2(j)], dpot(k), cjp, dtmp, tmp)
  !       call tab%get([eta1m, eta2(j)], dpot(k), cjm, dtmp, tmp)
  !       dcj1 = (cjp - cj) / dd
  !       dcj2 = (cjp - cjm) / (2 * dd)
  !       if (2*abs(dcj1 - dcj2) < abs(dcjdeta(1) - dcj2) .and. abs((dcjdeta(1) - dcj2) / dcj2) > 1e-3) then
  !         print "(A)", "dcjdeta(1): "
  !         print "(A,ES25.16E3)", "  eta1 = ", eta1(i)
  !         print "(A,ES25.16E3)", "  eta2 = ", eta2(j)
  !         print "(A,ES25.16E3)", "  dpot = ", dpot(k)
  !         print "(A,3ES25.16E3)", "  ", dcjdeta(1), dcj1, dcj2
  !         print *
  !       end if

  !       call tab%get([eta1(i), eta2p], dpot(k), cjp, dtmp, tmp)
  !       call tab%get([eta1(i), eta2m], dpot(k), cjm, dtmp, tmp)
  !       dcj1 = (cjp - cj) / dd
  !       dcj2 = (cjp - cjm) / (2 * dd)
  !       if (2*abs(dcj1 - dcj2) < abs(dcjdeta(2) - dcj2) .and. abs((dcjdeta(2) - dcj2) / dcj2) > 1e-3) then
  !         print "(A)", "dcjdeta(2): "
  !         print "(A,ES25.16E3)", "  eta1 = ", eta1(i)
  !         print "(A,ES25.16E3)", "  eta2 = ", eta2(j)
  !         print "(A,ES25.16E3)", "  dpot = ", dpot(k)
  !         print "(A,3ES25.16E3)", "  ", dcjdeta(2), dcj1, dcj2
  !         print *
  !       end if

  !       call tab%get([eta1(i), eta2(j)], dpotp, cjp, dtmp, tmp)
  !       call tab%get([eta1(i), eta2(j)], dpotm, cjm, dtmp, tmp)
  !       dcj1 = (cjp - cj) / dd
  !       dcj2 = (cjp - cjm) / (2 * dd)
  !       if (2*abs(dcj1 - dcj2) < abs(dcjddpot - dcj2) .and. abs((dcjddpot - dcj2) / dcj2) > 1e-3) then
  !         print "(A)", "dcjddpot: "
  !         print "(A,ES25.16E3)", "  eta1 = ", eta1(i)
  !         print "(A,ES25.16E3)", "  eta2 = ", eta2(j)
  !         print "(A,ES25.16E3)", "  dpot = ", dpot(k)
  !         print "(A,3ES25.16E3)", "  ", dcjddpot, dcj1, dcj2
  !         print *
  !       end if
  !     end do
  !   end do
  ! end do
  ! stop


  ! block
  !   real :: eta(2), dpot, j, djdeta(2), djddpot

  !   eta  = [-65.0, -17.0]
  !   dpot = 48.0001

  !   eta1p = eta(1) + dd
  !   eta1m = eta(1) - dd
  !   eta2p = eta(2) + dd
  !   eta2m = eta(2) - dd
  !   dpotp = dpot + dd
  !   dpotm = dpot - dd

  !   call tab%get(eta, dpot, j, djdeta, djddpot)
  !   print *, j

  !   call tab%get([eta1p, eta(2)], dpot, cjp, dtmp, tmp)
  !   call tab%get([eta1m, eta(2)], dpot, cjm, dtmp, tmp)
  !   dcj1 = (cjp - j) / dd
  !   dcj2 = (cjp - cjm) / (2 * dd)
  !   print *, djdeta(1), dcj1, dcj2

  !   call tab%get([eta(1), eta2p], dpot, cjp, dtmp, tmp)
  !   call tab%get([eta(1), eta2m], dpot, cjm, dtmp, tmp)
  !   dcj1 = (cjp - j) / dd
  !   dcj2 = (cjp - cjm) / (2 * dd)
  !   print *, djdeta(2), dcj1, dcj2

  !   call tab%get(eta, dpotp, cjp, dtmp, tmp)
  !   call tab%get(eta, dpotm, cjm, dtmp, tmp)
  !   dcj1 = (cjp - j) / dd
  !   dcj2 = (cjp - cjm) / (2 * dd)
  !   print *, djddpot, dcj1, dcj2

  !   stop
  ! end block

  ! dpot = 20.0

  ! eta = linspace(-100.0, 40.0, NE)
  ! open (newunit = funit, file = "eta.csv", status = "replace", action = "write")
  ! do i = 1, NE
  !   write (funit, "(ES25.16E3)") eta(i)
  ! end do
  ! close (funit)

  ! allocate (cc(NE,1))
  ! do i = 1, NE
  !   call tab%get([eta(i), 15.0], dpot, cc(i,1), dtmp, tmp)
  ! end do

  ! open (newunit = funit, file = "t", status = "replace", action = "write")
  ! do i = 1, NE
  !   write (funit, "(2ES25.16E3)") eta(i), cc(i,1)
  ! end do
  ! close (funit)

  ! eta1 = linspace(-18.0, -10.0, NE1)
  ! eta2 = linspace(13.0, 19.0, NE2)
  ! open (newunit = funit, file = "eta1.csv", status = "replace", action = "write")
  ! do i = 1, NE1
  !   write (funit, "(ES25.16E3)") eta1(i)
  ! end do
  ! close (funit)
  ! open (newunit = funit, file = "eta2.csv", status = "replace", action = "write")
  ! do i = 1, NE2
  !   write (funit, "(ES25.16E3)") eta2(i)
  ! end do
  ! close (funit)

  ! allocate (cc(NE1,NE2))
  ! do j = 1, NE2; do i = 1, NE1
  !   call tab%get([eta1(i), eta2(j)], dpot, cc(i,j), dtmp, tmp)
  ! end do; end do

  ! open (newunit = funit, file = "t", status = "replace", action = "write")
  ! do j = 1, NE2
  !   do i = 1, NE1
  !     write (funit, "(ES25.16E3)", advance = "no") cc(i,j)
  !   end do
  !   write (funit,*)
  ! end do
  ! close (funit)
end program
