m4_include(util/macro.f90.inc)

module degen_test_module_m

  use distributions_m, only: fermi_dirac_integral_1h, fermi_dirac_integral_m1h, inv_fermi_dirac_integral_1h
  use dual_m
  use device_params_m
  use math_m
  use normalization_m
  use newton_m
  use radau5_m

  implicit none

  integer, parameter :: ci = CR_ELEC
  integer, parameter :: N = 1001
  real               :: ch
  real               :: edos

contains

  subroutine degen_test_eval()
    integer            :: i, idens, funit_dj
    real               :: len, Efield0, Efield1, pot(2), mob
    real               :: j, djdpot(2), djdmob
    real, allocatable  :: Efield(:), dens(:,:), djddens(:,:,:)

    call init_normconst(300.0)

    ch   = CR_CHARGE(ci)
    edos = norm(3e19, "1/cm^3")
    len  = norm(0.1, "um")

    Efield0    = - 20 / len
    Efield1    =   20 / len
    Efield     = linspace(Efield0, Efield1, N)
    allocate (dens(2,4), djddens(2,0:4,N))
    dens(:,1)  = norm([1e18, 1e0 ], "1/cm^3")
    dens(:,2)  = norm([1e19, 1e18], "1/cm^3")
    dens(:,3)  = norm([1e20, 1e18], "1/cm^3")
    dens(:,4)  = norm([1e21, 1e18], "1/cm^3")
    mob        = norm(1430.0, "cm^2/V/s")

    ! $omp parallel do schedule(dynamic) private(i,idens,pot,j,djdpot,djdmob) shared(len,dens,mob,djddens)
    do i = 1, N
      print *, i
      pot = [1.0, -1.0] * 0.5 * len * Efield(i)
      call eval_sg(len, pot, dens(:,1), mob, j, djdpot, djddens(:,0,i), djdmob)
      do idens = 1, size(dens, 2)
        call eval_degen(len, pot, dens(:,idens), mob, j, djdpot, djddens(:,idens,i), djdmob)
      end do
    end do
    ! $omp end parallel do
    print *, "done"

    open (newunit = funit_dj, file = "../tex/SISPAD2023/dj.csv", status = "replace", action = "write")
    write (funit_dj, "(A)") "dphi j1_sg j2_sg j1_1e18 j2_1e18 j1_1e19 j2_1e19 j1_1e20 j2_1e20 j1_1e21 j2_1e21"
    do i = 1, N
      pot = [1.0, -1.0] * 0.5 * len * Efield(i)
      write (funit_dj, "(3ES24.16)", advance = "no") pot(2) - pot(1), djddens(:,0,i) * len / mob
      do idens = 1, size(dens, 2)
        write (funit_dj, "(2ES24.16)", advance = "no") djddens(:,idens,i) * len / mob
      end do
      write (funit_dj, *)
    end do
    close (funit_dj)
  end subroutine

  subroutine eval_sg(len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! Scharfetter-Gummel stabilization
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    real :: ber1, ber2, dber1, dber2

    ber1  = ber(ch * (pot(1) - pot(2)))
    ber2  = ber(ch * (pot(2) - pot(1)))
    dber1 = ch * dberdx(ch * (pot(1) - pot(2)))
    dber2 = ch * dberdx(ch * (pot(2) - pot(1)))

    j = - mob * (ber1 * dens(2) - ber2 * dens(1)) / len

    djdpot(1) = - mob * (dber1 * dens(2) + dber2 * dens(1)) / len
    djdpot(2) =   mob * (dber1 * dens(2) + dber2 * dens(1)) / len

    djddens(1)  =   mob * ber2 / len
    djddens(2)  = - mob * ber1 / len

    djdmob = -(ber1 * dens(2) - ber2 * dens(1)) / len
  end subroutine

  subroutine eval_degen_old(len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! Generalized Scharfetter-Gummel stabilization for degenerate case
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    integer            :: num_eval
    real               :: jsg, djsgdpot(2), djsgddens(2), djsgdmob
    real               :: jj0, jj, djjdp(3), nn(2)
    type(newton1D_opt) :: newt_opt
    type(ode_options)  :: ode_opt
    type(ode_result)   :: ode_res1, ode_res2

    ! initial guess
    call eval_sg(len, pot, dens, mob, jsg, djsgdpot, djsgddens, djsgdmob)
    jj0 = jsg * len / (mob * edos)
    nn  = dens / edos

    ! solve with newton iteration
    call ode_opt%init(1, atol = [minval(nn*1e-14)], rtol = [1e-10], max_rejected = 50)
    call newt_opt%init()
    call newton1D(newton_fun, [pot(2) - pot(1), nn(1), nn(2)], newt_opt, jj0, jj, djjdp)

    ! extract solution + derivatives
    j       = jj * mob * edos / len
    djdpot  = [-1.0, 1.0] * djjdp(1) * mob * edos / len
    djddens = djjdp(2:3) * mob / len
    djdmob  = jj * edos / len

  contains

    subroutine newton_fun(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (jj)
      real,              intent(in)  :: p(:)
        !! parameters (pot(2) - pot(1), nn(1), nn(2))
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: nnr(2), djj(2), dpot(2), dnn0(2)

      if (p(2) < p(3) * 1e-3) then
        ! solve ode from left to right
        call radau5(ode_fun, 0.0, 1.0, [1.0], [p(2)], [p(1), x], ode_opt, ode_res1)
        nnr( 1) = ode_res1%Usmp(    1,  1)
        djj( 1) = ode_res1%dUsmpdP( 1,2,1)
        dpot(1) = ode_res1%dUsmpdP( 1,1,1)
        dnn0(1) = ode_res1%dUsmpdU0(1,1,1)
        nnr( 2) = p(3)
        djj( 2) = 0
        dpot(2) = 0
        dnn0(2) = 1
      elseif (p(3) < p(2) * 1e-3) then
        ! solve ode from right to left
        call radau5(ode_fun, 1.0, 0.0, [0.0], [p(3)], [p(1), x], ode_opt, ode_res2)
        nnr( 1) = p(2)
        djj( 1) = 0
        dpot(1) = 0
        dnn0(1) = 1
        nnr( 2) = ode_res2%Usmp(    1,  1)
        djj( 2) = ode_res2%dUsmpdP( 1,2,1)
        dpot(2) = ode_res2%dUsmpdP( 1,1,1)
        dnn0(2) = ode_res2%dUsmpdU0(1,1,1)
      else
        ! solve ode from left to center and from right to center
        call radau5(ode_fun, 0.0, 0.5, [0.5], [p(2)], [p(1), x], ode_opt, ode_res1)
        call radau5(ode_fun, 1.0, 0.5, [0.5], [p(3)], [p(1), x], ode_opt, ode_res2)
        nnr( 1) = ode_res1%Usmp(    1,  1)
        djj( 1) = ode_res1%dUsmpdP( 1,2,1)
        dpot(1) = ode_res1%dUsmpdP( 1,1,1)
        dnn0(1) = ode_res1%dUsmpdU0(1,1,1)
        nnr( 2) = ode_res2%Usmp(    1,  1)
        djj( 2) = ode_res2%dUsmpdP( 1,2,1)
        dpot(2) = ode_res2%dUsmpdP( 1,1,1)
        dnn0(2) = ode_res2%dUsmpdU0(1,1,1)
      end if

      ! residual
      f = nnr(1) - nnr(2)
      if (present(dfdx)) dfdx = djj(1) - djj(2)
      if (present(dfdp)) dfdp = [dpot(1) - dpot(2), dnn0(1), -dnn0(2)]
    end subroutine

    subroutine ode_fun(x, U, P, f, dfdU, dfdP)
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (nn)
      real,           intent(in)  :: P(:)
        !! parameters (pot(2) - pot(1), jj)
      real, optional, intent(out) :: f(:)
        !! output dnn/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt nn
      real, optional, intent(out) :: dfdP(:,:)
        !! output derivative of f wrt P

      real :: eta, deta, Fm1h, dFm12, alpha, dalpha, fsg, dfsg

      m4_ignore(x)

      ! get "degeneracy-factor" (~1 for small densities; < 1 for large densities)
      call inv_fermi_dirac_integral_1h(U(1), eta, deta)
      call fermi_dirac_integral_m1h(eta, Fm1h, dFm12)
      alpha = Fm1h / U(1)
      dalpha = (dFm12 * deta - alpha) / U(1)

      ! Scharfetter-Gummel
      fsg  = - ch * P(1) * U(1) - P(2)
      dfsg = - ch * P(1)

      ! scale Scharfetter-Gummel by degeneracy factor
      if (present(f)) then
        f(1) = alpha * fsg
      end if
      if (present(dfdU)) then
        dfdU(1,1) = dalpha * fsg + alpha * dfsg
      end if
      if (present(dfdp)) then
        dfdp(1,1) = - alpha * ch * U(1)
        dfdp(1,2) = - alpha
      end if
    end subroutine

  end subroutine

  subroutine eval_degen(len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! Generalized Scharfetter-Gummel stabilization for degenerate case
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    integer            :: num_eval
    real               :: jsg, djsgdpot(2), djsgddens(2), djsgdmob
    real               :: jj0, jj, djjdp(3), nn(2)
    type(newton1D_opt) :: newt_opt
    type(ode_options)  :: ode_opt
    type(ode_result)   :: ode_res1, ode_res2

    ! initial guess
    call eval_sg(len, pot, dens, mob, jsg, djsgdpot, djsgddens, djsgdmob)
    jj0 = jsg * len / (mob * edos)
    nn  = dens / edos

    ! solve with newton iteration
    call ode_opt%init(1, atol = [1e-10], rtol = [1e-10], max_rejected = 50)
    call newt_opt%init()
    call newton1D(newton_fun, [nn(1), nn(2), - ch * (pot(2) - pot(1))], newt_opt, jj0, jj, djjdp)

print "(A,ES24.16)", "newton"
print *
! stop

    ! extract solution + derivatives
    j       = jj * mob * edos / len
    djdpot  = [ch, -ch] * djjdp(3) * mob * edos / len
    djddens = djjdp(1:2) * mob / len
    djdmob  = jj * edos / len

  contains

    subroutine newton_fun(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (jj)
      real,              intent(in)  :: p(:)
        !! parameters (nn(1), nn(2), pot(2) - pot(1))
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: nn1(2), dnn1(2,3)

      if (p(1) < p(2) * 1e-3) then
        ! solve ode from left to right
        call solve_ode(0.0, 1.0, p(1), x, p(3), nn1(1), dnn1(1,:))
        nn1(2)      = p(2)
        dnn1(2,1)   = 1.0
        dnn1(2,2:3) = 0.0
      elseif (p(1) > p(2) * 1e-3) then
        ! solve ode from right to left
        call solve_ode(1.0, 0.0, p(2), x, p(3), nn1(2), dnn1(2,:))
        nn1(1)      = p(1)
        dnn1(1,1)   = 1.0
        dnn1(1,2:3) = 0.0
      else
        ! solve ode from left to center and from right to center
        call solve_ode(0.0, 0.5, p(1), x, p(3), nn1(1), dnn1(1,:))
        call solve_ode(1.0, 0.5, p(2), x, p(3), nn1(2), dnn1(2,:))
      end if

      ! residual
      f = nn1(1) - nn1(2)
      if (present(dfdx)) dfdx = dnn1(1,2) - dnn1(2,2)
      if (present(dfdp)) dfdp = [dnn1(1,1), -dnn1(2,1), dnn1(1,3) - dnn1(2,3)]
    end subroutine

    subroutine solve_ode(x0, x1, nn0, jj, dpot, nn1, dnn1)
      real, intent(in)  :: x0
        !! initial x, normalized to edge length
      real, intent(in)  :: x1
        !! final x, normalized to edge length
      real, intent(in)  :: nn0
        !! initial normalized density
      real, intent(in)  :: jj
        !! normalized current density (constant along edge)
      real, intent(in)  :: dpot
        !! normalized potential difference
      real, intent(out) :: nn1
        !! output final normalized density
      real, intent(out) :: dnn1(:)
        !! output derivatives of nn1 [dnn1dnn0, dnn1djj, dnn1ddpot]

      real             :: dx, sx
      type(ode_result) :: ode_res
      type(dual_3)     :: beta, nnsg, nn1_dual

      ! always start at x = 0 and go in positive direction
      if (x1 >= x0) then
        dx = x1 - x0
        sx = 1.0
      else
        dx = x0 - x1
        sx = -1.0
      end if
block
  integer, parameter :: nsmp = 201
  integer :: ii
  real, allocatable :: xsmp(:)

  xsmp = linspace(0.0, 1.0, nsmp)

  call radau5(ode_fun, 0.0, 1.0, xsmp, [0.0], [3.3333333333333333E-02, -6.6236643183010870E+02, -1.8719999999999999E+01], ode_opt, ode_res)

  do ii = 1, nsmp
    print "(2ES24.16)", xsmp(ii), ode_res%Usmp(1,ii)
  end do
  stop
end block

      ! solve ode, start from beta = 1, flip sign of jj and dpot if necessary
      call radau5(ode_fun, 0.0, dx, [dx], [0.0], [nn0, sx*jj, sx*dpot], ode_opt, ode_res)

      ! final beta (first and only sample)
      beta%x     = ode_res%Usmp(1,1)
      beta%dx(1) = ode_res%dUsmpdP(1,1,1)
      beta%dx(2) = ode_res%dUsmpdP(1,2,1) * sx
      beta%dx(3) = ode_res%dUsmpdP(1,3,1) * sx

      ! get final density
      call get_nnsg(dx, nn0, sx * jj, sx * dpot, nnsg%x, nnsg%dx)
      nnsg%dx(2) = nnsg%dx(2) * sx
      nnsg%dx(3) = nnsg%dx(3) * sx
      nn1_dual = exp(beta) * nnsg

      ! extract value + derivatives
      nn1     = nn1_dual%x
      dnn1(1) = nn1_dual%dx(1)
      dnn1(2) = nn1_dual%dx(2)
      dnn1(3) = nn1_dual%dx(3)
    end subroutine

    subroutine ode_fun(x, U, P, f, dfdU, dfdP)
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (b)
      real,           intent(in)  :: P(:)
        !! parameters (nn0, jj, dpot)
      real, optional, intent(out) :: f(:)
        !! output db/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt b
      real, optional, intent(out) :: dfdP(:,:)
        !! output derivative of f wrt P

      real         :: deta, dFm12
      type(dual_4) :: alpha, beta, dpot, eta, ff, Fm1h, jj, nn, nnsg

      call beta%init(U(1), i = 1)
      call jj%init(  P(2), i = 3)
      call dpot%init(P(3), i = 4)

      ! get Scharfetter-Gummel density
      call get_nnsg(x, P(1), P(2), P(3), nnsg%x, nnsg%dx(2:4))
      nnsg%dx(1) = 0

      ! get alpha ("degeneracy"-factor)
      nn = exp(beta) * nnsg
      call inv_fermi_dirac_integral_1h(nn%x, eta%x, deta)
      eta%dx = deta * nn%dx
      call fermi_dirac_integral_m1h(eta%x, Fm1h%x, dFm12)
      Fm1h%dx = dFm12 * eta%dx
      alpha = Fm1h / nn

      ! dbeta/dx
      ff = (alpha - 1.0) * dpot - (alpha * exp(-beta) - 1.0) * jj / nnsg

      ! extract
      if (present(f   )) f(1)        = ff%x
      if (present(dfdU)) dfdU(1,1)   = ff%dx(1)
      if (present(dfdP)) dfdP(1,1:3) = ff%dx(2:4)
    end subroutine

    subroutine get_nnsg(x, nn0, jj, dpot, nn1, dnn1)
      !! get Scharfetter-Gummel density along edge
      real, intent(in)  :: x
        !! normalized x coordinate in [0, 1]
      real, intent(in)  :: nn0
        !! normalized density at x = 0
      real, intent(in)  :: jj
        !! normalized current density
      real, intent(in)  :: dpot
        !! normalized potential difference
      real, intent(out) :: nn1
        !! output density at x
      real, intent(out) :: dnn1(:)
        !! output derivatives of nn1 wrt [nn0, jj, dpot]

      real :: e, em1

      e   = exp(  dpot * x)
      em1 = expm1(dpot * x)
      nn1 = nn0 * e - jj / dpot * em1

      dnn1(1) = e
      dnn1(2) = - em1 / dpot
      dnn1(3) = (nn0 - jj / dpot) * x * e + jj / dpot**2 * em1
    end subroutine

  end subroutine

end module
