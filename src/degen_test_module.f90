m4_include(util/macro.f90.inc)

module degen_test_module_m

  use distributions_m
  use dual_m
  use error_m
  use ieee_arithmetic
  use math_m
  use radau5_m

  implicit none

  real, parameter :: alpha = sqrt(0.125)
  real, parameter :: etaF = -16.0

contains

  subroutine degen_test()
    integer, parameter :: Neta2 = 101

    integer           :: ieta2
    real              :: dpot, eta(2)
    real, allocatable :: eta2(:), j(:), djdeta(:,:), djddpot(:)

    allocate (eta2(Neta2), j(Neta2), djdeta(2,Neta2), djddpot(Neta2))
    eta2 = linspace(-100.0, 100.0, Neta2)

    dpot   = 50.0
    eta(1) = -10.0

    do ieta2 = 1, Neta2
      eta(2) = eta2(ieta2)

      call get_current(eta, dpot, j(ieta2), djdeta(:,ieta2), djddpot(ieta2))

      ! call fermi_dirac_integral_1h(eta2(ieta2), j(ieta2), djddpot(ieta2))
      ! j(ieta2) = dpot * j(ieta2)

      print "(2ES24.16)", eta2(ieta2), j(ieta2)
    end do

  end subroutine

  function get_current_sg(n, dpot) result(j)
    real, intent(in)  :: n(2)
    real, intent(in)  :: dpot
    real              :: j

    j = ber(dpot) * n(1) - ber(-dpot) * n(2)
  end function

  subroutine get_current(eta, dpot, j, djdeta, djddpot)
    real, intent(in)  :: eta(2)
    real, intent(in)  :: dpot
    real, intent(out) :: j
    real, intent(out) :: djdeta(2)
    real, intent(out) :: djddpot

    integer           :: dir
    logical           :: small_eta
    real              :: deta, jmin, jmax, jsgn, n(2), dndeta(2), t, nc
    type(ode_options) :: opt

    deta = eta(2) - eta(1)
    call fermi_dirac_integral_1h(eta(1), n(1), dndeta(1))
    call fermi_dirac_integral_1h(eta(2), n(2), dndeta(2))

    if (all(eta <= etaF)) then
      ! special case: 1/F12(eta) ~ exp(-eta) + sqrt(1/8) -> Bernoulli iteration
      j = get_current_sg(n, dpot)
      small_eta = .true.
      call get_current_newton(j, djdeta, djddpot)
    else
      ! get jmin, jmax by slope (detadx must be equal to deta for some x in [0, 1])
      jmin = min(n(1), n(2)) * abs(dpot - deta)
      jmax = max(n(1), n(2)) * abs(dpot - deta)
      jsgn = sign(1.0, dpot - deta)
      if (jsgn < 0) then
        ! correct sign and swap jmin, jmax
        t    = jmin
        jmin = - jmax
        jmax = - t
      end if

      ! determine shooting direction and reduce j range further if possible
      if (dpot < 0) then
        if (eta(2) > eta(1)) then
          dir  = 1
          jmax = min(jmax, dpot * n(2))
        elseif (eta(2) >= eta(1) + dpot) then
          dir  = 0
          jmin = max(jmin, dpot * n(2))
          jmax = min(jmax, 0.0)
        else
          dir  = -1
          call fermi_dirac_integral_1h(eta(1) + 0.5 * dpot, nc, t)
          jmin = max(jmin, 0.0)
          jmax = min(jmax, abs(dpot - deta) * nc)
        end if
      else
        if (eta(2) < eta(1)) then
          dir  = -1
          jmin = max(jmin, dpot * n(1))
        elseif (eta(2) <= eta(1) + dpot) then
          dir  = 0
          jmin = max(jmin, 0.0)
          jmax = min(jmax, dpot * n(1))
        else
          dir  = 1
          call fermi_dirac_integral_1h(eta(2) - 0.5 * dpot, nc, t)
          jmin = max(jmin, - abs(dpot - deta) * nc)
          jmax = min(jmax, 0.0)
        end if
      end if

      t = 1e-10 * (max(abs(eta(1)), abs(eta(2)), abs(dpot)) + 1)
      if ((abs(deta - dpot) < t) .or. (abs(deta) < t)) then
        ! eta changes approximately linear
        j         = 0.5 * (n(1) + n(2)) * (dpot - deta)
        djddpot   = 0.5 * (n(1) + n(2))
        djdeta(1) = 0.5 * (dndeta(1) * (dpot - deta) + (n(1) + n(2)))
        djdeta(2) = 0.5 * (dndeta(2) * (dpot - deta) - (n(1) + n(2)))
      else
        if (dir == 0) then
          j = jmin
          if (dpot > 0) j = jmax
        else
          j = get_current_sg(n, dpot)
        end if

        call get_current_newton(j, djdeta, djddpot)
      end if
    end if

  contains

    subroutine get_current_newton(j, djdeta, djddpot)
      real,    intent(inout) :: j
      real,    intent(out)   :: djdeta(2)
      real,    intent(out)   :: djddpot

      integer :: it
      logical :: status
      real    :: atol, rtol, dj, err, err0, f, dfdj, dfdeta(2), dfddpot, fmin, fmax, jmin0, jmax0, smin, smax, s

      ! clear
      djdeta  = 0
      djddpot = 0

      if (.not. small_eta) then
        ! init ode options
        call opt%init(1, atol = [1e-16], rtol = [1e-14], max_rejected = 200)

        ! lower bound
        jmin0 = - huge(1.0)
        do while (.true.)
          call residual(jmin, status, f = fmin)
          if (status) exit
          jmin  = jmin + max(abs(jmin) * 1e-15, 1e-100)
          jmin0 = jmin
        end do
        if (fmin == 0) then
          j = jmin
          goto 100
        end if

        ! upper bound
        jmax0 = huge(1.0)
        do while (.true.)
          call residual(jmax, status, f = fmax)
          if (status) exit
          jmax  = jmax - max(abs(jmax) * 1e-15, 1e-100)
          jmax0 = jmax
        end do
        if (fmax == 0) then
          j = jmax
          goto 100
        end if

        ! widen bounds if necessary
        smin = sign(1.0, fmin)
        smax = sign(1.0, fmax)
        do while (smin == smax)
          if (abs(fmin) < abs(fmax)) then
            if (jmin == jmin0) then
              j = jmin
              goto 100
            end if

            j = jmin - max(abs(jmin) * 1e-12, 1e-100)
            do while (.true.)
              call residual(j, status, f = f)
              if (status) exit
              j     = j + max(abs(j) * 1e-15, 1e-100)
              jmin0 = j
            end do
            if (f == 0) goto 100

            jmax = jmin
            fmax = fmin
            jmin = j
            fmin = f
            smin = sign(1.0, fmin)
          else
            if (jmax == jmax0) then
              j = jmax
              goto 100
            end if

            j = jmax + max(abs(jmax) * 1e-10, 1e-100)
            do while (.true.)
              call residual(j, status, f = f)
              if (status) exit
              j     = j - max(abs(j) * 1e-15, 1e-100)
              jmax0 = j
            end do
            if (f == 0) goto 100

            jmin = jmax
            fmin = fmax
            jmax = j
            fmax = f
            smax = sign(1.0, fmax)
          end if
          j = 0.5 * (jmin + jmax)
        end do
      end if

      ! tolerances
      atol = max(1e-16 * abs(j), 1e-100)
      rtol = 1e-14
      err0 = huge(1.0)
      if (.not. small_eta) then
        err  = 0.5 * (jmax - jmin)
      else
        err = 0.5 * err0
      end if

      ! newton iteration
      it = 0
      do while ((err > atol) .and. (err > abs(j) * rtol))
        it = it + 1

        ! bisection
        if (.not. small_eta) then
          if ((j < jmin) .or. (j > jmax) .or. (err0 <= err)) then
            j = 0.5 * (jmin + jmax)
          end if
        end if

        ! evaluate residual
        call residual(j, status, f = f, dfdj = dfdj)
        if (.not. status) call program_error("could not evaluate residual, even though jmin <= j <= jmax")
        if (f == 0) goto 100

        ! calculate newton update and new error
        dj   = f / dfdj
        err0 = err
        err  = abs(dj)

        ! update bounds
        if (.not. small_eta) then
          s = sign(1.0, f)
          if (s * smax > 0) then
            jmax = j
            fmax = f
          else
            jmin = j
            fmin = f
          end if
        end if

        ! update solution
        j = j - dj

        ! print "(I6,2ES24.16)", it, j, err

        ! exit if close to solution
        if (.not. small_eta) then
          if (((jmax - jmin) < 0.5 * abs(jmin + jmax) * rtol) .or. (0.5 * (jmax - jmin) < atol)) then
            j = 0.5 * (jmax + jmin)
            goto 100
          end if
        end if
      end do

      ! calculate derivatives and return
      100 call residual(j, status, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      if (.not. status) call program_error("could not evaluate residual at solution")
      djdeta  = - dfdeta / dfdj
      djddpot = - dfddpot / dfdj
    end subroutine

    subroutine residual(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      if (small_eta) then
        call residual_small_eta(j, status, f = f, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      else
        call residual_shooting(j, status, f = f, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      end if
    end subroutine

    subroutine residual_small_eta(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      real :: arg, B1, B2, dB1, dB2, e(2)

      e = exp(eta)
      status = ieee_is_finite(j)

      ! shifted bernoulli argument
      arg = dpot - alpha * j
      B1  = ber(   -arg)
      B2  = ber(    arg)
      dB1 = dberdx(-arg)
      dB2 = dberdx( arg)

      if (present(f)) then
        f = j - B1 * e(1) + B2 * e(2)
      end if
      if (present(dfdj)) then
        dfdj = 1 + alpha * (dB1 * e(1) + dB2 * e(2))
      end if
      if (present(dfdeta)) then
        dfdeta(1) = - B1 * e(1)
        dfdeta(2) =   B2 * e(2)
      end if
      if (present(dfddpot)) then
        dfddpot = - dB1 * e(1) + dB2 * e(2)
      end if
    end subroutine

    subroutine residual_shooting(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      real :: e, e2, dedeta, de2de, dedj, de2dj, deddpot, de2ddpot, xknee

      status = .true.

      if (dir > 0) then
        ! solve ode from left to right
        call solve_ode(0.0, 1.0, eta(1), j, status, e, dedeta, dedj, deddpot)
        if (.not. status) return

        if (present(f)) then
          f = e - eta(2)
        end if
        if (present(dfdj)) then
          dfdj = dedj
        end if
        if (present(dfdeta)) then
          dfdeta(1) = dedeta
          dfdeta(2) = - 1.0
        end if
        if (present(dfddpot)) then
          dfddpot = deddpot
        end if
      elseif (dir < 0) then
        ! solve ode from right to left
        call solve_ode(1.0, 0.0, eta(2), j, status, e, dedeta, dedj, deddpot)
        if (.not. status) return

        if (present(f)) then
          f = e - eta(1)
        end if
        if (present(dfdj)) then
          dfdj = dedj
        end if
        if (present(dfdeta)) then
          dfdeta(1) = - 1.0
          dfdeta(2) = dedeta
        end if
        if (present(dfddpot)) then
          dfddpot = deddpot
        end if
      else
        xknee = 0.5
        if (dpot < 0) then
          xknee = (eta(2) - eta(1)) / dpot
          call solve_ode(0.0, xknee, eta(1), j, status, e, dedeta, dedj, deddpot)
          if (.not. status) return
          call solve_ode(xknee, 1.0, e, j, status, e2, de2de, de2dj, de2ddpot)
          if (.not. status) return

          if (present(f)) then
            f = e2 - eta(2)
          end if
          if (present(dfdj)) then
            dfdj = de2dj + de2de * dedj
          end if
          if (present(dfdeta)) then
            dfdeta(1) = de2de * dedeta
            dfdeta(2) = - 1.0
          end if
          if (present(dfddpot)) then
            dfddpot = de2ddpot + de2de * deddpot
          end if
        elseif (dpot > 0) then
          xknee = 1 - (eta(2) - eta(1)) / dpot
          call solve_ode(1.0, xknee, eta(2), j, status, e, dedeta, dedj, deddpot)
          if (.not. status) return
          call solve_ode(xknee, 0.0, e, j, status, e2, de2de, de2dj, de2ddpot)
          if (.not. status) return

          if (present(f)) then
            f = e2 - eta(1)
          end if
          if (present(dfdj)) then
            dfdj = de2dj + de2de * dedj
          end if
          if (present(dfdeta)) then
            dfdeta(1) = - 1.0
            dfdeta(2) = de2de * dedeta
          end if
          if (present(dfddpot)) then
            dfddpot = de2ddpot + de2de * deddpot
          end if
        end if

      end if
    end subroutine

    recursive subroutine solve_ode(x0, x1, eta0, j, status, eta1, deta1deta0, deta1dj, deta1ddpot)
      real,    intent(in)  :: x0
      real,    intent(in)  :: x1
      real,    intent(in)  :: eta0
      real,    intent(in)  :: j
      logical, intent(out) :: status
      real,    intent(out) :: eta1
      real,    intent(out) :: deta1deta0
      real,    intent(out) :: deta1dj
      real,    intent(out) :: deta1ddpot

      real             :: deta1dx(1)
      type(dual_3)     :: dl_B, dl_dpot, dl_eta0, dl_eta1, dl_j, dl_t, dl_t1, dl_t2, dl_xF
      type(ode_result) :: result

      status = .true.

      if (eta0 < etaF) then
        ! use dual numbers to avoid most of the manual derivative calculations
        call dl_eta0%init(eta0, i = 1)
        call dl_dpot%init(dpot, i = 2)
        call dl_j%init(      j, i = 3)
        dl_t = dl_dpot - alpha * dl_j

        ! find x where eta crosses etaF
        dl_t1 = dl_t * exp(   etaF) - dl_j
        dl_t2 = dl_t * exp(dl_eta0) - dl_j
        if (((dl_t1%x > 0) .and. (dl_t2%x > 0)) .or. ((dl_t1%x < 0) .and. (dl_t2%x < 0))) then
          if (abs(dl_t%x) < 1e-6) then
            dl_xF = x0 +      (exp(    dl_eta0) - exp(    etaF)) /      dl_j     &
              &   + dl_t    * (exp(2 * dl_eta0) - exp(2 * etaF)) / (2 * dl_j**2) &
              &   + dl_t**2 * (exp(3 * dl_eta0) - exp(3 * etaF)) / (3 * dl_j**3)
          else
            dl_xF = x0 + log(dl_t1 / dl_t2) / dl_t
          end if

          if (((x1 > x0) .and. (dl_xF%x >= x0) .and. (dl_xF%x < x1)) .or. &
            & ((x1 < x0) .and. (dl_xF%x <= x0) .and. (dl_xF%x > x1))) then
            ! go from (xF,etaF) to (x1,eta1) using ode solver
            call solve_ode(dl_xF%x, x1, etaF, j, status, eta1, deta1deta0, deta1ddpot, deta1dj)
            if (.not. status) return
            call detadx(x1, [eta1], [dpot, j], status, f = deta1dx)
            if (.not. status) return
            deta1deta0 =            - deta1dx(1) * dl_xF%dx(1)
            deta1ddpot = deta1ddpot - deta1dx(1) * dl_xF%dx(2)
            deta1dj    = deta1dj    - deta1dx(1) * dl_xF%dx(3)
          else
            ! stay below -16
            dl_B%x  = ber(   dl_t%x * (x1 - x0))
            dl_B%dx = dberdx(dl_t%x * (x1 - x0)) * dl_t%dx
            dl_eta1 = dl_eta0 + log1p((dl_t - dl_j * exp(-dl_eta0)) * (x1 - x0) / dl_B)

            eta1       = dl_eta1%x
            deta1deta0 = dl_eta1%dx(1)
            deta1ddpot = dl_eta1%dx(2)
            deta1dj    = dl_eta1%dx(3)
          end if
        else
          ! stay constant (avoid numerical issues)
          eta1       = eta0
          deta1deta0 = 1
          deta1ddpot = 0
          deta1dj    = 0
        end if
      else
        ! go from (x0,eta0) to (x1,eta1) using ode solver
        call radau5(detadx, x0, x1, [x1], [eta0], [dpot, j], opt, status, result)
        if (.not. status) return
        eta1       = result%Usmp(    1,  1)
        deta1deta0 = result%dUsmpdU0(1,1,1)
        deta1ddpot = result%dUsmpdP( 1,1,1)
        deta1dj    = result%dUsmpdP( 1,2,1)
      end if
    end subroutine

    subroutine detadx(x, U, p, status, f, dfdU, dfdp)
      !! ode right-hand side
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (eta)
      real,           intent(in)  :: p(:)
        !! parameters (dpot, j)
      logical,        intent(out) :: status
        !! success/fail
      real, optional, intent(out) :: f(:)
        !! output deta/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt eta
      real, optional, intent(out) :: dfdp(:,:)
        !! output derivative of f wrt P

      real :: F12, dF12

      m4_ignore(x)

      call fermi_dirac_integral_1h(U(1), F12, dF12)
      status = ieee_is_finite(F12)

      if (present(f)) then
        f(1) = p(1) - p(2) / F12
      end if
      if (present(dfdU)) then
        dfdU(1,1) = p(2) / F12**2 * dF12
      end if
      if (present(dfdp)) then
        dfdp(1,1) = 1.0
        dfdp(1,2) = - 1.0 / F12
      end if
    end subroutine

  end subroutine

end module
