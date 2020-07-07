#include "../macro.f90.inc"

module radau5_m
  use lapack95
  use ode_m
  implicit none

  private
  public :: ode_options, ode_result, radau5

  ! radau5 IIa parameters
  real, parameter, private :: C(3) = [ &
    (4.0 - sqrt(6.0)) / 10.0,          &
    (4.0 + sqrt(6.0)) / 10.0,          &
    1.0                                &
  ]
  real, parameter, private :: A(3,3) = reshape([ &
    (  88.0 -   7.0 * sqrt(6.0)) /  360.0,       &
    ( 296.0 + 169.0 * sqrt(6.0)) / 1800.0,       &
    (  16.0 -         sqrt(6.0)) /   36.0,       &
    ( 296.0 - 169.0 * sqrt(6.0)) / 1800.0,       &
    (  88.0 +   7.0 * sqrt(6.0)) /  360.0,       &
    (  16.0 +         sqrt(6.0)) /   36.0,       &
    (-  2.0 +   3.0 * sqrt(6.0)) /  225.0,       &
    (-  2.0 -   3.0 * sqrt(6.0)) /  225.0,       &
    (   1.0                    ) /    9.0        &
  ], [ 3, 3 ])

  real, parameter, private :: G0   = 2.74888829595676675854321047154372e-01
  real, parameter, private :: E(3) = [      &
      G0 * (-13.0 - 7.0 * sqrt(6.0)) / 3.0, &
      G0 * (-13.0 + 7.0 * sqrt(6.0)) / 3.0, &
      G0 * (- 1.0                  ) / 3.0  &
  ]

contains

  subroutine radau5(fun, x0, x1, xsmp, U0, P, opt, res)
    !! radau5 ode solver
    procedure(ode_fun)             :: fun
      !! pointer to function to integrate
    real,              intent(in)  :: x0
      !! initial x
    real,              intent(in)  :: x1
      !! final x
    real,              intent(in)  :: xsmp(:)
      !! x sample points
    real,              intent(in)  :: U0(:)
      !! initial state
    real,              intent(in)  :: P(:)
      !! parameters
    type(ode_options), intent(in)  :: opt
      !! solver options
    type(ode_result),  intent(out) :: res
      !! output result object

    ! call base solver with radau5 kernel (3 stages)
    call ode_solve(radau5_kernel, 3, fun, x0, x1, xsmp, U0, P, opt, res)
  end subroutine

  subroutine radau5_kernel(fun, xold, x, dxk, Uk, dUkdQ, fk, dfkdUk, dfkdP, polyk, &
    P, opt, dxn, Un, dUndQ, fn, dfndUn, dfndP, polyn, dpolyndQ, err, status)
    !! kernel for radau5 ode solver
    procedure(ode_fun)               :: fun
      !! function to integrate
    real,              intent(in)    :: xold
      !! initial x position of previous step
    real,              intent(in)    :: x
      !! initial x position of step
    real,              intent(in)    :: dxk
      !! stepsize
    real,              intent(in)    :: Uk(:)
      !! state at beginning of step
    real,              intent(in)    :: dUkdQ(:,:)
      !! derivatives of Uk wrt Q = [U0; P]
    real,              intent(in)    :: fk(:)
      !! dUk/dx
    real,              intent(in)    :: dfkdUk(:,:)
      !! derivatives of fk wrt Uk
    real,              intent(in)    :: dfkdP(:,:)
      !! derivatives of fk wrt P
    real,              intent(in)    :: polyk(:,:)
      !! interpolation polynomial (from previous step)
    real,              intent(in)    :: P(:)
      !! parameters
    type(ode_options), intent(in)    :: opt
      !! options
    real,              intent(out)   :: dxn
      !! output new stepsize
    real,              intent(out)   :: Un(:)
      !! output state at end of step
    real,              intent(out)   :: dUndQ(:,:)
      !! output derivatives of Un wrt Q = [U0; P]
    real,              intent(out)   :: fn(:)
      !! output dUn/dx
    real,              intent(out)   :: dfndUn(:,:)
      !! output derivatives of fn wrt Un
    real,              intent(out)   :: dfndP(:,:)
      !! output derivatives of fn wrt P
    real,              intent(out)   :: polyn(:,:)
      !! output interpolation polynomial for this step
    real,              intent(out)   :: dpolyndQ(:,:,:)
      !! output derivatives of polyn wrt Q = [U0; P]
    real,              intent(inout) :: err
      !! output scalar error estimate
    logical,           intent(inout) :: status
      !! output success (true) or fail (false)

    ! local variables
    integer :: nU, nP

    IGNORE(dfkdP)

    ! get system size, number of parameters and number of stages
    nU = size(Uk)
    nP = size(P)

    block
      integer :: it
      real    :: newt_err, newt_err0, err_k
      real    :: f(nU,3), dfdz(nU,nU,3), dfdP(nU,nP,3)
      real    :: z(nU*3), h(nU*3,1), dhdz(nU*3,nU*3)
      real    :: dhdQ(nU*3,nU+nP), dzdQ(nU*3,nU+nP)

      ! reset output
      dxn      = 0
      Un       = 0
      dUndQ    = 0
      fn       = 0
      dfndUn   = 0
      polyn    = 0
      dpolyndQ = 0

      ! initial guess for z (use interpolation polynomial from last step)
      call initial_guess(nU, xold, x, dxk, polyk, z)

      ! newton iteration
      newt_err0 = 2e99
      newt_err  = 1e99
      it = 0
      do while (newt_err > opt%newton_rtol)
        ! count iterations
        it = it + 1

        ! fail if too many iterations
        if (it > opt%newton_max_it) then
          ! try half step size next time
          dxn = 0.5 * dxk

          ! fail
          status = .false.
          return
        end if

        ! evaluate function and derivatives at all stages
        call eval_f(nU, fun, x, dxk, Uk, z, P, f, dfdz)

        ! get h and dhdz
        call eval_h(nU, dxk, f, dfdz, z, h(:,1), dhdz)

        ! solve
        call gesv(dhdz, h)

        ! get error
        newt_err0 = newt_err
        newt_err  = maxval(abs(h(:,1)) / (abs(z) + opt%newton_atol / opt%newton_rtol))

        ! update z
        z = z - h(:,1)
      end do

      ! get U at beginning of next step
      Un = Uk + z((2*nU+1):(3*nU))

      ! get error estimate
      call error_estimate(nU, fun, x, dxk, Uk, fk, dfkdUk, P, opt, z, Un, status, err_k)

      ! step size prediction
      call step_size_prediction(x, x - xold, dxk, err, err_k, opt, it, dxn)

      ! fail if err > 1
      if (err_k > 1.0) then
        status = .false.
      else
        ! success
        status = .true.
        err    = err_k

        ! derivatives
        call eval_f(nU, fun, x, dxk, Uk, z, P, f, dfdz, dfdP = dfdP)
        call eval_h(nU, dxk, f, dfdz, z, h(:,1), dhdz)
        call eval_dhdQ(nU, nP, dxk, dfdz, dfdP, dUkdQ, dhdQ)

        ! calculate dzdQ = - (dhdz)^(-1) * dhdQ
        call gesv(dhdz, dhdQ)
        dzdQ = - dhdQ

        ! get dUndQ
        dUndQ = dUkdQ + dzdQ((nU*2+1):(nU*3),:)

        ! set fn, dfndUn and dfndP
        fn     = f(:,3)
        dfndUn = dfdz(:,:,3)
        dfndP  = dfdP(:,:,3)

        ! new polynomial interpolation
        call interpolate(dxk, z(1:nU), z((nU+1):(2*nU)), z((2*nU+1):(3*nU)), &
                         dzdQ(1:nU,:), dzdQ((nU+1):(2*nU),:), dzdQ((2*nU+1):(3*nU),:), polyn, dpolyndQ)
      end if
    end block
  end subroutine

  subroutine initial_guess(nU, xold, x, dxk, polyk, z)
    integer, intent(in)  :: nU
      !! system size
    real,    intent(in)  :: xold
      !! initial x position of previous step
    real,    intent(in)  :: x
      !! initial x position of step
    real,    intent(in)  :: dxk
      !! stepsize
    real,    intent(in)  :: polyk(:,:)
      !! interpolation polynomial (from previous step)
    real,    intent(out) :: z(:)
      !! output initial guess for z

    ! local variables
    integer :: i, j, i0, i1

    ! use interpolation polynomial from previous step to extrapolate z(x) = U(x) - U_k
    i1 = 0
    do i = 1, 3
      i0 = i1 + 1
      i1 = i1 + nU

      z(i0:i1) = 0
      do j = 1, 3
        z(i0:i1) = z(i0:i1) + polyk(:,j) * (x + C(i) * dxk - xold)**j
      end do
    end do
  end subroutine

  subroutine eval_f(nU, fun, x, dxk, Uk, z, P, f, dfdz, dfdP)
    integer,        intent(in)  :: nU
      !! system size
    procedure(ode_fun)          :: fun
      !! function to integrate
    real,           intent(in)  :: x
      !! beginning of step
    real,           intent(in)  :: dxk
      !! step size
    real,           intent(in)  :: Uk(:)
      !! state at beginning of step
    real,           intent(in)  :: z(:)
      !! delta states
    real,           intent(in)  :: P(:)
      !! parameters
    real,           intent(out) :: f(:,:)
      !! output function values
    real,           intent(out) :: dfdz(:,:,:)
      !! output derivatives of f wrt z
    real, optional, intent(out) :: dfdP(:,:,:)
      !! output derivatives of f wrt parameters

    ! local variables
    integer :: i, i0, i1

    i1 = 0
    do i = 1, 3
      i0 = i1 + 1
      i1 = i1 + nU
      if (present(dfdP)) then
        call fun(x + C(i) * dxk, Uk + z(i0:i1), P, f = f(:,i), dfdU = dfdz(:,:,i), dfdP = dfdP(:,:,i))
      else
        call fun(x + C(i) * dxk, Uk + z(i0:i1), P, f = f(:,i), dfdU = dfdz(:,:,i))
      end if
    end do
  end subroutine

  subroutine eval_h(nU, dxk, f, dfdz, z, h, dhdz)
    integer, intent(in)  :: nU
      !! system size
    real,    intent(in)  :: dxk
      !! step size
    real,    intent(in)  :: f(:,:)
      !! function values for all stages
    real,    intent(in)  :: dfdz(:,:,:)
      !! derivatives of f wrt z
    real,    intent(in)  :: z(:)
      !! delta states
    real,    intent(out) :: h(:)
      !! output residuals
    real,    intent(out) :: dhdz(:,:)
      !! output derivatives of h wrt z

    ! local variables
    integer :: i, j, i0, i1, j0, j1

    dhdz = 0

    i1 = 0
    do i = 1, 3
      i0 = i1 + 1
      i1 = i1 + nU

      h(i0:i1) = z(i0:i1)
      do j = i0, i1
        dhdz(j,j) = 1
      end do

      j1 = 0
      do j = 1, 3
        j0 = j1 + 1
        j1 = j1 + nU

        h(i0:i1) = h(i0:i1) - dxk * A(i,j) * f(:,j)
        dhdz(i0:i1,j0:j1) = dhdz(i0:i1,j0:j1) - dxk * A(i,j) * dfdz(:,:,j)
      end do
    end do
  end subroutine

  subroutine eval_dhdQ(nU, nP, dxk, dfdz, dfdP, dUkdQ, dhdQ)
    integer, intent(in)  :: nU
      !! system size
    integer, intent(in)  :: nP
      !! number of parameters
    real,    intent(in)  :: dxk
      !! step size
    real,    intent(in)  :: dfdz(:,:,:)
      !! derivatives of f wrt z (= dfdUk)
    real,    intent(in)  :: dfdP(:,:,:)
      !! derivatives of f wrt parameters
    real,    intent(in)  :: dUkdQ(:,:)
      !! derivatives of Uk wrt Q = [U0; P]
    real,    intent(out) :: dhdQ(3*nU,nU+nP)
      !! output derivatives of h wrt Q = [U0; P]

    ! local variables
    integer :: i, j, i0, i1
    real    :: dhdUk(nU,nU), dhdP(nU,nP)

    i1 = 0
    do i = 1, 3
      i0 = i1 + 1
      i1 = i1 + nU

      dhdUk = 0
      dhdP  = 0
      do j = 1, 3
        dhdUk = dhdUk - dxk * A(i,j) * dfdz(:,:,j)
        dhdP  = dhdP  - dxk * A(i,j) * dfdP(:,:,j)
      end do

      dhdQ(i0:i1,:) = matmul(dhdUk, dUkdQ)
      dhdQ(i0:i1,(nU+1):(nU+nP)) = dhdQ(i0:i1,(nU+1):(nU+nP)) + dhdP
    end do
  end subroutine

  subroutine error_estimate(nU, fun, x, dxk, Uk, fk, dfkdUk, P, opt, z, Un, status, err)
    integer,           intent(in)  :: nU
      !! system size
    procedure(ode_fun)             :: fun
      !! function to integrate
    real,              intent(in)  :: x
      !! initial x position of step
    real,              intent(in)  :: dxk
      !! stepsize
    real,              intent(in)  :: Uk(:)
      !! state at beginning of step
    real,              intent(in)  :: fk(:)
      !! dUk/dx
    real,              intent(in)  :: dfkdUk(:,:)
      !! derivatives of fk wrt Uk
    real,              intent(in)  :: P(:)
      !! parameters
    type(ode_options), intent(in)  :: opt
      !! options
    real,              intent(in)  :: z(:)
      !! delta states
    real,              intent(in)  :: Un(:)
      !! state at end of step
    logical,           intent(in)  :: status
      !! old status
    real,              intent(out) :: err
      !! output scalar error estimate

    ! local variables
    integer :: i, i0, i1, ipiv(nU)
    real    :: eU(nU,1), eU_tmp(nU), eUmat(nU,nU)

    ! delta U
    eU(:,1) = G0 * dxk * fk
    i1 = 0
    do i = 1, 3
      i0 = i1 + 1
      i1 = i1 + nU

      eU(:,1) = eU(:,1) + E(i) * z(i0:i1)
    end do

    ! matrix to refine error estimate
    eUmat = - dxk * G0 * dfkdUk
    do i = 1, nU
      eUmat(i,i) = eUmat(i,i) + 1.0
    end do

    ! solve
    call gesv(eUmat, eU, ipiv = ipiv)

    ! if last step was rejected, refine error estimate again
    if (.not. status) then
      call fun(x, Uk + eU(:,1), P, f = eU_tmp)
      eU(:,1) = G0 * dxk * eU_tmp
      i1 = 0
      do i = 1, 3
        i0 = i1 + 1
        i1 = i1 + nU

        eU(:,1) = eU(:,1) + E(i) * z(i0:i1)
      end do

      ! solve again
      call getrs(eUmat, ipiv, eU)
    end if

    ! get scalar error estimate
    err = 0
    do i = 1, nU
      err = err + (eU(i,1) / (opt%atol(i) + max(abs(Uk(i)), abs(Un(i))) * opt%rtol(i)))**2
    end do
    err = sqrt(err / nU)
  end subroutine

  subroutine interpolate(dxk, z1, z2, z3, dz1dQ, dz2dQ, dz3dQ, polyn, dpolyndQ)
    real, intent(in)  :: dxk
      !! step size
    real, intent(in)  :: z1(:)
      !! delta state at 1st stage
    real, intent(in)  :: z2(:)
      !! delta state at 2nd stage
    real, intent(in)  :: z3(:)
      !! delta state at 3rd stage
    real, intent(in)  :: dz1dQ(:,:)
      !! derivatives of z1 wrt Q = [U0; P]
    real, intent(in)  :: dz2dQ(:,:)
      !! derivatives of z2 wrt Q = [U0; P]
    real, intent(in)  :: dz3dQ(:,:)
      !! derivatives of z3 wrt Q = [U0; P]
    real, intent(out) :: polyn(:,:)
      !! output interpolation polynomial for this step
    real, intent(out) :: dpolyndQ(:,:,:)
      !! output derivatives of polyn wrt Q = [U0; P]

    polyn(:,1) =  (13*(z1 + z2) +    z3 +  7*sqrt(6.0)*(z1 - z2)) / (3 * dxk   )
    polyn(:,2) = -(23*(z1 + z2) +  8*z3 + 22*sqrt(6.0)*(z1 - z2)) / (3 * dxk**2)
    polyn(:,3) =  (10*(z1 + z2) + 10*z3 + 15*sqrt(6.0)*(z1 - z2)) / (3 * dxk**3)

    dpolyndQ(:,:,1) =  (13*(dz1dQ + dz2dQ) +    dz3dQ +  7*sqrt(6.0)*(dz1dQ - dz2dQ)) / (3 * dxk   )
    dpolyndQ(:,:,2) = -(23*(dz1dQ + dz2dQ) +  8*dz3dQ + 22*sqrt(6.0)*(dz1dQ - dz2dQ)) / (3 * dxk**2)
    dpolyndQ(:,:,3) =  (10*(dz1dQ + dz2dQ) + 10*dz3dQ + 15*sqrt(6.0)*(dz1dQ - dz2dQ)) / (3 * dxk**3)
  end subroutine

  subroutine step_size_prediction(x, dx_old, dxk, err_old, err, opt, it, dxn)
    real,              intent(in)  :: x
      !! initial x position of step
    real,              intent(in)  :: dx_old
      !! old step size from step before
    real,              intent(in)  :: dxk
      !! step size
    real,              intent(in)  :: err_old
      !! old error estimate from step before
    real,              intent(in)  :: err
      !! error estimate
    type(ode_options), intent(in)  :: opt
      !! ode solver options
    integer,           intent(in)  :: it
      !! number of Newton iterations needed
    real,              intent(out) :: dxn
      !! output new step size

    ! local variables
    real :: fac, dxn1, dxn2

    ! get safety factor
    fac = 0.9 * (2 * opt%newton_max_it + 1.0) / (2 * opt%newton_max_it + it)

    ! estimate new step size by two different methods
    dxn1 = fac * dxk * err**(-0.25)
    dxn2 = dxn1 * dxk / dx_old * (err_old / err)**(0.25)

    ! take minimum of both results
    dxn = min(dxn1, dxn2)

    if (abs(dxn) < opt%min_rx * abs(x+dxk)) then
      dxn = sign(opt%min_rx * abs(x+dxk), dxn)
    end if
  end subroutine

end module