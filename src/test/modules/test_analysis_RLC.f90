module test_analysis_RLC_m

  use capacitor_m,     only: capacitor
  use circuit_m,       only: circuit
  use eigenvalues_m,   only: eigenvalues
  use ground_m,        only: ground
  use inductor_m,      only: inductor
  use input_src_m,     only: const_src
  use math_m,          only: PI, linspace
  use normalization_m, only: norm, init_normconst
  use resistor_m,      only: resistor
  use small_signal_m,  only: small_signal
  use steady_state_m,  only: steady_state
  use string_m,        only: new_string
  use test_case_m,     only: test_case
  use transient_m,     only: transient, TRANS_TRAPZ, TRANS_TRBDF2
  use vsource_m,       only: vsource

  implicit none

  private
  public test_analysis_RLC

  real, parameter :: RTOL = 1e-12
  real, parameter :: ATOL = 1e-10

contains

  subroutine test_analysis_RLC
    type(test_case) :: tc

    type(circuit)   :: crt
    type(ground)    :: GND
    type(resistor)  :: R
    type(inductor)  :: L
    type(capacitor) :: C
    type(vsource)   :: V1

    complex, allocatable :: s(:), Y_exp(:), eig_exp(:)
    integer              :: i, Nt
    real                 :: Vinp, Rval, Lval, Cval, wres, alpha, xi, tmax
    real,    allocatable :: t_exp(:), I_exp(:)
    type(const_src)      :: src
    type(eigenvalues)    :: ev
    type(small_signal)   :: sm
    type(steady_state)   :: ss
    type(transient)      :: tr

    call init_normconst(300.0)

    Vinp  = norm(2.0, "V")
    Rval  = norm(0.4, "Ohm")
    Lval  = norm(1e-9, "H")
    Cval  = norm(2e-9, "F")
    wres  = 1 / sqrt(Lval*Cval)
    alpha = Rval / 2.0 / Lval
    xi    = Rval / 2.0 * sqrt(Cval/Lval)

    call crt%init("circuit")
    call GND%init(crt, "GND", "N0")
    call V1%init(crt, "V1", "N0", "N1")
    call R%init(crt, "R", "N1", "N2", Rval)
    call L%init(crt, "L", "N2", "N3", Lval)
    call C%init(crt, "C", "N3", "N0", Cval)
    call crt%init_final()

    call tc%init("analysis RLC")

    !! ----------------------- steady state: const input voltage Vinp --------------------------------------------------
    call src%init([Vinp])
    call ss%init(crt%sys)
    ! linear => maximal 2 iterations
    call ss%set_newton_params(max_it = 2)
    call ss%init_output([new_string("N0.V"), new_string("N1.V"), new_string("N2.V"), new_string("N3.V")], "ss_volt.fbs")
    call ss%run(input = src)

    call tc%assert_eq(crt%nodes%n, 4, "RLC circuit: number of nodes")
    ! expected result: all voltage drops over the capacitor
    do i = 1, 4
      if (crt%nodes%d(i)%p%volt%name == "N0.V") then
        call tc%assert_eq(crt%nodes%d(i)%p%volt%get(), [0.0], RTOL, ATOL, "steady-state: ground node")
      else
        call tc%assert_eq(crt%nodes%d(i)%p%volt%get(), [Vinp], RTOL, ATOL, "steady-state: other nodes")
      end if
    end do

    !! ----------------------- small signal ----------------------------------------------------------------------------
    call sm%init(crt%sys)
    call sm%init_output([new_string("V1.I1"), new_string("V1.V")], "sm_admittance.fbs")
    s = (0.0, 1.0) * linspace(0.0, 2*wres, 5)
    call sm%run(s)

    ! expected admittance parameter of RLC series circuit
    Y_exp = s / (Lval * (s**2 + s * Rval / Lval + 1 / (Lval * Cval)))
    call tc%assert_eq(sm%get_scalar("V1.I1", "V1.V"), Y_exp, RTOL, ATOL, "small-signal: admittance parameters")
    call tc%assert_eq(sm%get_scalar("L.I1",  "V1.V"), Y_exp, RTOL, ATOL, "small-signal: admittance parameters")
    call tc%assert_eq(sm%get_scalar("C.I1",  "V1.V"), Y_exp, RTOL, ATOL, "small-signal: admittance parameters")
    ! expected small-signal parameter of input voltage: 1.0
    Y_exp = (1.0, 0.0)
    call tc%assert_eq(sm%get_scalar("V1.V",  "V1.V"), Y_exp, RTOL, ATOL, "small-signal: admittance parameters")
    ! expected small-signal parameter of GND node: 0.0
    Y_exp = (0.0, 0.0)
    call tc%assert_eq(sm%get_scalar("N0.V",  "V1.V"), Y_exp, RTOL, ATOL, "small-signal: admittance parameters")

    !! ----------------------- eigenvalues -----------------------------------------------------------------------------
    call ev%run_dense(crt%sys)
    eig_exp = [-alpha + sqrt(cmplx(alpha**2-wres**2)), -alpha - sqrt(cmplx(alpha**2-wres**2))]
    call tc%assert_eq(ev%s, eig_exp, RTOL, ATOL, "eigenvalues: poles of admittance")

    !! ----------------------- transient: voltage step from 0.0 to Vinp at t=0 -----------------------------------------
    tmax = 3 / (wres/2.0/PI) ! simulate three resonance frequency periods
    Nt = 101 ! number of time points
    call ss%init(crt%sys)
    call ss%set_newton_params(max_it = 2)

    ! expected result
    allocate (t_exp(0 : Nt-1))
    allocate (I_exp(0 : Nt-1))
    t_exp = linspace(0.0, tmax, Nt)
    if (xi <= 1-1e-15) then
      I_exp = Vinp / Lval / sqrt(wres**2 - alpha**2) * exp(-alpha*t_exp) * sin(sqrt(wres**2 - alpha**2)*t_exp)
    else if (abs(xi-1) < 1e-15) then
      I_exp = Vinp / Lval * t_exp * exp(-alpha*t_exp)
    else
      I_exp = Vinp / Lval / (2*wres*sqrt(xi**2-1)) * (exp(-wres*(xi-sqrt(xi**2-1))*t_exp) - exp(-wres*(xi+sqrt(xi**2-1))*t_exp))
    end if
    ! trapezoidal rule
    call src%init([0.0])
    call ss%run(input = src)
    call src%init([Vinp])
    call crt%nodes%d(2)%p%volt%set([Vinp])
    call crt%nodes%d(3)%p%volt%set([Vinp])
    call tr%init(crt%sys)
    call tr%set_ode_params(method = TRANS_TRAPZ)
    call tr%set_newton_params(max_it = 2)
    call tr%set_var_params("GND.I1", atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("V1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("V1.I1",  atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("L.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("C.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("N0.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N2.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N3.V",   atol = Vinp      * 1e-10)
    call tr%init_output([new_string("V1.I1")], "trapz.fbs")
    call tr%run([0.0, tmax], n_steps = [Nt-1], input = src)
    call tc%assert_eq(tr%t, t_exp, RTOL, ATOL, "trapezoidal rule: time points")
    i = crt%sys%search_main_var("V1.I1")
    call tc%assert_eq(tr%x(i,:), I_exp, 1e-2, maxval(abs(I_exp))*1e-2, "trapezoidal rule: transient current")

    ! tr-bdf2: fixed time step
    call src%init([0.0])
    call ss%run(input = src)
    call src%init([Vinp])
    call crt%nodes%d(2)%p%volt%set([Vinp])
    call crt%nodes%d(3)%p%volt%set([Vinp])
    call tr%init(crt%sys)
    call tr%set_ode_params(method = TRANS_TRBDF2)
    call tr%set_newton_params(max_it = 2)
    call tr%set_var_params("GND.I1", atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("V1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("V1.I1",  atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("L.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("C.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("N0.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N2.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N3.V",   atol = Vinp      * 1e-10)
    call tr%init_output([new_string("V1.I1")], "trbdf2_fixed.fbs")
    call tr%init_cache("trbdf2_cache.fbs")
    call tr%run([0.0, tmax], n_steps = [Nt-1], input = src)
    call tc%assert_eq(tr%t, t_exp, RTOL, ATOL, "tr-bdf2: time points")
    i = crt%sys%search_main_var("V1.I1")
    call tc%assert_eq(tr%x(i,:), I_exp, 1e-2, maxval(abs(I_exp))*1e-2, "tr-bdf2: transient current")

    ! tr-bdf2: start from result in cachefile
    call tr%init(crt%sys)
    call tr%set_ode_params(method = TRANS_TRBDF2)
    call tr%set_newton_params(max_it = 2)
    call tr%set_var_params("GND.I1", atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("V1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("V1.I1",  atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("L.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("C.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("N0.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N1.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N2.V",   atol = Vinp      * 1e-10)
    call tr%set_var_params("N3.V",   atol = Vinp      * 1e-10)
    call tr%init_output([new_string("V1.I1")], "trbdf2_readin.fbs")
    ! start a simulation at 1/4 of the previous simulation
    call tr%run([t_exp(Nt/4), tmax], n_steps = [Nt-1-Nt/4], input = src, start_file = "trbdf2_cache.fbs", start_i = Nt/4)
    call tc%assert_eq(tr%t, t_exp(Nt/4:Nt-1), RTOL, ATOL, "tr-bdf2 readin: time points")
    i = crt%sys%search_main_var("V1.I1")
    call tc%assert_eq(tr%x(i,:), I_exp(Nt/4:Nt-1), 1e-2, maxval(abs(I_exp))*1e-2, "tr-bdf2 readin: transient current")

    ! tr-bdf2: adaptive time step
    call src%init([0.0])
    call ss%run(input = src)
    call src%init([Vinp])
    call crt%nodes%d(2)%p%volt%set([Vinp])
    call crt%nodes%d(3)%p%volt%set([Vinp])
    call tr%init(crt%sys)
    call tr%set_ode_params(method = TRANS_TRBDF2, adaptive = .true., eabs = maxval(I_exp) * 1e-4)
    call tr%set_newton_params(max_it = 2)
    call tr%set_var_params("GND.I1", atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("V1.V",   atol = Vinp      * 1e-10, eabs = Vinp * 1e-4)
    call tr%set_var_params("V1.I1",  atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("L.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("C.I1",   atol = Vinp/Rval * 1e-10)
    call tr%set_var_params("N0.V",   atol = Vinp      * 1e-10, eabs = Vinp * 1e-4)
    call tr%set_var_params("N1.V",   atol = Vinp      * 1e-10, eabs = Vinp * 1e-4)
    call tr%set_var_params("N2.V",   atol = Vinp      * 1e-10, eabs = Vinp * 1e-4)
    call tr%set_var_params("N3.V",   atol = Vinp      * 1e-10, eabs = Vinp * 1e-4)
    call tr%init_output([new_string("V1.I1")], "trbdf2_adapt.fbs")
    call tr%run([0.0, tmax], dt0 = tmax / 300, input = src)
    if (xi <= 1-1e-15) then
      I_exp = Vinp / Lval / sqrt(wres**2 - alpha**2) * exp(-alpha*tr%t) * sin(sqrt(wres**2 - alpha**2)*tr%t)
    else if (abs(xi-1) < 1e-15) then
      I_exp = Vinp / Lval * tr%t * exp(-alpha*tr%t)
    else
      I_exp = Vinp / Lval / (2*wres*sqrt(xi**2-1)) * (exp(-wres*(xi-sqrt(xi**2-1))*tr%t) - exp(-wres*(xi+sqrt(xi**2-1))*tr%t))
    end if
    i = crt%sys%search_main_var("V1.I1")
    call tc%assert_eq(tr%x(i,:), I_exp, 1e-2, maxval(abs(I_exp))*1e-2, "adaptive tr-bdf2: transient current")

    call tc%finish()

    ! destruct circuit
    call GND%destruct()
    call V1%destruct()
    call R%destruct()
    call L%destruct()
    call C%destruct()
    call crt%destruct()
  end subroutine

end module