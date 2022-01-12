program example

  use esystem_m,              only: esystem
  use example_steady_state_m
  use example_contact_m,      only: contacts
  use input_src_m,            only: const_src, harmonic_src
  use normalization_m,        only: norm, denorm
  use transient_m,            only: backward_euler, trapezoidal, bdf2, mbdf2, tr_bdf2
  use small_signal_m,         only: small_signal
  use math_m,                 only: PI, logspace

  implicit none

  real                 :: c(2,0:1), s(2,1)
  real                 :: t_0, t_1, delta_t
  type(const_src)      :: input_cons
  type(harmonic_src)   :: input_harm

  type(backward_euler) :: be_full
  type(trapezoidal)    :: tr_full
  type(bdf2)           :: bdf2_full
  type(mbdf2)          :: mbdf2_full
  type(tr_bdf2)        :: tr_bdf2_full

  complex, allocatable :: omega(:)
  real,    allocatable :: x0(:), freq(:)
  integer              :: i, ifreq, nfreq, iounit
  type(small_signal)   :: ac

  call init_configuration("src/example/example.inp")
  call init_dd()
  call init_full()
  call init_nlpe()

  ! init constant source
  call input_cons%init([(contacts(i)%volt%x, i = 1, size(contacts))])

  ! init sine input source (voltages)
  c(:,0) = [(contacts(i)%volt%x, i = 1, size(contacts))]
  c(:,1) = 0
  s(:,1) = norm([0.0, 0.01], "V")
  call input_harm%init(norm(1.0, "THz"), c, s)

  ! solve steady-state
  call solve_nlpe()
  call solve_gummel()
  call solve_full_newton(input_cons)
  call output()
  x0 = sys_full%get_x()

  ! small-signal
  nfreq  = 401
  freq  = norm(logspace(1e9, 1e20, nfreq), "Hz")
  omega = (0.0, 1.0) * 2 * PI * freq
  call ac%run_analysis(sys_full, omega)
  open (newunit = iounit, file = "small_signal.csv", action = "write")
  do ifreq = 1, nfreq
    call ac%select_real(2, ifreq)
    write (iounit, '(2ES24.16)', advance = "no") denorm(freq(ifreq), "Hz"), denorm(contacts(2)%curr%get(), "A")
    call ac%select_imag(2, ifreq)
    write (iounit, '(1ES24.16)') denorm(contacts(2)%curr%get(), "A")
  end do
  close (iounit)

  ! harmonic balance and transient
  t_0     = norm(0.0,   "s")
  t_1     = 6 / input_harm%freq
  delta_t = (t_1 - t_0) / 500

  print *, "Harmonic Balance"
  call sys_full%set_x(x0)
  call solve_harmonic_balance(input_harm, t_0, t_1, (t_1 - t_0)*1e-3)

  print *, "Backward-Euler"
  call sys_full%set_x(x0)
  call solve_transient_full_newton(be_full, input_harm, t_0, t_1, delta_t)

  print *, "Trapezoidal"
  call sys_full%set_x(x0)
  call solve_transient_full_newton(tr_full, input_harm, t_0, t_1, delta_t)

  print *, "BDF2"
  call sys_full%set_x(x0)
  call solve_transient_full_newton(bdf2_full, input_harm, t_0, t_1, delta_t)

  print *, "MBDF2"
  call sys_full%set_x(x0)
  call solve_transient_full_newton(mbdf2_full, input_harm, t_0, t_1, delta_t)

  print *, "TR-BDF2"
  call sys_full%set_x(x0)
  call solve_transient_full_newton(tr_bdf2_full, input_harm, t_0, t_1, delta_t)
end program
