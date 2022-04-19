module test_analysis_m

  use test_case_m,    only: test_case
  use esystem_m,      only: esystem
  use grid_data_m,    only: grid_data0_real
  use jacobian_m,     only: jacobian, jacobian_ptr
  use math_m,         only: PI, linspace, logspace
  use res_equation_m, only: res_equation
  use small_signal_m, only: small_signal
  use stencil_m,      only: dirichlet_stencil
  use variable_m,     only: variable_real

  implicit none

  private
  public test_analysis

  type, extends(variable_real) :: phi
    !! pendulum angle
    real, pointer :: x => null()
  contains
    procedure :: init => phi_init
  end type

  type, extends(variable_real) :: omega
    !! pendulum angular velocity
    real, pointer :: x => null()
  contains
    procedure :: init => omega_init
  end type

  type, extends(res_equation) :: dt_phi
    !! d/dt phi - omega = 0
    type(omega), pointer :: om => null()
  contains
    procedure :: init => dt_phi_init
    procedure :: eval => dt_phi_eval
  end type

  type, extends(res_equation) :: dt_omega
    type(phi),   pointer :: ph  => null()
    type(omega), pointer :: om  => null()
    type(omega), pointer :: om0 => null()

    real :: alpha
      !! friction coefficient
    real :: Omega_res
      !! approx. resonant frequency
    real :: cpl
      !! coupling to rod

    type(jacobian), pointer :: df_ph  => null()
    type(jacobian), pointer :: df_om  => null()
    type(jacobian), pointer :: dft_om => null()
    type(jacobian), pointer :: df_om0 => null()
  contains
    procedure :: init => dt_omega_init
    procedure :: eval => dt_omega_eval
  end type

  type, extends(res_equation) :: dt_omega0
    type(omega), pointer :: om0    => null()
    type(omega), pointer :: om(:)  => null()
    type(omega), pointer :: torque => null()

    real, allocatable :: cpl(:)
      !! coupling to pendulums
    real              :: cpl_torque
      !! coupling to torque

    type(jacobian),     pointer     :: df_om0
    type(jacobian),     pointer     :: dft_om0
    type(jacobian_ptr), allocatable :: df_om(:)
    type(jacobian),     pointer     :: df_torque
  contains
    procedure :: init => dt_omega0_init
    procedure :: eval => dt_omega0_eval
  end type

contains

  subroutine phi_init(this, name)
    class(phi),   intent(out) :: this
    character(*), intent(in)  :: name

    type(grid_data0_real), pointer :: p

    call this%variable_init(name, "1")

    ! get pointer to data
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

  subroutine omega_init(this, name)
    class(omega), intent(out) :: this
    character(*), intent(in)  :: name

    type(grid_data0_real), pointer :: p

    call this%variable_init(name, "1/s")

    ! get pointer to data
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

  subroutine dt_phi_init(this, ph, om)
    class(dt_phi),       intent(out) :: this
    type(phi),           intent(in)  :: ph
    type(omega), target, intent(in)  :: om

    integer :: idx(0)
    type(jacobian), pointer :: jaco


    ! init base
    call this%equation_init("dt_"//ph%name)

    ! save pointers to variables
    this%om => om

    ! init residuals (main variable phi)
    call this%init_f(ph)

    ! depend on omega
    jaco => this%init_jaco_f(this%depend(om), const = .true.)
    call jaco%set(idx, idx, -1.0)

    ! depend on phi
    jaco => this%init_jaco_f(this%depend(ph), dtime = .true.)
    call jaco%set(idx, idx, 1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dt_phi_eval(this)
    class(dt_phi),  intent(inout) :: this

    call this%f%set([-this%om%x])
  end subroutine

  subroutine dt_omega_init(this, ph, om, om0, alpha, Omega_res, cpl)
    class(dt_omega),     intent(out) :: this
    type(phi),   target, intent(in)  :: ph
    type(omega), target, intent(in)  :: om
    type(omega), target, intent(in)  :: om0
    real,                intent(in)  :: alpha
    real,                intent(in)  :: Omega_res
    real,                intent(in)  :: cpl

    integer :: idep, idx1(0), idx2(0)

    ! init base
    call this%equation_init("dt_"//om%name)

    ! set members
    this%ph        => ph
    this%om        => om
    this%om0       => om0
    this%alpha     =  alpha
    this%Omega_res =  Omega_res
    this%cpl       =  cpl

    ! init residuals (main variable omega)
    call this%init_f(om)

    ! depend on phi
    idep = this%depend(ph)
    this%df_ph => this%init_jaco_f(idep)

    ! depend on omega
    idep = this%depend(om)
    this%df_om => this%init_jaco_f(idep, const = .true.)
    call this%df_om%set(idx1, idx2, alpha + cpl)

    ! time derivative
    this%dft_om => this%init_jaco_f(idep, dtime = .true.)
    call this%dft_om%set(idx1, idx2, 1.0)

    ! depend on omega0
    idep = this%depend(om0)
    this%df_om0 => this%init_jaco_f(idep, const = .true.)
    call this%df_om0%set(idx1, idx2, - cpl)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dt_omega_eval(this)
    class(dt_omega), intent(inout) :: this

    integer :: idx1(0), idx2(0)

    ! d/dt omega + w0^2 * sin(ph) + alpha * om == 0
    call this%f%set([this%Omega_res**2 * sin(this%ph%x) + this%alpha * this%om%x - this%cpl * (this%om0%x - this%om%x)])
    call this%df_ph%set(idx1, idx2, this%Omega_res**2 * cos(this%ph%x))
  end subroutine

  subroutine dt_omega0_init(this, om0, om, torque, cpl, cpl_torque)
    class(dt_omega0),    intent(out) :: this
    type(omega), target, intent(in)  :: om0
    type(omega), target, intent(in)  :: om(:)
    type(omega), target, intent(in)  :: torque
    real,                intent(in)  :: cpl(:)
    real,                intent(in)  :: cpl_torque

    integer :: i, idep, idx1(0), idx2(0)

    ! init base
    call this%equation_init("dt_"//om0%name)

    ! set members
    this%om0        => om0
    this%om         => om
    this%torque     => torque
    this%cpl        =  cpl
    this%cpl_torque = cpl_torque

    ! init residuals
    call this%init_f(om0)

    ! depend on omega0
    idep = this%depend(om0)
    this%df_om0 => this%init_jaco_f(idep, const = .true.)
    call this%df_om0%set(idx1, idx2, sum(cpl) + cpl_torque)

    ! time derivative
    this%dft_om0 => this%init_jaco_f(idep, dtime = .true.)
    call this%dft_om0%set(idx1, idx2, 1.0)

    ! depend on omega
    allocate (this%df_om(size(om)))
    do i = 1, size(om)
      idep = this%depend(om(i))
      this%df_om(i)%p => this%init_jaco_f(idep, const = .true.)
      call this%df_om(i)%p%set(idx1, idx2, - cpl(i))
    end do

    ! depend on torque
    idep = this%depend(torque)
    this%df_torque => this%init_jaco_f(idep, const = .true.)
    call this%df_torque%set(idx1, idx2, - cpl_torque)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dt_omega0_eval(this)
    class(dt_omega0), intent(inout) :: this

    integer :: i
    real    :: f(1)

    f = sum(this%cpl) * this%om0%x
    do i = 1, size(this%om)
      f = f - this%cpl(i) * this%om(i)%x
    end do
    f = f - this%cpl_torque * this%torque%x
    call this%f%set(f)
  end subroutine

  subroutine test_analysis()
    integer, parameter :: Np = 5
    integer            :: i
    character(16)      :: name
    real               :: om_res(Np), alpha(Np), cpl(Np), cpl_torque
    type(phi)          :: ph(1:Np)
    type(omega)        :: om(0:Np), torque
    type(dt_phi)       :: dt_ph(1:Np)
    type(dt_omega)     :: dt_om(1:Np)
    type(dt_omega0)    :: dt_om0
    type(esystem)      :: es
    type(test_case)    :: tc

    call tc%init("analysis")

    alpha      = 0.1
    cpl        = 5.0
    cpl_torque = 1.0
    om_res     = 2 * PI * linspace(1.0, 5.0, Np)

    call om(0)%init("omega_0")
    call torque%init("torque")
    do i = 1, Np
      write (name, "(A,I0)") "phi_", i
      call ph(i)%init(trim(name))
      write (name, "(A,I0)") "omega_", i
      call om(i)%init(trim(name))
      call dt_ph(i)%init(ph(i), om(i))
      call dt_om(i)%init(ph(i), om(i), om(0), alpha(i), om_res(i), cpl(i))
    end do
    call dt_om0%init(om(0), om(1:Np), torque, cpl, cpl_torque)

    ! forced oscillation
    call es%init("es")
    do i = 1, Np
      call es%add_equation(dt_ph(i))
      call es%add_equation(dt_om(i))
    end do
    call es%add_equation(dt_om0)
    call es%provide(torque, input=.true.)
    call es%init_final()

    ! test small-signal
    block
      integer, parameter   :: Ns = 501
      real,    allocatable :: f(:)
      complex, allocatable :: s(:)
      type(small_signal)   :: ac

      f = logspace(0.1, 10.0, Ns)
      s = (0.0, 1.0) * 2 * PI * f
      call ac%run_analysis(es, s)

      do i = 1, Ns
        call ac%select_abs(1, i)
        ! print "(2ES24.16)", f(i), om(0)%x
      end do
    end block

    call tc%finish()
  end subroutine

end module
