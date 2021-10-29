#include "../util/macro.f90.inc"

module steady_state_m

  use error_m,     only: assert_failed, program_error
  use esystem_m,   only: esystem
  use gmres_m,     only: gmres_options
  use input_src_m, only: input_src
  use matrix_m,    only: matrix_real, sparse_real
  use newton_m,    only: newton, newton_opt

  implicit none

  private
  public steady_state

  type steady_state
    !! steady-state analyzer

    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system
    real,      allocatable :: x(:,:)
      !! steady-state data
  contains
    procedure :: run    => steady_state_run
    procedure :: select => steady_state_select
  end type

  abstract interface
    subroutine gummel()
    end subroutine
  end interface

contains

  recursive subroutine steady_state_run(this, sys, nopt, gopt, input, t_input, gum)
    !! perform steady-state analysis for one or multiple sets of input parameters
    class(steady_state),           intent(out)   :: this
    type(esystem),       target,   intent(inout) :: sys
      !! equation system
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)
    class(input_src),    optional, intent(in)    :: input
      !! optional input source
    real,                optional, intent(in)    :: t_input(:)
      !! perform quasi-stationary simulations for these time-points (default: [0.0])
    procedure(gummel),   optional                :: gum
      !! optional procedure to run before newton iteration

    integer                   :: nt, i
    real                      :: p(0)
    real, allocatable         :: t(:)
    type(newton_opt)          :: nopt_
    type(sparse_real), target :: df, dfp

    ! save pointer to equation system for later use in select routine
    this%sys => sys

    ! parse optional newton options
    if (present(nopt)) then
      ASSERT(size(nopt%atol) == sys%n)
      nopt_ = nopt
    else
      call nopt_%init(sys%n)
    end if

    ! input
    if (present(input)) then
      ASSERT(input%n == sys%ninput)
      if (present(t_input)) then
        t = t_input
      else
        t = [0.0]
      end if
      nt = size(t)
    else
      ASSERT(sys%ninput == 0)
      ASSERT(.not. present(t_input))
      nt = 1
    end if

    ! allocate memory to store results
    allocate (this%x(sys%n, nt))

    ! solve steady-state for each input time
    do i = 1, nt
      if (present(input)) call sys%set_input(input%get(t(i)))
      if (present(gum))   call gum()
      call newton(fun, p, nopt_, sys%get_x(), this%x(:,i), gmres_opt = gopt)

      ! release factorization
      if (nopt_%it_solver) then
        call dfp%destruct()
      else
        call df%destruct()
      end if

      ! save results in esystem variables
      call sys%set_x(this%x(:,i))
    end do

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      IGNORE(p)
      IGNORE(dfdp)

      ! set variables
      call sys%set_x(x)

      ! evaluate system to compute residuals, jacobian, and possible preconditioner
      if (nopt_%it_solver) then
        call sys%eval(f = f, df = df, dfp = dfp)
      else
        call sys%eval(f = f, df = df)
      end if

      ! output jacobian pointers
      if (present(dfdx     )) dfdx      => df
      if (present(dfdx_prec)) dfdx_prec => dfp
    end subroutine
  end subroutine

  subroutine steady_state_select(this, i)
    !! save steady-state results (overwrite esystem variables)
    class(steady_state), intent(in) :: this
    integer,             intent(in) :: i
      !! input set index

    call this%sys%set_x(this%x(:,i))
  end subroutine

end module
