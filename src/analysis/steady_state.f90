#include "../util/macro.f90.inc"

module steady_state_m

  use error_m,   only: assert_failed, program_error
  use esystem_m, only: esystem
  use gmres_m,   only: gmres_options
  use matrix_m,  only: matrix_real, sparse_real
  use newton_m,  only: newton, newton_opt

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

  recursive subroutine steady_state_run(this, sys, nopt, gopt, input, gum)
    !! perform steady-state analysis for one or multiple sets of input parameters
    class(steady_state),           intent(out)   :: this
    type(esystem),       target,   intent(inout) :: sys
      !! equation system
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)
    real,                optional, intent(in)    :: input(:,:)
      !! optional input parameters, dimension: sys%ninput, number of input sets
    procedure(gummel),   optional                :: gum
      !! optional procedure to run before newton iteration

    integer                   :: n_inputs, i
    real                      :: p(0)
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

    ! check input parameters
    if (present(input)) then
      ASSERT(size(input, 1) == sys%ninput)
      n_inputs = size(input, 2)
    else
      ASSERT(sys%ninput == 0)
      n_inputs = 1
    end if

    ! allocate memory to store results
    allocate (this%x(sys%n, n_inputs))

    ! solve steady-state for each set of input parameters
    do i = 1, n_inputs
      if (present(input)) call sys%set_input(input(:,i))
      if (present(gum))   call gum()
      call newton(fun, p, nopt_, sys%get_x(), this%x(:,i), gmres_opt = gopt)

      ! save results in esystem variables
      call sys%set_x(this%x(:,i))
    end do

    ! release factorizations
    if (nopt_%it_solver) then
      call dfp%destruct()
    else
      call df%destruct()
    end if

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
      ASSERT(size(p) == 0)
      IGNORE(p)
      ASSERT(present(dfdx))
      ASSERT(.not. present(dfdp))
      IGNORE(dfdp)

      ! save input variable
      call sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (nopt_%it_solver) then
        ASSERT(present(dfdx_prec))
        call sys%eval(f = f, df = df, dfp = dfp)
        dfdx_prec => dfp
      else
        ASSERT(.not. present(dfdx_prec))
        call sys%eval(f = f, df = df)
      end if
      dfdx => df
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
