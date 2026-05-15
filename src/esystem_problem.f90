module esystem_problem_m
  !! Adapter that lets jacobian_validator_m operate on an esystem.
  !!
  !! Wraps a pointer to an esystem and exposes its residual + Jacobian in the
  !! flat (x, f, J) form the validators expect. Use it to check any cybear
  !! Jacobian (poisson, continuity, schottky_bc, dd) against numerical
  !! references: finite differences, Taylor slope, JVP+Richardson.
  !!
  !! Does NOT extend jacobian_problem_cmplx because the DD residual is not
  !! analytic in general (mobility max(), Fermi-Dirac, sign-conditioned
  !! tunneling). The complex-step validator therefore is not applicable.

  use jacobian_validator_m, only: jacobian_problem
  use esystem_m,            only: esystem
  use dense_m,              only: dense_real

  implicit none

  private
  public :: esystem_problem

  type, extends(jacobian_problem) :: esystem_problem
    type(esystem), pointer :: sys => null()
  contains
    procedure :: eval_real => esp_eval_real
    procedure :: eval_jac  => esp_eval_jac
  end type

contains

  subroutine esp_eval_real(this, x, f)
    class(esystem_problem), intent(in)  :: this
    real,                   intent(in)  :: x(:)
    real,                   intent(out) :: f(:)

    call this%sys%set_x(x)
    call this%sys%eval(f = f)
  end subroutine

  subroutine esp_eval_jac(this, x, J)
    class(esystem_problem), intent(in)  :: this
    real,                   intent(in)  :: x(:)
    real,                   intent(out) :: J(:,:)

    type(dense_real) :: D

    call this%sys%set_x(x)
    call this%sys%eval()
    call this%sys%get_df(D)
    J = D%d
  end subroutine

end module
