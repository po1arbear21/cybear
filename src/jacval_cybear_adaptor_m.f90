module jacval_cybear_adaptor_m
  !! Phase 2 integration: wrap cybear's DD residual + analytical Jacobian as a
  !! jacobian_problem_cmplx so the snapshot dumper / validator suite can act on
  !! the real device residual.
  !!
  !! Status: SKELETON. Two prerequisites must be completed before the
  !! eval_cmplx procedure can be filled in:
  !!
  !!   1. Audit cybear's residual code (src/dd.f90, src/poisson.f90,
  !!      src/continuity.f90, src/schottky.f90, src/mobility.f90, ...) for
  !!      non-holomorphic operations:
  !!          abs, max, min          --> wrap to branch on real(.)
  !!          sqrt near zero         --> guard against imaginary branch cut
  !!          comparison branches    --> branch on real(.)
  !!          piecewise mobility     --> branch on real(.)
  !!          inner iterative solves --> verify the converged result is
  !!                                     holomorphic in the outer state
  !!      Enumerate each occurrence; create a small csd_intrinsics.F90 module
  !!      with the wrappers.
  !!
  !!   2. Add a CYBEAR_CSD preprocessor switch that retypes the residual's
  !!      real(wp) state arrays as complex(cwp). Compile the residual twice:
  !!      libcybear_real.a (default) and libcybear_csd.a (with -DCYBEAR_CSD).
  !!      The snapshot dumper links against both: eval_real uses the real
  !!      build, eval_cmplx uses the CSD build.
  !!
  !! Once those are in place the body of eval_cmplx below becomes:
  !!
  !!      call cybear_csd_residual(z, fz, this%device, this%bias_V, this%flags)
  !!
  !! and analogously eval_real / eval_jac call the real-arithmetic entry
  !! points. The cybear-side device/state machinery is *not* part of this
  !! adaptor; this module is the thin boundary between the abstract
  !! jacobian_problem_cmplx interface and whatever shape the cybear residual
  !! takes after the audit.

  use jacobian_validator_m, only: jacobian_problem_cmplx
  use jacval_snapshot_m,    only: physics_flags_t

  implicit none

  private
  public :: cybear_dd_problem

  type, extends(jacobian_problem_cmplx) :: cybear_dd_problem
    !! Wraps a cybear DD device residual evaluation at fixed bias.
    !!
    !! In a complete Phase-2 implementation this type owns (or borrows) the
    !! cybear device handle (mesh, materials, doping, boundary conditions)
    !! and routes eval_real / eval_jac / eval_cmplx through cybear's
    !! existing residual + assembly routines.
    real :: bias_V              = 0.0
    type(physics_flags_t) :: flags
    ! TODO Phase 2: device handle, snapshot of converged state, etc.
  contains
    procedure :: eval_real  => cybear_eval_real
    procedure :: eval_jac   => cybear_eval_jac
    procedure :: eval_cmplx => cybear_eval_cmplx
  end type

contains

  subroutine cybear_eval_real(this, x, f)
    class(cybear_dd_problem), intent(in)  :: this
    real,                     intent(in)  :: x(:)
    real,                     intent(out) :: f(:)
    ! TODO Phase 2: dispatch to cybear's real-arithmetic residual
    !   call cybear_residual_real(x, f, this%bias_V, this%flags, ...)
    f = 0.0
    error stop "cybear_eval_real: not yet implemented; see jacval.md Phase 2"
  end subroutine

  subroutine cybear_eval_jac(this, x, J)
    class(cybear_dd_problem), intent(in)  :: this
    real,                     intent(in)  :: x(:)
    real,                     intent(out) :: J(:,:)
    ! TODO Phase 2: dispatch to cybear's analytical-Jacobian assembly
    !   call cybear_jacobian_assemble(x, J, this%bias_V, this%flags, ...)
    J = 0.0
    error stop "cybear_eval_jac: not yet implemented; see jacval.md Phase 2"
  end subroutine

  subroutine cybear_eval_cmplx(this, z, fz)
    class(cybear_dd_problem), intent(in)  :: this
    complex,                  intent(in)  :: z(:)
    complex,                  intent(out) :: fz(:)
    ! TODO Phase 2: dispatch to cybear's complex-arithmetic residual
    !   call cybear_residual_cmplx(z, fz, this%bias_V, this%flags, ...)
    fz = (0.0, 0.0)
    error stop "cybear_eval_cmplx: not yet implemented; see jacval.md Phase 2"
  end subroutine

end module
