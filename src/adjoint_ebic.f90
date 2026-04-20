module adjoint_ebic_m
  !! Adjoint-method EBIC kernel.
  !!
  !! Solves A^T phi = c once per bias point, then evaluates
  !! I_EBIC(r_beam) = phi^T s(r_beam) per beam position. Reuses existing
  !! Ramo-Shockley current DOFs (c as unit indicator) and the continuity
  !! beam-generation coupling (s via sys_full residual at u* with G switched on).
  !!
  !! Usage:
  !!   call ss%run(...)                          ! forward DC solve at dark state u*
  !!   call refactorize_at_converged(dev, ss%solver)
  !!   call build_collection_functional(dev, k, c)
  !!   call solve_adjoint(ss%solver, c, phi)
  !!   do i_beam = 1, n_beam
  !!     call build_adjoint_source(dev, r_beam(i_beam), s)
  !!     I_EBIC(i_beam) = evaluate_ebic(phi, s)
  !!   end do

  use block_m,       only: block_real
  use device_m,      only: device
  use error_m,       only: program_error
  use solver_base_m, only: solver_real

  implicit none

  private
  public build_collection_functional
  public build_adjoint_source
  public refactorize_at_converged
  public solve_adjoint
  public evaluate_ebic
  public extract_phi_for_var

contains

  subroutine build_collection_functional(dev, contact_id, c)
    !! c_k = unit indicator vector at the sys_full DOF index of dev%curr(k)%x.
    !!
    !! Uses the Ramo-Shockley terminal-current DOF already present in sys_full
    !! (the "currents" vselector main variable, one block with nct values in
    !! contact-id order). No hand-assembly of dI/du needed.
    class(device), target, intent(inout) :: dev
    integer,               intent(in)    :: contact_id
    real, allocatable,     intent(out)   :: c(:)

    integer :: iimvar, ibl, dof

    if (contact_id < 1 .or. contact_id > dev%par%nct) then
      call program_error("build_collection_functional: contact_id out of range")
    end if

    iimvar = dev%sys_full%search_main_var("currents")
    ibl    = dev%sys_full%res2block(iimvar)%d(1)

    if (dev%sys_full%i1(ibl) - dev%sys_full%i0(ibl) + 1 /= dev%par%nct) then
      call program_error("build_collection_functional: currents block size mismatch")
    end if

    allocate (c(dev%sys_full%n), source = 0.0)
    dof = dev%sys_full%i0(ibl) + contact_id - 1
    c(dof) = 1.0
  end subroutine

  subroutine refactorize_at_converged(dev, solver)
    !! Re-evaluate Jacobian at the current state and re-factorize the solver.
    !!
    !! The Newton loop in steady_state stores the LU at x_{final-1} (the state
    !! before the converged update). For strict adjoint correctness we want
    !! LU(J(u*)), so one extra eval + factorize pins the LU to the converged
    !! state. Negligible cost vs. the full Newton loop.
    class(device),       target, intent(inout) :: dev
    class(solver_real),          intent(inout) :: solver

    type(block_real), pointer :: dfdx, dfdx_precon

    call dev%sys_full%eval(df = dfdx, dfp = dfdx_precon)
    call solver%factorize(dfdx, P = dfdx_precon)
  end subroutine

  subroutine build_adjoint_source(dev, r_beam, s)
    !! Linearized source s = V * G(r; r_beam), scattered into the continuity
    !! rows of sys_full, with zero elsewhere.
    !!
    !! Implementation: set beam_pos, then call sys_full%eval. At the converged
    !! dark state u*, F(u*; G=0) = 0; the generation term is the only driver
    !! of non-zero F since continuity has dF/dG = -V. So F(u*; G) = -V*G on
    !! continuity rows, ~0 elsewhere (within Newton residual tolerance).
    !! Then s = -F.
    !!
    !! Relies on u* being the current esystem state (i.e. called right after
    !! Newton convergence without other set_x calls in between).
    class(device), target, intent(inout) :: dev
    real,                  intent(in)    :: r_beam
    real, allocatable,     intent(out)   :: s(:)

    real, allocatable :: f(:)

    allocate (s(dev%sys_full%n))
    allocate (f(dev%sys_full%n))

    dev%beam_pos%x = r_beam
    call dev%sys_full%eval(f = f)
    s = -f
  end subroutine

  subroutine solve_adjoint(solver, c, phi)
    !! Transpose-solve A^T phi = c using the existing factorization.
    !! PARDISO/MUMPS both support this via a flag; no re-factorization.
    class(solver_real), intent(inout) :: solver
    real,               intent(in)    :: c(:)
    real, allocatable,  intent(out)   :: phi(:)

    allocate (phi(size(c)))
    call solver%solve(c, phi, trans = 'T')
  end subroutine

  function evaluate_ebic(phi, s) result(I_EBIC)
    !! I_EBIC = phi * s (scalar dot product in sys_full flat indexing).
    !! phi's Poisson rows are naturally skipped because s vanishes there.
    real, intent(in) :: phi(:), s(:)
    real             :: I_EBIC

    if (size(phi) /= size(s)) then
      call program_error("evaluate_ebic: phi and s size mismatch")
    end if
    I_EBIC = dot_product(phi, s)
  end function

  subroutine extract_phi_for_var(dev, phi, var_name, tab_idx, phi_slice)
    !! Extract the portion of phi corresponding to main variable var_name,
    !! tab tab_idx. For "ndens"/"pdens" the tabs are:
    !!   tab 1       = interior vertices (transport_vct(0))
    !!   tab 2..nct+1 = contact vertices (transport_vct(1..nct))
    !! The interior slice is what's needed for the collection-probability map.
    class(device), target, intent(inout) :: dev
    real,                  intent(in)    :: phi(:)
    character(*),          intent(in)    :: var_name
    integer,               intent(in)    :: tab_idx
    real, allocatable,     intent(out)   :: phi_slice(:)

    integer :: iimvar, ibl

    iimvar = dev%sys_full%search_main_var(var_name)
    if (tab_idx < 1 .or. tab_idx > size(dev%sys_full%res2block(iimvar)%d)) then
      call program_error("extract_phi_for_var: tab_idx out of range")
    end if

    ibl      = dev%sys_full%res2block(iimvar)%d(tab_idx)
    phi_slice = phi(dev%sys_full%i0(ibl):dev%sys_full%i1(ibl))
  end subroutine

end module
