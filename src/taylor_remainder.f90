module taylor_remainder_m
  !! Taylor remainder test for Jacobian validation
  !!
  !! Verifies that the analytical Jacobian J(x) is correct by checking:
  !!   T0(eps) = ||F(x0 + eps*dx) - F(x0)||               ~ O(eps)   (always)
  !!   T1(eps) = ||F(x0 + eps*dx) - F(x0) - eps*J*dx||    ~ O(eps^2) (iff J correct)
  !!
  !! On a log-log plot: T0 has slope 1, T1 has slope 2 if Jacobian is correct.

  use esystem_m, only: esystem
  use block_m,   only: block_real

  implicit none

  private
  public :: taylor_test, taylor_test_core, taylor_test_per_block

  abstract interface
    subroutine eval_f_sub(x, f)
      !! evaluate residual F(x)
      real, intent(in)  :: x(:)
      real, intent(out) :: f(:)
    end subroutine
  end interface

contains

  subroutine taylor_test_core(n, x0, dx, f0, Jdx, eval_f, neps, name)
    !! Core Taylor remainder test — generic version with callback
    integer, intent(in) :: n
      !! system size
    real, intent(in) :: x0(n)
      !! linearization point
    real, intent(in) :: dx(n)
      !! perturbation direction
    real, intent(in) :: f0(n)
      !! residual at x0: F(x0)
    real, intent(in) :: Jdx(n)
      !! Jacobian-vector product: J(x0) * dx
    procedure(eval_f_sub) :: eval_f
      !! callback to evaluate F(x)
    integer, intent(in), optional :: neps
      !! number of eps decades (default: 10)
    character(*), intent(in), optional :: name
      !! custom header name (default: "Taylor Remainder Test")

    integer :: k, nk
    real, allocatable :: f1(:), diff(:)
    real :: eps, T0, T1, T0_prev, T1_prev, order0, order1

    nk = 10; if (present(neps)) nk = neps

    allocate(f1(n), diff(n))

    ! header
    print "(A)", ""
    if (present(name)) then
      print "(2X,A)", name
    else
      print "(A)", "  Taylor Remainder Test"
    end if
    print "(A)", "  ====================="
    print "(A,I0)", "  n = ", n
    print "(A)", ""
    print "(A)", "    eps          T0             T1           order0  order1"
    print "(A)", "    ----------   -----------    -----------  ------  ------"

    T0_prev = 0.0
    T1_prev = 0.0

    do k = 1, nk
      eps = 10.0**(-k)

      ! evaluate F(x0 + eps*dx)
      call eval_f(x0 + eps * dx, f1)

      ! T0 = ||F(x0+eps*dx) - F(x0)||
      diff = f1 - f0
      T0 = sqrt(dot_product(diff, diff))

      ! T1 = ||F(x0+eps*dx) - F(x0) - eps*J*dx||
      diff = f1 - f0 - eps * Jdx
      T1 = sqrt(dot_product(diff, diff))

      if (k > 1 .and. T0_prev > 0.0 .and. T1_prev > 0.0 &
        & .and. T0 > 0.0 .and. T1 > 0.0) then
        order0 = log10(T0_prev / T0)
        order1 = log10(T1_prev / T1)
        print "(4X, ES12.2, 2ES14.4, 2F8.2)", eps, T0, T1, order0, order1
      else
        print "(4X, ES12.2, 2ES14.4, A)", eps, T0, T1, "       -       -"
      end if

      T0_prev = T0
      T1_prev = T1
    end do

    print "(A)", ""

    deallocate(f1, diff)
  end subroutine


  subroutine taylor_test_per_block(sys, neps, seed)
    !! Run Taylor test for each variable block separately.
    !! Perturbs one block at a time and shows the full convergence table
    !! (same output as taylor_test) so the user can judge correctness.
    class(esystem), intent(inout) :: sys
      !! equation system to test
    integer, intent(in), optional :: neps
      !! number of eps decades (default: 10)
    integer, intent(in), optional :: seed
      !! RNG seed (default: 42)

    integer :: n, ibl, iseed, isz, bsize
    real, allocatable :: x0(:), f0(:), dx(:), Jdx(:)
    type(block_real), pointer :: df
    character(:), allocatable :: vname, tname
    character(80) :: header

    n = sys%n
    iseed = 42; if (present(seed)) iseed = seed

    allocate(x0(n), f0(n), dx(n), Jdx(n))

    ! save state and evaluate F(x0), J(x0)
    x0 = sys%get_x()
    call sys%eval(f = f0, df = df)

    print "(A)", ""
    print "(A)", "  Taylor Test Per Block"
    print "(A,I0,A,I0,A)", "  n = ", n, ",  nbl = ", sys%nbl, " blocks"

    do ibl = 1, sys%nbl
      ! restore state and re-evaluate J(x0) (eval overwrites internal df)
      call sys%set_x(x0)
      call sys%eval(f = f0, df = df)

      ! block size
      bsize = sys%i1(ibl) - sys%i0(ibl) + 1

      ! get variable name and table name
      associate (nd => sys%g%nodes%d(sys%g%imvar%d(sys%block2res(1,ibl)))%p)
        vname = nd%v%name
        tname = nd%v%tab(sys%block2res(2,ibl))%p%name
      end associate
      write(header, "(A,I0,A,A,A,A,A,I0,A)") &
        "Block ", ibl, ": ", trim(vname), " / ", trim(tname), &
        " (", bsize, " vars)"

      ! skip empty blocks
      if (bsize <= 0) then
        print "(A)", ""
        print "(2X,A,A)", trim(header), " — skipped (empty)"
        cycle
      end if

      ! perturbation: random in this block only, zero elsewhere
      dx = 0.0
      call random_seed(size = isz)
      block
        integer, allocatable :: seed_arr(:)
        allocate(seed_arr(isz))
        seed_arr = iseed + ibl
        call random_seed(put = seed_arr)
      end block
      call random_number(dx(sys%i0(ibl):sys%i1(ibl)))
      dx(sys%i0(ibl):sys%i1(ibl)) = dx(sys%i0(ibl):sys%i1(ibl)) - 0.5

      ! compute J(x0) * dx for this perturbation
      Jdx = 0.0
      call df%mul_vec(dx, Jdx)

      ! run full convergence table for this block
      call taylor_test_core(n, x0, dx, f0, Jdx, sys_eval_f, neps, name = trim(header))
    end do

    ! restore original state
    call sys%set_x(x0)
    deallocate(x0, f0, dx, Jdx)

  contains

    subroutine sys_eval_f(x, f)
      real, intent(in)  :: x(:)
      real, intent(out) :: f(:)
      call sys%set_x(x)
      call sys%eval(f = f)
    end subroutine

  end subroutine


  subroutine taylor_test(sys, neps, seed)
    !! Convenience wrapper: run Taylor test on an equation system
    class(esystem), intent(inout) :: sys
      !! equation system to test (must be initialized and have a valid state)
    integer, intent(in), optional :: neps
      !! number of eps decades to test (default: 10)
    integer, intent(in), optional :: seed
      !! RNG seed for reproducibility (default: 42)

    integer :: n, iseed, isz
    real, allocatable :: x0(:), f0(:), dx(:), Jdx(:)
    type(block_real), pointer :: df

    n = sys%n
    iseed = 42; if (present(seed)) iseed = seed

    allocate(x0(n), f0(n), dx(n), Jdx(n))

    ! save current state and evaluate F(x0), J(x0)
    x0 = sys%get_x()
    call sys%eval(f = f0, df = df)

    ! generate reproducible random perturbation
    call random_seed(size = isz)
    block
      integer, allocatable :: seed_arr(:)
      allocate(seed_arr(isz))
      seed_arr = iseed
      call random_seed(put = seed_arr)
    end block
    call random_number(dx)
    dx = dx - 0.5

    ! compute J(x0) * dx
    Jdx = 0.0
    call df%mul_vec(dx, Jdx)

    ! run core test
    call taylor_test_core(n, x0, dx, f0, Jdx, sys_eval_f, neps)

    ! restore original state
    call sys%set_x(x0)

    deallocate(x0, f0, dx, Jdx)

  contains

    subroutine sys_eval_f(x, f)
      real, intent(in)  :: x(:)
      real, intent(out) :: f(:)
      call sys%set_x(x)
      call sys%eval(f = f)
    end subroutine

  end subroutine

end module
