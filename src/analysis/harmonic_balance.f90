m4_include(../util/macro.f90.inc)

module harmonic_balance_m

  use error_m,          only: assert_failed, program_error
  use esystem_m,        only: esystem
  use gmres_m,          only: gmres_options
  use high_precision_m, only: hp_real, real_to_hp, hp_to_real, operator (+)
  use input_src_m,      only: periodic_src
  use blas95,           only: gemv, ger
  use math_m,           only: PI, linspace
  use matrix_m,         only: matrix_real, matrix_add, matrix_convert, sparse_real, block_real
  use normalization_m,  only: denorm
  use newton_m,         only: newton, newton_opt

  implicit none

  type harmonic_balance
    type(esystem), pointer :: sys => null()
      !! pointer to corresponding equation system

    integer :: nH
      !! number of harmonics

    real, allocatable :: freq(:)
      !! frequencies

    real, allocatable :: x(:,:,:)
      !! coefficients (1:sys%n, 0:2*nH, 1:nfreq)
  contains
    procedure :: run             => harmonic_balance_run
    procedure :: select_harmonic => harmonic_balance_select_harmonic
    procedure :: select_time     => harmonic_balance_select_time
  end type

contains

  subroutine harmonic_balance_run(this, sys, nH, freq, input, nopt, nt)
    !! performa harmonic balance analysis
    class(harmonic_balance),       intent(out)   :: this
    type(esystem),       target,   intent(inout) :: sys
      !! equation system
    integer,                       intent(in)    :: nH
      !! number of harmonics
    real,                          intent(in)    :: freq(:)
      !! frequencies to analyze
    class(periodic_src),           intent(in)    :: input
      !! periodic input source (input%freq is ignored)
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver (size must be sys%n)
    integer,             optional, intent(in)    :: nt
      !! number of time steps for numerical integration (default: 0 => adaptive integration)

    integer                          :: ifreq, nfreq, it, it0, it1, nt0, nt1, k, l
    real                             :: dum(0)
    real, allocatable                :: t(:), w(:,:), sol(:), sol0(:), xi(:), fi(:), ff(:,:), ff_old(:,:)
    type(newton_opt)                 :: nopt_
    class(periodic_src), allocatable :: input_

    type(sparse_real), target              :: dft, dfi, dg
    type(sparse_real), target, allocatable :: dff(:)
    type(sparse_real), pointer             :: ptr, df0, dfc(:), dfs(:)
    type(block_real)                       :: dg_block

    ! check for errors
    m4_assert(nH > 0)
    m4_assert(input%n == sys%ninput)

    ! parse optional newton options
    if (present(nopt)) then
      m4_assert(size(nopt%atol) == sys%n*(1+2*nH))
      m4_assert(.not. nopt%it_solver)
      nopt_ = nopt
    else
      call nopt_%init(sys%n*(1+2*nH))
    end if

    ! number of time points (normalized to 1/f)
    nt0 = 8
    nt1 = 128
    if (present(nt)) then
      m4_assert(nt > 0)
      nt0 = nt
      nt1 = nt
    end if

    ! normalize input source frequency
    input_ = input
    input_%freq = 1.0

    ! set members
    this%sys  => sys
    this%nH   =  nH
    this%freq =  freq
    nfreq     =  size(freq)
    allocate (this%x(1:sys%n,0:2*nH,1:nfreq), source = 0.0)

    ! generate time points suitable for recursive refinement
    allocate (t(nt1), source = 0.0)
    it1 = 1
    do while (it1 <= nt1/2)
      it0 = it1 + 1
      it1 = 2 * it1
      t(it0:it1) = linspace(1.0/it1, 1.0 - 1.0/it1, it1/2)
    end do

    ! precompute fourier weights (jacobians need weights up to two times the highest frequency)
    allocate (w(0:4*nH,nt1), source = 0.0)
    w(0,:) = 1
    do k = 1, 2*nH
      w(2*k-1,:) = 2 * cos(2*PI*k*t)
      w(2*k  ,:) = 2 * sin(2*PI*k*t)
    end do

    ! allocate temporary memory
    allocate (sol(sys%n*(1+2*nH)), sol0(sys%n*(1+2*nH)), xi(sys%n), fi(sys%n), ff(sys%n,0:2*nH), ff_old(sys%n,0:2*nH), source = 0.0)
    allocate (dff(0:4*nH))
    df0 => dff(0)
    dfc => dff(1:4*nH-1:2)
    dfs => dff(2:4*nH  :2)
    do k = 0, 4*nH
      call dff(k)%init(sys%n)
    end do
    call dg_block%init([(sys%n, k = 0, 2*nH)])

    ! get time derivative matrix
    call sys%get_dft(dft)

    ! load steady-state as starting point for newton iteration
    sol0(1:sys%n) = sys%get_x()

    ! solve for all frequencies
    do ifreq = 1, nfreq
      print *
      print "(A,ES24.16)", "harmonic balance: freq = ", denorm(freq(ifreq), "Hz")

      ! solve system by newton iteration
      call newton(fun, dum, nopt_, sol0, sol)

      ! release factorization
      call dg%destruct()

      ! save solution and update newton starting point for next frequency
      this%x(:,:,ifreq) = reshape(sol, [sys%n, 1+2*nH])
      sol0 = sol
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

      integer :: i0, i1, j0, j1

      m4_ignore(p)
      m4_ignore(dfdx_prec)
      m4_ignore(dfdp)

      ! set fourier coefficients for f and df
      call eval_fourier(reshape(x, [sys%n, 1+2*nH]))

      ! get residuals of large system
      if (present(f)) then
        f = reshape(ff, [sys%n*(1+2*nH)])
        j1 = sys%n
        do k = 1, nH
          i0 = j1 + 1
          i1 = j1 + sys%n
          j0 = i1 + 1
          j1 = i1 + sys%n
          call dft%mul_vec( 2*PI*k*freq(ifreq)*x(j0:j1), f(i0:i1), fact_y = 1.0)
          call dft%mul_vec(-2*PI*k*freq(ifreq)*x(i0:i1), f(j0:j1), fact_y = 1.0)
        end do
      end if

      ! get jacobian of large system
      if (present(dfdx)) then
        call dg_block%set(1, 1, df0)
        do l = 1, nH
          call dg_block%set(1, 2*l  , dfc(l), fact = 0.5)
          call dg_block%set(1, 2*l+1, dfs(l), fact = 0.5)
        end do
        do k = 1, nH
          call dg_block%set(2*k  , 1, dfc(k))
          call dg_block%set(2*k+1, 1, dfs(k))
          do l = 1, nH
            call dg_block%set(2*k, 2*l, dfc(k+l), fact = 0.5)
            call dg_block%get(2*k, 2*l, ptr)
            if (k == l) then
              call matrix_add(df0, ptr)
            else
              call matrix_add(dfc(abs(k-l)), ptr, fact = 0.5)
            end if

            call dg_block%set(2*k, 2*l+1, dfs(k+l), fact = 0.5)
            call dg_block%get(2*k, 2*l+1, ptr)
            if (k == l) then
              call matrix_add(dft, ptr, fact = 2*PI*k*freq(ifreq))
            else
              call matrix_add(dfs(abs(k-l)), ptr, fact = -0.5*sign(1.0, real(k-l)))
            end if

            call dg_block%set(2*k+1, 2*l, dfs(k+l), fact = 0.5)
            call dg_block%get(2*k+1, 2*l, ptr)
            if (k == l) then
              call matrix_add(dft, ptr, fact = -2*PI*k*freq(ifreq))
            else
              call matrix_add(dfs(abs(k-l)), ptr, fact = 0.5*sign(1.0, real(k-l)))
            end if

            call dg_block%set(2*k+1, 2*l+1, dfc(k+l), fact = -0.5)
            call dg_block%get(2*k+1, 2*l+1, ptr)
            if (k == l) then
              call matrix_add(df0, ptr)
            else
              call matrix_add(dfc(abs(k-l)), ptr, fact = 0.5)
            end if
          end do
        end do
        call matrix_convert(dg_block, dg)
        dfdx => dg
      end if
    end subroutine

    subroutine eval_fourier(xx)
      real, intent(in) :: xx(:,0:)

      real, parameter :: ATOL = 1e-15, RTOL = 1e-12
      logical         :: cond

      ! reset
      ff = 0
      do k = 0, 4*nH
        call dff(k)%reset()
      end do

      it1 = 1
      do while (it1 <= nt1/2)
        if (it1 < nt0) then
          ! start integral estimate with nt0 points
          it0 = 1
          it1 = nt0
        else
          ! exit if accurate enough
          cond = .true.
          do k = 0, 2*nH
            cond = cond .and. (maxval(abs(ff_old(:,k) - ff(:,k)) / (abs(ff(:,k)) + ATOL / RTOL)) <= RTOL)
          end do
          if (cond) exit

          ! prepare for time-step halving
          ff_old = ff
          ff = 0.5 * ff
          do k = 0, 4*nH
            call dff(k)%scale(0.5)
          end do

          ! go to next set of timepoints
          it0 = it1 + 1
          it1 = 2 * it1
        end if

        do it = it0, it1
          ! get variables for time t(it)
          xi = xx(:,0)
          call gemv(xx(:,1:2*nH), w(1:2*nH,it), xi, alpha = 0.5, beta = 1.0)

          ! evaluate equation system
          call sys%set_x(xi)
          call sys%set_input(input_%get(t(it)))
          call sys%eval(f = fi, df = dfi)

          ! update residual fourier coefficients
          call ger(ff, fi, w(0:2*nH,it), alpha = 1.0 / it1)

          ! update jacobian fourier coefficients
          do k = 0, 4*nH
            call matrix_add(dfi, dff(k), fact = w(k,it) / it1)
          end do
        end do
      end do
    end subroutine

  end subroutine

  subroutine harmonic_balance_select_harmonic(this, k, ifreq, cosine, sine)
    !! select and save specific harmonic (overwrite equation system variables)
    class(harmonic_balance), intent(inout) :: this
    integer,                 intent(in)    :: k
      !! harmonic number (k = 0..nH; 0 not allowed for cosine/sine)
    integer, optional,       intent(in)    :: ifreq
      !! frequency index (default: 1)
    logical, optional,       intent(in)    :: cosine
      !! select cosine harmonic (default: false)
    logical, optional,       intent(in)    :: sine
      !! select sine harmonic (default: false)

    integer :: ifreq_, n, l
    logical :: cosine_, sine_

    ! optional arguments
    ifreq_ = 1
    if (present(ifreq)) ifreq_ = ifreq
    cosine_ = .false.
    if (present(cosine)) cosine_ = cosine
    sine_ = .false.
    if (present(sine)) sine_ = sine

    ! error checks
    m4_assert((k >= 0) .and. (k <= this%nH))
    m4_assert((ifreq_ >= 1) .and. (ifreq_ <= size(this%freq)))
    m4_assert(.not. (sine_ .and. cosine_))
    m4_assert((k > 0) .or. .not. (sine_ .or. cosine_))
    m4_assert((k == 0) .or. cosine_ .or. sine_)

    ! select harmonic and save in esystem
    n = this%sys%n
    l = 0
    if (cosine_) l = 2 * k - 1
    if (  sine_) l = 2 * k
    call this%sys%set_x(this%x(:,l,ifreq_))
  end subroutine

  subroutine harmonic_balance_select_time(this, t, ifreq)
    !! calculate and save variables from harmonics (overwrite equation system variables)
    class(harmonic_balance), intent(inout) :: this
    real,                    intent(in)    :: t
      !! time
    integer, optional,       intent(in)    :: ifreq
      !! frequency index (default: 1)

    integer           :: ifreq_, k
    real              :: t_
    real, allocatable :: x(:)

    ! optional arguments
    ifreq_ = 1
    if (present(ifreq)) ifreq_ = ifreq

    ! error checks
    m4_assert((ifreq_ >= 1) .and. (ifreq_ <= size(this%freq)))

    ! normalize time to one period
    t_ = t * this%freq(ifreq_)

    ! get variables at time t
    x = this%x(:,0,ifreq_)
    do k = 1, this%nH
      x = x + this%x(:,2*k-1,ifreq_) * cos(2*PI*k*t_) + this%x(:,2*k,ifreq_) * sin(2*PI*k*t_)
    end do

    ! save variables in esystem
    call this%sys%set_x(x)
  end subroutine

end module
