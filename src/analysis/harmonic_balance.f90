#include "../util/macro.f90.inc"

module harmonic_balance_m

  use error_m,          only: assert_failed, program_error
  use esystem_m,        only: esystem
  use gmres_m,          only: gmres_options
  use high_precision_m, only: hp_real, real_to_hp, hp_to_real, operator (+)
  use input_src_m,      only: periodic_src
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

    real, allocatable :: x(:,:)
      !! coefficients (sys%n*(1+2*nH), nfreq)
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
    class(periodic_src),           intent(inout) :: input
      !! periodic input source (input%freq may get changed)
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    integer,             optional, intent(in)    :: nt
      !! number of time steps for numerical integration (default: 128)

    integer                        :: i, i0, i1, j, j0, j1, k, l, nfreq, nt_
    real                           :: dum(0)
    real,              allocatable :: t(:), x0(:), ff(:)
    type(hp_real),     allocatable :: xx(:), fc(:,:), fs(:,:)
    type(newton_opt)               :: nopt_
    type(sparse_real)              :: df, dft
    type(sparse_real), target      :: dg
    type(sparse_real), allocatable :: dfc(:), dfs(:)
    type(sparse_real), pointer     :: ptr
    type(block_real)               :: dg_blk

    ! error checking
    ASSERT(nH > 0)
    ASSERT(input%n == sys%ninput)

    ! parse optional newton options
    if (present(nopt)) then
      ASSERT(size(nopt%atol) == sys%n*(1+2*nH))
      ASSERT(.not. nopt%it_solver)
      nopt_ = nopt
    else
      call nopt_%init(sys%n*(1+2*nH))
    end if

    ! time points (normalized to 1/f)
    nt_ = 128
    if (present(nt)) nt_ = nt
    ASSERT(nt_ > 0)
    t = linspace(0.5/nt_, 1.0 - 0.5/nt_, nt_)

    ! normalize input source frequency
    input%freq = 1.0

    ! set members
    this%sys  => sys
    this%nH   =  nH
    this%freq =  freq
    nfreq     =  size(freq)
    allocate (this%x(sys%n*(1+2*nH),nfreq), source = 0.0)

    ! allocate temporary memory
    allocate (x0(sys%n*(1+2*nH)), ff(sys%n), source = 0.0)
    allocate (xx(sys%n), fc(sys%n,0:nH), fs(sys%n,1:nH))
    allocate (dfc(0:2*nH), dfs(0:2*nH))
    call dfc(0)%init(sys%n)
    do k = 1, 2*nH
      call dfc(k)%init(sys%n)
      call dfs(k)%init(sys%n)
    end do
    call dg_blk%init([(sys%n, i = 1, 1+2*nH)])

    ! get time derivative matrix
    call sys%get_dft(dft)

    ! load steady-state as starting point for newton iteration
    x0(1:sys%n) = sys%get_x()

    ! solve for all frequencies
    do i = 1, nfreq
      print *
      print "(A,ES24.16)", "harmonic balance: freq = ", denorm(freq(i), "Hz")

      ! solve system by newton iteration
      call newton(fun, dum, nopt_, x0, this%x(:,i))

      ! release factorization
      call dg%destruct()

      ! update newton starting point for next frequency (should reduce number of iterations)
      x0 = this%x(:,i)
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

      IGNORE(p)
      IGNORE(dfdx_prec)
      IGNORE(dfdp)

      ! evaluate fourier coefficients of f and df
      fc = real_to_hp(0.0)
      fs = real_to_hp(0.0)
      call dfc(0)%reset()
      do k = 1, 2*nH
        call dfc(k)%reset()
        call dfs(k)%reset()
      end do
      do j = 1, nt_
        ! get variables for time t(j)
        xx = real_to_hp(x(1:sys%n))
        j1 = sys%n
        do k = 1, nH
          i0 = j1 + 1
          i1 = j1 + sys%n
          j0 = i1 + 1
          j1 = i1 + sys%n
          xx = (xx + x(i0:i1) * cos(2*PI*k*t(j))) + x(j0:j1) * sin(2*PI*k*t(j))
        end do

        ! evaluate system
        call sys%set_x(hp_to_real(xx))
        call sys%set_input(input%get(t(j)))
        call sys%eval(f = ff, df = df)

        ! integrate residuals
        fc(:,0) = fc(:,0) + ff / nt_
        do k = 1, nH
          fc(:,k) = fc(:,k) + (2*cos(2*PI*k*t(j))/nt_) * ff
          fs(:,k) = fs(:,k) + (2*sin(2*PI*k*t(j))/nt_) * ff
        end do

        ! integrate jacobians
        call matrix_add(df, dfc(0), fact = 1.0/nt_)
        do k = 1, 2*nH
          call matrix_add(df, dfc(k), fact = 2*cos(2*PI*k*t(j))/nt_)
          call matrix_add(df, dfs(k), fact = 2*sin(2*PI*k*t(j))/nt_)
        end do
      end do

      ! get residuals of large system
      if (present(f)) then
        f(1:sys%n) = hp_to_real(fc(:,0))
        j1 = sys%n
        do k = 1, nH
          i0 = j1 + 1
          i1 = j1 + sys%n
          j0 = i1 + 1
          j1 = i1 + sys%n
          f(i0:i1) = hp_to_real(fc(:,k))
          f(j0:j1) = hp_to_real(fs(:,k))
          call dft%mul_vec( 2*PI*k*freq(i)*x(j0:j1), f(i0:i1), fact_y = 1.0)
          call dft%mul_vec(-2*PI*k*freq(i)*x(i0:i1), f(j0:j1), fact_y = 1.0)
        end do
      end if

      ! get jacobian of large system
      if (present(dfdx)) then
        call dg_blk%set(1, 1, dfc(0))
        do l = 1, nH
          call dg_blk%set(1, 2*l  , dfc(l), fact = 0.5)
          call dg_blk%set(1, 2*l+1, dfs(l), fact = 0.5)
        end do
        do k = 1, nH
          call dg_blk%set(2*k  , 1, dfc(k))
          call dg_blk%set(2*k+1, 1, dfs(k))
          do l = 1, nH
            call dg_blk%set(2*k, 2*l, dfc(k+l), fact = 0.5)
            call dg_blk%get(2*k, 2*l, ptr)
            if (k == l) then
              call matrix_add(dfc(0), ptr)
            else
              call matrix_add(dfc(abs(k-l)), ptr, fact = 0.5)
            end if

            call dg_blk%set(2*k, 2*l+1, dfs(k+l), fact = 0.5)
            call dg_blk%get(2*k, 2*l+1, ptr)
            if (k == l) then
              call matrix_add(dft, ptr, fact = 2*PI*k*freq(i))
            else
              call matrix_add(dfs(abs(k-l)), ptr, fact = -0.5*sign(1.0, real(k-l)))
            end if

            call dg_blk%set(2*k+1, 2*l, dfs(k+l), fact = 0.5)
            call dg_blk%get(2*k+1, 2*l, ptr)
            if (k == l) then
              call matrix_add(dft, ptr, fact = -2*PI*k*freq(i))
            else
              call matrix_add(dfs(abs(k-l)), ptr, fact = 0.5*sign(1.0, real(k-l)))
            end if

            call dg_blk%set(2*k+1, 2*l+1, dfc(k+l), fact = -0.5)
            call dg_blk%get(2*k+1, 2*l+1, ptr)
            if (k == l) then
              call matrix_add(dfc(0), ptr)
            else
              call matrix_add(dfc(abs(k-l)), ptr, fact = 0.5)
            end if
          end do
        end do
        call matrix_convert(dg_blk, dg)
        dfdx => dg
      end if
    end subroutine

  end subroutine

  subroutine harmonic_balance_select_harmonic(this, k, ifreq, cosine, sine)
    !! select and save specific harmonic (overwrite equation system variables)
    class(harmonic_balance), intent(inout) :: this
    integer,                 intent(in)    :: k
      !! harmonic number (k = 0..nH; 0 not allowed for sine)
    integer, optional,       intent(in)    :: ifreq
      !! frequency index (default: 1)
    logical, optional,       intent(in)    :: cosine
      !! select cosine harmonic (default: false)
    logical, optional,       intent(in)    :: sine
      !! select sine harmonic (default: false)

    integer :: ifreq_, n, i0, i1
    logical :: cosine_, sine_

    ! optional arguments
    ifreq_ = 1
    if (present(ifreq)) ifreq_ = ifreq
    cosine_ = .false.
    if (present(cosine)) cosine_ = cosine
    sine_ = .false.
    if (present(sine)) sine_ = sine

    ! error checks
    ASSERT((k >= 0) .and. (k <= this%nH))
    ASSERT((ifreq_ >= 1) .and. (ifreq_ <= size(this%freq)))
    ASSERT(.not. (sine_ .and. cosine_))
    ASSERT((k > 0) .or. (.not. sine_))
    ASSERT((k == 0) .or. cosine_ .or. sine_)

    ! select harmonic and save in esystem
    n = this%sys%n
    if (k == 0) then
      ! constant
      call this%sys%set_x(this%x(1:n,ifreq_))
    else
      if (cosine_) then
        i0 = (2 * k - 1) * n + 1
        i1 = (2 * k    ) * n
        call this%sys%set_x(this%x(i0:i1,ifreq_))
      else ! sine_
        i0 = (2 * k    ) * n + 1
        i1 = (2 * k + 1) * n
        call this%sys%set_x(this%x(i0:i1,ifreq_))
      endif
    end if
  end subroutine

  subroutine harmonic_balance_select_time(this, t, ifreq)
    !! calculate and save variables from harmonics (overwrite equation system variables)
    class(harmonic_balance), intent(inout) :: this
    real,                    intent(in)    :: t
      !! time
    integer, optional,       intent(in)    :: ifreq
      !! frequency index (default: 1)

    integer           :: ifreq_, n, i0, i1, j0, j1, k
    real              :: t_
    real, allocatable :: x(:)

    ! optional arguments
    ifreq_ = 1
    if (present(ifreq)) ifreq_ = ifreq

    ! error checks
    ASSERT((ifreq_ >= 1) .and. (ifreq_ <= size(this%freq)))

    ! normalize time to one period
    t_ = t * this%freq(ifreq_)

    ! get variables at time t
    n = this%sys%n
    x = this%x(1:n,ifreq_)
    j1 = n
    do k = 1, this%nH
      i0 = j1 + 1
      i1 = j1 + n
      j0 = i1 + 1
      j1 = i1 + n
      x = x + this%x(i0:i1,ifreq_) * cos(2*PI*k*t_) + this%x(j0:j1,ifreq_) * sin(2*PI*k*t_)
    end do

    ! save variables in esystem
    call this%sys%set_x(x)
  end subroutine

end module
