m4_include(util/macro.f90.inc)

module distribution_table_m

  use bin_search_m,    only: bin_search, BS_LESS
  use error_m,         only: assert_failed, program_error
  use ieee_arithmetic, only: ieee_value, IEEE_POSITIVE_INF
  use math_m,          only: linspace
  use omp_lib,         only: omp_get_thread_num, omp_get_num_threads
  use quad_m,          only: quad
  use qsort_m,         only: qsort
  use util_m,          only: hash
  use vector_m,        only: vector_int, vector_log, vector_real

  implicit none

  type distribution_table
    !! lookup table for cumulative distribution function

    character(:), allocatable :: name
      !! name for loading and saving
    integer                   :: kmax
      !! maximum derivative stored

    real,    allocatable :: eta(:)
      !! chemical potential
    real,    allocatable :: val(:,:)
      !! values and derivatives (0:kmax+1 x num_entries)
    logical, allocatable :: lg(:,:)
      !! logarithmic interpolation? per interval (0:kmax+1 x num_entries)
  contains
    procedure :: init   => distribution_table_init
    procedure :: gen    => distribution_table_gen
    procedure :: get    => distribution_table_get
    procedure :: inv    => distribution_table_inv
    procedure :: output => distribution_table_output
    procedure :: load   => distribution_table_load
    procedure :: save   => distribution_table_save
  end type

  interface
    function density_of_states(t) result(Z)
      !! get density of states
      real(kind=16), intent(in) :: t
        !! energy relative to band edge (in units of k_B T)
      real(kind=16)             :: Z
        !! return density of states
    end function

    function distribution_density(u, k) result(f)
      !! k-th derivative of distribution density (e.g. fermi-dirac or maxwell-boltzmann)
      real(kind=16),    intent(in) :: u
        !! energy relative to chemical potential/fermi level (in units of k_B T)
      integer,          intent(in) :: k
        !! k-th derivative (possible values from 0 to kmax)
      real(kind=16)                :: f
        !! return k-th distribution density
    end function
  end interface

  logical :: distribution_table_disable_save = .false.

contains

  subroutine distribution_table_init(this, name, dos, dist, eta_min, eta_max, kmax)
    !! initialize distribution table
    class(distribution_table), intent(out) :: this
    character(*),              intent(in)  :: name
      !! name for loading and saving
    procedure(density_of_states)           :: dos
    procedure(distribution_density)        :: dist
    real,                      intent(in)  :: eta_min
      !! minimal supported eta
    real,                      intent(in)  :: eta_max
      !! maximal supported eta
    integer,                   intent(in)  :: kmax
      !! maximal derivative

    integer, parameter :: N0 = 1024+1
    real,    parameter :: RTOL = 1e-13, ATOL = 1e-16

    integer              :: i, i1, i2, i3, j1, j2, j3, k, k1, k2, k3, ithread, nthreads, n
    integer, allocatable :: nth(:), perm(:)
    logical              :: status, lg1(0:kmax)
    logical, allocatable :: lg(:,:)
    real                 :: eta1, eta2, eta3, deta, val1(0:kmax), val2(0:kmax), err1(0:kmax), err2(0:kmax), sgn
    real,    allocatable :: eta(:), val(:,:)
    type(vector_int)     :: stack
    type(vector_log)     :: vlg
    type(vector_real)    :: veta, vval

    this%name = name
    this%kmax = kmax

    ! try to load from file
    call this%load("/tmp", eta_min, eta_max, kmax, status)
    if (status) return

    ! initial coarse grid
    allocate (eta(N0), val(0:kmax+1,N0), lg(0:kmax,N0))
    eta = linspace(eta_min, eta_max, N0)

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,j1,j2,j3,k,k1,k2,k3,ithread,lg1,eta1,eta2,eta3,deta,val1,val2,err1,err2,sgn,stack,vlg,veta,vval) &
    !$omp shared(this,kmax,nthreads,n,nth,eta,val,lg)

    ithread = omp_get_thread_num() + 1

    ! allocate number of table entries per thread
    !$omp single
    nthreads = omp_get_num_threads()
    allocate (nth(0:nthreads + 1), source = 0)
    !$omp end single

    ! generate entries for coarse grid
    !$omp do schedule(dynamic)
    do i = 1, N0
      call this%gen(dos, dist, eta(i), val(:,i))
      lg(:,i) = .false.
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    call vlg%init( N0 * (kmax + 1), c = 4 * N0 * (kmax + 1), x = reshape(lg,  [N0 * (kmax + 1)]))
    call veta%init(N0,              c = 4 * N0,              x = eta)
    call vval%init(N0 * (kmax + 2), c = 4 * N0 * (kmax + 2), x = reshape(val, [N0 * (kmax + 2)]))

    ! interval stack
    call stack%init(0, c = 16)

    ! refinement
    !$omp do schedule(dynamic) reduction(.or.:lg)
    do i = 1, N0 - 1
      ! add interval to stack
      call stack%push(i)
      call stack%push(i+1)

      ! refine intervals on stack
      do while (stack%n > 0)
        ! get interval from stack
        i1      = stack%d(stack%n - 1)
        i2      = stack%d(stack%n    )
        stack%n = stack%n - 2

        eta1 = veta%d(i1)
        eta2 = veta%d(i2)
        eta3 = 0.5 * (eta1 + eta2)
        deta = eta2 - eta1

        ! save midpoint
        call veta%push(eta3)
        i3 = veta%n

        ! data indices
        j1 = (i1 - 1) * (kmax + 2) + 1
        j2 = (i2 - 1) * (kmax + 2) + 1
        j3 = (i3 - 1) * (kmax + 2) + 1
        k1 = i1 * (kmax + 2)
        k2 = i2 * (kmax + 2)
        k3 = i3 * (kmax + 2)

        ! make room for new point
        call vval%resize(k3)

        ! generate data for new point
        call this%gen(dos, dist, eta3, vval%d(j3:k3))

        ! direct interpolation from existing data
        val1(0:kmax) = 0.5 * (vval%d(j1:k1-1) + vval%d(j2:k2-1)) + 0.125 * deta * (vval%d(j1+1:k1) - vval%d(j2+1:k2))
        err1 = abs(val1 - vval%d(j3:k3-1)) / (abs(vval%d(j3:k3-1)) + ATOL)

        ! logarithmic interpolation from existing data
        lg1 = .false.
        do k = 0, kmax
          ! values must have the same sign and not be zero
          if (vval%d(j1+k) * vval%d(j2+k) <= 0) cycle
          sgn = sign(1.0, vval%d(j1+k))

          ! logarithmic interpolation
          val2(k) = sgn * sqrt(vval%d(j1+k) * vval%d(j2+k)) * exp(0.125 * deta * (vval%d(j1+k+1) / vval%d(j1+k) - vval%d(j2+k+1) / vval%d(j2+k)))

          ! decide if direct or logarithmic interpolation is better
          err2(k) = abs(val2(k) - vval%d(j3+k)) / (abs(vval%d(j3+k)) + ATOL)
          if (err2(k) < err1(k)) lg1(k) = .true.
        end do

        ! lg for interval i3 --- i2
        call vlg%push(lg1)

        ! check if refinement is necessary
        if (any(merge(err2, err1, lg1) > RTOL)) then
          ! add 2 new intervals to stack
          call stack%push([i1, i3])
          call stack%push([i3, i2])
        else
          ! lg for interval i1 -- i3
          vlg%d(((i1 - 1) * (kmax + 1) + 1):(i1 * (kmax + 1))) = lg1
        end if
      end do

      ! store lg for coarse grid
      lg(:,i) = vlg%d(((i - 1) * (kmax + 1) + 1):(i * (kmax + 1)))

      print "(I4,A,I4)", i, " / ", N0 - 1
    end do
    !$omp end do nowait

    ! save number of points generated by this thread
    nth(ithread + 1) = veta%n - N0

    !$omp barrier

    !$omp single

    ! count elements
    nth(0) = 1
    nth(1) = N0
    do i = 1, nthreads + 1
      nth(i) = nth(i) + nth(i-1)
    end do
    n = nth(nthreads + 1) - 1

    ! allocate global memory + fill in coarse grid
    allocate (this%eta(n), this%val(0:kmax+1,n), this%lg(0:kmax,n))
    this%eta(         1:N0) = eta
    this%val(0:kmax+1,1:N0) = val
    this%lg( 0:kmax,  1:N0) = lg

    !$omp end single

    ! copy local values to global memory
    j1 = N0 * (kmax + 2) + 1
    k1 = veta%n * (kmax + 2)
    this%eta(  nth(ithread):nth(ithread+1)-1) = veta%d(N0+1:veta%n)
    this%val(:,nth(ithread):nth(ithread+1)-1) = reshape(vval%d(j1:k1), [kmax+2, veta%n - N0])
    this%lg( :,nth(ithread):nth(ithread+1)-1) = reshape( vlg%d(j1-N0:vlg%n), [kmax+1, veta%n - N0])

    ! cleanup thread-local memory
    call stack%destruct()
    call vlg%destruct()
    call veta%destruct()
    call vval%destruct()

    !$omp end parallel

    ! sort by eta
    allocate (perm(size(this%eta)))
    call qsort(this%eta, perm = perm)
    this%val = this%val(:,perm)
    this%lg  = this%lg( :,perm)

    ! save to file
    call this%save("/tmp")
  end subroutine

  subroutine distribution_table_gen(this, dos, dist, eta, val)
    !! generate single entry
    class(distribution_table), intent(in)  :: this
    procedure(density_of_states)           :: dos
    procedure(distribution_density)        :: dist
    real,                      intent(in)  :: eta
      !! chemical potential
    real,                      intent(out) :: val(0:)
      !! output value + derivatives

    integer       :: k
    real(kind=16) :: inf16, p16(0), dFda16, dFdb16, dFdp16(0), val16

    inf16 = ieee_value(1.0_16, IEEE_POSITIVE_INF)

    do k = 0, this%kmax + 1
      call quad(integrand, -real(eta,kind=16), INF16, p16, val16, dFda16, dFdb16, dFdp16, rtol = 1e-14_16, max_levels = 20)
      val(k) = real(val16)
    end do

  contains

    subroutine integrand(t, p, g, dgdt, dgdp)
      real(kind=16), intent(in)  :: t
      real(kind=16), intent(in)  :: p(:)
      real(kind=16), intent(out) :: g
      real(kind=16), intent(out) :: dgdt
      real(kind=16), intent(out) :: dgdp(:)

      m4_ignore(p)
      dgdt = 0
      m4_ignore(dgdp)

      g = dos(t + eta) * dist(t, k)
      if (mod(k, 2) == 1) g = -g
    end subroutine

  end subroutine

  subroutine distribution_table_get(this, eta, k, val, dvaldeta)
    !! get k-th derivative of cumulative distribution function from table
    class(distribution_table), intent(in)  :: this
    real,                      intent(in)  :: eta
      !! chemical potential
    integer,                   intent(in)  :: k
      !! derivative, must be between 0 and kmax
    real,                      intent(out) :: val
      !! output k-th derivative of cumulative distribution function
    real,                      intent(out) :: dvaldeta
      !! output derivative of val wrt eta

    integer :: i
    real    :: eta1, eta2, deta, val1, dval1, val2, dval2, t, h00, h10, h01, h11, g00, g10, g01, g11

    m4_assert((eta >= this%eta(1)) .and. (eta <= this%eta(size(this%eta))))
    m4_assert((k >= 0) .and. (k <= this%kmax))

    ! find interval
    i = bin_search(this%eta, eta, mode = BS_LESS)
    if (i == ubound(this%eta, 1)) i = i - 1

    eta1 = this%eta(i)
    eta2 = this%eta(i+1)
    deta = eta2 - eta1
    t    = (eta - eta1) / deta

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    g00 = 6 * t * (t - 1)
    g10 = (3 * t - 1) * (t - 1)
    g01 = - 6 * t * (t - 1)
    g11 = t * (3 * t - 2)

    val1  = this%val(k,  i  )
    dval1 = this%val(k+1,i  )
    val2  = this%val(k,  i+1)
    dval2 = this%val(k+1,i+1)

    if (this%lg(k,i)) then
      ! logarithmic interpolation
      val      =  h00 * log(abs(val1)) + h10 * deta * dval1 / val1 + h01 * log(abs(val2)) + h11 * deta * dval2 / val2
      dvaldeta = (g00 * log(abs(val1)) + g10 * deta * dval1 / val1 + g01 * log(abs(val2)) + g11 * deta * dval2 / val2) / deta

      val      = sign(1.0, val1) * exp(val)
      dvaldeta = val * dvaldeta
    else
      ! direct interpolation
      val      =  h00 * val1 + h10 * deta * dval1 + h01 * val2 + h11 * deta * dval2
      dvaldeta = (g00 * val1 + g10 * deta * dval1 + g01 * val2 + g11 * deta * dval2) / deta
    end if
  end subroutine

  subroutine distribution_table_inv(this, F, eta, detadF)
    !! get inverse of distribution function (only for k = 0)
    class(distribution_table), intent(in)  :: this
    real,                      intent(in)  :: F
      !! value of distribution function
    real,                      intent(out) :: eta
      !! output corresponding eta
    real,                      intent(out) :: detadF
    !! output derivative of eta wrt F

    real,    parameter :: ATOL   = 5e-13
    integer, parameter :: MAX_IT = 10

    integer :: i, it
    real    :: eta1, eta2, deta, F1, dF1, F2, dF2, t, t_min, t_max, t_old, res, dresdt, dt, err

    m4_assert((F >= this%val(0,1)) .and. (F <= this%val(0,size(this%eta))))

    ! find interval
    i = bin_search(this%val(0,:), F, mode = BS_LESS)

    eta1 = this%eta(i)
    eta2 = this%eta(i+1)
    deta = eta2 - eta1

    F1   = this%val(0,i  )
    dF1  = this%val(1,i  )
    F2   = this%val(0,i+1)
    dF2  = this%val(1,i+1)

    if (i == ubound(this%val, 2)) i = i - 1

    ! start value
    if (this%lg(0,i)) then
      t = log(F / F1) / log(F2 / F1)
    else
      t = (F - F1) / (F2 - F1)
    end if

    ! bounds
    t_min = 0
    t_max = 1

    ! Newton iteration to get t
    err = huge(1.0)
    it  = 0
    do while (err > ATOL)
      it = it + 1
      if (it > MAX_IT) then
        print "(A,ES25.16E3)", "F = ", F
        call program_error("Newton did not converge")
      end if

      ! evaluate resiudal and get Newton update
      call residual(t, res, dresdt)
      dt = - res / dresdt
      err = abs(dt)

      ! update bounds
      if (dt > 0) then
        t_min = t
      else
        t_max = t
      end if

      ! update solution
      t_old = t
      t     = t + dt

      ! bisection
      if ((t < t_min) .or. (t > t_max) .or. ((t_old == t_min) .and. (t == t_max))) then
        t = 0.5 * (t_min + t_max)
        err = min(err, t_max - t_min)
      end if
    end do

    ! get eta and derivative
    call residual(t, res, dresdt)
    eta    = eta1 + deta * t
    detadF = deta / dresdt

  contains

    subroutine residual(t, res, dresdt)
      real, intent(in)  :: t
      real, intent(out) :: res
      real, intent(out) :: dresdt

      real :: h00, h10, h01, h11, g00, g10, g01, g11

      h00 = (1 + 2 * t) * (1 - t)**2
      h10 = t * (1 - t)**2
      h01 = t**2 * (3 - 2 * t)
      h11 = t**2 * (t - 1)

      g00 = 6 * t * (t - 1)
      g10 = (3 * t - 1) * (t - 1)
      g01 = - 6 * t * (t - 1)
      g11 = t * (3 * t - 2)

      if (this%lg(0,i)) then
        ! logarithmic interpolation
        res    = h00 * log(F1) + h10 * deta * dF1 / F1 + h01 * log(F2) + h11 * deta * dF2 / F2
        dresdt = g00 * log(F1) + g10 * deta * dF1 / F1 + g01 * log(F2) + g11 * deta * dF2 / F2

        res    = exp(res)
        dresdt = res * dresdt
        res    = res - F
      else
        ! direct interpolation
        res    = h00 * F1 + h10 * deta * dF1 + h01 * F2 + h11 * deta * dF2 - F
        dresdt = g00 * F1 + g10 * deta * dF1 + g01 * F2 + g11 * deta * dF2
      end if
    end subroutine

  end subroutine

  subroutine distribution_table_output(this, fname, eta1, eta2, k, nsamples)
    !! output distribution table to file
    class(distribution_table), intent(in) :: this
    character(*),              intent(in) :: fname
      !! filename
    real,                      intent(in) :: eta1
      !! lower bound
    real,                      intent(in) :: eta2
      !! upper bound
    integer,                   intent(in) :: k
      !! output k-th derivative
    integer,                   intent(in) :: nsamples
      !! number of points

    integer           :: i, funit
    real              :: val, dval
    real, allocatable :: eta(:)

    m4_assert((k >= 0) .and. (k <= this%kmax))

    eta = linspace(eta1, eta2, nsamples)

    open (newunit = funit, file = fname, status = "replace", action = "write")
    do i = 1, nsamples
      call this%get(eta(i), k, val, dval)
      write (funit, "(2ES25.16E3)") eta(i), val
    end do
    close (funit)
  end subroutine

  subroutine distribution_table_load(this, dir, eta_min, eta_max, kmax, status)
    !! load distribution table from file
    class(distribution_table), intent(inout) :: this
    character(*),              intent(in)    :: dir
      !! directory
    real,                      intent(in)    :: eta_min
    real,                      intent(in)    :: eta_max
    integer,                   intent(in)    :: kmax
    logical,                   intent(out)   :: status
      !! success (true) or fail (false)

    character(256) :: fname
    integer        :: h, funit, num_entries

    status = .false.
    if (distribution_table_disable_save) return

    ! filename
    write (fname, "(2ES25.16E3,I0)") eta_min, eta_max, kmax
    h = hash(fname)
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    ! load
    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")
    read (funit) num_entries
    allocate (this%eta(num_entries), this%val(0:kmax+1,num_entries), this%lg(0:kmax,num_entries))
    read (funit) this%eta
    read (funit) this%val
    read (funit) this%lg
    close (funit)
    status = .true.
  end subroutine

  subroutine distribution_table_save(this, dir)
    !! save distribution table to file
    class(distribution_table), intent(in) :: this
    character(*),              intent(in) :: dir
      !! directory

    character(256) :: fname
    integer        :: h, funit, num_entries

    if (distribution_table_disable_save) return

    ! filename
    num_entries = size(this%eta)
    write (fname, "(2ES25.16E3,I0)") this%eta(1), this%eta(num_entries), this%kmax
    h = hash(fname)
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) num_entries
    write (funit) this%eta
    write (funit) this%val
    write (funit) this%lg
    close (funit)
  end subroutine

end module
