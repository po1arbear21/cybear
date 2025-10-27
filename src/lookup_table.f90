m4_include(util/macro.f90.inc)

module lookup_table_m

  use, intrinsic :: iso_fortran_env, only: real128
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF, ieee_is_finite

  use bin_search_m,    only: bin_search, BS_LESS
  use error_m,         only: assert_failed, program_error
  use math_m,          only: linspace, expm1, log1p
  use omp_lib,         only: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
  use quad_m,          only: quad
  use qsort_m,         only: qsort
  use util_m,          only: hash
  use vector_m,        only: vector_int, vector_log, vector_real128

  implicit none

  private
  public :: par_int_table, def_int_table

  interface
    subroutine integrand(t, x, f, dfdt, dfdx)
      !! integrand
      import real128
      real(real128), intent(in)  :: t
        !! inner argument
      real(real128), intent(in)  :: x(:)
        !! parameters
      real(real128), intent(out) :: f
        !! output integrand
      real(real128), intent(out) :: dfdt
        !! output derivative of integrand wrt t
      real(real128), intent(out) :: dfdx(:)
        !! output derivative of integrand wrt x
    end subroutine
  end interface

  type par_int_table
    !! lookup table for parameter-dependent integrals
    !! q(x) = integral_{a}^{b} f(t,x) dt

    character(:), allocatable :: name
      !! name for loading and saving

    real(real128), allocatable :: x(:)
      !! sample points
    real(real128), allocatable :: q(:)
      !! function values at sample points
    real(real128), allocatable :: d(:)
      !! function derivatives at sample points
    logical,       allocatable :: l(:)
      !! logarithmic interpolation flags
  contains
    procedure :: init   => par_int_table_init
    procedure :: get    => par_int_table_get
    procedure :: inv    => par_int_table_inv
    procedure :: output => par_int_table_output
    procedure :: load   => par_int_table_load
    procedure :: save   => par_int_table_save
  end type

  type def_int_table
    !! lookup table for definite integrals
    !! q(a,b) = integral_{a}^{b} f(x) dx

    character(:), allocatable :: name
      !! name for loading and saving

    integer                    :: ng
      !! number of gauss-legendre nodes
    real(real128), allocatable :: xg(:)
      !! gauss-legendre nodes
    real(real128), allocatable :: wg(:)
      !! gauss-legendre weights

    real(real128), allocatable :: p(:)
      !! parameters passed to integrand

    real(real128), allocatable :: x(:)
      !! x nodes
    real(real128), allocatable :: q(:)
      !! q(iq) = integral_{x(i1)}^{x(i2)} f(t) dt
    integer,       allocatable :: i(:)
      !! indices of x interval midpoint (i(iq) = i3)
  contains
    procedure :: init => def_int_table_init
    procedure :: get  => def_int_table_get
    procedure :: load => def_int_table_load
    procedure :: save => def_int_table_save

    procedure, private :: gauss => def_int_table_gauss
  end type

  real(real128), parameter :: xg1(1) = [            &
    0.000000000000000000000000000000000000_real128  &
  ]
  real(real128), parameter :: wg1(1) = [            &
    1.000000000000000000000000000000000000_real128  &
  ]
  real(real128), parameter :: xg2(2) = [            &
   -0.577350269189625764509148780501957456_real128, &
    0.577350269189625764509148780501957456_real128  &
  ]
  real(real128), parameter :: wg2(2) = [            &
    1.000000000000000000000000000000000000_real128, &
    1.000000000000000000000000000000000000_real128  &
  ]
  real(real128), parameter :: xg3(3) = [            &
   -0.774596669241483377035853079956479922_real128, &
    0.000000000000000000000000000000000000_real128, &
    0.774596669241483377035853079956479922_real128  &
  ]
  real(real128), parameter :: wg3(3) = [            &
    0.555555555555555555555555555555555556_real128, &
    0.888888888888888888888888888888888889_real128, &
    0.555555555555555555555555555555555556_real128  &
  ]
  real(real128), parameter :: xg4(4) = [            &
   -0.861136311594052575223946488892809505_real128, &
   -0.339981043584856264802665759103244687_real128, &
    0.339981043584856264802665759103244687_real128, &
    0.861136311594052575223946488892809505_real128  &
  ]
  real(real128), parameter :: wg4(4) = [            &
    0.347854845137453857373063949221999407_real128, &
    0.652145154862546142626936050778000593_real128, &
    0.652145154862546142626936050778000593_real128, &
    0.347854845137453857373063949221999407_real128  &
  ]
  real(real128), parameter :: xg5(5) = [            &
   -0.906179845938663992797626878299392965_real128, &
   -0.538469310105683091036314420700208805_real128, &
    0.000000000000000000000000000000000000_real128, &
    0.538469310105683091036314420700208805_real128, &
    0.906179845938663992797626878299392965_real128  &
  ]
  real(real128), parameter :: wg5(5) = [            &
    0.236926885056189087514264040719917363_real128, &
    0.478628670499366468041291514835638193_real128, &
    0.568888888888888888888888888888888889_real128, &
    0.478628670499366468041291514835638193_real128, &
    0.236926885056189087514264040719917363_real128  &
  ]

  integer, parameter :: QUAD_MIN_LVL = 10
  integer, parameter :: QUAD_MAX_LVL = 20

contains

  subroutine par_int_table_init(this, dir, name, xmin, xmax, f, a, b, rtol, atol, xtol)
    class(par_int_table), intent(out) :: this
    character(*),         intent(in)  :: dir
      !! directory
    character(*),         intent(in)  :: name
      !! name for loading and saving
    real(real128),        intent(in)  :: xmin
      !! minimal supported x
    real(real128),        intent(in)  :: xmax
      !! maximal supported x
    procedure(integrand)              :: f
      !! integrand
    real(real128),        intent(in)  :: a
      !! lower integration boundary
    real(real128),        intent(in)  :: b
      !! upper integration boundary
    real(real128),        intent(in)  :: rtol
      !! relative error tolerance
    real(real128),        intent(in)  :: atol
      !! absolute error tolerance
    real(real128),        intent(in)  :: xtol
      !! minimum interval size for x

    integer, parameter :: MIN_LVL = 10

    integer              :: lvl, max_lvl, i, i1, i2, ithread, n, n0, nthreads
    integer, allocatable :: ns(:), nt(:), perm(:)
    logical              :: status, l1
    real(real128)        :: dx, x, q, d, q1, q2, e1, e2, sgn
    type(vector_real128) :: gx, tx, gq, tq, gd, td
    type(vector_log)     :: gl, tl
    type(vector_int)     :: gs, ts

    ! try to load from file
    this%name = name
    call this%load(dir, xmin, xmax, a, b, rtol, atol, xtol, status)
    if (status) return

    print "(2A)", "Generating table ", name

    max_lvl = max(ceiling(log((xmax - xmin) / xtol) / log(2.0)), MIN_LVL)

    ! init global vectors
    call gx%init(0, c = 128)
    call gq%init(0, c = 128)
    call gd%init(0, c = 128)
    call gl%init(0, c = 128)
    call gs%init(0, c = 128)

    nthreads = omp_get_max_threads()
    allocate (ns(0:nthreads), nt(0:nthreads), source = 0)

    ! 0-th level
    call gx%push(xmin)
    call gl%push(.false.)
    call gen(xmin, q, d)
    call gq%push(q)
    call gd%push(d)
    call gx%push(xmax)
    call gl%push(.false.)
    call gen(xmax, q, d)
    call gq%push(q)
    call gd%push(d)
    call gs%push([1, 2])

    do lvl = 1, max_lvl
      n  = gs%n / 2
      n0 = gx%n
      ns = 0
      nt = 0

      print "(I2,A,I2,A,I0)", lvl, " / ", max_lvl, ": ", n

      !$omp parallel default(none) &
      !$omp private(i,i1,i2,ithread,l1,dx,x,q,d,q1,q2,e1,e2,sgn,tx,tq,td,tl,ts) &
      !$omp shared(this,rtol,atol,gx,gq,gd,gl,gs,lvl,n,n0,nthreads,ns,nt)

      ithread = omp_get_thread_num() + 1

      ! init thread-local memory
      call tx%init(0, c = 128)
      call tq%init(0, c = 128)
      call td%init(0, c = 128)
      call tl%init(0, c = 128)
      call ts%init(0, c = 128)

      !$omp do schedule(dynamic)
      do i = 1, n
        ! get interval from stack
        i1 = gs%d(i*2-1)
        i2 = gs%d(i*2  )

        ! interval size and midpoint
        dx = gx%d(i2) - gx%d(i1)
        x  = 0.5 * (gx%d(i1) + gx%d(i2))

        ! add interval midpoint
        call gen(x, q, d)
        call tx%push(x)
        call tq%push(q)
        call td%push(d)

        ! direct interpolation from existing data
        q1 = 0.5 * (gq%d(i1) + gq%d(i2)) + 0.125 * dx * (gd%d(i1) - gd%d(i2))
        e1 = abs(q1 - q) / (abs(q) + atol)

        ! logarithmic interpolation from existing data
        l1 = .false.
        if (gq%d(i1) * gq%d(i2) > 0) then
          sgn = sign(1.0_real128, gq%d(i1))

          ! logarithmic interpolation
          q2 = sgn * sqrt(gq%d(i1) * gq%d(i2)) * exp(0.125 * dx * (gd%d(i1) / gq%d(i1) - gd%d(i2) / gq%d(i2)))

          ! decide if direct or logarithmic interpolation is better
          e2 = abs(q2 - q) / (abs(q) + atol)
          if (e2 < e1) l1 = .true.
        end if

        ! l for interval x --- x2
        call tl%push(l1)

        ! check if refinement is necessary
        if ((lvl < MIN_LVL) .or. (merge(e2, e1, l1) > rtol)) then
          ! add 2 new intervals to stack
          call ts%push([i1, n0 + tx%n, n0 + tx%n, i2])
        else
          ! l for interval x1 --- x
          gl%d(i1) = l1
        end if
      end do
      !$omp end do nowait

      ! store number of generated points and stack intervals in global memory
      ns(ithread) = ts%n
      nt(ithread) = tx%n

      !$omp barrier

      !$omp single

      ! accumulate numbers to get start indices in global memory
      do i = 2, nthreads
        ns(i) = ns(i-1) + ns(i)
        nt(i) = nt(i-1) + nt(i)
      end do

      ! resize global memory
      call gx%resize(n0 + nt(nthreads))
      call gq%resize(n0 + nt(nthreads))
      call gd%resize(n0 + nt(nthreads))
      call gl%resize(n0 + nt(nthreads))
      call gs%resize(     ns(nthreads))

      !$omp end single

      ! convert indices on stack to global indices
      do i = 1, ts%n / 2
        if (ts%d(2*i - 1) > n0) ts%d(2*i - 1) = ts%d(2*i - 1) + nt(ithread-1)
        if (ts%d(2*i    ) > n0) ts%d(2*i    ) = ts%d(2*i    ) + nt(ithread-1)
      end do

      ! copy thread-local values to global memory
      gx%d((n0 + nt(ithread-1) + 1):(n0 + nt(ithread))) = tx%d(1:tx%n)
      gq%d((n0 + nt(ithread-1) + 1):(n0 + nt(ithread))) = tq%d(1:tq%n)
      gd%d((n0 + nt(ithread-1) + 1):(n0 + nt(ithread))) = td%d(1:td%n)
      gl%d((n0 + nt(ithread-1) + 1):(n0 + nt(ithread))) = tl%d(1:tl%n)
      gs%d((ns(ithread-1) + 1):ns(ithread)) = ts%d(1:ts%n)

      ! cleanup thread-local memory
      call tx%destruct()
      call tq%destruct()
      call td%destruct()
      call tl%destruct()
      call ts%destruct()

      !$omp end parallel

      ! check if refinement done
      if (gs%n == 0) exit
    end do

    if (gs%n > 0) print *, "Warning: some intervals are imprecise"

    this%x = gx%to_array()
    this%q = gq%to_array()
    this%d = gd%to_array()
    this%l = gl%to_array()

    ! sort by x
    allocate (perm(size(this%x)))
    call qsort(this%x, perm = perm)
    this%q = this%q(perm)
    this%d = this%d(perm)
    this%l = this%l(perm)

    ! save
    call this%save(dir, xmin, xmax, a, b, rtol, atol, xtol)

  contains

    subroutine gen(x, q, d)
      real(real128), intent(in)  :: x
      real(real128), intent(out) :: q
      real(real128), intent(out) :: d

      real(real128) :: p(1), dqda, dqdb, dqdp(1)

      p(1) = x
      if ((x > a) .and. (x < b)) then
        call quad(f, a, b, p, q, dqda, dqdb, dqdp, x0=x, rtol=rtol, min_levels=QUAD_MIN_LVL, max_levels=QUAD_MAX_LVL)
      else
        call quad(f, a, b, p, q, dqda, dqdb, dqdp, rtol=rtol, min_levels=QUAD_MIN_LVL, max_levels=QUAD_MAX_LVL)
      end if
      d = dqdp(1)
    end subroutine

  end subroutine

  subroutine par_int_table_get(this, x, q, dqdx)
    !! get function value and derivative
    class(par_int_table), intent(in)  :: this
    real(real128),        intent(in)  :: x
      !! argument
    real(real128),        intent(out) :: q
      !! output function value
    real(real128),        intent(out) :: dqdx
      !! output derivative of q wrt x

    integer       :: i
    real(real128) :: x1, x2, dx, q1, d1, q2, d2, t, h00, h10, h01, h11, d00, d10, d01, d11

    if ((x < this%x(1)) .or. (x > this%x(size(this%x)))) then
      print "(3ES25.16E3)", x, this%x(1), this%x(size(this%x))
      call program_error("x out of bounds")
    end if

    ! find interval
    i = bin_search(this%x, x, mode = BS_LESS)
    if (i == ubound(this%x, 1)) i = i - 1

    x1 = this%x(i  )
    x2 = this%x(i+1)
    dx = x2 - x1
    t  = (x - x1) / dx

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    d00 = 6 * t * (t - 1)
    d10 = (3 * t - 1) * (t - 1)
    d01 = - 6 * t * (t - 1)
    d11 = t * (3 * t - 2)

    q1 = this%q(i  )
    d1 = this%d(i  )
    q2 = this%q(i+1)
    d2 = this%d(i+1)

    if (this%l(i)) then
      ! logarithmic interpolation
      q    =  h00 * log(abs(q1)) + h10 * dx * d1 / q1 + h01 * log(abs(q2)) + h11 * dx * d2 / q2
      dqdx = (d00 * log(abs(q1)) + d10 * dx * d1 / q1 + d01 * log(abs(q2)) + d11 * dx * d2 / q2) / dx

      q    = sign(1.0_real128, q1) * exp(q)
      dqdx = q * dqdx
    else
      ! direct interpolation
      q    =  h00 * q1 + h10 * dx * d1 + h01 * q2 + h11 * dx * d2
      dqdx = (d00 * q1 + d10 * dx * d1 + d01 * q2 + d11 * dx * d2) / dx
    end if
  end subroutine

  subroutine par_int_table_inv(this, q, x, dxdq)
    !! get inverse of function (must be monotone)
    class(par_int_table), intent(in)  :: this
    real(real128),        intent(in)  :: q
      !! function value
    real(real128),        intent(out) :: x
      !! output argument
    real(real128),        intent(out) :: dxdq
      !! output derivative of x wrt q

    real(real128), parameter :: ATOL   = 1e-30_real128
    integer,       parameter :: MAX_IT = 10

    integer       :: i, it
    logical       :: l
    real(real128) :: x1, x2, dx, q1, d1, q2, d2, t, tmin, tmax, told, res, dresdt, dt, err

    if ((q < this%q(1)) .or. (q > this%q(size(this%x)))) then
      print "(A,ES25.16E3)", "q    = ", q
      print "(A,ES25.16E3)", "qmin = ", this%q(1)
      print "(A,ES25.16E3)", "qmax = ", this%q(size(this%x))
      call program_error("q is out of range")
    end if

    ! find interval
    i = bin_search(this%q, q, mode = BS_LESS)
    if (i == ubound(this%q, 1)) i = i - 1

    ! load interval values
    x1 = this%x(i)
    x2 = this%x(i+1)
    dx = x2 - x1
    q1 = this%q(i  )
    d1 = this%d(i  )
    q2 = this%q(i+1)
    d2 = this%d(i+1)
    l  = this%l(i)

    ! start estimate
    if (l) then
      t = log(q / q1) / log(q2 / q1)
    else
      t = (q - q1) / (q2 - q1)
    end if

    ! search bounds
    tmin = 0
    tmax = 1

    ! Newton iteration to get t
    err = huge(1.0_real128)
    it  = 0
    do while (err > ATOL)
      it = it + 1
      if (it > MAX_IT) then
        print "(A,ES25.16E3)", "q = ", q
        call program_error("Newton did not converge")
      end if

      ! evaluate resiudal and get Newton update
      call residual(t, res, dresdt)
      dt = - res / dresdt
      err = abs(dt)

      ! update bounds
      if (dt > 0) then
        tmin = t
      else
        tmax = t
      end if

      ! update solution
      told = t
      t    = t + dt

      ! bisection
      if ((t < tmin) .or. (t > tmax) .or. ((told == tmin) .and. (t == tmax))) then
        t   = 0.5 * (tmin + tmax)
        err = min(err, tmax - tmin)
      end if
    end do

    ! get x and derivative
    call residual(t, res, dresdt)
    x    = x1 + dx * t
    dxdq = dx / dresdt

  contains

    subroutine residual(t, res, dresdt)
      real(real128), intent(in)  :: t
      real(real128), intent(out) :: res
      real(real128), intent(out) :: dresdt

      real(real128) :: h00, h10, h01, h11, d00, d10, d01, d11

      h00 = (1 + 2 * t) * (1 - t)**2
      h10 = t * (1 - t)**2
      h01 = t**2 * (3 - 2 * t)
      h11 = t**2 * (t - 1)

      d00 = 6 * t * (t - 1)
      d10 = (3 * t - 1) * (t - 1)
      d01 = - 6 * t * (t - 1)
      d11 = t * (3 * t - 2)

      if (l) then
        ! logarithmic interpolation
        res    = h00 * log(q1) + h10 * dx * d1 / q1 + h01 * log(q2) + h11 * dx * d2 / q2
        dresdt = d00 * log(q1) + d10 * dx * d1 / q1 + d01 * log(q2) + d11 * dx * d2 / q2

        res    = exp(res)
        dresdt = res * dresdt
        res    = res - q
      else
        ! direct interpolation
        res    = h00 * q1 + h10 * dx * d1 + h01 * q2 + h11 * dx * d2 - q
        dresdt = d00 * q1 + d10 * dx * d1 + d01 * q2 + d11 * dx * d2
      end if
    end subroutine

  end subroutine

  subroutine par_int_table_output(this, fname, x1, x2, nsamples)
    !! output parametric integral table to file
    class(par_int_table), intent(in) :: this
    character(*),         intent(in) :: fname
      !! filename
    real(real128),        intent(in) :: x1
      !! lower bound
    real(real128),        intent(in) :: x2
      !! upper bound
    integer,              intent(in) :: nsamples
      !! number of points

    integer                    :: i, funit
    real(real128)              :: q, dqdx
    real(real128), allocatable :: x(:)

    x = linspace(x1, x2, nsamples)

    open (newunit = funit, file = fname, status = "replace", action = "write")
    do i = 1, nsamples
      call this%get(x(i), q, dqdx)
      write (funit, "(3ES25.16E3)") x(i), q, dqdx
    end do
    close (funit)
  end subroutine

  subroutine par_int_table_load(this, dir, xmin, xmax, a, b, rtol, atol, xtol, status)
    !! load parametric integral table from file
    class(par_int_table), intent(inout) :: this
    character(*),         intent(in)    :: dir
      !! directory
    real(real128),        intent(in)    :: xmin
      !! minimal supported x
    real(real128),        intent(in)    :: xmax
      !! maximal supported x
    real(real128),        intent(in)    :: a
      !! lower integration boundary
    real(real128),        intent(in)    :: b
      !! upper integration boundary
    real(real128),        intent(in)    :: rtol
      !! relative error tolerance
    real(real128),        intent(in)    :: atol
      !! absolute error tolerance
    real(real128),        intent(in)    :: xtol
      !! minimum interval size for x
    logical,              intent(out)   :: status
      !! success (true) or fail (false)

    character(512) :: fname
    integer        :: h, funit, n

    status = .false.

    ! filename
    write (fname, "(7ES41.32E3)") xmin, xmax, a, b, rtol, atol, xtol
    h = hash(fname) ! hash parameters to get unique name
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    ! load
    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")
    read (funit) n
    allocate (this%x(n), this%q(n), this%d(n), this%l(n))
    read (funit) this%x
    read (funit) this%q
    read (funit) this%d
    read (funit) this%l
    close (funit)
    status = .true.
  end subroutine

  subroutine par_int_table_save(this, dir, xmin, xmax, a, b, rtol, atol, xtol)
    !! save parametric integral table to file
    class(par_int_table), intent(in) :: this
    character(*),         intent(in) :: dir
      !! directory
    real(real128),        intent(in) :: xmin
      !! minimal supported x
    real(real128),        intent(in) :: xmax
      !! maximal supported x
    real(real128),        intent(in) :: a
      !! lower integration boundary
    real(real128),        intent(in) :: b
      !! upper integration boundary
    real(real128),        intent(in) :: rtol
      !! relative error tolerance
    real(real128),        intent(in) :: atol
      !! absolute error tolerance
    real(real128),        intent(in) :: xtol
      !! minimum interval size for x

    character(512) :: fname
    integer        :: h, funit, num_entries

    ! filename
    num_entries = size(this%x)
    write (fname, "(7ES41.32E3)") xmin, xmax, a, b, rtol, atol, xtol
    h = hash(fname)
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) num_entries
    write (funit) this%x
    write (funit) this%q
    write (funit) this%d
    write (funit) this%l
    close (funit)
  end subroutine

  subroutine def_int_table_init(this, dir, name, xmin, xmax, ng, intg, p, rtol, atol, xtol)
    class(def_int_table), intent(out) :: this
    character(*),         intent(in)  :: dir
      !! directory
    character(*),         intent(in)  :: name
      !! name for loading and saving
    real(real128),        intent(in)  :: xmin
      !! minimal supported x
    real(real128),        intent(in)  :: xmax
      !! maximal supported x
    integer,              intent(in)  :: ng
      !! number of gauss-legendre points for integration on lowest level
    procedure(integrand)              :: intg
      !! integrand
    real(real128),        intent(in)  :: p(:)
      !! parameters for integrand
    real(real128),        intent(in)  :: rtol
      !! relative error tolerance
    real(real128),        intent(in)  :: atol
      !! absolute error tolerance
    real(real128),        intent(in)  :: xtol
      !! minimum interval size for x

    integer, parameter :: MIN_LVL = 10

    integer              :: lvl, max_lvl, i, i1, i2, i3, iq, iq0, ithread, n, nx0, nq0, nthreads
    integer, allocatable :: ns(:), nt(:)
    logical              :: status
    real(real128)        :: dx, x, d1, d2, q, q1, q2, dqdp(size(p)), err
    type(vector_real128) :: gx, tx, gq, tq
    type(vector_int)     :: gs, ts, gi

    ! set members
    this%name =  name
    this%ng   = ng
    this%p    =  p

    ! set gauss-legendre quadrature nodes and weights
    select case (ng)
    case (1)
      this%xg = xg1
      this%wg = wg1
    case (2)
      this%xg = xg2
      this%wg = wg2
    case (3)
      this%xg = xg3
      this%wg = wg3
    case (4)
      this%xg = xg4
      this%wg = wg4
    case (5)
      this%xg = xg5
      this%wg = wg5
    case default
      call program_error("ng must be between 1 and 5")
    end select

    ! try to load table from file
    call this%load(dir, xmin, xmax, rtol, atol, xtol, status)
    if (status) return

    print "(2A)", "Generating table ", name

    max_lvl = max(ceiling(log((xmax - xmin) / xtol) / log(2.0)), MIN_LVL)

    ! init global vectors
    call gx%init(0, c = 128)
    call gq%init(0, c = 128)
    call gs%init(0, c = 128)
    call gi%init(0, c = 128)

    nthreads = omp_get_max_threads()
    allocate (ns(0:nthreads), nt(0:nthreads), source = 0)

    ! 0-th level
    call gx%push([xmin, xmax])
    call this%gauss(intg, xmin, xmax, q, d1, d2)
    call gq%push(q)
    call gs%push([1, 2, 1])
    call gi%push([0])

    do lvl = 1, max_lvl
      n   = gs%n / 3
      nx0 = gx%n
      nq0 = gq%n
      ns  = 0
      nt  = 0

      print "(I2,A,I2,A,I0)", lvl, " / ", max_lvl, ": ", n

      !$omp parallel default(none) &
      !$omp private(i,i1,i2,iq,ithread,dx,x,d1,d2,q,q1,q2,dqdp,err,tx,tq,ts) &
      !$omp shared(this,xmin,xmax,rtol,atol,p,gx,gq,gs,gi,lvl,n,nx0,nq0,nthreads,ns,nt)

      ithread = omp_get_thread_num() + 1

      ! init thread-local memory
      call tx%init(0, c = 128)
      call tq%init(0, c = 128)
      call ts%init(0, c = 128)

      !$omp do schedule(dynamic)
      do i = 1, n
        ! get interval from stack
        i1 = gs%d(i*3-2)
        i2 = gs%d(i*3-1)
        iq = gs%d(i*3  )

        ! interval size and midpoint
        dx = gx%d(i2) - gx%d(i1)
        x  = 0.5 * (gx%d(i1) + gx%d(i2))
        call tx%push(x)

        ! left sub-interval
        call this%gauss(intg, gx%d(i1), x, q1, d1, d2)
        call tq%push(q1)

        ! right sub-interval
        call this%gauss(intg, x, gx%d(i2), q2, d1, d2)
        call tq%push(q2)

        ! refinement: compare sum over the two sub-intervals to previous value for complete interval
        q   = q1 + q2
        err = abs(q - gq%d(iq)) / (abs(q) + atol)
        if ((lvl < MIN_LVL) .or. (err > rtol)) call ts%push([i1, nx0 + tx%n, tq%n-1, nx0 + tx%n, i2, tq%n])

        ! update structure
        gs%d(i*3-2) = ithread ! save thread number temporarily in gs
        gi%d(iq) = nx0 + tx%n
      end do
      !$omp end do nowait

      ! store number of generated points and stack intervals in global memory
      ns(ithread) = ts%n
      nt(ithread) = tx%n ! tq%n == tx%n * 2

      !$omp barrier

      !$omp single

      ! accummulate numbers to get start indices in global memory
      do i = 2, nthreads
        ns(i) = ns(i-1) + ns(i)
        nt(i) = nt(i-1) + nt(i)
      end do

      ! convert interval pointers to global indices, use gs to get iq and ithread
      do i = 1, n
        iq = gs%d(i*3)
        gi%d(iq) = gi%d(iq) + nt(gs%d(i*3-2)-1)
      end do

      ! resize global memory
      call gx%resize(nx0 +     nt(nthreads))
      call gq%resize(nq0 + 2 * nt(nthreads))
      call gs%resize(          ns(nthreads))
      call gi%resize(gi%n + 2*n)
      gi%d(gi%n-2*n+1:gi%n) = 0

      !$omp end single

      ! convert indices on stack to global indices
      do i = 1, ts%n / 3
        if (ts%d(3*i - 2) > nx0) ts%d(3*i - 2) = ts%d(3*i - 2) + nt(ithread-1)
        if (ts%d(3*i - 1) > nx0) ts%d(3*i - 1) = ts%d(3*i - 1) + nt(ithread-1)
        ts%d(3*i) = ts%d(3*i) + nq0 + 2 * nt(ithread-1)
      end do

      ! copy thread-local values to global memory
      gx%d((nx0 +   nt(ithread-1) + 1):(nx0 +   nt(ithread))) = tx%d(1:tx%n)
      gq%d((nq0 + 2*nt(ithread-1) + 1):(nq0 + 2*nt(ithread))) = tq%d(1:tq%n)
      gs%d((ns(ithread-1) + 1):ns(ithread)) = ts%d(1:ts%n)

      ! cleanup thread-local memory
      call tx%destruct()
      call tq%destruct()
      call ts%destruct()

      !$omp end parallel

      ! check if refinement done
      if (gs%n == 0) exit
    end do

    if (gs%n > 0) print *, "Warning: some intervals are imprecise"

    ! delete q values of upper levels
    call gs%resize(0)
    call gs%push(1)
    do while (gs%n > 0)
      iq = gs%back()
      call gs%pop()

      ! check if leaf node
      i3 = gi%d(iq)
      if (i3 == 0) cycle

      ! not a leaf node: clear value
      gq%d(iq) = 0

      ! add children to stack
      call gs%push([2*(i3-2)+1, 2*(i3-2)])
    end do

    ! recalculate upper levels by summing over lower levels
    call gs%push(1)
    do while (gs%n > 0)
      ! go to left-most leaf of subtree
      do
        iq = gs%back()
        i3 = gi%d(iq)
        if (i3 == 0) exit ! found leaf

        ! add left sub-tree to stack
        call gs%push(2*(i3-2))
      end do

      iq0 = iq
      do
        ! remove node from stack
        iq = iq0
        call gs%pop()

        ! root node
        if (gs%n == 0) exit

        ! update q value of direct parent
        iq0 = gs%back()
        gq%d(iq0) = gq%d(iq0) + gq%d(iq)

        ! check if we come from left sub-tree
        i3 = gi%d(iq0)
        if (iq == 2*(i3-2)) then
          ! add right sub-tree to stack and go back to top loop
          call gs%push(2*(i3-2)+1)
          exit
        end if
      end do
    end do

    this%x = gx%to_array()
    this%q = gq%to_array()
    this%i = gi%to_array()

    ! save
    call this%save(dir, xmin, xmax, rtol, atol, xtol)
  end subroutine

  recursive subroutine def_int_table_get(this, intg, a, b, q, dqda, dqdb)
    !! get value of integral
    class(def_int_table), intent(in)  :: this
    procedure(integrand)              :: intg
      !! integrand
    real(real128),        intent(in)  :: a
      !! lower integration bound
    real(real128),        intent(in)  :: b
      !! upper integration bound
    real(real128),        intent(out) :: q
      !! output value of integral
    real(real128),        intent(out) :: dqda
      !! output derivative of q wrt a
    real(real128),        intent(out) :: dqdb
      !! output derivative of q wrt b

    integer       :: i1, i2, i3, j1, j2, j3, iq
    real(real128) :: q1, dum

    if (b < a) then
      call this%get(intg, b, a, q, dqdb, dqda)
      q    = - q
      dqda = - dqda
      dqdb = - dqdb
      return
    end if

    ! reset output
    q    = 0
    dqda = 0
    dqdb = 0

    ! start at root node
    i1 = 1
    i2 = 2
    iq = 1

    if ((a < this%x(i1)) .or. (b > this%x(i2))) then
      print "(A,ES41.32E3)", "a    = ", a
      print "(A,ES41.32E3)", "b    = ", b
      print "(A,ES41.32E3)", "xmin = ", this%x(i1)
      print "(A,ES41.32E3)", "xmin = ", this%x(i2)
      call program_error("out of bounds")
    end if

    ! find last common ancestor
    do
      i3 = this%i(iq) ! x-index of interval midpoint

      if (i3 == 0) then
        ! both integration boundaries lie within the same unsplit interval => use gauss-legendre
        call this%gauss(intg, a, b, q, dqda, dqdb)
        return
      end if

      if ((a <= this%x(i3)) .and. (b <= this%x(i3))) then
        ! both integration boundaries lie within the left sub-interval
        i2 = i3
        iq = 2 * (i3 - 2) ! go to left sub-tree
      elseif ((a >= this%x(i3)) .and. (b >= this%x(i3))) then
        ! both integration boundaries lie within the right sub-interval
        i1 = i3
        iq = 2 * (i3 - 2) + 1 ! go to right sub-tree
      else ! (a < this%x(i3)) .and. (b > this%x(i3))
        ! found last common ancestor, save i1, i2 and i3
        j1 = i1
        j2 = i2
        j3 = i3
        exit
      end if
    end do

    ! integrate from lower bound to midpoint
    i1 = j1
    i2 = j3
    iq = 2 * (j3 - 2)
    do
      i3 = this%i(iq)

      if (i3 == 0) then
        ! unsplit interval => use gauss-legendre
        call this%gauss(intg, a, this%x(i2), q1, dqda, dum)
        q = q + q1
        exit
      end if

      if (a <= this%x(i3)) then
        ! go to left sub-tree
        i2 = i3
        iq = 2 * (i3 - 2)

        ! add complete integral over right sub-tree
        q = q + this%q(iq + 1)
      else
        ! go to right sub-tree
        i1 = i3
        iq = 2 * (i3 - 2) + 1
      end if
    end do

    ! integrate from midpoint to upper bound
    i1 = j3
    i2 = j2
    iq = 2 * (j3 - 2) + 1
    do
      i3 = this%i(iq)

      if (i3 == 0) then
        ! unsplit interval => use gauss-legendre
        call this%gauss(intg, this%x(i1), b, q1, dum, dqdb)
        q = q + q1
        exit
      end if

      if (b < this%x(i3)) then
        ! go to left sub-tree
        i2 = i3
        iq = 2 * (i3 - 2)
      else
        ! go to right sub-tree
        i1 = i3
        iq = 2 * (i3 - 2) + 1

        ! add complete integral over left sub-tree
        q = q + this%q(iq - 1)
      end if
    end do
  end subroutine

  subroutine def_int_table_load(this, dir, xmin, xmax, rtol, atol, xtol, status)
    !! load definite integral table from file
    class(def_int_table), intent(inout) :: this
    character(*),         intent(in)    :: dir
      !! directory
    real(real128),        intent(in)    :: xmin
      !! minimal supported x
    real(real128),        intent(in)    :: xmax
      !! maximal supported x
    real(real128),        intent(in)    :: rtol
      !! relative error tolerance
    real(real128),        intent(in)    :: atol
      !! absolute error tolerance
    real(real128),        intent(in)    :: xtol
      !! minimum interval size for x
    logical,              intent(out)   :: status
      !! success (true) or fail (false)

    character(:), allocatable :: tmp
    character(32)             :: fmt
    character(256)            :: fname
    integer                   :: n, h, funit, nx, nq, ni

    status = .false.

    ! filename
    n = 5 + size(this%p)
    write (fmt, "(A,I0,A)") "(", n, "ES41.32E3,I1)"
    n = n * 41 + 1
    allocate (character(len=n) :: tmp)
    write (tmp, fmt) xmin, xmax, rtol, atol, xtol, this%p, this%ng
    h = hash(tmp) ! hash parameters to get unique name
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    ! load
    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")
    read (funit) nx
    read (funit) nq
    read (funit) ni
    allocate (this%x(nx), this%q(nq), this%i(ni))
    read (funit) this%x
    read (funit) this%q
    read (funit) this%i
    close (funit)
    status = .true.
  end subroutine

  subroutine def_int_table_save(this, dir, xmin, xmax, rtol, atol, xtol)
    !! save definite integral table to file
    class(def_int_table), intent(in) :: this
    character(*),         intent(in) :: dir
      !! directory
    real(real128),        intent(in) :: xmin
      !! minimal supported x
    real(real128),        intent(in) :: xmax
      !! maximal supported x
    real(real128),        intent(in) :: rtol
      !! relative error tolerance
    real(real128),        intent(in) :: atol
      !! absolute error tolerance
    real(real128),        intent(in) :: xtol
      !! minimum interval size for x

    character(:), allocatable :: tmp
    character(32)             :: fmt
    character(256)            :: fname
    integer                   :: n, h, funit

    ! filename
    n = 5 + size(this%p)
    write (fmt, "(A,I0,A)") "(", n, "ES41.32E3,I1)"
    n = n * 41 + 1
    allocate (character(len=n) :: tmp)
    write (tmp, fmt) xmin, xmax, rtol, atol, xtol, this%p, this%ng
    h = hash(tmp) ! hash parameters to get unique name
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    ! load
    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) size(this%x)
    write (funit) size(this%q)
    write (funit) size(this%i)
    write (funit) this%x
    write (funit) this%q
    write (funit) this%i
    close (funit)
  end subroutine

  subroutine def_int_table_gauss(this, intg, a, b, q, dqda, dqdb)
    !! Gauss-Legendre quadrature
    class(def_int_table), intent(in)  :: this
    procedure(integrand)              :: intg
      !! integrand
    real(real128),        intent(in)  :: a
      !! lower integration boundary
    real(real128),        intent(in)  :: b
      !! upper integration boundary
    real(real128),        intent(out) :: q
      !! output value of integral
    real(real128),        intent(out) :: dqda
      !! output derivative of q wrt a
    real(real128),        intent(out) :: dqdb
      !! output derivative of q wrt b

    integer       :: i
    real(real128) :: x, w, f, dfdx, dfdp(size(this%p))

    q    = 0
    dqda = 0
    dqdb = 0
    do i = 1, size(this%xg)
      x = 0.5 * (a + b) + 0.5 * (b - a) * this%xg(i)
      w =                 0.5 * (b - a) * this%wg(i)

      call intg(x, this%p, f, dfdx, dfdp)
      q    = q + w * f
      dqda = dqda + 0.5 * (- this%wg(i) * f + w * dfdx * (1 - this%xg(i)))
      dqdb = dqdb + 0.5 * (  this%wg(i) * f + w * dfdx * (1 + this%xg(i)))
    end do
  end subroutine

end module
