m4_include(util/macro.f90.inc)

module lookup_table_m

  use, intrinsic :: iso_fortran_env, only: real128
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_POSITIVE_INF

  use bin_search_m,    only: bin_search, BS_LESS
  use error_m,         only: assert_failed, program_error
  use math_m,          only: linspace
  use omp_lib,         only: omp_get_thread_num, omp_get_num_threads
  use quad_m,          only: quad
  use qsort_m,         only: qsort
  use util_m,          only: hash
  use vector_m,        only: vector_int, vector_log, vector_real

  implicit none

  private
  public :: lookup_table

  type lookup_table
    !! lookup table for a function defined by an integral
    !! g(x) = integral_{a(x)}^{b(x)} f(t,x) dt

    character(:), allocatable :: name
      !! name for loading and saving

    real,    allocatable :: x(:)
      !! sample points
    real,    allocatable :: g(:)
      !! function values at sample points
    real,    allocatable :: dg(:)
      !! function derivatives at sample points
    logical, allocatable :: lg(:)
  contains
    procedure :: init   => lookup_table_init
    procedure :: gen    => lookup_table_gen
    generic   :: get    => lookup_table_get_D, &
      &                    lookup_table_get_Q
    procedure :: inv    => lookup_table_inv
    procedure :: output => lookup_table_output
    procedure :: load   => lookup_table_load
    procedure :: save   => lookup_table_save

    procedure, private :: lookup_table_get_D, lookup_table_get_Q
  end type

  interface
    subroutine bound(x, b, dbdx)
      !! integration boundary
      import real128
      real(real128), intent(in)  :: x
        !! argument
      real(real128), intent(out) :: b
        !! output boundary
      real(real128), intent(out) :: dbdx
        !! output derivative of boundary
    end subroutine

    subroutine integrand(t, x, f, dfdt, dfdx)
      !! integrand
      import real128
      real(real128), intent(in)  :: t
        !! inner argument
      real(real128), intent(in)  :: x(:)
        !! argument
      real(real128), intent(out) :: f
        !! output integrand
      real(real128), intent(out) :: dfdt
        !! output derivative of integrand wrt t
      real(real128), intent(out) :: dfdx(:)
        !! output derivative of integrand wrt x
    end subroutine
  end interface

contains

  subroutine lookup_table_init(this, dir, name, xmin, xmax, a, b, f)
    !! initialize lookup table
    class(lookup_table), intent(out) :: this
    character(*),        intent(in)  :: dir
      !! directory
    character(*),        intent(in)  :: name
      !! name for loading and saving
    real,                intent(in)  :: xmin
      !! minimal supported x
    real,                intent(in)  :: xmax
      !! maximal supported x
    procedure(bound)                 :: a
      !! lower integration boundary
    procedure(bound)                 :: b
      !! upper integration boundary
    procedure(integrand)             :: f
      !! integrand

    integer, parameter :: N0 = 1024 + 1
    real,    parameter :: RTOL = 1e-13, ATOL = 1e-16, XTOL = 1e-8

    integer              :: i, i1, i2, i3, ithread, nthreads, n
    integer, allocatable :: ntab(:), perm(:)
    logical              :: status, lg1
    logical, allocatable :: lg(:)
    real                 :: x1, x2, x3, dx, g1, g2, e1, e2, sgn
    real,    allocatable :: x(:), g(:), dg(:)
    type(vector_int)     :: stack
    type(vector_log)     :: vlg
    type(vector_real)    :: vx, vg, vdg

    this%name = name

    ! try to load from file
    call this%load(dir, xmin, xmax, status)
    if (status) return

    ! init coarse grid
    allocate (x(N0), g(N0), dg(N0), source = 0.0)
    allocate (lg(N0), source = .false.)
    x = linspace(xmin, xmax, N0)

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,ithread,vx,vg,vdg,vlg,stack,x1,x2,x3,dx,lg1,g1,g2,e1,e2,sgn) &
    !$omp shared(this,nthreads,n,ntab,g,dg,lg,x,xmin,xmax)

    ithread = omp_get_thread_num() + 1

    ! allocate number of table entries per thread
    !$omp single
    nthreads = omp_get_num_threads()
    allocate (ntab(0:nthreads + 1), source = 0)
    !$omp end single

    ! generate entries for coarse grid
    !$omp do schedule(dynamic)
    do i = 1, N0
      call this%gen(a, b, f, x(i), g(i), dg(i))
      print "(I4,A,I4)", i, " / ", N0
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    call vx%init( N0, c = 4 * N0, x =     x)
    call vg%init( N0, c = 4 * N0, x =  g)
    call vdg%init(N0, c = 4 * N0, x = dg)
    call vlg%init(N0, c = 4 * N0, x =    lg)

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

        x1 = vx%d(i1)
        x2 = vx%d(i2)
        x3 = 0.5 * (x1 + x2)
        dx = x2 - x1

        ! save midpoint
        call vx%push(x3)
        i3 = vx%n

        ! make room for new point
        call vg%resize(i3)
        call vdg%resize(i3)

        ! generate data for new point
        call this%gen(a, b, f, x3, vg%d(i3), vdg%d(i3))

        ! direct interpolation from existing data
        g1 = 0.5 * (vg%d(i1) + vg%d(i2)) + 0.125 * dx * (vdg%d(i1) - vdg%d(i2))
        e1 = abs(g1 - vg%d(i3)) / (abs(vg%d(i3)) + ATOL)

        ! logarithmic interpolation from existing data
        lg1 = .false.
        if (vg%d(i1) * vg%d(i2) > 0) then
          sgn = sign(1.0, vg%d(i1))

          ! logarithmic interpolation
          g2 = sgn * sqrt(vg%d(i1) * vg%d(i2)) * exp(0.125 * dx * (vdg%d(i1) / vg%d(i1) - vdg%d(i2) / vg%d(i2)))

          ! decide if direct or logarithmic interpolation is better
          e2 = abs(g2 - vg%d(i3)) / (abs(vg%d(i3)) + ATOL)
          if (e2 < e1) lg1 = .true.
        end if

        ! lg for interval i3 --- i2
        call vlg%push(lg1)

        ! check if refinement is necessary
        if ((merge(e2, e1, lg1) > RTOL) .and. (dx > XTOL)) then
          ! add 2 new intervals to stack
          call stack%push([i1, i3])
          call stack%push([i3, i2])
        else
          ! lg for interval i1 -- i3
          vlg%d(i1) = lg1
        end if
      end do

      ! store lg for coarse grid
      lg(i) = vlg%d(i)

      ! print interval number
      print "(I4,A,I4)", i, " / ", N0 - 1
    end do
    !$omp end do nowait

    ! save number of points generated by this thread
    ntab(ithread + 1) = vx%n - N0

    !$omp barrier

    !$omp single

    ! count elements
    ntab(0) = 1
    ntab(1) = N0
    do i = 1, nthreads + 1
      ntab(i) = ntab(i) + ntab(i-1)
    end do
    n = ntab(nthreads + 1) - 1

    ! allocate global memory + fill in coarse grid
    allocate (this%x(n), this%g(n), this%dg(n), this%lg(n))
    this%x(1:N0)  = x
    this%g(1:N0)  = g
    this%dg(1:N0) = dg
    this%lg(1:N0) = lg

    !$omp end single

    ! copy local values to global memory
    this%x( ntab(ithread):ntab(ithread+1)-1) = vx%d( N0+1:vx%n )
    this%g( ntab(ithread):ntab(ithread+1)-1) = vg%d( N0+1:vg%n )
    this%dg(ntab(ithread):ntab(ithread+1)-1) = vdg%d(N0+1:vdg%n)
    this%lg(ntab(ithread):ntab(ithread+1)-1) = vlg%d(N0+1:vlg%n)

    ! cleanup thread-local memory
    call stack%destruct()
    call vlg%destruct()
    call vx%destruct()
    call vg%destruct()
    call vdg%destruct()

    !$omp end parallel

    ! sort by x
    allocate (perm(size(this%x)))
    call qsort(this%x, perm = perm)
    this%g  = this%g( perm)
    this%dg = this%dg(perm)
    this%lg = this%lg(perm)

    ! save to file
    call this%save(dir)
  end subroutine

  subroutine lookup_table_gen(this, a, b, f, x, g, dg)
    !! generate single entry
    class(lookup_table), intent(in)  :: this
    procedure(bound)                 :: a
      !! lower integration boundary
    procedure(bound)                 :: b
      !! upper integration boundary
    procedure(integrand)             :: f
      !! integrand
    real,                intent(in)  :: x
      !! argument
    real,                intent(out) :: g
      !! value of the integral
    real,                intent(out) :: dg
      !! derivative of the integral

    real(real128), parameter :: RTOL = 1e-14_16

    real(real128) :: a16, da16, b16, db16, p(1), g16, dgda16, dgdb16, dgdp16(1)

    m4_ignore(this)

    p(1) = real(x, kind = real128)

    ! integration bounds
    call a(p(1), a16, da16)
    call b(p(1), b16, db16)

    ! integral
    call quad(f, a16, b16, p, g16, dgda16, dgdb16, dgdp16, rtol = RTOL, max_levels = 20)

    g  = real(g16)
    dg = real(dgda16 * da16 + dgdb16 * db16 + dgdp16(1))
  end subroutine

  subroutine lookup_table_get_D(this, x, g, dgdx)
    !! get function value and derivative in double precision
    class(lookup_table), intent(in)  :: this
    real,                intent(in)  :: x
      !! argument
    real,                intent(out) :: g
      !! output function value
    real,                intent(out) :: dgdx
      !! output derivative of g wrt x

    integer :: i
    real    :: x1, x2, dx, g1, dg1, g2, dg2, t, h00, h10, h01, h11, dh00, dh10, dh01, dh11

    if ((x < this%x(1)) .or. (x > this%x(size(this%x)))) then
      print "(3ES25.16E3)", x, this%x(1), this%x(size(this%x))
      call program_error("x out of bounds")
    end if

    ! find interval
    i = bin_search(this%x, x, mode = BS_LESS)
    if (i == ubound(this%x, 1)) i = i - 1

    x1 = this%x(i)
    x2 = this%x(i+1)
    dx = x2 - x1
    t  = (x - x1) / dx

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    dh00 = 6 * t * (t - 1)
    dh10 = (3 * t - 1) * (t - 1)
    dh01 = - 6 * t * (t - 1)
    dh11 = t * (3 * t - 2)

    g1  = this%g( i  )
    dg1 = this%dg(i  )
    g2  = this%g( i+1)
    dg2 = this%dg(i+1)

    if (this%lg(i)) then
      ! logarithmic interpolation
      g    =   h00 * log(abs(g1)) +  h10 * dx * dg1 / g1 +  h01 * log(abs(g2)) +  h11 * dx * dg2 / g2
      dgdx = (dh00 * log(abs(g1)) + dh10 * dx * dg1 / g1 + dh01 * log(abs(g2)) + dh11 * dx * dg2 / g2) / dx

      g    = sign(1.0, g1) * exp(g)
      dgdx = g * dgdx
    else
      ! direct interpolation
      g    =   h00 * g1 +  h10 * dx * dg1 +  h01 * g2 + h11 * dx *  dg2
      dgdx = (dh00 * g1 + dh10 * dx * dg1 + dh01 * g2 + dh11 * dx * dg2) / dx
    end if
  end subroutine

  subroutine lookup_table_get_Q(this, x, g, dgdx)
    !! get function value and derivative in quad precision
    class(lookup_table), intent(in)  :: this
    real(real128),       intent(in)  :: x
      !! argument
    real(real128),       intent(out) :: g
      !! output function value
    real(real128),       intent(out) :: dgdx
      !! output derivative of g wrt x

    integer       :: i
    real(real128) :: x1, x2, dx, g1, dg1, g2, dg2, t, h00, h10, h01, h11, dh00, dh10, dh01, dh11

    if ((x < this%x(1)) .or. (x > this%x(size(this%x)))) then
      print "(3ES25.16E3)", x, this%x(1), this%x(size(this%x))
      call program_error("x out of bounds")
    end if

    ! find interval
    i = bin_search(this%x, real(x), mode = BS_LESS)
    if (i == ubound(this%x, 1)) i = i - 1

    x1 = this%x(i)
    x2 = this%x(i+1)
    dx = x2 - x1
    t  = (x - x1) / dx

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    dh00 = 6 * t * (t - 1)
    dh10 = (3 * t - 1) * (t - 1)
    dh01 = - 6 * t * (t - 1)
    dh11 = t * (3 * t - 2)

    g1  = this%g( i  )
    dg1 = this%dg(i  )
    g2  = this%g( i+1)
    dg2 = this%dg(i+1)

    if (this%lg(i)) then
      ! logarithmic interpolation
      g    =   h00 * log(abs(g1)) +  h10 * dx * dg1 / g1 +  h01 * log(abs(g2)) +  h11 * dx * dg2 / g2
      dgdx = (dh00 * log(abs(g1)) + dh10 * dx * dg1 / g1 + dh01 * log(abs(g2)) + dh11 * dx * dg2 / g2) / dx

      g    = sign(1.0_16, g1) * exp(g)
      dgdx = g * dgdx
    else
      ! direct interpolation
      g    =   h00 * g1 +  h10 * dx * dg1 +  h01 * g2 + h11 * dx *  dg2
      dgdx = (dh00 * g1 + dh10 * dx * dg1 + dh01 * g2 + dh11 * dx * dg2) / dx
    end if
  end subroutine

  subroutine lookup_table_inv(this, g, x, dxdg)
    !! get inverse of function (must be monotone)
    class(lookup_table), intent(in)  :: this
    real,                intent(in)  :: g
      !! function value
    real,                intent(out) :: x
      !! output argument
    real,                intent(out) :: dxdg
      !! output derivative of x wrt g

    real,    parameter :: ATOL   = 5e-13
    integer, parameter :: MAX_IT = 10

    integer :: i, it
    logical :: lg
    real    :: x1, x2, dx, g1, dg1, g2, dg2, t, tmin, tmax, told, res, dresdt, dt, err

    if ((g < this%g(1)) .or. (g > this%g(size(this%x)))) then
      print "(A,ES25.16E3)", "g    = ", g
      print "(A,ES25.16E3)", "gmin = ", this%g(1)
      print "(A,ES25.16E3)", "gmax = ", this%g(size(this%x))
      call program_error("g is out of range")
    end if

    ! find interval
    i = bin_search(this%g, g, mode = BS_LESS)
    if (i == ubound(this%g, 1)) i = i - 1

    ! load interval values
    x1  = this%x(i)
    x2  = this%x(i+1)
    dx  = x2 - x1
    g1  = this%g( i  )
    dg1 = this%dg(i  )
    g2  = this%g( i+1)
    dg2 = this%dg(i+1)
    lg  = this%lg(i)

    ! start estimate
    if (lg) then
      t = log(g / g1) / log(g2 / g1)
    else
      t = (g - g1) / (g2 - g1)
    end if

    ! search bounds
    tmin = 0
    tmax = 1

    ! Newton iteration to get t
    err = huge(1.0)
    it  = 0
    do while (err > ATOL)
      it = it + 1
      if (it > MAX_IT) then
        print "(A,ES25.16E3)", "g = ", g
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
    dxdg = dx / dresdt

  contains

    subroutine residual(t, res, dresdt)
      real, intent(in)  :: t
      real, intent(out) :: res
      real, intent(out) :: dresdt

      real :: h00, h10, h01, h11, dh00, dh10, dh01, dh11

      h00 = (1 + 2 * t) * (1 - t)**2
      h10 = t * (1 - t)**2
      h01 = t**2 * (3 - 2 * t)
      h11 = t**2 * (t - 1)

      dh00 = 6 * t * (t - 1)
      dh10 = (3 * t - 1) * (t - 1)
      dh01 = - 6 * t * (t - 1)
      dh11 = t * (3 * t - 2)

      if (lg) then
        ! logarithmic interpolation
        res    =  h00 * log(g1) +  h10 * dx * dg1 / g1 +  h01 * log(g2) +  h11 * dx * dg2 / g2
        dresdt = dh00 * log(g1) + dh10 * dx * dg1 / g1 + dh01 * log(g2) + dh11 * dx * dg2 / g2

        res    = exp(res)
        dresdt = res * dresdt
        res    = res - g
      else
        ! direct interpolation
        res    =  h00 * g1 +  h10 * dx * dg1 +  h01 * g2 +  h11 * dx * dg2 - g
        dresdt = dh00 * g1 + dh10 * dx * dg1 + dh01 * g2 + dh11 * dx * dg2
      end if
    end subroutine

  end subroutine

  subroutine lookup_table_output(this, fname, x1, x2, nsamples)
    !! output lookup table to file
    class(lookup_table), intent(in) :: this
    character(*),              intent(in) :: fname
      !! filename
    real,                      intent(in) :: x1
      !! lower bound
    real,                      intent(in) :: x2
      !! upper bound
    integer,                   intent(in) :: nsamples
      !! number of points

    integer           :: i, funit
    real              :: g, dg
    real, allocatable :: x(:)

    x = linspace(x1, x2, nsamples)

    open (newunit = funit, file = fname, status = "replace", action = "write")
    ! write (funit, "(A)") "x g dgdx"
    do i = 1, nsamples
      call this%get(x(i), g, dg)
      write (funit, "(3ES25.16E3)") x(i), g, dg
    end do
    close (funit)
  end subroutine

  subroutine lookup_table_load(this, dir, xmin, xmax, status)
    !! load lookup table from file
    class(lookup_table), intent(inout) :: this
    character(*),        intent(in)    :: dir
      !! directory
    real,                intent(in)    :: xmin
      !! minimal supported x
    real,                intent(in)    :: xmax
      !! maximal supported x
    logical,             intent(out)   :: status
      !! success (true) or fail (false)

    character(256) :: fname
    integer        :: h, funit, num_entries

    status = .false.

    ! filename
    write (fname, "(2ES25.16E3)") xmin, xmax
    h = hash(fname) ! hash parameters to get unique name
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    ! load
    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")
    read (funit) num_entries
    allocate (this%x(num_entries), this%g(num_entries), this%dg(num_entries), this%lg(num_entries))
    read (funit) this%x
    read (funit) this%g
    read (funit) this%dg
    read (funit) this%lg
    close (funit)
    status = .true.
  end subroutine

  subroutine lookup_table_save(this, dir)
    !! save lookup table to file
    class(lookup_table), intent(in) :: this
    character(*),        intent(in) :: dir
      !! directory

    character(256) :: fname
    integer        :: h, funit, num_entries

    ! filename
    num_entries = size(this%x)
    write (fname, "(2ES25.16E3)") this%x(1), this%x(num_entries)
    h = hash(fname)
    write (fname, "(4A,Z0.8,A)") dir, "/", this%name, "_", h, ".bin"

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) num_entries
    write (funit) this%x
    write (funit) this%g
    write (funit) this%dg
    write (funit) this%lg
    close (funit)
  end subroutine

end module
