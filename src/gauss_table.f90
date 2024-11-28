m4_include(util/macro.f90.inc)

module gauss_table_m

  use bin_search_m, only: bin_search, BS_LESS
  use error_m,      only: assert_failed, program_error
  use gauss_m,      only: gauss, gauss_legendre, gauss_laguerre
  use math_m,       only: linspace, logspace, expm1, ber, dberdx
  use mpfr_m,       only: add, div, exp, fma, mpfr, mpfr_cleanup, mul, neg, sub
  use omp_lib,      only: omp_get_thread_num, omp_get_num_threads
  use qsort_m,      only: qsort
  use vector_m,     only: vector_int, vector_real

  implicit none

  type gauss_table
    !! gauss quadrature for int_0^1 exp(a*x) f(x) dx

    integer :: n
      !! number of quadrature nodes
    logical :: log_a
      !! use logarithmically spaced a

    real, allocatable :: a(:)
    real, allocatable :: x(:,:)
    real, allocatable :: dx(:,:)
    real, allocatable :: w(:,:)
    real, allocatable :: dw(:,:)

    real, allocatable :: xleg(:)
    real, allocatable :: wleg(:)
    real, allocatable :: dxleg(:)
    real, allocatable :: dwleg(:)
    real, allocatable :: xlag(:)
    real, allocatable :: wlag(:)
  contains
    procedure :: init     => gauss_table_init
    procedure :: init_leg => gauss_table_init_leg
    procedure :: gen      => gauss_table_gen
    procedure :: get      => gauss_table_get
    procedure :: interp   => gauss_table_interp
    procedure :: output   => gauss_table_output
    procedure :: load     => gauss_table_load
    procedure :: save     => gauss_table_save
  end type

  type mpfr_tmp
    type(mpfr), allocatable :: fact(:)
    type(mpfr), allocatable :: ak(:)
    type(mpfr)              :: e
    type(mpfr)              :: t, u, v
  contains
    procedure :: init     => mpfr_tmp_init
    procedure :: destruct => mpfr_tmp_destruct
  end type

  integer, parameter :: NTAYLOR = 10

contains

  subroutine gauss_table_init(this, n, log_a)
    !! initialize gauss quadrature lookup table
    class(gauss_table), intent(out) :: this
    integer,            intent(in)  :: n
      !! number of gauss quadrature nodes
    logical, optional,  intent(in)  :: log_a
      !! use logarithmically spaced a (default: false)

    integer, parameter :: N0   = 1001
    real,    parameter :: RTOL = 1e-14, ATOL = 1e-16

    integer              :: i, i1, i2, i3, j1, j2, j3, k1, k2, k3, ithread, nthreads, ntot
    integer, allocatable :: nt(:), perm(:)
    logical              :: status
    real,    allocatable :: a(:), x(:,:), x1(:), dx(:,:), dx1(:), w(:,:), w1(:), dw(:,:), dw1(:)
    type(gauss)          :: gs
    type(mpfr_tmp)       :: tmp
    type(vector_real)    :: va, vx, vdx, vw, vdw
    type(vector_int)     :: vi

    this%n = n
    this%log_a = .false.
    if (present(log_a)) this%log_a = log_a

    ! try to load from file
    call this%load("/tmp", status)
    if (status) return

    ! initial coarse grid
    allocate (a(N0), x(n,N0), x1(n), dx(n,N0), dx1(n), w(n,N0), w1(n), dw(n,N0), dw1(n))
    if (this%log_a) then
      a = logspace(1e-8, 1e2, N0)
    else
      a = linspace(0.0, 100.0, N0)
    end if

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,j1,j2,j3,k1,k2,k3,ithread,x1,dx1,w1,dw1,gs,tmp,va,vx,vdx,vw,vdw,vi) &
    !$omp shared(this,n,nthreads,ntot,nt,a,x,dx,w,dw)

    nthreads = omp_get_num_threads()
    ithread  = omp_get_thread_num() + 1

    ! init mpfr
    call gs%init(n, NP = 1)
    call tmp%init(n)

    !$omp single
    ! init Gauss-Legendre (not needed for linearly spaced a)
    if (this%log_a) call this%init_leg(gs, n)

    ! init Gauss-Laguerre
    allocate (this%xlag(n), this%wlag(n))
    call gauss_laguerre(this%xlag, this%wlag)

    ! number of table entries per thread
    allocate (nt(0:nthreads + 1), source = 0)
    !$omp end single

    ! generate entries for coarse grid
    !$omp do schedule(dynamic)
    do i = 1, N0
      call this%gen(gs, tmp, a(i), x(:,i), dx(:,i), w(:,i), dw(:,i))
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    call va%init( N0,     c = 4 * N0,     x = a)
    call vx%init( N0 * n, c = 4 * N0 * n, x = reshape( x, [N0 * n]))
    call vdx%init(N0 * n, c = 4 * N0 * n, x = reshape(dx, [N0 * n]))
    call vw%init( N0 * n, c = 4 * N0 * n, x = reshape( w, [N0 * n]))
    call vdw%init(N0 * n, c = 4 * N0 * n, x = reshape(dw, [N0 * n]))

    ! interval stack
    call vi%init(0, c = 16)

    ! refinement
    !$omp do schedule(dynamic)
    do i = 1, N0 - 1
      ! add interval
      call vi%push(i)
      call vi%push(i+1)

      do while (vi%n > 0)
        ! get interval from stack
        i1   = vi%d(vi%n-1)
        i2   = vi%d(vi%n  )
        vi%n = vi%n - 2

        ! midpoint
        if (this%log_a) then
          call va%push(exp(0.5 * (log(va%d(i1) * va%d(i2)))))
        else
          call va%push(0.5 * (va%d(i1) + va%d(i2)))
        end if
        i3 = va%n

        ! data indices
        j1 = (i1 - 1) * n + 1
        j2 = (i2 - 1) * n + 1
        j3 = (i3 - 1) * n + 1
        k1 = i1 * n
        k2 = i2 * n
        k3 = i3 * n

        ! make room for new point
        call vx%resize( k3)
        call vdx%resize(k3)
        call vw%resize( k3)
        call vdw%resize(k3)

        ! generate data for new point and additionally interpolate it from existing data
        call this%gen(gs, tmp, va%d(i3), vx%d(j3:k3), vdx%d(j3:k3), vw%d(j3:k3), vdw%d(j3:k3))
        call this%interp(va%d(i1), vx%d(j1:k1), vdx%d(j1:k1), va%d(i2), vx%d(j2:k2), vdx%d(j2:k2), va%d(i3), x1, dx1)
        call this%interp(va%d(i1), vw%d(j1:k1), vdw%d(j1:k1), va%d(i2), vw%d(j2:k2), vdw%d(j2:k2), va%d(i3), w1, dw1)

        ! check if refinement is necessary
        if (any(abs(x1 - vx%d(j3:k3)) / (abs(vx%d(j3:k3)) + ATOL) > RTOL) .or. &
            any(abs(w1 - vw%d(j3:k3)) / (    vw%d(j3:k3)  + ATOL) > RTOL)) then
          ! add 2 new intervals to stack
          call vi%push(i1)
          call vi%push(i3)
          call vi%push(i3)
          call vi%push(i2)
        end if
      end do

      print "(I4,A,I4)", i, " / ", 1000
    end do
    !$omp end do

    ! save number of points generated by this thread
    nt(ithread+1) = va%n - N0

    !$omp barrier

    !$omp single

    ! count elements
    nt(0) = 1
    nt(1) = N0
    do i = 1, nthreads + 1
      nt(i) = nt(i) + nt(i-1)
    end do
    ntot = nt(nthreads + 1) - 1

    ! allocate global memory + fill in coarse grid
    allocate (this%a(ntot), this%x(n,ntot), this%dx(n,ntot), this%w(n,ntot), this%dw(n,ntot))
    this%a(   1:N0) = a
    this%x( :,1:N0) = x
    this%dx(:,1:N0) = dx
    this%w( :,1:N0) = w
    this%dw(:,1:N0) = dw

    !$omp end single

    ! copy local values to global memory
    j1 = N0 * n + 1
    k1 = va%n * n
    this%a(   nt(ithread):nt(ithread+1)-1) = va%d(N0+1:va%n)
    this%x( :,nt(ithread):nt(ithread+1)-1) = reshape( vx%d(j1:k1), [n, va%n - N0])
    this%dx(:,nt(ithread):nt(ithread+1)-1) = reshape(vdx%d(j1:k1), [n, va%n - N0])
    this%w( :,nt(ithread):nt(ithread+1)-1) = reshape( vw%d(j1:k1), [n, va%n - N0])
    this%dw(:,nt(ithread):nt(ithread+1)-1) = reshape(vdw%d(j1:k1), [n, va%n - N0])

    ! cleanup thread-local memory
    call va%destruct()
    call vx%destruct()
    call vdx%destruct()
    call vw%destruct()
    call vdw%destruct()
    call gs%destruct()
    call tmp%destruct()

    !$omp end parallel

    ! sort by a
    allocate (perm(size(this%a)))
    call qsort(this%a, perm = perm)
    this%x  = this%x(:, perm)
    this%dx = this%dx(:,perm)
    this%w  = this%w(:, perm)
    this%dw = this%dw(:,perm)

    ! save to file
    call this%save("/tmp")
  end subroutine

  subroutine gauss_table_init_leg(this, gs, n)
    !! initialize Gauss-Legendre + derivatives
    class(gauss_table), intent(inout) :: this
    type(gauss),        intent(inout) :: gs
      !! gauss quadrature generator
    integer,            intent(in)    :: n
      !! number of quadrature nodes

    real :: s, dsdp(1), dx(n,1), dw(n,1)

    allocate (this%xleg(n), this%dxleg(n), this%wleg(n), this%dwleg(n))
    call gs%generate(moms, [0.0], this%xleg, this%wleg, s, dx, dw, dsdp)
    this%dxleg = dx(:,1)
    this%dwleg = dw(:,1)

  contains

    subroutine moms(p, s, dsdp)
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p

      integer :: k

      m4_ignore(p)

      ! s(k) = 1 / (k + 1); ds(k) = 1 / (k + 2)
      call s(0)%set(1)
      do k = 1, 2 * n - 1
        call s(k)%set(k + 1)
        call div(s(k), 1, s(k))
        call dsdp(k-1,1)%set(s(k))
      end do
      call dsdp(2*n-1,1)%set(2 * n + 1)
      call div(dsdp(2*n-1,1), 1, dsdp(2*n-1,1))
    end subroutine

  end subroutine

  subroutine gauss_table_gen(this, gs, tmp, a, x, dx, w, dw)
    !! generate entry
    class(gauss_table), intent(in)    :: this
    type(gauss),        intent(inout) :: gs
      !! gauss quadrature generator
    type(mpfr_tmp),     intent(inout) :: tmp
      !! temporary mpfr values
    real,               intent(in)    :: a
      !! weight function parameter
    real,               intent(out)   :: x(:)
      !! output gauss node positions
    real,               intent(out)   :: dx(:)
      !! output derivatives of x wrt a
    real,               intent(out)   :: w(:)
      !! output gauss weights
    real,               intent(out)   :: dw(:)
      !! output derivatives of w wrt a

    integer :: k, l
    real    :: dx_(size(dx),1), dw_(size(dw),1), s, ds(1)

    m4_ignore(this)

    ! a**k
    call tmp%ak(0)%init_set(1)
    do k = 1, ubound(tmp%ak, 1)
      call mul(tmp%ak(k), tmp%ak(k-1), a)
    end do

    ! exp(a)
    call tmp%e%set(a)
    call exp(tmp%e, tmp%e)

    call gs%generate(moms, [a], x, w, s, dx_, dw_, ds)
    dx = dx_(:,1)
    dw = dw_(:,1)

  contains

    subroutine moms(p, s, dsdp)
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p

      m4_ignore(p) ! p(1) = a

      if (abs(a) > 1e-2) then
        do k = 0, gs%m
          call tmp%t%set(0)
          call tmp%u%set(0)
          do l = 0, k
            call div(tmp%v, tmp%fact(k), tmp%fact(l))
            if (mod(l, 2) == 1) call neg(tmp%v, tmp%v)
            call fma(tmp%t, tmp%v, tmp%ak(l), tmp%t)
            if (l > 0) then
              call mul(tmp%v, tmp%v, l)
              call fma(tmp%u, tmp%v, tmp%ak(l-1), tmp%u)
            end if
          end do
          call mul(tmp%t, tmp%t, tmp%e)
          call mul(tmp%u, tmp%u, tmp%e)
          call sub(s(k), tmp%t, tmp%fact(k))
          call div(s(k), s(k), tmp%ak(k+1))
          if (mod(k, 2) == 1) call neg(s(k), s(k))
          call add(dsdp(k,1), tmp%t, tmp%u)
          call div(dsdp(k,1), dsdp(k,1), tmp%ak(k+1))
          if (mod(k, 2) == 1) call neg(dsdp(k,1), dsdp(k,1))
          call mul(tmp%t, s(k), k + 1)
          call div(tmp%t, tmp%t, tmp%ak(1))
          call sub(dsdp(k,1), dsdp(k,1), tmp%t)
        end do
      else
        ! taylor: s(k) = sum_{l=0}^NTAYLOR a^l / ((l + 1)! + k * l!)
        do k = 0, gs%m
          call s(k)%set(0)
          call dsdp(k,1)%set(0)
          do l = 0, NTAYLOR
            call mul(tmp%t, tmp%fact(l), k)
            call add(tmp%t, tmp%t, tmp%fact(l+1))
            call div(tmp%t, 1.0, tmp%t)
            call fma(s(k), tmp%ak(l), tmp%t, s(k))
            if (l > 0) then
              call mul(tmp%u, tmp%t, l)
              call fma(dsdp(k,1), tmp%ak(l-1), tmp%u, dsdp(k,1))
            end if
          end do
        end do
      end if
    end subroutine

  end subroutine

  subroutine gauss_table_get(this, a, x, dx, w, dw, scale)
    class(gauss_table), intent(in)    :: this
    real,               intent(in)    :: a
      !! weight function parameter
    real,               intent(out)   :: x(:)
      !! output gauss node positions
    real,               intent(out)   :: dx(:)
      !! output derivatives of x wrt a
    real,               intent(out)   :: w(:)
      !! output gauss weigmts
    real,               intent(out)   :: dw(:)
      !! output derivatives of w wrt a
    logical, optional,  intent(in)    :: scale
      !! scale weights? (default: true)

    integer :: i
    logical :: scale_
    real    :: abs_a, B, dB, e(size(x)), s(size(x)), ds(size(x))

    abs_a = abs(a)

    if (abs_a < this%a(1)) then
      ! Gauss-Legendre
      x  = this%xleg + this%dxleg * abs_a
      dx = this%dxleg
      w  = this%wleg + this%dwleg * abs_a
      dw = this%dwleg
    elseif (abs_a > this%a(size(this%a))) then
      ! Gauss-Laguerrre
      x  = 1 - this%xlag(size(x):1:-1) / abs_a
      dx =     this%xlag(size(x):1:-1) / abs_a**2
      w  =   - this%wlag(size(x):1:-1) / expm1(-abs_a)
      dw =   - this%wlag(size(x):1:-1) * exp(-abs_a) / expm1(-abs_a)**2
    else
      ! interpolate from table
      i = bin_search(this%a, abs_a, mode = BS_LESS)
      if (i == ubound(this%a,1)) i = i - 1
      call this%interp(this%a(i), this%x(:,i), this%dx(:,i), this%a(i+1), this%x(:,i+1), this%dx(:,i+1), abs_a, x, dx)
      call this%interp(this%a(i), this%w(:,i), this%dw(:,i), this%a(i+1), this%w(:,i+1), this%dw(:,i+1), abs_a, w, dw)
    end if

    ! use symmetry for negative a
    if (a < 0) then
      x  = 1 -  x(size(x):1:-1)
      dx =     dx(size(x):1:-1)
      w  =      w(size(w):1:-1)
      dw =   - dw(size(w):1:-1)
    end if

    ! weight scaling
    scale_ = .true.
    if (present(scale)) then
      scale_ = scale
    end if
    if (scale_) then
      if (a <= 1) then
        B  = ber(a)
        dB = dberdx(a)
        e  = exp(- a * x)
        s  = e / B
        ds = - e * (x + a * dx + dB / B) / B
      else
        ! avoid overflow for a >> 0
        B  = ber(-a)
        dB = - dberdx(-a)
        e  = exp(- a * (x - 1))
        s  = e / B
        ds = - e * (x - 1 + a * dx + dB / B) / B
      end if
      dw = dw * s + w * ds
      w  =  w * s
    end if
  end subroutine

  subroutine gauss_table_interp(this, a1, f1, df1, a2, f2, df2, a, f, df)
    !! Hermite interpolation
    class(gauss_table), intent(in)  :: this
    real,               intent(in)  :: a1
      !! left point
    real,               intent(in)  :: f1(:)
      !! left function values
    real,               intent(in)  :: df1(:)
      !! left derivatives
    real,               intent(in)  :: a2
      !! right point
    real,               intent(in)  :: f2(:)
      !! right function values
    real,               intent(in)  :: df2(:)
      !! right derivatives
    real,               intent(in)  :: a
      !! interpolation point
    real,               intent(out) :: f(:)
      !! output interpolated function values at a
    real,               intent(out) :: df(:)
      !! output interpolated derivatives at a

    real :: e, e1, e2, de, da, t
    real :: h00, h10, h01, h11, g00, g10, g01, g11

    m4_ignore(this)

    if (this%log_a) then
      e  = log(a)
      e1 = log(a1)
      e2 = log(a2)
      de = e2 - e1
      t  = (e - e1) / de
    else
      da = a2 - a1
      t  = (a - a1) / da
    end if

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    g00 = 6 * t * (t - 1)
    g10 = (3 * t - 1) * (t - 1)
    g01 = - 6 * t * (t - 1)
    g11 = t * (3 * t - 2)

    if (this%log_a) then
      f  =  h00 * f1 + h10 * de * a1 * df1 + h01 * f2 + h11 * de * a2 * df2
      df = (g00 * f1 + g10 * de * a1 * df1 + g01 * f2 + g11 * de * a2 * df2) / (de * a)
    else
      f  =  h00 * f1 + h10 * da * df1 + h01 * f2 + h11 * da * df2
      df = (g00 * f1 + g10 * da * df1 + g01 * f2 + g11 * da * df2) / da
    end if
  end subroutine

  subroutine gauss_table_output(this, fname_x, fname_w, amax, nsamples)
    !! output gauss table to file
    class(gauss_table), intent(in) :: this
    character(*),       intent(in) :: fname_x
    character(*),       intent(in) :: fname_w
    real,               intent(in) :: amax
    integer,            intent(in) :: nsamples

    character(32)     :: fmt
    integer           :: i, funit_x, funit_w
    real, allocatable :: a(:), x(:), dx(:), w(:), dw(:)

    allocate (a(nsamples), x(this%n), dx(this%n), w(this%n), dw(this%n))

    a = linspace(-amax, amax, nsamples)

    write (fmt, "(A,I0,A)") "(", 1 + this%n, "ES25.16E3)"

    open (newunit = funit_x, file = fname_x, status = "replace", action = "write")
    open (newunit = funit_w, file = fname_w, status = "replace", action = "write")
    write (funit_x, "(A)", advance = "no") "a"
    do i = 1, this%n
      write (funit_x, "(A,I0)", advance = "no") " x", i
    end do
    write (funit_x, *)
    write (funit_w, "(A)", advance = "no") "a"
    do i = 1, this%n
      write (funit_w, "(A,I0)", advance = "no") " w", i
    end do
    write (funit_w, *)
    do i = 1, size(a)
      call this%get(a(i), x, dx, w, dw, scale = .false.)
      write (funit_x, fmt) a(i), x
      write (funit_w, fmt) a(i), w
    end do
    close (funit_x)
    close (funit_w)
  end subroutine

  subroutine gauss_table_load(this, dir, status)
    !! load gauss table from file
    class(gauss_table), intent(inout) :: this
    character(*),       intent(in)    :: dir
      !! directory
    logical,            intent(out)   :: status
      !! success (true) or fail (false)

    character(256) :: fname
    integer        :: funit, ntot

    ! filename
    write (fname, "(2A,I0,L1,A)") dir, "/gtab_", this%n, this%log_a, ".bin"

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")

    read (funit) ntot
    allocate (this%a(ntot), this%x(this%n,ntot), this%dx(this%n,ntot), this%w(this%n,ntot), this%dw(this%n,ntot), this%xleg(this%n), this%wleg(this%n), this%dxleg(this%n), this%dwleg(this%n), this%xlag(this%n), this%wlag(this%n))
    read (funit) this%a
    read (funit) this%x
    read (funit) this%dx
    read (funit) this%w
    read (funit) this%dw
    if (this%log_a) then
      read (funit) this%xleg
      read (funit) this%wleg
      read (funit) this%dxleg
      read (funit) this%dwleg
    end if
    read (funit) this%xlag
    read (funit) this%wlag

    close (funit)
    status = .true.
  end subroutine

  subroutine gauss_table_save(this, dir)
    !! save gauss table to file
    class(gauss_table), intent(in) :: this
    character(*),       intent(in) :: dir
      !! directory

    character(256) :: fname
    integer        :: funit, ntot

    ! filename
    write (fname, "(2A,I0,L1,A)") dir, "/gtab_", this%n, this%log_a, ".bin"

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")

    ntot = size(this%a)
    print *, "ntot = ", ntot

    write (funit) ntot
    write (funit) this%a
    write (funit) this%x
    write (funit) this%dx
    write (funit) this%w
    write (funit) this%dw
    if (this%log_a) then
      write (funit) this%xleg
      write (funit) this%wleg
      write (funit) this%dxleg
      write (funit) this%dwleg
    end if
    write (funit) this%xlag
    write (funit) this%wlag

    close (funit)
  end subroutine

  subroutine mpfr_tmp_init(this, n)
    !! initialize temporary mpfr values
    class(mpfr_tmp), intent(out) :: this
    integer,         intent(in)  :: n
      !! number of gauss points

    integer :: k, m

    m = 2 * n - 1

    ! factorial
    allocate (this%fact(0:max(m,NTAYLOR+1)))
    call this%fact(0)%init_set(1)
    do k = 1, ubound(this%fact, 1)
      call this%fact(k)%init_set(this%fact(k-1))
      call mul(this%fact(k), this%fact(k), k)
    end do

    ! a**k (not initialized)
    allocate (this%ak(0:max(m+1,NTAYLOR)))
    do k = 0, ubound(this%ak, 1)
      call this%ak(k)%init()
    end do

    ! exp(a) (not initialize)
    call this%e%init()

    ! temp
    call this%t%init()
    call this%u%init()
    call this%v%init()
  end subroutine

  subroutine mpfr_tmp_destruct(this)
    !! free temporary memory
    class(mpfr_tmp), intent(inout) :: this

    integer :: k

    do k = lbound(this%fact, 1), ubound(this%fact, 1)
      call this%fact(k)%destruct()
    end do
    do k = lbound(this%ak, 1), ubound(this%ak, 1)
      call this%ak(k)%destruct()
    end do
    call this%e%destruct()
    call this%t%destruct()
    call this%u%destruct()
    call this%v%destruct()
    call mpfr_cleanup()
  end subroutine

end module
