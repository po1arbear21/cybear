m4_include(macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module gauss_m
  !! Calculate Gauss quadrature nodes and weights using Golub-Welsh algorithm
  !! requires library MPFR

  use blas95,   only: gemm
  use error_m,  only: assert_failed, program_error
  use lapack95, only: stev
  use mpfr_m,   only: add, div, fma, fmma, fmms, fms, gamma_mpfr, mpfr, mpfr_startup, mpfr_cleanup, mul, neg, sqr, sqrt_mpfr, sub

  implicit none

  private
  public gauss, gauss_legendre, gauss_laguerre, gauss_hermite

  type gauss
    !! gauss quadrature generator for custom weights

    integer                 :: n
      !! number of quadrature points
    integer                 :: m
      !! highest moment (M = 2 * N - 1)
    integer                 :: np
      !! number of parameters

    type(mpfr), allocatable :: a(:), dadp(:,:), b(:), dbdp(:,:), c(:), dcdp(:,:), d(:), dddp(:,:)
    type(mpfr), allocatable :: s(:), dsdp(:,:), sg(:,:), dsgdp(:,:,:)
    type(mpfr)              :: t
    type(mpfr)              :: u
    type(mpfr), allocatable :: v(:)
    type(mpfr)              :: x
  contains
    procedure :: init     => gauss_init
    procedure :: destruct => gauss_destruct
    procedure :: generate => gauss_generate
  end type

  interface
    subroutine moments(p, s, dsdp)
      !! moments for custom weight function
      import mpfr
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p
    end subroutine
  end interface

contains

  subroutine gauss_init(this, n, NP, prec)
    !! initialize memory
    class(gauss),      intent(out) :: this
    integer,           intent(in)  :: n
      !! number of quadrature points
    integer, optional, intent(in)  :: np
      !! number of parameters (default: 0)
    integer, optional, intent(in)  :: prec
      !! MPFR precision (default: 1024)

    integer :: i, ip, j, prec_

    this%n  = n
    this%m  = 2 * n - 1
    this%np = 0
    if (present(np)) this%np = np

    prec_ = 1024
    if (present(prec)) prec_ = prec
    call mpfr_startup(prec = prec_)

    ! allocate memory
    allocate (this%a(n), this%b(n-1), this%c(n), this%d(n-1), this%s(0:this%m), this%sg(n,n+1), this%v(n))
    allocate (this%dadp(n,this%np), this%dbdp(n-1,this%np), this%dcdp(n,this%np), this%dddp(n-1,this%np))
    allocate (this%dsdp(0:this%m,this%np), this%dsgdp(n,n+1,this%np))
    do i = 1, n
      call this%a(i)%init()
      call this%c(i)%init()
      call this%v(i)%init()
      if (i < n) then
        call this%b(i)%init()
        call this%d(i)%init()
      end if
    end do
    do i = 0, this%m
      call this%s(i)%init()
    end do
    do j = 1, n + 1
      do i = 1, n
        call this%sg(i,j)%init()
      end do
    end do
    call this%t%init()
    call this%u%init()
    call this%x%init()

    do ip = 1, this%np
      do i = 1, n
        call this%dadp(i,ip)%init()
        call this%dcdp(i,ip)%init()
        if (i < n) then
          call this%dbdp(i,ip)%init()
          call this%dddp(i,ip)%init()
        end if
      end do
      do i = 0, this%m
        call this%dsdp(i,ip)%init()
      end do
      do j = 1, n + 1
        do i = 1, n
          call this%dsgdp(i,j,ip)%init()
        end do
      end do
    end do
  end subroutine

  subroutine gauss_destruct(this)
    !! cleanup memory
    class(gauss), intent(inout) :: this

    integer :: i, j, ip

    do i = 1, this%n
      call this%a(i)%destruct()
      call this%c(i)%destruct()
      call this%v(i)%destruct()
      if (i < this%n) then
        call this%b(i)%destruct()
        call this%d(i)%destruct()
      end if
    end do
    do i = 0, this%m
      call this%s(i)%destruct()
    end do
    do j = 1, this%n + 1
      do i = 1, this%n
        call this%sg(i,j)%destruct()
      end do
    end do
    call this%t%destruct()
    call this%u%destruct()
    call this%x%destruct()
    do ip = 1, this%np
      do i = 1, this%n
        call this%dadp(i,ip)%destruct()
        call this%dcdp(i,ip)%destruct()
        if (i < this%n) then
          call this%dbdp(i,ip)%destruct()
          call this%dddp(i,ip)%destruct()
        end if
      end do
      do i = 0, this%m
        call this%dsdp(i,ip)%destruct()
      end do
      do j = 1, this%n + 1
        do i = 1, this%n
          call this%dsgdp(i,j,ip)%destruct()
        end do
      end do
    end do

    call mpfr_cleanup()
  end subroutine

  subroutine gauss_generate(this, momfun, p, x, w, s, dxdp, dwdp, dsdp)
    !! generate gauss nodes and weights
    class(gauss), intent(inout) :: this
    procedure(moments)          :: momfun
      !! calculate moments of weight function
    real,         intent(in)    :: p(:)
      !! parameters
    real,         intent(out)   :: x(:)
      !! output gauss quadrature nodes
    real,         intent(out)   :: w(:)
      !! output unscaled gauss quadrature weights (sum(w) == 1)
    real,         intent(out)   :: s
      !! output scaling factor for weights (0-th moment of weight function)
    real,         intent(out)   :: dxdp(:,:)
      !! output derivatives of x wrt p
    real,         intent(out)   :: dwdp(:,:)
      !! output derivatives of w wrt p
    real,         intent(out)   :: dsdp(:)
      !! output derivatives of s wrt p

    integer           :: i, ip, j, k, n, np, m
    real, allocatable :: b(:), v(:,:), dadp(:), dbdp(:), dAdpv(:,:), vdAdpv(:,:)

    n  = this%N
    np = this%NP
    m  = this%M

    ! check dimensions
    m4_assert(size(p) == np)
    m4_assert(size(x) == n)
    m4_assert(size(w) == n)
    m4_assert(size(dxdp, 1) == n)
    m4_assert(size(dxdp, 2) == np)
    m4_assert(size(dwdp, 1) == n)
    m4_assert(size(dwdp, 2) == np)
    m4_assert(size(dsdp, 1) == np)

    ! get moments
    if (.not. allocated(this%dsdp)) allocate (this%dsdp(0:n,1:np))
    call momfun(p, this%s, this%dsdp)

    ! extract 0-th moment
    s = this%s(0)%to_real()
    do ip = 1, np
      dsdp(ip) = this%dsdp(0,ip)%to_real()
    end do

    ! simple case n=1: x = s1/s0; w = 1
    if (n == 1) then
      call div(this%x, this%s(1), this%s(0))
      x(1) = this%x%to_real()
      w(1) = 1.0

      do ip = 1, np
        call mul(this%t, this%dsdp(0,ip), this%x)
        call sub(this%t, this%dsdp(1,ip), this%t)
        call div(this%t, this%t, this%s(0))
        dxdp(1,ip) = this%t%to_real()
        dwdp(1,ip) = 0
      end do

      return
    end if

    ! calculate recurrence relation
    do i = 1, n
      do j = 1, n + 1
        ! reset sg, dsgdp
        call this%sg(i,j)%set(0)
        do ip = 1, np
          call this%dsgdp(i,j,ip)%set(0)
        end do

        do k = 1, i - 1
          ! sg(i,j) += sg(k,i) * sg(k,j) / sg(k,k)
          call mul(this%t, this%sg(k,i), this%sg(k,j))
          call div(this%t, this%t, this%sg(k,k))
          call add(this%sg(i,j), this%sg(i,j), this%t)

          ! dsgdp(i,j,:) += (sg(k,j) * dsgdp(k,i,:) + sg(k,i) * dsgdp(k,j,:) - sg(k,i) * sg(k,j) / sg(k,k) * dsgdp(k,k,:)) / sg(k,k)
          do ip = 1, np
            call fmms(this%u, this%sg(k,i), this%dsgdp(k,j,ip), this%t, this%dsgdp(k,k,ip))
            call fma(this%u, this%sg(k,j), this%dsgdp(k,i,ip), this%u)
            call div(this%u, this%u, this%sg(k,k))
            call add(this%dsgdp(i,j,ip), this%dsgdp(i,j,ip), this%u)
          end do
        end do

        ! sg(i,j) = s(i+j-2) - sg(i,j)
        call sub(this%sg(i,j), this%s(i+j-2), this%sg(i,j))
        do ip = 1, np
          call sub(this%dsgdp(i,j,ip), this%dsdp(i+j-2,ip), this%dsgdp(i,j,ip))
        end do
      end do
    end do

    ! matrix diagonals
    do i = 1, n
      ! t = 1 / sg(i,i)
      call div(this%t, 1.0, this%sg(i,i))

      ! c(i) = sg(i,i+1) / sg(i,i)
      call mul(this%c(i), this%sg(i,i+1), this%t)

      ! dcdp(i,:) = (dsgdp(i,i+1,:) - c(i) * dsgdp(i,i,:)) / sg(i,i)
      do ip = 1, np
        call mul(this%dcdp(i,ip), this%c(i), this%dsgdp(i,i,ip))
        call sub(this%dcdp(i,ip), this%dsgdp(i,i+1,ip), this%dcdp(i,ip))
        call mul(this%dcdp(i,ip), this%dcdp(i,ip), this%t)
      end do

      ! a(i) = c(i) - c(i-1); c(0) = 0
      if (i == 1) then
        call this%a(i)%set(this%c(i))
        do ip = 1, np
          call this%dadp(i,ip)%set(this%dcdp(i,ip))
        end do
      else
        call sub(this%a(i), this%c(i), this%c(i-1))
        do ip = 1, np
          call sub(this%dadp(i,ip), this%dcdp(i,ip), this%dcdp(i-1,ip))
        end do
      end if

      if (i < n) then
        ! d(i) = sg(i+1,i+1) / sg(i,i)
        call mul(this%d(i), this%sg(i+1,i+1), this%t)

        ! dddp(i,:) = (dsgdp(i+1,i+1,:) - d(i) * dsgdp(i,i,:)) / sg(i,i)
        do ip = 1, np
          call mul(this%dddp(i,ip), this%d(i), this%dsgdp(i,i,ip))
          call sub(this%dddp(i,ip), this%dsgdp(i+1,i+1,ip), this%dddp(i,ip))
          call mul(this%dddp(i,ip), this%dddp(i,ip), this%t)
        end do

        ! b(i) = sqrt(d(i))
        call sqrt_mpfr(this%b(i), this%d(i))

        ! dbdp(i,:) = dddp(i,:) / (2 * sqrt(d(i)))
        if (np > 0) then
          call add(this%u, this%b(i), this%b(i))
          call div(this%u, 1.0, this%u) ! u = 1 / (2 * sqrt(d(i)))
          do ip = 1, np
            call mul(this%dbdp(i,ip), this%dddp(i,ip), this%u)
          end do
        end if
      end if
    end do

    ! solve approximate matrix using LAPACK
    allocate (b(n-1), v(n,n))
    do i = 1, n - 1
      x(i) = this%a(i)%to_real()
      b(i) = this%b(i)%to_real()
    end do
    x(n) = this%a(n)%to_real()
    call stev(x, b, v)

    ! refine using single Rayleigh quotient iteration step for each eigenvalue
    do j = 1, n
      ! load estimate
      call this%x%set(x(j))
      do i = 1, n
        call this%v(i)%set(v(i,j))
      end do

      ! forward substitution
      call sub(this%t, this%a(1), this%x)    ! t = a(1) - x
      call div(this%t, 1.0, this%t)          ! t = 1 / (a(1) - x)
      call mul(this%c(1), this%b(1), this%t) ! c(1) = b(1) / (a(1) - x)
      call mul(this%v(1), this%v(1), this%t) ! v(1) = v(1) / (a(1) - x)
      do i = 2, n - 1
        call sub(this%t, this%a(i), this%x)                   ! t = a(i) - x
        call fms(this%t, this%b(i-1), this%c(i-1), this%t)    ! t = b(i-1) * c(i-1) - (a(i) - x)
        call div(this%t, 1.0, this%t)                         ! t = 1 / (b(i-1) * c(i-1) - (a(i) - x))
        call mul(this%c(i), this%b(i), this%t)                ! c(i) = b(i) / (b(i-1) * c(i-1) - (a(i) - x))
        call neg(this%c(i), this%c(i))                        ! c(i) = b(i) / (a(i) - x - b(i-1) * c(i-1))
        call fms(this%u, this%b(i-1), this%v(i-1), this%v(i)) ! u = b(i-1) * v(i-1) - v(i)
        call mul(this%v(i), this%u, this%t)                   ! v(i) = (v(i) - b(i-1) * v(i-1)) / (a(i) - x - b(i-1) * c(i-1))
      end do
      call sub(this%t, this%a(n), this%x)                   ! t = a(n) - x
      call fms(this%t, this%b(n-1), this%c(n-1), this%t)    ! t = b(n-1) * c(n-1) - (a(n) - x)
      call fms(this%u, this%b(n-1), this%v(n-1), this%v(n)) ! u = b(n-1) * v(n-1) - v(n)
      call div(this%v(n), this%u, this%t)                   ! v(n) = (v(n) - b(n-1) * v(n-1)) / (a(n) - x - b(n-1) * c(n-1))

      ! back substitution: v(i) -= c(i) * v(i+1)
      do i = n-1, 1, -1
        call fms(this%v(i), this%c(i), this%v(i+1), this%v(i))
        call neg(this%v(i), this%v(i))
      end do

      ! normalization
      call this%t%set(0)
      do i = 1, n
        call sqr(this%u, this%v(i))            ! u = v(i)**2
        call add(this%t, this%t, this%u)       ! t += v(i)**2
      end do
      call sqrt_mpfr(this%t, this%t)           ! t = sqrt(v(1)**2 + v(2)**2 + ...)
      call div(this%t, 1.0, this%t)            ! t = 1 / sqrt(v(1)**2 + v(2)**2 + ...)
      do i = 1, n
        call mul(this%v(i), this%v(i), this%t) ! v(i) /= sqrt(v(1)**2 + v(2)**2 + ...)
      end do

      ! improve eigenvalue estimate
      call fmma(this%t, this%a(1), this%v(1), this%b(1), this%v(2))       ! t = a(1) * v(1) + b(1) * v(2)
      call mul(this%x, this%v(1), this%t)                                 ! x = v(1) * (a(1) * v(1) + b(1) * v(2))
      do i = 2, n - 1
        call fmma(this%t, this%b(i-1), this%v(i-1), this%a(i), this%v(i)) ! t = b(i-1) * v(i-1) + a(i) * v(i)
        call fma(this%t, this%b(i), this%v(i+1), this%t)                  ! t = b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1)
        call fma(this%x, this%v(i), this%t, this%x)                       ! x += v(i) * (b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1))
      end do
      call fmma(this%t, this%b(n-1), this%v(n-1), this%a(n), this%v(n))   ! t = b(n-1) * v(n-1) + a(n) * v(n)
      call fma(this%x, this%v(n), this%t, this%x)                         ! x += v(n) * (b(n-1) * v(n-1) + a(n) * v(n))

      ! convert to real
      x(j) = this%x%to_real()
      if (np > 0) then
        do i = 1, n
          v(i,j) = this%v(i)%to_real()
        end do
      end if

      ! extract unscaled weight w = v(1)**2
      call sqr(this%t, this%v(1))
      ! call mul(this%t, this%t, this%s(0))
      w(j) = this%t%to_real()
    end do

    ! derivatives
    if (np > 0) then
      allocate (dadp(n), dbdp(n-1), dAdpv(n,n), vdAdpv(n,n))
      do ip = 1, np
        ! convert dAdp to real
        do i = 1, n-1
          dadp(i) = this%dadp(i,ip)%to_real()
          dbdp(i) = this%dbdp(i,ip)%to_real()
        end do
        dadp(n) = this%dadp(n,ip)%to_real()

        ! dAdp * v
        do i = 1, n
          dAdpv(:,    i) =                  dadp * v(:,    i)
          dAdpv(1:n-1,i) = dAdpv(1:n-1,i) + dbdp * v(2:n,  i)
          dAdpv(2:n,  i) = dAdpv(2:n,  i) + dbdp * v(1:n-1,i)
        end do

        ! v' * dAdp * v
        call gemm(v, dAdpv, vdAdpv, transA = "T")

        do i = 1, n
          dxdp(i,ip) = vdAdpv(i,i)
          dwdp(i,ip) = 0
          do j = 1, n
            if (j == i) cycle
            dwdp(i,ip) = dwdp(i,ip) + vdAdpv(j,i) / (x(i) - x(j)) * v(1,j)
          end do
        end do
        ! dwdp(:,ip) = v(1,:) * (2 * dwdp(:,ip) * this%s(0)%to_real() + v(1,:) * this%dsdp(0,ip)%to_real())
        dwdp(:,ip) = 2 * v(1,:) * dwdp(:,ip)
      end do
    end if
  end subroutine

  subroutine gauss_legendre(x, w)
    !! get gauss-legendre nodes and weights
    real, intent(out) :: x(:)
      !! gauss nodes
    real, intent(out) :: w(:)
      !! gauss weights

    integer     :: n
    real        :: s, dum1(0), dum2(size(x),0), dum3(size(x),0), dum4(0)
    type(gauss) :: gs

    n = size(x)
    m4_assert(n > 0)
    m4_assert(size(w) == n)

    call gs%init(n)
    call gs%generate(moms, dum1, x, w, s, dum2, dum3, dum4)
    w = w * s
    call gs%destruct()

  contains

    subroutine moms(p, s, dsdp)
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p

      integer :: i, k

      m4_ignore(p)
      m4_ignore(dsdp)

      do i = 0, n - 1
        k = 2 * i

        ! s(k) = 2 / (k + 1)
        call s(k)%set(k + 1)
        call div(s(k), 2, s(k))

        ! s(k+1) = 0
        call s(k+1)%set(0)
      end do
    end subroutine

  end subroutine

  subroutine gauss_laguerre(x, w)
    !! get gauss-laguerre nodes and weights
    real, intent(out) :: x(:)
      !! gauss nodes
    real, intent(out) :: w(:)
      !! gauss weights

    integer     :: n
    real        :: s, dum1(0), dum2(size(x),0), dum3(size(x),0), dum4(0)
    type(gauss) :: gs

    n = size(x)
    m4_assert(n > 0)
    m4_assert(size(w) == n)

    call gs%init(n)
    call gs%generate(moms, dum1, x, w, s, dum2, dum3, dum4)
    w = w * s
    call gs%destruct()

  contains

    subroutine moms(p, s, dsdp)
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p

      integer :: i

      m4_ignore(p)
      m4_ignore(dsdp)

      call s(0)%set(1)
      call s(1)%set(1)
      do i = 2, 2*n-1
        call mul(s(i), s(i-1), i) ! s(i) = i!
      end do
    end subroutine

  end subroutine

  subroutine gauss_hermite(x, w)
    !! get gauss-hermite nodes and weights
    real, intent(out) :: x(:)
      !! gauss nodes
    real, intent(out) :: w(:)
      !! gauss weights

    integer     :: n
    real        :: s, dum1(0), dum2(size(x),0), dum3(size(x),0), dum4(0)
    type(gauss) :: gs

    n = size(x)
    m4_assert(n > 0)
    m4_assert(size(w) == n)

    call gs%init(n)
    call gs%generate(moms, dum1, x, w, s, dum2, dum3, dum4)
    w = w * s
    call gs%destruct()

  contains

    subroutine moms(p, s, dsdp)
      real,       intent(in)    :: p(:)
        !! parameters
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
      type(mpfr), intent(inout) :: dsdp(0:,:)
        !! derivatives of s wrt p

      integer :: i, k

      m4_ignore(p)
      m4_ignore(dsdp)

      do i = 0, n - 1
        k = 2 * i
        call s(k)%set(0.5 * (k + 1)) ! s(k  ) = (k+1)/2
        call gamma_mpfr(s(k), s(k))  ! s(k  ) = gamma((k+1)/2)
        call s(k+1)%set(0)           ! s(k+1) = 0
      end do
    end subroutine

  end subroutine

end module

m4_divert(0)
