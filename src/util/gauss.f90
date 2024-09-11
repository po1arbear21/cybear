m4_include(macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module gauss_m
  !! Calculate Gauss quadrature nodes and weights using Golub-Welsh algorithm
  !! requires library MPFR

  use error_m,  only: assert_failed, program_error
  use lapack95, only: stev
  use mpfr_m,   only: add, div, fma, fmma, fms, gamma_mpfr, mpfr, mpfr_startup, mpfr_cleanup, mul, neg, sqr, sqrt_mpfr, sub

  implicit none

  private
  public gauss

  interface
    subroutine moments(s)
      !! calculate moments for custom weight function
      import mpfr
      type(mpfr), intent(inout) :: s(0:)
        !! moments, already initialized
    end subroutine
  end interface

contains

  subroutine gauss(weight, x, w, prec, momfun)
    character(*),                 intent(in)  :: weight
      !! weight function ("legendre", "laguerre", "hermite" or "custom")
    real,                         intent(out) :: x(:)
      !! output gauss quadrature nodes
    real,                         intent(out) :: w(:)
      !! output gauss quadrature weights
    integer,            optional, intent(in)  :: prec
      !! precision for MPFR calculations (default: 1024)
    procedure(moments), optional              :: momfun
      !! calculate moments for weight = "custom"

    integer                 :: i, k, n, m, prec_
    type(mpfr), allocatable :: s(:), a(:), b(:)

    prec_ = 1024
    if (present(prec)) prec_ = prec
    call mpfr_startup(prec = prec_)

    n = size(x)
    m = 2 * n - 1

    ! allocate memory
    allocate (s(0:m), a(n), b(n-1))
    do i = 0, m
      call s(i)%init()
    end do
    do i = 1, n
      call a(i)%init()
      if (i < n) call b(i)%init()
    end do

    ! set moments
    select case (weight)
    case ("legendre")
      do i = 0, n-1
        k = 2 * i
        call s(k)%set(k + 1)    ! s(k  ) = k + 1
        call div(s(k), 2, s(k)) ! s(k  ) = 2 / (k + 1)
        call s(k+1)%set(0)      ! s(k+1) = 0
      end do
    case ("laguerre")
      call s(0)%set(1)
      call s(1)%set(1)
      do i = 2, m
        call mul(s(i), s(i-1), i) ! s(i) = i!
      end do
    case ("hermite")
      do i = 0, n-1
        k = 2 * i
        call s(k)%set(0.5 * (k + 1)) ! s(k  ) = (k+1)/2
        call gamma_mpfr(s(k), s(k))  ! s(k  ) = gamma((k+1)/2)
        call s(k+1)%set(0)           ! s(k+1) = 0
      end do
    case ("custom")
      m4_assert(present(momfun))
      call momfun(s)
    case default
      call program_error("Unknown weight function '"//weight//"'")
    end select

    ! calculate recurrence relation
    call recurrence(s, a, b)

    ! solve eigenvalue problem to get gauss nodes and weights
    call eigenvalues(s, a, b, x, w)

    ! cleanup
    do i = 0, m
      call s(i)%destruct()
    end do
    do i = 1, n
      call a(i)%destruct()
      if (i < n) call b(i)%destruct()
    end do
    call mpfr_cleanup()
  end subroutine

  subroutine recurrence(s, a, b)
    type(mpfr), intent(in)    :: s(0:)
    type(mpfr), intent(inout) :: a(:)
    type(mpfr), intent(inout) :: b(:)

    integer                 :: i, j, k, n, m
    type(mpfr)              :: t, u
    type(mpfr), allocatable :: sg(:,:)

    m = ubound(s, 1)
    n = size(a)

    ! allocate memory
    allocate (sg(n,n+1))
    call t%init()
    call u%init()
    do i = 1, n
      do j = 1, n + 1
        call sg(i,j)%init()
      end do
    end do

    ! recurrence relation
    do i = 1, n
      do j = 1, n + 1
        call sg(i,j)%set(0)
        do k = 1, i - 1
          call mul(t, sg(k,i), sg(k,j))            ! t = sg(k,i) * sg(k,j)
          call div(t, t, sg(k,k))                  ! t = sg(k,i) * sg(k,j) / sg(k,k)
          call add(sg(i,j), sg(i,j), t)            ! sg(i,j) = sg(i,j) + sg(k,i) * sg(k,j) / sg(k,k)
        end do
        call sub(sg(i,j), s(i+j-2), sg(i,j))       ! sg(i,j) = s(i+j-2) - sg(i,j)
      end do
    end do
    call t%set(0)
    do i = 1, n
      call div(u, sg(i,i+1), sg(i,i))              ! u = c_i = sg(i,i+1) / sg(i,i)
      call add(t, t,  u)                           ! t = c_i - c_{i-1}
      call a(i)%set(t)                             ! a(i) = c_i - c_{i-1}
      call neg(t, u)                               ! t = - c_i
    end do
    do i = 1, n-1
      call div(t, sg(i+1,i+1), sg(i,i))            ! t = c_i = sg(i+1,i+1) / sg(i,i)
      call sqrt_mpfr(b(i), t)                      ! b_i = sqrt(c_i)
    end do

    ! cleanup
    call t%destruct()
    call u%destruct()
    do i = 1, n
      do j = 1, n + 1
        call sg(i,j)%destruct()
      end do
    end do
  end subroutine

  subroutine eigenvalues(s, a, b, x, w)
    type(mpfr), intent(in)  :: s(0:)
    type(mpfr), intent(in)  :: a(:)
    type(mpfr), intent(in)  :: b(:)
    real,       intent(out) :: x(:)
    real,       intent(out) :: w(:)

    integer                 :: i, j, m, n
    real,       allocatable :: d(:), v(:,:)
    type(mpfr)              :: t, u, xx
    type(mpfr), allocatable :: c(:), vv(:)

    m = ubound(s, 1)
    n = size(a)

    ! allocate memory
    allocate (c(n), d(n-1), v(n,n), vv(n))
    call t%init()
    call u%init()
    call xx%init()
    do i = 1, n
      call c(i)%init()
      call vv(i)%init()
    end do

    ! solve approximate matrix using LAPACK
    do i = 1, n-1
      x(i) = a(i)%to_real()
      d(i) = b(i)%to_real()
    end do
    x(n) = a(n)%to_real()
    call stev(x, d, v)

    ! refine using a single Rayleigh quotient iteration step for each eigenvalue
    do j = 1, n
      ! load estimate for eigenvalue and vector
      call xx%set(x(j))
      do i = 1, n
        call vv(i)%set(v(i,j))
      end do

      ! forward substitution
      call sub(t, a(1), xx)                        ! t = a(1) - x
      call div(t, 1.0, t)                          ! t = 1 / (a(1) - x)
      call mul(c(1), b(1), t)                      ! c(1) = b(1) / (a(1) - x)
      call mul(vv(1), vv(1), t)                    ! v(1) = v(1) / (a(1) - x)
      do i = 2, n-1
        call sub(t, a(i), xx)                      ! t = a(i) - x
        call fms(t, b(i-1), c(i-1), t)             ! t = b(i-1) * c(i-1) - (a(i) - x)
        call div(t, 1.0, t)                        ! t = 1 / (b(i-1) * c(i-1) - (a(i) - x))
        call mul(c(i), b(i), t)                    ! c(i) = b(i) / (b(i-1) * c(i-1) - (a(i) - x))
        call neg(c(i), c(i))                       ! c(i) = b(i) / (a(i) - x - b(i-1) * c(i-1))
        call fms(u, b(i-1), vv(i-1), vv(i))        ! u = b(i-1) * v(i-1) - v(i)
        call mul(vv(i), u, t)                      ! v(i) = (v(i) - b(i-1) * v(i-1)) / (a(i) - x - b(i-1) * c(i-1))
      end do
      call sub(t, a(n), xx)                        ! t = a(n) - x
      call fms(t, b(n-1), c(n-1), t)               ! t = b(n-1) * c(n-1) - (a(n) - x)
      call fms(u, b(n-1), vv(n-1), vv(n))          ! u = b(n-1) * v(n-1) - v(n)
      call div(vv(n), u, t)                        ! v(n) = (v(n) - b(n-1) * v(n-1)) / (a(n) - x - b(n-1) * c(n-1))

      ! back substitution
      do i = n-1, 1, -1
        call fms(vv(i), c(i), vv(i+1), vv(i))      ! v(i) = c(i) * v(i+1) - v(i)
        call neg(vv(i), vv(i))                     ! v(i) = v(i) - c(i) * v(i+1)
      end do

      ! normalization
      call t%set(0.0)
      do i = 1, n
        call sqr(u, vv(i))                         ! u  = v(i)**2
        call add(t, t, u)                          ! t += v(i)**2
      end do
      call sqrt_mpfr(t, t)                         ! t = sqrt(v(1)**2 + v(2)**2 + ...)
      do i = 1, n
        call div(vv(i), vv(i), t)                  ! v(i) /= sqrt(v(1)**2 + v(2)**2 + ....)
      end do

      ! improve eigenvalue estimate
      call fmma(t, a(1), vv(1), b(1), vv(2))       ! t = a(1) * v(1) + b(1) * v(2)
      call mul(xx, vv(1), t)                       ! x = v(1) * (a(1) * v(1) + b(1) * v(2))
      do i = 2, n-1
        call fmma(t, b(i-1), vv(i-1), a(i), vv(i)) ! t = b(i-1) * v(i-1) + a(i) * v(i)
        call fma(t, b(i), vv(i+1), t)              ! t = b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1)
        call fma(xx, vv(i), t, xx)                 ! x = x + v(i) * (b(i-1) * v(i-1) + a(i) * v(i) + b(i) * v(i+1))
      end do
      call fmma(t, b(n-1), vv(n-1), a(n), vv(n))   ! t = b(n-1) * v(n-1) + a(n) * v(n)
      call fma(xx, vv(n), t, xx)                   ! x = x + v(n) * (b(n-1) * v(n-1) + a(n) * v(n))

      ! extract node and weight w = v(1)**2 * s(0)
      x(j) = xx%to_real()
      call sqr(t, vv(1))
      call mul(t, t, s(0))
      w(j) = t%to_real()
    end do

    ! cleanup
    call t%destruct()
    call u%destruct()
    call xx%destruct()
    do i = 1, n
      call c(i)%destruct()
      call vv(i)%destruct()
    end do
  end subroutine

end module

m4_divert(0)
