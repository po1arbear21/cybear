#include "macro.f90.inc"

module poly_m
  use math_m
  implicit none

  private
  public poly1D
  public poly2D
  public poly3D

  type poly1D
    !! 1D interpolation polynomial (in Newton base)

    integer           :: n
      !! rank of polynomial (maximal exponent)

    real, allocatable :: x(:)
      !! interpolation nodes (0:n)
    real, allocatable :: f(:)
      !! interpolation values (0:n)

    real, allocatable :: a(:)
      !! newton base coefficients (0:n)
    real, allocatable :: dadf(:,:)
      !! derivatives of a wrt f (0:n,0:n)
  contains
    procedure :: init   => poly1D_init
    procedure :: interp => poly1D_interp
    procedure :: poly1D_eval_scalar
    procedure :: poly1D_eval_array
    generic   :: eval   => poly1D_eval_scalar, poly1D_eval_array
    procedure :: poly1D_grad_scalar
    procedure :: poly1D_grad_array
    generic   :: grad   => poly1D_grad_scalar, poly1D_grad_array
  end type

  type poly2D
    !! 2D interpolation polynomial (in newton base; nodes are grid points)

    integer :: nx
      !! x rank of polynomial (maximal exponent)
    integer :: ny
      !! y rank of polynomial (maximal exponent)

    real, allocatable :: x(:)
      !! x coords of interpolation nodes (0:nx)
    real, allocatable :: y(:)
      !! y coords of interpolation nodes (0:ny)
    real, allocatable :: f(:,:)
      !! interpolation values (0:nx,0:ny)

    real, allocatable :: a(:,:)
      !! newton base coefficients (0:nx,0:ny)
    real, allocatable :: dadf(:,:,:,:)
      !! derivatives of a wrt f (0:nx,0:ny,0:nx,0:ny)
  contains
    procedure :: init   => poly2D_init
    procedure :: interp => poly2D_interp
    procedure :: poly2D_eval_scalar
    procedure :: poly2D_eval_array
    generic   :: eval   => poly2D_eval_scalar, poly2D_eval_array
    procedure :: poly2D_grad_scalar
    procedure :: poly2D_grad_array
    generic   :: grad   => poly2D_grad_scalar, poly2D_grad_array
  end type

  type poly3D
    !! 3D interpolation polynomial (in newton base; nodes are grid points)

    integer :: nx
      !! x rank of polynomial (maximal exponent)
    integer :: ny
      !! y rank of polynomial (maximal exponent)
    integer :: nz
      !! z rank of polynomial (maximal exponent)

    real, allocatable :: x(:)
      !! x coords of interpolation nodes (0:nx)
    real, allocatable :: y(:)
      !! y coords of interpolation nodes (0:ny)
    real, allocatable :: z(:)
      !! z coords of interpolation nodes (0:nz)
    real, allocatable :: f(:,:,:)
      !! interpolation values (0:nx,0:ny,0:nz)

    real, allocatable :: a(:,:,:)
      !! newton base coefficients (0:nx,0:ny,0:nz)
    real, allocatable :: dadf(:,:,:,:,:,:)
      !! derivatives of a wrt f (0:nx,0:ny,0:nz,0:nx,0:ny,0:nz)
  contains
    procedure :: init   => poly3D_init
    procedure :: interp => poly3D_interp
    procedure :: poly3D_eval_scalar
    procedure :: poly3D_eval_array
    generic   :: eval   => poly3D_eval_scalar, poly3D_eval_array
    procedure :: poly3D_grad_scalar
    procedure :: poly3D_grad_array
    generic   :: grad   => poly3D_grad_scalar, poly3D_grad_array
  end type

contains

  subroutine poly1D_init(this, n)
    !! initialize 1D interpolation polynomial
    class(poly1D), intent(out) :: this
    integer,       intent(in)  :: n
      !! rank

    this%n = n
    allocate (this%x(0:n))
    allocate (this%f(0:n))
    allocate (this%a(0:n))
    allocate (this%dadf(0:n,0:n))
  end subroutine

  subroutine poly1D_interp(this, x, f)
    !! calculate interpolation coefficients
    class(poly1D), intent(inout) :: this
    real,          intent(in)    :: x(0:)
      !! interpolation nodes (0:n)
    real,          intent(in)    :: f(0:)
      !! interpolation values (0:n)

    ! local variables
    integer :: i, j, k0, k1

    ! save nodes and values
    this%x = x
    this%f = f

    ! initialize coefficients
    this%a = f
    this%dadf = eye_real(this%n+1)

    ! divided difference algorithm
    do i = 1, this%n
      do j = this%n, i, -1
        this%a(j) = (this%a(j) - this%a(j-1)) / (x(j) - x(j-i))

        ! derivatives (lower band matrix that grows each step, lower triangular matrix once i == n)
        k0 = max(j-i,0)
        k1 = j
        this%dadf(j,k0:k1) = (this%dadf(j,k0:k1) - this%dadf(j-1,k0:k1)) / (x(j) - x(j-i))
      end do
    end do
  end subroutine

  subroutine poly1D_eval_scalar(this, x, p, dpdf)
    !! evaluate 1D interpolation polynomial at one point (using Horner's scheme)
    class(poly1D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! sampling point
    real,           intent(out) :: p
      !! output sampling value
    real, optional, intent(out) :: dpdf(0:)
      !! optional: output derivatives of p wrt this%f

    ! local variables
    integer :: i

    ! Horner's scheme
    p = this%a(this%n)
    do i = this%n-1, 0, -1
      p = p * (x - this%x(i)) + this%a(i)
    end do

    ! derivatives
    if (present(dpdf)) then
      dpdf = this%dadf(this%n,:)
      do i = this%n-1, 0, -1
        dpdf = dpdf * (x - this%x(i))
        dpdf(0:i) = dpdf(0:i) + this%dadf(i,0:i) ! use triangular structure of dadf
      end do
    end if
  end subroutine

  subroutine poly1D_eval_array(this, x, p, dpdf)
    !! evaluate 1D interpolation polynomial for multiple points
    class(poly1D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! sampling points
    real,           intent(out) :: p(:)
      !! output sampling values
    real, optional, intent(out) :: dpdf(:,0:)
      !! optional output derivatives of p wrt this%f

    ! local variables
    integer :: i

    ! evaluate all points
    if (present(dpdf)) then
      do i = 1, size(x)
        call poly1D_eval_scalar(this, x(i), p(i), dpdf = dpdf(i,:))
      end do
    else
      do i = 1, size(x)
        call poly1D_eval_scalar(this, x(i), p(i))
      end do
    end if
  end subroutine

  subroutine poly1D_grad_scalar(this, x, gradp, dgradpdf)
    !! evaluate gradient of 1D interpolation polynomial at one point (using Horner's scheme)
    class(poly1D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! sampling point
    real,           intent(out) :: gradp
      !! optional: output gradient of p (derivative of p wrt x)
    real, optional, intent(out) :: dgradpdf(0:)
      !! optional: output derivatives of gradp wrt this%f

    ! local variables
    integer :: i
    real    :: p, dpdf(0:this%n)

    ! Horner's scheme
    gradp = 0
    p     = this%a(this%n)
    do i = this%n-1, 0, -1
      gradp = p + gradp * (x - this%x(i))
      p     = p * (x - this%x(i)) + this%a(i)
    end do

    ! derivatives
    if (present(dgradpdf)) then
      dgradpdf = 0
      dpdf     = this%dadf(this%n,:)
      do i = this%n-1, 0, -1
        dgradpdf  = dpdf + dgradpdf * (x - this%x(i))
        dpdf      = dpdf * (x - this%x(i))
        dpdf(0:i) = dpdf(0:i) + this%dadf(i,0:i) ! use triangular structure of dadf
      end do
    end if
  end subroutine

  subroutine poly1D_grad_array(this, x, gradp, dgradpdf)
    !! evaluate gradient of 1D interpolation polynomial for multiple points
    class(poly1D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! sampling point
    real,           intent(out) :: gradp(:)
      !! optional: output gradient of p at x (derivative of p wrt x)
    real, optional, intent(out) :: dgradpdf(:,0:)
      !! optional: output derivatives of gradp wrt this%f

    ! local variables
    integer :: i

    ! evaluate gradient for all points
    if (present(dgradpdf)) then
      do i = 1, size(x)
        call poly1D_grad_scalar(this, x(i), gradp(i), dgradpdf = dgradpdf(i,:))
      end do
    else
      do i = 1, size(x)
        call poly1D_grad_scalar(this, x(i), gradp(i))
      end do
    end if
  end subroutine

  subroutine poly2D_init(this, nx, ny)
    !! initialize 2D interpolation polynomial
    class(poly2D), intent(out) :: this
    integer,       intent(in)  :: nx
      !! x rank
    integer,       intent(in)  :: ny
      !! y rank

    this%nx = nx
    this%ny = ny
    allocate (this%x(0:nx))
    allocate (this%y(0:ny))
    allocate (this%f(0:nx,0:ny))
    allocate (this%a(0:nx,0:ny))
    allocate (this%dadf(0:nx,0:ny,0:nx,0:ny))
  end subroutine

  subroutine poly2D_interp(this, x, y, f)
    !! calculate interpolation coefficients
    class(poly2D), intent(inout) :: this
    real,          intent(in)    :: x(0:)
      !! x coords of interpolation nodes (0:nx)
    real,          intent(in)    :: y(0:)
      !! y coords of interpolation nodes (0:ny)
    real,          intent(in)    :: f(0:,0:)
      !! interpolation values (0:nx, 0:ny)

    ! local variables
    integer :: i, j, k0, k1

    ! save nodes and values
    this%x = x
    this%y = y
    this%f = f

    ! initialize coefficients
    this%a = f
    this%dadf = 0
    do j = 0, this%ny; do i = 0, this%nx
      this%dadf(i,j,i,j) = 1.0
    end do; end do

    ! divided difference algorithm in x direction
    do i = 1, this%nx
      do j = this%nx, i, -1
        this%a(j,:) = (this%a(j,:) - this%a(j-1,:)) / (x(j) - x(j-i))

        ! derivatives
        k0 = max(j-i,0)
        k1 = j
        this%dadf(j,:,k0:k1,:) = (this%dadf(j,:,k0:k1,:) - this%dadf(j-1,:,k0:k1,:)) / (x(j) - x(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%ny
      do j = this%ny, i, -1
        this%a(:,j) = (this%a(:,j) - this%a(:,j-1)) / (y(j) - y(j-i))

        ! derivatives
        k0 = max(j-i,0)
        k1 = j
        this%dadf(:,j,:,k0:k1) = (this%dadf(:,j,:,k0:k1) - this%dadf(:,j-1,:,k0:k1)) / (y(j) - y(j-i))
      end do
    end do
  end subroutine

  subroutine poly2D_eval_scalar(this, x, y, p, dpdf)
    !! evaluate 2D interpolation polynomial at one point (using Horner's scheme)
    class(poly2D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! x coord of sampling point
    real,           intent(in)  :: y
      !! y coord of sampling point
    real,           intent(out) :: p
      !! output sampling value
    real, optional, intent(out) :: dpdf(0:,0:)
      !! optional output derivatives of p wrt this%f (0:nx,0:ny)

    ! local variables
    integer :: i, j
    real    :: px, dpxdf(0:this%nx,0:this%ny)

    associate (nx => this%nx, ny => this%ny, a => this%a, dadf => this%dadf)
      ! Horner's scheme in 2D
      p = a(nx,ny)
      do i = nx-1, 0, -1
        p = p * (x - this%x(i)) + a(i,ny)
      end do
      do j = ny-1, 0, -1
        px = a(nx,j)
        do i = nx-1, 0, -1
          px = px * (x - this%x(i)) + a(i,j)
        end do
        p = p * (y - this%y(j)) + px
      end do

      ! derivatives
      if (present(dpdf)) then
        dpdf = dadf(nx,ny,:,:)
        do i = nx-1, 0, -1
          dpdf = dpdf * (x - this%x(i))
          dpdf(0:i,:) = dpdf(0:i,:) + dadf(i,ny,0:i,:)
        end do
        do j = ny-1, 0, -1
          dpxdf(:,0:j) = dadf(nx,j,:,0:j)
          do i = nx-1, 0, -1
            dpxdf( : ,0:j) = dpxdf( : ,0:j) * (x - this%x(i))
            dpxdf(0:i,0:j) = dpxdf(0:i,0:j) + dadf(i,j,0:i,0:j)
          end do
          dpdf = dpdf * (y - this%y(j))
          dpdf(:,0:j) = dpdf(:,0:j) + dpxdf(:,0:j)
        end do
      end if
    end associate
  end subroutine

  subroutine poly2D_eval_array(this, x, y, p, dpdf)
    !! evaluate 2D interpolation polynomial for points in a grid
    class(poly2D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! x coords of sampling points
    real,           intent(in)  :: y(:)
      !! y coords of sampling points
    real,           intent(out) :: p(:,:)
      !! output sampling value
    real, optional, intent(out) :: dpdf(:,:,0:,0:)
      !! optional output derivatives of p wrt this%f (1:size(x),1:size(y),0:nx,0:ny)

    ! local variables
    integer :: i, j

    if (present(dpdf)) then
      do j = 1, size(y); do i = 1, size(x)
        call poly2D_eval_scalar(this, x(i), y(j), p(i,j), dpdf(i,j,:,:))
      end do; end do
    else
      do j = 1, size(y); do i = 1, size(x)
        call poly2D_eval_scalar(this, x(i), y(j), p(i,j))
      end do; end do
    end if
  end subroutine

  subroutine poly2D_grad_scalar(this, x, y, gradp, dgradpdf)
    !! evaluate gradient of 2D interpolation polynomial at one point (using Horner's scheme)
    class(poly2D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! x coord of sampling point
    real,           intent(in)  :: y
      !! y coord of sampling point
    real,           intent(out) :: gradp(:)
      !! output gradient of p (2)
    real, optional, intent(out) :: dgradpdf(:,0:,0:)
      !! optional: output derivatives of gradp wrt this%f (2,0:nx,0:ny)

    associate (nx => this%nx, ny => this%ny, a => this%a, dadf => this%dadf)
      block
        integer :: i, j
        real    :: p, dpdf(0:nx,0:ny), px, dpxdf(0:nx,0:ny), gradpx, dgradpxdf(0:nx,0:ny)
        ! Horner's scheme in 2D
        gradp = 0
        p     = a(nx,ny)
        do i = nx-1, 0, -1
          gradp(1) = p + gradp(1) * (x - this%x(i))
          p = p * (x - this%x(i)) + a(i,ny)
        end do
        do j = ny-1, 0, -1
          gradpx = 0
          px     = a(nx,j)
          do i = nx-1, 0, -1
            gradpx = px + gradpx * (x - this%x(i))
            px     = px * (x - this%x(i)) + a(i,j)
          end do
          gradp(1) =     gradp(1) * (y - this%y(j)) + gradpx
          gradp(2) = p + gradp(2) * (y - this%y(j))
          p        = p * (y - this%y(j)) + px
        end do

        ! derivatives
        if (present(dgradpdf)) then
          dgradpdf = 0
          dpdf     = dadf(nx,ny,:,:)
          do i = nx-1, 0, -1
            dgradpdf(1,:,:) = dpdf + dgradpdf(1,:,:) * (x - this%x(i))
            dpdf            = dpdf * (x - this%x(i))
            dpdf(0:i,:)     = dpdf(0:i,:) + dadf(i,ny,0:i,:)
          end do
          do j = ny-1, 0, -1
            dgradpxdf(:,0:j) = 0
            dpxdf(    :,0:j) = dadf(nx,j,:,0:j)
            do i = nx-1, 0, -1
              dgradpxdf(: ,0:j) = dpxdf( : ,0:j) + dgradpxdf(:,0:j) * (x - this%x(i))
              dpxdf(    : ,0:j) = dpxdf( : ,0:j) * (x - this%x(i))
              dpxdf(   0:i,0:j) = dpxdf(0:i,0:j) + dadf(i,j,0:i,0:j)
            end do
            dgradpdf(1,:, : ) =        dgradpdf(1,:,:) * (y - this%y(j))
            dgradpdf(1,:,0:j) = dgradpdf(1,:,0:j) + dgradpxdf(:,0:j)
            dgradpdf(2,:, : ) = dpdf + dgradpdf(2,:,:) * (y - this%y(j))
            dpdf              = dpdf * (y - this%y(j))
            dpdf(:,0:j)       = dpdf(:,0:j) + dpxdf(:,0:j)
          end do
        end if
      end block
    end associate
  end subroutine

  subroutine poly2D_grad_array(this, x, y, gradp, dgradpdf)
    !! evaluate gradient of 2D interpolation polynomial for multiple points
    class(poly2D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! x coord of sampling point
    real,           intent(in)  :: y(:)
      !! y coord of sampling point
    real,           intent(out) :: gradp(:,:,:)
      !! output gradient of p (2 x size(x) x size(y))
    real, optional, intent(out) :: dgradpdf(:,:,:,0:,0:)
      !! optional: output derivatives of gradp wrt this%f (2,size(x),size(y),0:nx,0:ny)

    ! local variables
    integer :: i, j

    if (present(dgradpdf)) then
      do j = 1, size(y); do i = 1, size(x)
        call poly2D_grad_scalar(this, x(i), y(j), gradp(:,i,j), dgradpdf(:,i,j,:,:))
      end do; end do
    else
      do j = 1, size(y); do i = 1, size(x)
        call poly2D_grad_scalar(this, x(i), y(j), gradp(:,i,j))
      end do; end do
    end if
  end subroutine

  subroutine poly3D_init(this, nx, ny, nz)
    !! initialize 3D interpolation polynomial
    class(poly3D), intent(out) :: this
    integer,       intent(in)  :: nx
      !! x rank
    integer,       intent(in)  :: ny
      !! y rank
    integer,       intent(in)  :: nz
      !! z rank

    this%nx = nx
    this%ny = ny
    this%nz = nz
    allocate (this%x(0:nx))
    allocate (this%y(0:ny))
    allocate (this%z(0:nz))
    allocate (this%f(0:nx,0:ny,0:nz))
    allocate (this%a(0:nx,0:ny,0:nz))
    allocate (this%dadf(0:nx,0:ny,0:nz,0:nx,0:ny,0:nz))
  end subroutine

  subroutine poly3D_interp(this, x, y, z, f)
    !! calculate interpolation coefficients
    class(poly3D), intent(inout) :: this
    real,          intent(in)    :: x(0:)
      !! x coords of interpolation nodes (0:nx)
    real,          intent(in)    :: y(0:)
      !! y coords of interpolation nodes (0:ny)
    real,          intent(in)    :: z(0:)
      !! z coords of interpolation nodes (0:nz)
    real,          intent(in)    :: f(0:,0:,0:)
      !! interpolation values (0:nx,0:ny,0:nz)

    ! local variables
    integer :: i, j, k, k0, k1

    ! save nodes and values
    this%x = x
    this%y = y
    this%z = z
    this%f = f

    ! initialize coefficients
    this%a = f
    this%dadf = 0
    do k = 0, this%nz; do j = 0, this%ny; do i = 0, this%nx
      this%dadf(i,j,k,i,j,k) = 1.0
    end do; end do; end do

    ! divided difference algorithm in x direction
    do i = 1, this%nx
      do j = this%nx, i, -1
        this%a(j,:,:) = (this%a(j,:,:) - this%a(j-1,:,:)) / (x(j) - x(j-i))

        ! derivatives
        k0 = max(j-i,0)
        k1 = j
        this%dadf(j,:,:,k0:k1,:,:) = (this%dadf(j,:,:,k0:k1,:,:) - this%dadf(j-1,:,:,k0:k1,:,:)) / (x(j) - x(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%ny
      do j = this%ny, i, -1
        this%a(:,j,:) = (this%a(:,j,:) - this%a(:,j-1,:)) / (y(j) - y(j-i))

        ! derivatives
        k0 = max(j-i,0)
        k1 = j
        this%dadf(:,j,:,:,k0:k1,:) = (this%dadf(:,j,:,:,k0:k1,:) - this%dadf(:,j-1,:,:,k0:k1,:)) / (y(j) - y(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%nz
      do j = this%nz, i, -1
        this%a(:,:,j) = (this%a(:,:,j) - this%a(:,:,j-1)) / (z(j) - z(j-i))

        ! derivatives
        k0 = max(j-i,0)
        k1 = j
        this%dadf(:,:,j,:,:,k0:k1) = (this%dadf(:,:,j,:,:,k0:k1) - this%dadf(:,:,j-1,:,:,k0:k1)) / (z(j) - z(j-i))
      end do
    end do
  end subroutine

  subroutine poly3D_eval_scalar(this, x, y, z, p, dpdf)
    !! evaluate 3D interpolation polynomial at one point (using Horner's scheme)
    class(poly3D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! x coord of sampling point
    real,           intent(in)  :: y
      !! y coord of sampling point
    real,           intent(in)  :: z
      !! z coord of sampling point
    real,           intent(out) :: p
      !! output sampling value
    real, optional, intent(out) :: dpdf(0:,0:,0:)
      !! optional output derivatives of p wrt this%f (0:nx,0:ny,0:nz)

    associate (nx => this%nx, ny => this%ny, nz => this%nz, a => this%a, dadf => this%dadf)
      block
        integer :: i, j, k
        real    :: px, dpxdf(0:nx,0:ny,0:nz), pxy, dpxydf(0:nx,0:ny,0:nz)

        ! Horner's scheme in 3D
        p = a(nx,ny,nz)
        do i = nx-1, 0, -1
          p = p * (x - this%x(i)) + a(i,ny,nz)
        end do
        do j = ny-1, 0, -1
          px = a(nx,j,nz)
          do i = nx-1, 0, -1
            px = px * (x - this%x(i)) + a(i,j,nz)
          end do
          p = p * (y - this%y(j)) + px
        end do
        do k = nz-1, 0, -1
          pxy = a(nx,ny,k)
          do i = nx-1, 0, -1
            pxy = pxy * (x - this%x(i)) + a(i,ny,k)
          end do
          do j = ny-1, 0, -1
            px = a(nx,j,k)
            do i = nx-1, 0, -1
              px = px * (x - this%x(i)) + a(i,j,k)
            end do
            pxy = pxy * (y - this%y(j)) + px
          end do
          p = p * (z - this%z(k)) + pxy
        end do

        ! derivatives
        if (present(dpdf)) then
          ! j = ny; k = nz
          dpdf = this%dadf(nx,ny,nz,:,:,:)
          do i = nx-1, 0, -1
            dpdf          = dpdf * (x - this%x(i))
            dpdf(0:i,:,:) = dpdf(0:i,:,:) + dadf(i,ny,nz,0:i,:,:)
          end do

          ! k = nz
          do j = ny-1, 0, -1
            dpxdf(:,0:j,:) = dadf(nx,j,nz,:,0:j,:)
            do i = nx-1, 0, -1
              dpxdf( : ,0:j,:) = dpxdf( : ,0:j,:) * (x - this%x(i))
              dpxdf(0:i,0:j,:) = dpxdf(0:i,0:j,:) + dadf(i,j,nz,0:i,0:j,:)
            end do
            dpdf = dpdf * (y - this%y(j))
            dpdf(:,0:j,:) = dpdf(:,0:j,:) + dpxdf(:,0:j,:)
          end do

          do k = nz-1, 0, -1
            ! j = ny
            dpxydf(:,:,0:k) = dadf(nx,ny,k,:,:,0:k)
            do i = nx-1, 0, -1
              dpxydf( : ,:,0:k) = dpxydf( : ,:,0:k) * (x - this%x(i))
              dpxydf(0:i,:,0:k) = dpxydf(0:i,:,0:k) + dadf(i,ny,k,0:i,:,0:k)
            end do

            do j = ny-1, 0, -1
              dpxdf(:,0:j,0:k) = dadf(nx,j,k,:,0:j,0:k)
              do i = nx-1, 0, -1
                dpxdf( : ,0:j,0:k) = dpxdf(:,0:j,0:k) * (x - this%x(i))
                dpxdf(0:i,0:j,0:k) = dpxdf(0:i,0:j,0:k) + dadf(i,j,k,0:i,0:j,0:k)
              end do
              dpxydf(:, : ,0:k) = dpxydf(:,:,0:k) * (y - this%y(j))
              dpxydf(:,0:j,0:k) = dpxydf(:,0:j,0:k) + dpxdf(:,0:j,0:k)
            end do
            dpdf = dpdf * (z - this%z(k))
            dpdf(:,:,0:k) = dpdf(:,:,0:k) + dpxydf(:,:,0:k)
          end do
        end if
      end block
    end associate
  end subroutine

  subroutine poly3D_eval_array(this, x, y, z, p, dpdf)
    !! evaluate 3D interpolation polynomial for points in a grid
    class(poly3D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! x coords of sampling points
    real,           intent(in)  :: y(:)
      !! y coords of sampling points
    real,           intent(in)  :: z(:)
      !! z coords of sampling points
    real,           intent(out) :: p(:,:,:)
      !! output sampling value
    real, optional, intent(out) :: dpdf(:,:,:,0:,0:,0:)
      !! optional output derivatives of p wrt this%f (size(x),size(y),size(z),0:nx,0:ny,0:nz)

    ! local variables
    integer :: i, j, k

    if (present(dpdf)) then
      do k = 1, size(z); do j = 1, size(y); do i = 1, size(x)
        call poly3D_eval_scalar(this, x(i), y(j), z(k), p(i,j,k), dpdf(i,j,k,:,:,:))
      end do; end do; end do
    else
      do k = 1, size(z); do j = 1, size(y); do i = 1, size(x)
        call poly3D_eval_scalar(this, x(i), y(j), z(k), p(i,j,k))
      end do; end do; end do
    end if
  end subroutine

  subroutine poly3D_grad_scalar(this, x, y, z, gradp, dgradpdf)
    !! evaluate gradient of 3D interpolation polynomial at one point (using Horner's scheme)
    class(poly3D),  intent(in)  :: this
    real,           intent(in)  :: x
      !! x coord of sampling point
    real,           intent(in)  :: y
      !! y coord of sampling point
    real,           intent(in)  :: z
      !! z coord of sampling point
    real,           intent(out) :: gradp(:)
      !! output gradient of p (3)
    real, optional, intent(out) :: dgradpdf(:,0:,0:,0:)
      !! optional: output derivatives of gradp wrt this%f (3,0:nx,0:ny,0:nz)

    associate (nx => this%nx, ny => this%ny, nz => this%nz, a => this%a, dadf => this%dadf)
      block
        integer :: i, j, k
        real    :: p, dpdf(0:nx,0:ny,0:nz), px, dpxdf(0:nx,0:ny,0:nz), pxy, dpxydf(0:nx,0:ny,0:nz)
        real    :: gradpx, dgradpxdf(0:nx,0:ny,0:nz), gradpxy(2), dgradpxydf(2,0:nx,0:ny,0:nz)

        ! Horner's scheme in 3D
        gradp = 0
        p     = a(nx,ny,nz)
        do i = nx-1, 0, -1
          gradp(1) = p + gradp(1) * (x - this%x(i))
          p        = p * (x - this%x(i)) + a(i,ny,nz)
        end do
        do j = ny-1, 0, -1
          gradpx = 0
          px     = a(nx,j,nz)
          do i = nx-1, 0, -1
            gradpx = px + gradpx * (x - this%x(i))
            px     = px * (x - this%x(i)) + a(i,j,nz)
          end do
          gradp(1) =     gradp(1) * (y - this%y(j)) + gradpx
          gradp(2) = p + gradp(2) * (y - this%y(j))
          p        = p * (y - this%y(j)) + px
        end do
        do k = nz-1, 0, -1
          gradpxy = 0
          pxy     = a(nx,ny,k)
          do i = nx-1, 0, -1
            gradpxy(1) = pxy + gradpxy(1) * (x - this%x(i))
            pxy        = pxy * (x - this%x(i)) + a(i,ny,k)
          end do
          do j = ny-1, 0, -1
            gradpx = 0
            px     = a(nx,j,k)
            do i = nx-1, 0, -1
              gradpx = px + gradpx * (x - this%x(i))
              px     = px * (x - this%x(i)) + a(i,j,k)
            end do
            gradpxy(1) =       gradpxy(1) * (y - this%y(j)) + gradpx
            gradpxy(2) = pxy + gradpxy(2) * (y - this%y(j))
            pxy        = pxy * (y - this%y(j)) + px
          end do
          gradp(1) =     gradp(1) * (z - this%z(k)) + gradpxy(1)
          gradp(2) =     gradp(2) * (z - this%z(k)) + gradpxy(2)
          gradp(3) = p + gradp(3) * (z - this%z(k))
          p        = p * (z - this%z(k)) + pxy
        end do

        ! derivatives
        if (present(dgradpdf)) then
          ! j = ny; k = nz
          dgradpdf = 0
          dpdf     = this%dadf(nx,ny,nz,:,:,:)
          do i = nx-1, 0, -1
            dgradpdf(1,:,:,:) = dpdf + dgradpdf(1,:,:,:) * (x - this%x(i))
            dpdf              = dpdf * (x - this%x(i))
            dpdf(0:i,:,:)     = dpdf(0:i,:,:) + dadf(i,ny,nz,0:i,:,:)
          end do

          ! k = nz
          do j = ny-1, 0, -1
            dgradpxdf(:,0:j,:) = 0
            dpxdf(    :,0:j,:) = dadf(nx,j,nz,:,0:j,:)
            do i = nx-1, 0, -1
              dgradpxdf(: ,0:j,:) = dpxdf(:,0:j,:) + dgradpxdf(:,0:j,:) * (x - this%x(i))
              dpxdf(    : ,0:j,:) = dpxdf(:,0:j,:) * (x - this%x(i))
              dpxdf(   0:i,0:j,:) = dpxdf(0:i,0:j,:) + dadf(i,j,nz,0:i,0:j,:)
            end do
            dgradpdf(1,:, : ,:) =        dgradpdf(1,:, : ,:) * (y - this%y(j))
            dgradpdf(1,:,0:j,:) =        dgradpdf(1,:,0:j,:) + dgradpxdf(:,0:j,:)
            dgradpdf(2,:, : ,:) = dpdf + dgradpdf(2,:, : ,:) * (y - this%y(j))
            dpdf                = dpdf * (y - this%y(j))
            dpdf(      :,0:j,:) = dpdf(:,0:j,:) + dpxdf(:,0:j,:)
          end do

          do k = nz-1, 0, -1
            ! j = ny
            dgradpxydf(:,:,:,0:k) = 0
            dpxydf(      :,:,0:k) = dadf(nx,ny,k,:,:,0:k)
            do i = nx-1, 0, -1
              dgradpxydf(1, : ,:,0:k) = dpxydf(:,:,0:k) + dgradpxydf(1,:,:,0:k) * (x - this%x(i))
              dpxydf(       : ,:,0:k) = dpxydf(:,:,0:k) * (x - this%x(i))
              dpxydf(      0:i,:,0:k) = dpxydf(0:i,:,0:k) + dadf(i,ny,k,0:i,:,0:k)
            end do

            do j = ny-1, 0, -1
              dgradpxdf(:,0:j,0:k) = 0
              dpxdf(    :,0:j,0:k) = dadf(nx,j,k,:,0:j,0:k)
              do i = nx-1, 0, -1
                dgradpxdf( : ,0:j,0:k) = dpxdf(:,0:j,0:k) + dgradpxdf(:,0:j,0:k) * (x - this%x(i))
                dpxdf(     : ,0:j,0:k) = dpxdf(:,0:j,0:k) * (x - this%x(i))
                dpxdf(    0:i,0:j,0:k) = dpxdf(0:i,0:j,0:k) + dadf(i,j,k,0:i,0:j,0:k)
              end do
              dgradpxydf(1,:, : ,0:k) = dgradpxydf(1,:, : ,0:k) * (y - this%y(j))
              dgradpxydf(1,:,0:j,0:k) = dgradpxydf(1,:,0:j,0:k) + dgradpxdf(:,0:j,0:k)
              dgradpxydf(2,:, : ,0:k) = dpxydf(:,:,0:k) + dgradpxydf(2,:,:,0:k) * (y - this%y(j))
              dpxydf(      :, : ,0:k) = dpxydf(:,:,0:k) * (y - this%y(j))
              dpxydf(      :,0:j,0:k) = dpxydf(:,0:j,0:k) + dpxdf(:,0:j,0:k)
            end do
            dgradpdf(1,:,:, : ) = dgradpdf(1,:,:, : ) * (z - this%z(k))
            dgradpdf(1,:,:,0:k) = dgradpdf(1,:,:,0:k) + dgradpxydf(1,:,:,0:k)
            dgradpdf(2,:,:, : ) = dgradpdf(2,:,:, : ) * (z - this%z(k))
            dgradpdf(2,:,:,0:k) = dgradpdf(2,:,:,0:k) + dgradpxydf(2,:,:,0:k)
            dgradpdf(3,:,:, : ) = dpdf + dgradpdf(3,:,:,:) * (z - this%z(k))
            dpdf                = dpdf * (z - this%z(k))
            dpdf(      :,:,0:k) = dpdf(:,:,0:k) + dpxydf(:,:,0:k)
          end do
        end if
      end block
    end associate
  end subroutine

  subroutine poly3D_grad_array(this, x, y, z, gradp, dgradpdf)
    !! evaluate gradient of 3D interpolation polynomial for multiple points
    class(poly3D),  intent(in)  :: this
    real,           intent(in)  :: x(:)
      !! x coord of sampling point
    real,           intent(in)  :: y(:)
      !! y coord of sampling point
    real,           intent(in)  :: z(:)
      !! z coord of sampling point
    real,           intent(out) :: gradp(:,:,:,:)
      !! output gradient of p (3)
    real, optional, intent(out) :: dgradpdf(:,:,:,:,0:,0:,0:)
      !! optional: output derivatives of gradp wrt this%f (3,size(x),size(y),size(z),0:nx,0:ny,0:nz)

    ! local variables
    integer :: i, j, k

    if (present(dgradpdf)) then
      do k = 1, size(z); do j = 1, size(y); do i = 1, size(x)
        call poly3D_grad_scalar(this, x(i), y(j), z(k), gradp(:,i,j,k), dgradpdf(:,i,j,k,:,:,:))
      end do; end do; end do
    else
      do k = 1, size(z); do j = 1, size(y); do i = 1, size(x)
        call poly3D_grad_scalar(this, x(i), y(j), z(k), gradp(:,i,j,k))
      end do; end do; end do
    end if
  end subroutine

end module
