#include "assert.f90.inc"

module poly_m
  use math_m
  implicit none

  type poly1D
    !! 1D interpolation polynomial (in Newton base)

    integer           :: n
      !! rank of polynomial (maximal exponent)

    real, allocatable :: x(:)
      !! interpolation nodes
    real, allocatable :: f(:)
      !! interpolation values

    real, allocatable :: a(:)
      !! coefficients (newton base)
    real, allocatable :: dadf(:,:)
      !! derivatives of a wrt f
  contains
    procedure :: init   => poly1D_init
    procedure :: interp => poly1D_interp
    procedure :: poly1D_eval_scalar
    procedure :: poly1D_eval_array
    generic   :: eval   => poly1D_eval_scalar, poly1D_eval_array
  end type

  type poly2D
    !! 2D interpolation polynomial (in newton base; nodes are grid points)

    integer :: nx
      !! x rank of polynomial (maximal exponent)
    integer :: ny
      !! y rank of polynomial (maximal exponent)

    real, allocatable :: x(:)
      !! x coords of interpolation nodes
    real, allocatable :: y(:)
      !! y coords of interpolation nodes
    real, allocatable :: f(:,:)
      !! interpolation values

    real, allocatable :: a(:,:)
      !! coefficients (newton base)
    real, allocatable :: dadf(:,:,:,:)
      !! derivatives of a wrt f
  contains
    procedure :: init   => poly2D_init
    procedure :: interp => poly2D_interp
    procedure :: poly2D_eval_scalar
    procedure :: poly2D_eval_array
    generic   :: eval   => poly2D_eval_scalar, poly2D_eval_array
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
      !! x coords of interpolation nodes
    real, allocatable :: y(:)
      !! y coords of interpolation nodes
    real, allocatable :: z(:)
      !! z coords of interpolation nodes
    real, allocatable :: f(:,:,:)
      !! interpolation values

    real, allocatable :: a(:,:,:)
      !! coefficients (newton base)
    real, allocatable :: dadf(:,:,:,:,:,:)
      !! derivatives of a wrt f
  contains
    procedure :: init   => poly3D_init
    procedure :: interp => poly3D_interp
    procedure :: poly3D_eval_scalar
    procedure :: poly3D_eval_array
    generic   :: eval   => poly3D_eval_scalar, poly3D_eval_array
  end type

contains

  subroutine poly1D_init(this, n)
    !! initialize 1D interpolation polynomial
    class(poly1D), intent(out) :: this
    integer,       intent(in)  :: n
      !! rank

    this%n = n
    allocate (this%x(1:n+1))
    allocate (this%f(1:n+1))
    allocate (this%a(1:n+1))
    allocate (this%dadf(1:n+1,1:n+1))
  end subroutine

  subroutine poly1D_interp(this, x, f)
    !! calculate interpolation coefficients
    class(poly1D), intent(inout) :: this
    real,          intent(in)    :: x(this%n+1)
      !! interpolation nodes
    real,          intent(in)    :: f(this%n+1)
      !! interpolation values

    ! local variables
    integer :: i, j, k0, k1

    ! save nodes and values
    this%x = x
    this%f = f

    ! initialize coefficients
    this%a = f
! this%dadf = eye_real(this%n+1)

    ! divided difference algorithm
    do i = 1, this%n
      do j = this%n+1, i+1, -1
        this%a(j) = (this%a(j) - this%a(j-1)) / (x(j) - x(j-i))

        ! derivatives (lower band matrix that grows each step, lower triangular matrix once i == n)
        k0 = max(j-i,1)
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
    real, optional, intent(out) :: dpdf(:)
      !! optional output derivatives of p wrt this%f

    ! local variables
    integer :: i

    ! Horner's scheme
    p = this%a(this%n+1)
    do i = this%n, 1, -1
      p = this%a(i) + (x - this%x(i)) * p
    end do

    ! derivatives
    if (present(dpdf)) then
      dpdf = this%dadf(this%n+1,:)
      do i = this%n, 1, -1
        dpdf = (x - this%x(i)) * dpdf
        dpdf(1:i) = dpdf(1:i) + this%dadf(i,1:i) ! use triangular structure of dadf
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
    real, optional, intent(out) :: dpdf(:,:)
      !! optional output derivatives of p wrt this%f

    ! local variables
    integer :: i

    ! evaluate all points
    if (present(dpdf)) then
      do i = 1, size(x)
        call poly1D_eval_scalar(this, x(i), p(i), dpdf(i,:))
      end do
    else
      do i = 1, size(x)
        call poly1D_eval_scalar(this, x(i), p(i))
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
    allocate (this%x(1:nx+1))
    allocate (this%y(1:ny+1))
    allocate (this%f(1:nx+1,1:ny+1))
    allocate (this%a(1:nx+1,1:ny+1))
    allocate (this%dadf(1:nx+1,1:ny+1,1:nx+1,1:ny+1))
  end subroutine

  subroutine poly2D_interp(this, x, y, f)
    !! calculate interpolation coefficients
    class(poly2D), intent(inout) :: this
    real,          intent(in)    :: x(this%nx+1)
      !! x coords of interpolation nodes
    real,          intent(in)    :: y(this%ny+1)
      !! y coords of interpolation nodes
    real,          intent(in)    :: f(this%nx+1,this%ny+1)
      !! interpolation values

    ! local variables
    integer :: i, j, k0, k1

    ! save nodes and values
    this%x = x
    this%y = y
    this%f = f

    ! initialize coefficients
    this%a = f
    this%dadf = 0
    do j = 1, this%ny+1; do i = 1, this%nx+1
      this%dadf(i,j,i,j) = 1.0
    end do; end do

    ! divided difference algorithm in x direction
    do i = 1, this%nx
      do j = this%nx+1, i+1, -1
        this%a(j,:) = (this%a(j,:) - this%a(j-1,:)) / (x(j) - x(j-i))

        ! derivatives
        k0 = max(j-i,1)
        k1 = j
        this%dadf(j,:,k0:k1,:) = (this%dadf(j,:,k0:k1,:) - this%dadf(j-1,:,k0:k1,:)) / (x(j) - x(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%ny
      do j = this%ny+1, i+1, -1
        this%a(:,j) = (this%a(:,j) - this%a(:,j-1)) / (y(j) - y(j-i))

        ! derivatives
        k0 = max(j-i,1)
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
    real, optional, intent(out) :: dpdf(:,:)
      !! optional output derivatives of p wrt this%f

    associate (nx => this%nx, ny => this%ny, a => this%a, dadf => this%dadf)
      block
        integer :: i, j
        real    :: px, dpxdf(nx+1,ny+1)

        ! Horner's scheme in 2D
        p = a(nx+1,ny+1)
        do i = nx, 1, -1
          p = a(i,ny+1) + (x - this%x(i)) * p
        end do
        do j = ny, 1, -1
          px = a(nx+1,j)
          do i = nx, 1, -1
            px = a(i,j) + (x - this%x(i)) * px
          end do
          p = px + (y - this%y(j)) * p
        end do

        ! derivatives
        if (present(dpdf)) then
          dpdf = dadf(nx+1,ny+1,:,:)
          do i = nx, 1, -1
            dpdf = (x - this%x(i)) * dpdf
            dpdf(1:i,:) = dadf(i,ny+1,1:i,:) + dpdf(1:i,:)
          end do
          do j = ny, 1, -1
            dpxdf(:,1:j) = dadf(nx+1,j,:,1:j)
            do i = nx, 1, -1
              dpxdf(:,1:j) = (x - this%x(i)) * dpxdf(:,1:j)
              dpxdf(1:i,1:j) = dadf(i,j,1:i,1:j) + dpxdf(1:i,1:j)
            end do
            dpdf = (y - this%y(j)) * dpdf
            dpdf(:,1:j) = dpxdf(:,1:j) + dpdf(:,1:j)
          end do
        end if
      end block
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
    real, optional, intent(out) :: dpdf(:,:,:,:)
      !! optional output derivatives of p wrt this%f

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
    allocate (this%x(1:nx+1))
    allocate (this%y(1:ny+1))
    allocate (this%z(1:nz+1))
    allocate (this%f(1:nx+1,1:ny+1,1:nz+1))
    allocate (this%a(1:nx+1,1:ny+1,1:nz+1))
    allocate (this%dadf(1:nx+1,1:ny+1,1:nz+1,1:nx+1,1:ny+1,1:nz+1))
  end subroutine

  subroutine poly3D_interp(this, x, y, z, f)
    !! calculate interpolation coefficients
    class(poly3D), intent(inout) :: this
    real,          intent(in)    :: x(this%nx+1)
      !! x coords of interpolation nodes
    real,          intent(in)    :: y(this%ny+1)
      !! y coords of interpolation nodes
    real,          intent(in)    :: z(this%ny+1)
      !! z coords of interpolation nodes
    real,          intent(in)    :: f(this%nx+1,this%ny+1,this%nz+1)
      !! interpolation values

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
    do k = 1, this%nz+1; do j = 1, this%ny+1; do i = 1, this%nx+1
      this%dadf(i,j,k,i,j,k) = 1.0
    end do; end do; end do

    ! divided difference algorithm in x direction
    do i = 1, this%nx
      do j = this%nx+1, i+1, -1
        this%a(j,:,:) = (this%a(j,:,:) - this%a(j-1,:,:)) / (x(j) - x(j-i))

        ! derivatives
        k0 = max(j-i,1)
        k1 = j
        this%dadf(j,:,:,k0:k1,:,:) = (this%dadf(j,:,:,k0:k1,:,:) - this%dadf(j-1,:,:,k0:k1,:,:)) / (x(j) - x(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%ny
      do j = this%ny+1, i+1, -1
        this%a(:,j,:) = (this%a(:,j,:) - this%a(:,j-1,:)) / (y(j) - y(j-i))

        ! derivatives
        k0 = max(j-i,1)
        k1 = j
        this%dadf(:,j,:,:,k0:k1,:) = (this%dadf(:,j,:,:,k0:k1,:) - this%dadf(:,j-1,:,:,k0:k1,:)) / (y(j) - y(j-i))
      end do
    end do

    ! divided difference algorithm in y direction
    do i = 1, this%nz
      do j = this%nz+1, i+1, -1
        this%a(:,:,j) = (this%a(:,:,j) - this%a(:,:,j-1)) / (z(j) - z(j-i))

        ! derivatives
        k0 = max(j-i,1)
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
    real, optional, intent(out) :: dpdf(:,:,:)
      !! optional output derivatives of p wrt this%f

    associate (nx => this%nx, ny => this%ny, nz => this%nz, a => this%a, dadf => this%dadf)
      block
        integer :: i, j, k
        real    :: px, dpxdf(nx+1,ny+1,nz+1), pxy, dpxydf(nx+1,ny+1,nz+1)

        ! Horner's scheme in 3D
        p = a(nx+1,ny+1,nz+1)
        do i = nx, 1, -1
          p = a(i,ny+1,nz+1) + (x - this%x(i)) * p
        end do
        do j = ny, 1, -1
          px = a(nx+1,j,nz+1)
          do i = nx, 1, -1
            px = a(i,j,nz+1) + (x - this%x(i)) * px
          end do
          p = px + (y - this%y(j)) * p
        end do
        do k = nz, 1, -1
          pxy = a(nx+1,ny+1,k)
          do i = nx, 1, -1
            pxy = a(i,ny+1,k) + (x - this%x(i)) * pxy
          end do
          do j = ny, 1, -1
            px = a(nx+1,j,k)
            do i = nx, 1, -1
              px = a(i,j,k) + (x - this%x(i)) * px
            end do
            pxy = px + (y - this%y(j)) * pxy
          end do
          p = pxy + (z - this%z(k)) * p
        end do

        ! derivatives
        if (present(dpdf)) then
          ! j = ny + 1; k = nz + 1
          dpdf = this%dadf(nx+1,ny+1,nz+1,:,:,:)
          do i = nx, 1, -1
            dpdf = (x - this%x(i)) * dpdf
            dpdf(1:i,:,:) = dadf(i,ny+1,nz+1,1:i,:,:) + dpdf(1:i,:,:)
          end do

          ! k = nz + 1
          do j = ny, 1, -1
            dpxdf(:,1:j,:) = dadf(nx+1,j,nz+1,:,1:j,:)
            do i = nx, 1, -1
              dpxdf(:,1:j,:) = (x - this%x(i)) * dpxdf(:,1:j,:)
              dpxdf(1:i,1:j,:) = dadf(i,j,nz+1,1:i,1:j,:) + dpxdf(1:i,1:j,:)
            end do
            dpdf = (y - this%y(j)) * dpdf
            dpdf(:,1:j,:) = dpxdf(:,1:j,:) + dpdf(:,1:j,:)
          end do

          do k = nz, 1, -1
            ! j = ny + 1
            dpxydf(:,:,1:k) = dadf(nx+1,ny+1,k,:,:,1:k)
            do i = nx, 1, -1
              dpxydf(:,:,1:k) = (x - this%x(i)) * dpxydf(:,:,1:k)
              dpxydf(1:i,:,1:k) = dadf(i,ny+1,k,1:i,:,1:k) + dpxydf(1:i,:,1:k)
            end do

            do j = ny, 1, -1
              dpxdf(:,1:j,1:k) = dadf(nx+1,j,k,:,1:j,1:k)
              do i = nx, 1, -1
                dpxdf(:,1:j,1:k) = (x - this%x(i)) * dpxdf(:,1:j,1:k)
                dpxdf(1:i,1:j,1:k) = dadf(i,j,k,1:i,1:j,1:k) + dpxdf(1:i,1:j,1:k)
              end do
              dpxydf(:,:,1:k) = (y - this%y(j)) * dpxydf(:,:,1:k)
              dpxydf(:,1:j,1:k) = dpxdf(:,1:j,1:k) + dpxydf(:,1:j,1:k)
            end do
            dpdf = (z - this%z(k)) * dpdf
            dpdf(:,:,1:k) = dpxydf(:,:,1:k) + dpdf(:,:,1:k)
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
    real, optional, intent(out) :: dpdf(:,:,:,:,:,:)
      !! optional output derivatives of p wrt this%f

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

end module
