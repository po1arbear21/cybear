module test_poly_m
  use test_case_m
  use math_m
  use poly_m
  use util_m
  implicit none

contains

  subroutine test_poly()
    type(test_case) :: tc

    print "(1A)", "test_poly"
    call tc%init("poly")

    ! 1D quadratic
    block
      integer, parameter :: N = 3, Nxx = 20
      real               :: x(N), f(N), xx(Nxx), p(Nxx), dpdf(Nxx,N), pe(Nxx), dpdfe(Nxx,N), dadf(N), dbdf(N), dcdf(N)
      real,    parameter :: tol = 1e-13
      type(poly1D)       :: poly

      ! f(x) = x^2 + 2x + 3
      x = [-2.0, 1.0, 1.7]
      f = x**2 + 2*x + 3

      ! sampling points and expected values
      xx = linspace(-5.0, 5.0, Nxx)
      pe = xx**2 + 2*xx + 3

      call poly%init(N-1)
      call poly%interp(x, f)
      call poly%eval(xx, p, dpdf)

      dadf(1) =   (x(2) * x(3))/((x(1) - x(2))*(x(1) - x(3)))
      dadf(2) = - (x(1) * x(3))/((x(1) - x(2))*(x(2) - x(3)))
      dadf(3) =   (x(1) * x(2))/((x(1) - x(3))*(x(2) - x(3)))
      dbdf(1) = - (x(2) + x(3))/((x(1) - x(2))*(x(1) - x(3)))
      dbdf(2) =   (x(1) + x(3))/((x(1) - x(2))*(x(2) - x(3)))
      dbdf(3) = - (x(1) + x(2))/((x(1) - x(3))*(x(2) - x(3)))
      dcdf(1) =             1.0/((x(1) - x(2))*(x(1) - x(3)))
      dcdf(2) = -           1.0/((x(1) - x(2))*(x(2) - x(3)))
      dcdf(3) =             1.0/((x(1) - x(3))*(x(2) - x(3)))

      ! pe = a + b * xx + c * xx**2
      dpdfe(:,1) = dadf(1) + dbdf(1) * xx + dcdf(1) * xx**2
      dpdfe(:,2) = dadf(2) + dbdf(2) * xx + dcdf(2) * xx**2
      dpdfe(:,3) = dadf(3) + dbdf(3) * xx + dcdf(3) * xx**2

      call tc%assert_eq(pe, p, tol, "poly1D quadratic")
      call tc%assert_eq(dpdfe, dpdf, tol, "poly1D quadratic derivatives")
    end block

    ! 2D bilinear
    block
      integer, parameter :: Nx = 2, Ny = 2, Nxx = 4, Nyy = 4
      integer            :: i, j
      real               :: x(Nx), y(Ny), f(Nx,Ny), xx(Nxx), yy(Nyy)
      real               :: p(Nxx,Nyy), dpdf(Nxx,Nyy,Nx,Ny), pe(Nxx,Nyy), dpdfe(Nxx,Nyy,Nx,Ny)
      real,    parameter :: tol = 1e-13
      type(poly2D)       :: poly

      ! unit square
      x = [0.0, 1.0]
      y = [0.0, 1.0]
      f = reshape([0.7, 1.5, 9.2, 3.8], [2,2])

      ! sampling points and expected values
      xx = linspace(0.0, 1.0, Nxx)
      yy = linspace(0.0, 1.0, Nyy)
      do j = 1, Nyy; do i = 1, Nxx
        pe(i,j) = f(1,1) * (1.0 - xx(i)) * (1.0 - yy(j)) &
                + f(2,1) *        xx(i)  * (1.0 - yy(j)) &
                + f(1,2) * (1.0 - xx(i)) *        yy(j)  &
                + f(2,2) *        xx(i)  *        yy(j)
        dpdfe(i,j,1,1) = (1.0 - xx(i)) * (1.0 - yy(j))
        dpdfe(i,j,2,1) =        xx(i)  * (1.0 - yy(j))
        dpdfe(i,j,1,2) = (1.0 - xx(i)) *        yy(j)
        dpdfe(i,j,2,2) =        xx(i)  *        yy(j)
      end do; end do

      call poly%init(Nx-1, Ny-1)
      call poly%interp(x, y, f)
      call poly%eval(xx, yy, p, dpdf)

      call tc%assert_eq(pe, p, tol, "2D bilinear")

      do j = 1, Nyy; do i = 1, Nxx
        call tc%assert_eq(dpdfe(i,j,:,:), dpdf(i,j,:,:), tol, "2D bilinear derivatives "//int2str(i)//","//int2str(j))
      end do; end do
    end block

    ! 3D trilinear
    block
      integer, parameter :: Nx = 2, Ny = 2, Nz = 2, Nxx = 4, Nyy = 4, Nzz = 4
      integer            :: i, j, k
      real               :: x(Nx), y(Ny), z(Nz), f(Nx,Ny,Nz), xx(Nxx), yy(Nyy), zz(Nzz)
      real               :: p(Nxx,Nyy,Nzz), dpdf(Nxx,Nyy,Nzz,Nx,Ny,Nz), pe(Nxx,Nyy,Nzz), dpdfe(Nxx,Nyy,Nzz,Nx,Ny,Nz)
      real,    parameter :: tol = 1e-13
      type(poly3D)       :: poly

      ! unit square
      x = [0.0, 1.0]
      y = [0.0, 1.0]
      z = [0.0, 1.0]
      f = reshape([1.0, 1.2, 0.5, 0.8, 1.7, 1.4, 0.7, 0.9], [2,2,2])

      ! sampling points and expected values
      xx = linspace(0.0, 1.0, Nxx)
      yy = linspace(0.0, 1.0, Nyy)
      zz = linspace(0.0, 1.0, Nzz)
      do k = 1, Nzz; do j = 1, Nyy; do i = 1, Nxx
        pe(i,j,k) = f(1,1,1) * (1.0 - xx(i)) * (1.0 - yy(j)) * (1.0 - zz(k)) &
                  + f(2,1,1) *        xx(i)  * (1.0 - yy(j)) * (1.0 - zz(k)) &
                  + f(1,2,1) * (1.0 - xx(i)) *        yy(j)  * (1.0 - zz(k)) &
                  + f(2,2,1) *        xx(i)  *        yy(j)  * (1.0 - zz(k)) &
                  + f(1,1,2) * (1.0 - xx(i)) * (1.0 - yy(j)) *        zz(k)  &
                  + f(2,1,2) *        xx(i)  * (1.0 - yy(j)) *        zz(k)  &
                  + f(1,2,2) * (1.0 - xx(i)) *        yy(j)  *        zz(k)  &
                  + f(2,2,2) *        xx(i)  *        yy(j)  *        zz(k)
        dpdfe(i,j,k,1,1,1) = (1.0 - xx(i)) * (1.0 - yy(j)) * (1.0 - zz(k))
        dpdfe(i,j,k,2,1,1) =        xx(i)  * (1.0 - yy(j)) * (1.0 - zz(k))
        dpdfe(i,j,k,1,2,1) = (1.0 - xx(i)) *        yy(j)  * (1.0 - zz(k))
        dpdfe(i,j,k,2,2,1) =        xx(i)  *        yy(j)  * (1.0 - zz(k))
        dpdfe(i,j,k,1,1,2) = (1.0 - xx(i)) * (1.0 - yy(j)) *        zz(k)
        dpdfe(i,j,k,2,1,2) =        xx(i)  * (1.0 - yy(j)) *        zz(k)
        dpdfe(i,j,k,1,2,2) = (1.0 - xx(i)) *        yy(j)  *        zz(k)
        dpdfe(i,j,k,2,2,2) =        xx(i)  *        yy(j)  *        zz(k)
      end do; end do; end do

      call poly%init(Nx-1, Ny-1, Nz-1)
      call poly%interp(x, y, z, f)
      call poly%eval(xx, yy, zz, p, dpdf)

      call tc%assert_eq(pe, p, tol, "3D trilinear")

      do k = 1, Nzz; do j = 1, Nyy; do i = 1, Nxx
        call tc%assert_eq(dpdfe(i,j,k,:,:,:), dpdf(i,j,k,:,:,:), tol, "3D trilinear derivatives "//int2str(i)//","//int2str(j)//","//int2str(k))
      end do; end do; end do
    end block

    call tc%finish()
  end subroutine

end module
