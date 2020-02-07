#include "assert.f90.inc"

module polynomial_m
  use matrix_m
  implicit none

  type polynomial
    !! polynomial in multiple dimensions

    integer              :: ndim
      !! number of dimensions
    integer              :: ncoeff
      !! total number of coefficients
    integer, allocatable :: rank(:)
      !! rank for each dimension
    real,    allocatable :: p(:)
      !! coefficients (ncoeff)
  contains
    procedure :: init        => polynomial_init
    procedure :: interpolate => polynomial_interpolate
    procedure :: eval        => polynomial_eval
  end type

contains

  subroutine polynomial_init(this, ndim, rank)
    !! initialize polynomial
    class(polynomial), intent(out) :: this
    integer,           intent(in)  :: ndim
      !! number of dimensions
    integer,           intent(in)  :: rank(:)
      !! rank for each dimension

    this%ndim  = ndim
    this%ncoeff = product(rank + 1)
    allocate (this%rank(ndim), source = rank)
    allocate (this%p(this%ncoeff), source = 0.0)
  end subroutine

  subroutine polynomial_interpolate(this, r, f)
    !! interpolate polynomial through given points by solving Vandermonde system
    class(polynomial), intent(inout) :: this
    real,              intent(in)    :: r(:,:)
      !! interpolation points (ndim x ncoeff)
    real,              intent(in)    :: f(:)
      !! interpolation values

    ! local variables
    integer          :: row, col, dim, i(this%ndim)
    type(dense_real) :: mat

    ASSERT(size(r,2) == size(f))
    ASSERT(size(r,2) == this%ncoeff)

    ! initialize Vandermonde matrix
    call mat%init(this%ncoeff)
    do row = 1, this%ncoeff

      ! constant coefficient
      i = 0
      mat%d(row,1) = 1.0

      ! rest of coefficients
      do col = 2, this%ncoeff
        ! select next coefficient
        i(1) = i(1) + 1
        do dim = 1, this%ndim
          if (i(dim) > this%rank(dim)) then
            i(dim) = 0
            i(dim+1) = i(dim+1) + 1
          else
            exit
          end if
        end do

        ! set matrix entry
        mat%d(row,col) = 1.0
        do dim = 1, this%ndim
          mat%d(row,col) = mat%d(row,col) * r(dim,row)**i(dim)
        end do
      end do
    end do

    ! solve
    call mat%factorize()
    call mat%solve_vec(f, this%p)
  end subroutine

  subroutine polynomial_eval(this, r, f)
    !! evaluate polynomial for given point
    class(polynomial), intent(in)  :: this
    real,              intent(in)  :: r(:)
      !! point
    real,              intent(out) :: f
      !! output value of polynomial at point r

    ! local variables
    integer :: col, dim, i(this%ndim)
    real    :: fac

    ! constant coefficient
    i = 0
    f = this%p(1)

    ! rest of coefficients
    do col = 2, this%ncoeff
      ! select next coefficient
      i(1) = i(1) + 1
      do dim = 1, this%ndim
        if (i(dim) > this%rank(dim)) then
          i(dim) = 0
          i(dim+1) = i(dim+1) + 1
        else
          exit
        end if
      end do

      ! factor
      fac = 1.0
      do dim = 1, this%ndim
        fac = fac * r(dim)**i(dim)
      end do

      f = f + fac * this%p(col)
    end do
  end subroutine

end module