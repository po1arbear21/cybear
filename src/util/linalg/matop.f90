#include "../macro.f90.inc"

module matop_m
  use error_m
  use iso_c_binding,  only: c_loc, c_f_pointer
  use matrix_m

  implicit none

  private
  public matop_real
  public matop_cmplx
  public single_matop_real
  public single_matop_cmplx
  public chain_matop_real
  public chain_matop_cmplx
  public matop_c2r

#define T real
#define TT real
#include "matop_def.f90.inc"

#define T cmplx
#define TT complex
#include "matop_def.f90.inc"

  type, extends(matop_real) :: matop_c2r
    !! a complex matop which gets its input arguments in the exec calls as real arrays of double length.
    !! real and imaginary parts are included in the real array separately.
    !!
    !! use case: implementaion of zfgmres.

    class(matop_cmplx), pointer :: mop_c => null()
      !! complex matop
  contains
    procedure :: matop_c2r_init
    generic   :: init  => matop_c2r_init
    procedure :: exec1 => matop_c2r_exec1
    procedure :: exec2 => matop_c2r_exec2
  end type

contains

#define T real
#define TT real
#include "matop_imp.f90.inc"

#define T cmplx
#define TT complex
#include "matop_imp.f90.inc"

  subroutine matop_c2r_init(this, mop_c)
    ! initialzes by setting pointer to complex matop

    class(matop_c2r),   intent(out)        :: this
    class(matop_cmplx), intent(in), target :: mop_c
      !! underlying complex matop

    ! init base
    call this%init("", mop_c%nrows, mop_c%ncols)

    this%mop_c => mop_c
  end subroutine

  subroutine matop_c2r_exec1(this, x, y)
    !! compute mul_vec of input and saved mat_op.

    class(matop_c2r), intent(in)  :: this
    real, target,     intent(in)  :: x(:)
      !! input vector: real and imaginary parts.
    real, target,     intent(out) :: y(:)
      !! output vector: real and imaginary parts.

    complex, pointer :: x_ptr(:), y_ptr(:)

    ! create complex pointer to x, y
    call c_f_pointer(c_loc(x), x_ptr, shape=shape(x)/2)
    call c_f_pointer(c_loc(y), y_ptr, shape=shape(y)/2)

    call this%mop_c%exec(x_ptr, y_ptr)
  end subroutine

  subroutine matop_c2r_exec2(this, x, y)
    !! compute mul_mat of input and saved mat_op.

    class(matop_c2r), intent(in)  :: this
    real, target,     intent(in)  :: x(:,:)
    real, target,     intent(out) :: y(:,:)

    complex, pointer :: x_ptr(:,:), y_ptr(:,:)

    ! create complex pointer to x, y
    call c_f_pointer(c_loc(x), x_ptr, shape=[size(x, dim=1)/2, size(x, dim=2)])
    call c_f_pointer(c_loc(y), y_ptr, shape=[size(y, dim=1)/2, size(y, dim=2)])

    call this%mop_c%exec(x_ptr, y_ptr)
  end subroutine

end module
