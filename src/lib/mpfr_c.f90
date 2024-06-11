m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module mpfr_c_m

  use, intrinsic :: iso_c_binding

  implicit none

  ! integer kinds
  integer, parameter :: mpfr_prec_t       = c_long
  integer, parameter :: mpfr_exp_t        = c_long
  integer, parameter :: mpfr_sign_t       = c_int
  integer, parameter :: mpfr_real128      = c_long_double
  integer, parameter :: mpfr_flags_t      = c_int
  integer, parameter :: mpfr_rnd_t        = c_int
  integer, parameter :: mpfr_free_cache_t = c_int

  ! rounding modes
  enum, bind(C)
    enumerator :: MPFR_RNDN = 0,  & ! round to nearest, with ties to even
      &           MPFR_RNDZ,      & ! round toward zero
      &           MPFR_RNDU,      & ! round toward +Inf
      &           MPFR_RNDD,      & ! round toward -Inf
      &           MPFR_RNDA,      & ! round away from zero
      &           MPFR_RNDF,      & ! faithful rounding
      &           MPFR_RNDNA = -1   ! round to nearest, with ties away from zero (mpfr_round)
  end enum

  ! cache freeing policies
  enum, bind(C)
    enumerator :: MPFR_FREE_LOCAL_CACHE  = 1, &
      &           MPFR_FREE_GLOBAL_CACHE = 2
  end enum

  type, bind(C) :: mpfr_t
    !! multiprecision number (c type)
    integer(mpfr_prec_t) :: mpfr_prec = 0
    integer(mpfr_sign_t) :: mpfr_sign = 0
    integer(mpfr_exp_t)  :: mpfr_exp  = 0
    type(c_ptr)          :: mpfr_d    = c_null_ptr
  end type

  ! initialization
  interface

    subroutine mpfr_init2(x, prec) bind(C)
      import
      implicit none
      type(mpfr_t),                intent(out) :: x
      integer(mpfr_prec_t), value, intent(in)  :: prec
    end subroutine

    subroutine mpfr_clear(x) bind(C)
      import
      implicit none
      type(mpfr_t), intent(inout) :: x
    end subroutine

    subroutine mpfr_set_default_prec(prec) bind(C)
      import
      implicit none
      integer(mpfr_prec_t), value, intent(in) :: prec
    end subroutine

    pure function mpfr_get_default_prec() result(prec) bind(C)
      import
      integer(mpfr_prec_t) :: prec
    end function

    pure subroutine mpfr_set_prec(x, prec) bind(C)
      import
      type(mpfr_t),                intent(inout) :: x
      integer(mpfr_prec_t), value, intent(in)    :: prec
    end subroutine

    pure function mpfr_get_prec(x) result(prec) bind(C)
      import
      type(mpfr_t), intent(in) :: x
      integer(mpfr_prec_t)     :: prec
    end function

  end interface

  ! assignment
  interface

    function mpfr_set(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_si(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_flt(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(c_float),       value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_d(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(c_double),      value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_ld(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(c_long_double), value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_float128(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(mpfr_real128),  value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_si_2exp(rop, op, e, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_exp_t), value, intent(in)    :: e
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_set_str(rop, s, base, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      character(kind=c_char),     intent(in)    :: s(*)
      integer(c_int),      value, intent(in)    :: base
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_strtofr(rop, nptr, endptr, base, rnd) result(ret) bind(C)
      import
      type(mpfr_t),                    intent(inout) :: rop
      character(kind=c_char),          intent(in)    :: nptr(*)
      type(c_ptr),                     intent(in)    :: endptr(*)
      integer(c_int),         value,   intent(in)    :: base
      integer(mpfr_rnd_t),    value,   intent(in)    :: rnd
      integer(c_int)                                 :: ret
    end function

    pure subroutine mpfr_set_nan(x) bind(C)
      import
      type(mpfr_t), intent(inout) :: x
    end subroutine

    pure subroutine mpfr_set_inf(x, sign) bind(C)
      import
      type(mpfr_t),          intent(inout) :: x
      integer(c_int), value, intent(in)    :: sign
    end subroutine

    pure subroutine mpfr_set_zero(x, sign) bind(C)
      import
      type(mpfr_t),          intent(inout) :: x
      integer(c_int), value, intent(in)    :: sign
    end subroutine

    pure subroutine mpfr_swap(x, y) bind(C)
      import
      type(mpfr_t), intent(inout) :: x
      type(mpfr_t), intent(inout) :: y
    end subroutine

  end interface

  ! combined initialization and assignment
  interface

    function mpfr_init_set_str(x, s, base, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(out) :: x
      character(kind=c_char),     intent(in)  :: s(*)
      integer(c_int),      value, intent(in)  :: base
      integer(mpfr_rnd_t), value, intent(in)  :: rnd
      integer(c_int)                          :: ret
    end function

  end interface

  ! conversion functions
  interface

    pure function mpfr_get_flt(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      real(c_float)                          :: ret
    end function

    pure function mpfr_get_d(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      real(c_double)                         :: ret
    end function

    pure function mpfr_get_ld(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      real(c_long_double)                    :: ret
    end function

    pure function mpfr_get_float128(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      real(mpfr_real128)                     :: ret
    end function

    pure function mpfr_get_si(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_long)                        :: ret
    end function

    function mpfr_get_d_2exp(exp, op, rnd) result(ret) bind(C)
      import
      integer(c_long),            intent(out) :: exp
      type(mpfr_t),               intent(in)  :: op
      integer(mpfr_rnd_t), value, intent(in)  :: rnd
      real(c_double)                          :: ret
    end function

    function mpfr_get_ld_2exp(exp, op, rnd) result(ret) bind(C)
      import
      integer(c_long),            intent(out) :: exp
      type(mpfr_t),               intent(in)  :: op
      integer(mpfr_rnd_t), value, intent(in)  :: rnd
      real(c_long_double)                     :: ret
    end function

    function mpfr_frexp(exp, y, x, rnd) result(ret) bind(C)
      import
      integer(mpfr_exp_t),        intent(out)   :: exp
      type(mpfr_t),               intent(inout) :: y
      type(mpfr_t),               intent(in)    :: x
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    pure function mpfr_get_str_ndigits(b, p) result(ret) bind(C)
      import
      integer(c_int),       value, intent(in) :: b
      integer(mpfr_prec_t), value, intent(in) :: p
      integer(c_size_t)                       :: ret
    end function

    function mpfr_get_str(str, expptr, base, n, op, rnd) result(ret) bind(C)
      import
      type(c_ptr),         value, intent(in)  :: str
      integer(mpfr_exp_t),        intent(out) :: expptr
      integer(c_int),      value, intent(in)  :: base
      integer(c_size_t),   value, intent(in)  :: n
      type(mpfr_t),               intent(in)  :: op
      integer(mpfr_rnd_t), value, intent(in)  :: rnd
      type(c_ptr)                             :: ret
    end function

    subroutine mpfr_free_str(str) bind(C)
      import
      type(c_ptr), value, intent(in) :: str
    end subroutine

    pure function mpfr_fits_ulong_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_slong_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_uint_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_sint_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_ushort_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_sshort_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_uintmax_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

    pure function mpfr_fits_intmax_p(op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      integer(c_int)                         :: ret
    end function

  end interface

  ! arithmetic functions
  interface

    function mpfr_add(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_add_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_add_si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_add_d(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      real(c_double),      value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sub(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ui_sub(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sub_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_si_sub(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sub_si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_d_sub(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(c_double),      value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sub_d(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      real(c_double),      value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul_si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul_d(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      real(c_double),      value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sqr(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ui_div(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_si_div(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div_si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_d_div(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      real(c_double),      value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div_d(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      real(c_double),      value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sqrt(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sqrt_ui(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rec_sqrt(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_cbrt(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rootn_ui(rop, op, n, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(c_long),     value, intent(in)    :: n
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_root(rop, op, n, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(c_long),     value, intent(in)    :: n
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_neg(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_abs(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_dim(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul_2ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_mul_2si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div_2ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_div_2si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fac_ui(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fma(rop, op1, op2, op3, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      type(mpfr_t),               intent(in)    :: op3
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fms(rop, op1, op2, op3, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      type(mpfr_t),               intent(in)    :: op3
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fmma(rop, op1, op2, op3, op4, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      type(mpfr_t),               intent(in)    :: op3
      type(mpfr_t),               intent(in)    :: op4
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fmms(rop, op1, op2, op3, op4, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      type(mpfr_t),               intent(in)    :: op3
      type(mpfr_t),               intent(in)    :: op4
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_hypot(rop, x, y, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: x
      type(mpfr_t),               intent(in)    :: y
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sum(rop, tab, n, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: tab(*)
      integer(c_long),     value, intent(in)    :: n
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_dot(rop, a, b, n, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: a(*)
      type(mpfr_t),               intent(in)    :: b(*)
      integer(c_long),     value, intent(in)    :: n
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

  end interface

  ! comparison functions
  interface

    pure function mpfr_cmp(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_cmp_ui(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t),           intent(in) :: op1
      integer(c_long), value, intent(in) :: op2
      integer(c_int)                     :: ret
    end function

    pure function mpfr_cmp_si(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t),           intent(in) :: op1
      integer(c_long), value, intent(in) :: op2
      integer(c_int)                     :: ret
    end function

    pure function mpfr_cmp_d(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t),          intent(in) :: op1
      real(c_double), value, intent(in) :: op2
      integer(c_int)                    :: ret
    end function

    pure function mpfr_cmp_ld(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op1
      real(c_long_double), value, intent(in) :: op2
      integer(c_int)                         :: ret
    end function

    pure function mpfr_cmp_ui_2exp(op1, op2, e) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op1
      integer(c_long),     value, intent(in) :: op2
      integer(mpfr_exp_t), value, intent(in) :: e
      integer(c_int)                         :: ret
    end function

    pure function mpfr_cmp_si_2exp(op1, op2, e) result(ret) bind(C)
      import
      type(mpfr_t),               intent(in) :: op1
      integer(c_long),     value, intent(in) :: op2
      integer(mpfr_exp_t), value, intent(in) :: e
      integer(c_int)                         :: ret
    end function

    pure function mpfr_cmpabs(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_cmpabs_ui(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t),           intent(in) :: op1
      integer(c_long), value, intent(in) :: op2
      integer(c_int)                     :: ret
    end function

    pure function mpfr_nan_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    pure function mpfr_inf_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    pure function mpfr_number_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    pure function mpfr_zero_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    pure function mpfr_regular_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    pure function mpfr_greater_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_greaterequal_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_less_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_lessequal_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_equal_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_lessgreater_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_unordered_p(op1, op2) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op1
      type(mpfr_t), intent(in) :: op2
      integer(c_int)           :: ret
    end function

    pure function mpfr_total_order_p(x, y) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: x
      type(mpfr_t), intent(in) :: y
      integer(c_int)           :: ret
    end function

  end interface

  ! transcendental functions
  interface

    function mpfr_log(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_log_ui(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_log2(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_log10(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_log1p(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_exp(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_exp2(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_exp10(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_expm1(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_pow(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout)  :: rop
      type(mpfr_t),               intent(in)     :: op1
      type(mpfr_t),               intent(in)     :: op2
      integer(mpfr_rnd_t), value, intent(in)     :: rnd
      integer(c_int)                             :: ret
    end function

    function mpfr_pow_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_pow_si(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ui_pow_ui(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      integer(c_long),     value, intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ui_pow(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_cos(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sin(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_tan(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sin_cos(sop, cop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: sop
      type(mpfr_t),               intent(inout) :: cop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sec(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_csc(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_cot(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_acos(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_asin(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_atan(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_atan2(rop, y, x, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: y
      type(mpfr_t),               intent(in)    :: x
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_cosh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sinh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_tanh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sinh_cosh(sop, cop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: sop
      type(mpfr_t),               intent(inout) :: cop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_sech(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_csch(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_coth(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_acosh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_asinh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_atanh(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_eint(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_li2(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_gamma(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_gamma_inc(rop, op, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_lngamma(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_lgamma(rop, signp, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_int),             intent(out)   :: signp
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_digamma(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_beta(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_zeta(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_zeta_ui(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_erf(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_erfc(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_j0(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_j1(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_jn(rop, n, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: n
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_y0(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_y1(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_yn(rop, n, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(c_long),     value, intent(in)    :: n
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_agm(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ai(rop, x, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: x
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_const_log2(rop, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_const_pi(rop, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_const_euler(rop, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_const_catalan(rop, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

  end interface

  ! integer and remainder related functions
  interface

    function mpfr_rint(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_ceil(rop, op) result(ret) bind(C)
      import
      type(mpfr_t), intent(inout) :: rop
      type(mpfr_t), intent(in)    :: op
      integer(c_int)              :: ret
    end function

    function mpfr_floor(rop, op) result(ret) bind(C)
      import
      type(mpfr_t), intent(inout) :: rop
      type(mpfr_t), intent(in)    :: op
      integer(c_int)              :: ret
    end function

    function mpfr_round(rop, op) result(ret) bind(C)
      import
      type(mpfr_t), intent(inout) :: rop
      type(mpfr_t), intent(in)    :: op
      integer(c_int)              :: ret
    end function

    function mpfr_roundeven(rop, op) result(ret) bind(C)
      import
      type(mpfr_t), intent(inout) :: rop
      type(mpfr_t), intent(in)    :: op
      integer(c_int)              :: ret
    end function

    function mpfr_trunc(rop, op) result(ret) bind(C)
      import
      type(mpfr_t), intent(inout) :: rop
      type(mpfr_t), intent(in)    :: op
      integer(c_int)              :: ret
    end function

    function mpfr_rint_ceil(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rint_floor(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rint_round(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rint_roundeven(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_rint_trunc(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_frac(rop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_modf(iop, fop, op, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: iop
      type(mpfr_t),               intent(inout) :: fop
      type(mpfr_t),               intent(in)    :: op
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fmod(r, x, y, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: r
      type(mpfr_t),               intent(in)    :: x
      type(mpfr_t),               intent(in)    :: y
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_fmodquo(r, q, x, y, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: r
      integer(c_long),            intent(out)   :: q
      type(mpfr_t),               intent(in)    :: x
      type(mpfr_t),               intent(in)    :: y
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_remainder(r, x, y, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: r
      type(mpfr_t),               intent(in)    :: x
      type(mpfr_t),               intent(in)    :: y
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_remquo(r, q, x, y, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: r
      integer(c_long),            intent(out)   :: q
      type(mpfr_t),               intent(in)    :: x
      type(mpfr_t),               intent(in)    :: y
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    pure function mpfr_integer_p(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

  end interface

  ! rounding related functions
  interface

    subroutine mpfr_set_default_rounding_mode(rnd) bind(C)
      import
      integer(mpfr_rnd_t), value, intent(in) :: rnd
    end subroutine

    pure function mpfr_get_default_rounding_mode() result(ret) bind(C)
      import
      integer(mpfr_rnd_t) :: ret
    end function

    function mpfr_prec_round(x, prec, rnd) result(ret) bind(C)
      import
      type(mpfr_t),                intent(inout) :: x
      integer(mpfr_prec_t), value, intent(in)    :: prec
      integer(mpfr_rnd_t),  value, intent(in)    :: rnd
      integer(c_int)                             :: ret
    end function

    pure function mpfr_can_round(b, err, rnd1, rnd2, prec) result(ret) bind(C)
      import
      type(mpfr_t),                intent(in) :: b
      integer(mpfr_exp_t),  value, intent(in) :: err
      integer(mpfr_rnd_t),  value, intent(in) :: rnd1
      integer(mpfr_rnd_t),  value, intent(in) :: rnd2
      integer(mpfr_prec_t), value, intent(in) :: prec
      integer(c_int)                          :: ret
    end function

    pure function mpfr_min_prec(x) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: x
      integer(mpfr_prec_t)     :: ret
    end function

    pure function mpfr_print_rnd_mode(rnd) result(ret) bind(C)
      import
      integer(mpfr_rnd_t), value, intent(in) :: rnd
      type(c_ptr)                            :: ret
    end function

  end interface

  ! miscellaneous functions
  interface

    pure subroutine mpfr_nexttoward(x, y) bind(C)
      import
      type(mpfr_t), intent(inout) :: x
      type(mpfr_t), intent(in)    :: y
    end subroutine

    pure subroutine mpfr_nextabove(x) bind(C)
      import
      type(mpfr_t), intent(inout) :: x
    end subroutine

    pure subroutine mpfr_nextbelow(x) bind(C)
      import
      type(mpfr_t), intent(inout) :: x
    end subroutine

    function mpfr_min(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_max(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    pure function mpfr_get_exp(x) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: x
      integer(mpfr_exp_t)      :: ret
    end function

    function mpfr_set_exp(x, e) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: x
      integer(mpfr_exp_t), value, intent(in)    :: e
      integer(c_int)                            :: ret
    end function

    pure function mpfr_signbit(op) result(ret) bind(C)
      import
      type(mpfr_t), intent(in) :: op
      integer(c_int)           :: ret
    end function

    function mpfr_setsign(rop, op, s, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op
      integer(c_int),      value, intent(in)    :: s
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_copysign(rop, op1, op2, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: rop
      type(mpfr_t),               intent(in)    :: op1
      type(mpfr_t),               intent(in)    :: op2
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    pure function mpfr_get_version() result(ret) bind(C)
      import
      type(c_ptr) :: ret
    end function

    pure function mpfr_get_patches() result(ret) bind(C)
      import
      type(c_ptr) :: ret
    end function

    pure function mpfr_buildopt_tls_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_buildopt_float128_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_buildopt_decimal_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_buildopt_gmpinternals_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_buildopt_sharedcache_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_buildopt_tune_case() result(ret) bind(C)
      import
      type(c_ptr) :: ret
    end function

  end interface

  ! exception related functions
  interface

    pure function mpfr_get_emin() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    pure function mpfr_get_emax() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    function mpfr_set_emin(exp) result(ret) bind(C)
      import
      integer(mpfr_exp_t), value, intent(in) :: exp
      integer(c_int)                         :: ret
    end function

    function mpfr_set_emax(exp) result(ret) bind(C)
      import
      integer(mpfr_exp_t), value, intent(in) :: exp
      integer(c_int)                         :: ret
    end function

    pure function mpfr_get_emin_min() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    pure function mpfr_get_emin_max() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    pure function mpfr_get_emax_min() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    pure function mpfr_get_emax_max() result(ret) bind(C)
      import
      integer(mpfr_exp_t) :: ret
    end function

    function mpfr_check_range(x, t, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: x
      integer(c_int),      value, intent(in)    :: t
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    function mpfr_subnormalize(x, t, rnd) result(ret) bind(C)
      import
      type(mpfr_t),               intent(inout) :: x
      integer(c_int),      value, intent(in)    :: t
      integer(mpfr_rnd_t), value, intent(in)    :: rnd
      integer(c_int)                            :: ret
    end function

    subroutine mpfr_clear_underflow() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_overflow() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_divby0() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_nanflag() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_inexflag() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_erangeflag() bind(C)
      import
    end subroutine

    subroutine mpfr_clear_flags() bind(C)
      import
    end subroutine

    subroutine mpfr_set_underflow() bind(C)
      import
    end subroutine

    subroutine mpfr_set_overflow() bind(C)
      import
    end subroutine

    subroutine mpfr_set_divby0() bind(C)
      import
    end subroutine

    subroutine mpfr_set_nanflag() bind(C)
      import
    end subroutine

    subroutine mpfr_set_inexflag() bind(C)
      import
    end subroutine

    subroutine mpfr_set_erangeflag() bind(C)
      import
    end subroutine

    pure function mpfr_underflow_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_overflow_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_divby0_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_nanflag_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_inexflag_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    pure function mpfr_erangeflag_p() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

    subroutine mpfr_flags_clear(mask) bind(C)
      import
      integer(mpfr_flags_t), value, intent(in) :: mask
    end subroutine

    subroutine mpfr_flags_set(mask) bind(C)
      import
      integer(mpfr_flags_t), value, intent(in) :: mask
    end subroutine

    pure function mpfr_flags_test(mask) result(ret) bind(C)
      import
      integer(mpfr_flags_t), value, intent(in) :: mask
      integer(mpfr_flags_t)        :: ret
    end function

    pure function mpfr_flags_save() result(ret) bind(C)
      import
      integer(mpfr_flags_t) :: ret
    end function

    subroutine mpfr_flags_restore(flags, mask) bind(C)
      import
      integer(mpfr_flags_t), value, intent(in) :: flags
      integer(mpfr_flags_t), value, intent(in) :: mask
    end subroutine

  end interface

  ! memory handling functions
  interface

    subroutine mpfr_free_cache() bind(C)
      import
    end subroutine

    subroutine mpfr_free_cache2(way) bind(C)
      import
      integer(mpfr_free_cache_t), value, intent(in) :: way
    end subroutine

    subroutine mpfr_free_pool() bind(C)
      import
    end subroutine

    function mpfr_mp_memory_cleanup() result(ret) bind(C)
      import
      integer(c_int) :: ret
    end function

  end interface

end module

m4_divert(0)
