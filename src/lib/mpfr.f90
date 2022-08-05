m4_include(../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module mpfr_m

  use error_m,         only: program_error
  use iso_fortran_env, only: int64, real128
  use mpfr_c_m
  use string_m,        only: string
  use util_m,          only: strlen, f2cstring

  implicit none

  private

  public MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD, MPFR_RNDA, MPFR_RNDF, MPFR_RNDNA
  public MPFR_FREE_LOCAL_CACHE, MPFR_FREE_GLOBAL_CACHE
  public mpfr
  public mpfr_startup, mpfr_cleanup
  public get_default_rnd, set_default_rnd
  public get_default_prec, set_default_prec
  public add, sub, mul, sqr, div, pow, sqrt, neg, abs
  public log, log1p, exp, expm1
  public sin, asin, cos, acos, tan, atan
  public sinh, asinh, cosh, acosh, tanh, atanh

  type mpfr
    !! multiprecision number (fortran type)
    type(mpfr_t) :: m
      !! underlying c type
  contains
    procedure :: init     => mpfr_init
    procedure :: destruct => mpfr_destruct
    procedure :: get_prec => mpfr_get_precision

    generic   :: set      => set_mpfr, &
      &                      set_real, &
      &                      set_int,  &
      &                      set_string
    generic   :: init_set => init_set_mpfr, &
      &                      init_set_real, &
      &                      init_set_int,  &
      &                      init_set_string

    procedure :: to_real    => mpfr_to_real
    procedure :: to_real128 => mpfr_to_real128
    procedure :: to_int     => mpfr_to_int
    procedure :: to_int64   => mpfr_to_int64
    procedure :: to_string  => mpfr_to_string

    procedure, private :: set_mpfr        => mpfr_set_mpfr
    procedure, private :: set_real        => mpfr_set_real
    procedure, private :: set_int         => mpfr_set_int
    procedure, private :: set_string      => mpfr_set_string
    procedure, private :: init_set_mpfr   => mpfr_init_set_mpfr
    procedure, private :: init_set_real   => mpfr_init_set_real
    procedure, private :: init_set_int    => mpfr_init_set_int
    procedure, private :: init_set_string => mpfr_init_set_string
  end type

  integer(mpfr_rnd_t)  :: default_rnd
  integer(mpfr_prec_t) :: default_prec

  interface add
    module procedure :: add_mpfr_mpfr
    module procedure :: add_mpfr_real
    module procedure :: add_mpfr_int
  end interface

  interface sub
    module procedure :: sub_mpfr_mpfr
    module procedure :: sub_mpfr_real
    module procedure :: sub_real_mpfr
    module procedure :: sub_mpfr_int
    module procedure :: sub_int_mpfr
  end interface

  interface mul
    module procedure :: mul_mpfr_mpfr
    module procedure :: mul_mpfr_real
    module procedure :: mul_mpfr_int
  end interface

  interface sqr
    module procedure :: sqr_mpfr
  end interface

  interface div
    module procedure :: div_mpfr_mpfr
    module procedure :: div_mpfr_real
    module procedure :: div_real_mpfr
    module procedure :: div_mpfr_int
    module procedure :: div_int_mpfr
  end interface

  interface pow
    module procedure :: pow_mpfr_mpfr
    module procedure :: pow_mpfr_int
    module procedure :: pow_int_mpfr
  end interface

  interface sqrt
    module procedure :: sqrt_mpfr
  end interface

  interface neg
    module procedure :: neg_mpfr
  end interface

  interface abs
    module procedure :: abs_mpfr
  end interface

  interface log
    module procedure :: log_mpfr
  end interface

  interface log1p
    module procedure :: log1p_mpfr
  end interface

  interface exp
    module procedure :: exp_mpfr
  end interface

  interface expm1
    module procedure :: expm1_mpfr
  end interface

  interface sin
    module procedure :: sin_mpfr
  end interface

  interface asin
    module procedure :: asin_mpfr
  end interface

  interface cos
    module procedure :: cos_mpfr
  end interface

  interface acos
    module procedure :: acos_mpfr
  end interface

  interface tan
    module procedure :: tan_mpfr
  end interface

  interface atan
    module procedure :: atan_mpfr
  end interface

  interface sinh
    module procedure :: sinh_mpfr
  end interface

  interface asinh
    module procedure :: asinh_mpfr
  end interface

  interface cosh
    module procedure :: cosh_mpfr
  end interface

  interface acosh
    module procedure :: acosh_mpfr
  end interface

  interface tanh
    module procedure :: tanh_mpfr
  end interface

  interface atanh
    module procedure :: atanh_mpfr
  end interface


contains

  subroutine mpfr_startup(rnd, prec)
    !! startup MPFR by setting default rounding mode and precision (can safely be called multiple times)
    integer, optional, intent(in) :: rnd
      !! optional: default rounding mode (MPFR_RNDN)
    integer, optional, intent(in) :: prec
      !! optional: default precision (128)

    default_rnd  = int(MPFR_RNDN, kind = mpfr_rnd_t)
    if (present(rnd)) default_rnd = int(rnd, kind = mpfr_rnd_t)
    call mpfr_set_default_rounding_mode(default_rnd)

    default_prec = 128
    if (present(prec)) default_prec = int(prec, kind=mpfr_prec_t)
    call mpfr_set_default_prec(default_prec)
  end subroutine

  subroutine mpfr_cleanup()
    !! cleanup memory used internally
    integer(c_int) :: ret

    ret = mpfr_mp_memory_cleanup()
    if (ret /= 0) call program_error("MPFR internal error while cleaning up memory")
  end subroutine

  pure function get_default_rnd() result(rnd)
    !! get default rounding mode
    integer :: rnd
      !! return rounding mode

    rnd = int(default_rnd)
  end function

  subroutine set_default_rnd(rnd)
    !! set default rounding mode (not threadsafe)
    integer, intent(in) :: rnd
      !! new default rounding mode

    default_rnd = int(rnd, kind = mpfr_rnd_t)
    call mpfr_set_default_rounding_mode(default_rnd)
  end subroutine

  pure function get_default_prec() result(prec)
    !! get default precision
    integer :: prec
      !! return precision in bits

    prec = int(default_prec)
  end function

  subroutine set_default_prec(prec)
    !! set default precision in bits (not threadsafe)
    integer :: prec
      !! new default precision in bits

    default_prec = int(prec, kind = mpfr_prec_t)
    call mpfr_set_default_prec(default_prec)
  end subroutine


  subroutine mpfr_init(this, prec)
    !! initialize multiprecision number
    class(mpfr),       intent(out) :: this
    integer, optional, intent(in)  :: prec
      !! precision in bits

    integer(mpfr_prec_t) :: prec_

    prec_ = default_prec
    if (present(prec)) prec_ = int(prec, kind = mpfr_prec_t)

    call mpfr_init2(this%m, prec_)
  end subroutine

  subroutine mpfr_destruct(this)
    !! destruct multiprecision number (free memory)
    class(mpfr), intent(inout) :: this

    if (c_associated(this%m%mpfr_d)) call mpfr_clear(this%m)
  end subroutine

  elemental function mpfr_get_precision(this) result(prec)
    !! get precision of multiprecision number
    class(mpfr), intent(in) :: this
    integer                 :: prec
      !! return precision in bits

    prec = int(mpfr_get_prec(this%m))
  end function

  elemental function mpfr_to_real(this, rnd) result(r)
    !! convert multiprecision number to real
    class(mpfr),       intent(in) :: this
    integer, optional, intent(in) :: rnd
      !! optional rounding mode
    real                          :: r
      !! return real number

    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    r = mpfr_get_d(this%m, rnd_)
  end function

  elemental function mpfr_to_real128(this, rnd) result(r)
    !! convert multiprecision number to real128
    class(mpfr),       intent(in) :: this
    integer, optional, intent(in) :: rnd
      !! optional rounding mode
    real(real128)                 :: r
      !! return real number (quad precision)

    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    r = mpfr_get_float128(this%m, rnd_)
  end function

  elemental function mpfr_to_int(this, rnd) result(i)
    !! convert multiprecision number to integer
    class(mpfr),       intent(in) :: this
    integer, optional, intent(in) :: rnd
      !! optional rounding mode
    integer                       :: i
      !! return integer number of default kind

    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    i = int(mpfr_get_si(this%m, rnd_))
  end function

  elemental function mpfr_to_int64(this, rnd) result(i)
    !! convert multiprecision number to 64-bit integer
    class(mpfr),       intent(in) :: this
    integer, optional, intent(in) :: rnd
      !! optional rounding mode
    integer(int64)                :: i
      !! return 64-bit integer number

    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    i = int(mpfr_get_si(this%m, rnd_), kind = int64)
  end function

  function mpfr_to_string(this, base, n, rnd) result(s)
    !! convert multiprecision number to string
    class(mpfr),       intent(in) :: this
    integer, optional, intent(in) :: base
      !! optional base (default: 10)
    integer, optional, intent(in) :: n
      !! optional number of significant digits
    integer, optional, intent(in) :: rnd
      !! optional rounding mode
    character(:), allocatable     :: s
      !! return string

    integer(c_int)              :: base_
    integer(c_size_t)           :: n_
    integer(mpfr_rnd_t)         :: rnd_

    type(c_ptr)                 :: cp
    character, pointer          :: p(:)

    integer(mpfr_exp_t)         :: expptr
    character(digits(expptr)+1) :: expbuf
    integer                     :: explen, clen, i, i0

    base_ = 10
    if (present(base)) base_ = int(base, kind = c_int)

    n_ = 0
    if (present(n)) n_ = n

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ! convert to c string
    cp = mpfr_get_str(c_null_ptr, expptr, base_, n_, this%m, rnd_)

    ! convert to fortran pointer
    clen = int(strlen(cp))
    call c_f_pointer(cp, p, shape = [clen])

    ! convert exponent to string
    write (expbuf, "(SP,I0)") expptr - 1
    explen = len_trim(expbuf)

    ! allocate result string
    allocate (character(clen + explen + 2) :: s)

    ! output optional sign, first digit and decimal point
    if ((p(1) == "+") .or. (p(1) == "-")) then
      s(1:1) = p(1)
      s(2:2) = p(2)
      s(3:3) = "."
      i0 = 4
    else
      s(1:1) = p(1)
      s(2:2) = "."
      i0 = 3
    end if

    ! output rest of digits
    do i = i0, clen + 1
      s(i:i) = p(i-1)
    end do

    ! output exponent
    s(i:i) = "E"
    s(i+1:) = expbuf(1:explen)

    ! cleanup
    call mpfr_free_str(cp)
  end function

  subroutine mpfr_set_mpfr(this, mp, rnd)
    !! set multiprecision number to other mp number
    class(mpfr),       intent(inout) :: this
    type(mpfr),        intent(in)    :: mp
      !! other multiprecision number
    integer, optional, intent(in)    :: rnd
      !! optional rounding mode

    integer(mpfr_rnd_t) :: rnd_
    integer(c_int)      :: ret

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_set(this%m, mp%m, rnd_)
  end subroutine

  subroutine mpfr_set_real(this, r, rnd)
    !! set multiprecision number to real
    class(mpfr),       intent(inout) :: this
    real,              intent(in)    :: r
      !! real number
    integer, optional, intent(in)    :: rnd
      !! optional rounding mode

    integer(mpfr_rnd_t) :: rnd_
    integer(c_int)      :: ret

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_set_d(this%m, r, rnd_)
  end subroutine

  subroutine mpfr_set_int(this, i, rnd)
    !! set multiprecision number to integer
    class(mpfr),       intent(inout) :: this
    integer,           intent(in)    :: i
      !! integer number
    integer, optional, intent(in)    :: rnd
      !! optional rounding mode

    integer(mpfr_rnd_t) :: rnd_
    integer(c_int)      :: ret

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_set_si(this%m, int(i, kind = c_long), rnd_)
  end subroutine

  subroutine mpfr_set_string(this, s, rnd)
    !! set multiprecision number to string
    class(mpfr),       intent(inout) :: this
    character(*),      intent(in)    :: s
      !! string containing number
    integer, optional, intent(in)    :: rnd
      !! optional rounding mode

    integer(mpfr_rnd_t) :: rnd_
    integer(c_int)      :: ret

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_set_str(this%m, f2cstring(s), int(0, kind = c_int), rnd_)
  end subroutine

  subroutine mpfr_init_set_mpfr(this, mp, prec, rnd)
    !! initialize multiprecision number by other mp number
    class(mpfr),       intent(out) :: this
    type(mpfr),        intent(in)  :: mp
      !! other multiprecision number
    integer, optional, intent(in)  :: prec
      !! optional precision in bits (not provided: use mp%get_prec())
    integer, optional, intent(in)  :: rnd
      !! optional rounding mode

    integer :: prec_

    if (present(prec)) then
      prec_ = prec
    else
      prec_ = mp%get_prec()
    end if

    call this%init(prec = prec_)
    call this%set(mp, rnd = rnd)
  end subroutine

  subroutine mpfr_init_set_real(this, r, prec, rnd)
    !! initialize multiprecision number by real
    class(mpfr),       intent(out) :: this
    real,              intent(in)  :: r
      !! real number
    integer, optional, intent(in)  :: prec
      !! optional precision in bits
    integer, optional, intent(in)  :: rnd
      !! optional rounding mode

    call this%init(prec = prec)
    call this%set(r, rnd = rnd)
  end subroutine

  subroutine mpfr_init_set_int(this, i, prec, rnd)
    !! initialize multiprecision number by integer
    class(mpfr),       intent(out) :: this
    integer,           intent(in)  :: i
      !! integer number
    integer, optional, intent(in)  :: prec
      !! optional precision in bits
    integer, optional, intent(in)  :: rnd
      !! optional rounding mode

    call this%init(prec = prec)
    call this%set(i, rnd = rnd)
  end subroutine

  subroutine mpfr_init_set_string(this, s, rnd)
    !! initialize multiprecision number by string
    class(mpfr),       intent(out) :: this
    character(*),      intent(in)  :: s
      !! string containing number
    integer, optional, intent(in)  :: rnd
      !! optional rounding mode

    integer(mpfr_rnd_t) :: rnd_
    integer(c_int)      :: ret

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_init_set_str(this%m, f2cstring(s), int(0, kind = c_int), rnd_)
  end subroutine

  subroutine add_mpfr_mpfr(r, x, y, rnd)
    !! r = x + y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_add(r%m, x%m, y%m, rnd_)
  end subroutine

  subroutine add_mpfr_real(r, x, y, rnd)
    !! r = x + y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    real,              intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_add_d(r%m, x%m, y, rnd_)
  end subroutine

  subroutine add_mpfr_int(r, x, y, rnd)
    !! r = x + y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer,           intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_add_si(r%m, x%m, int(y, kind=c_long), rnd_)
  end subroutine

  subroutine sub_mpfr_mpfr(r, x, y, rnd)
    !! r = x - y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sub(r%m, x%m, y%m, rnd_)
  end subroutine

  subroutine sub_mpfr_real(r, x, y, rnd)
    !! r = x - y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    real,              intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sub_d(r%m, x%m, y, rnd_)
  end subroutine

  subroutine sub_real_mpfr(r, x, y, rnd)
    !! r = x - y
    type(mpfr),        intent(inout) :: r
    real,              intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_d_sub(r%m, x, y%m, rnd_)
  end subroutine

  subroutine sub_mpfr_int(r, x, y, rnd)
    !! r = x - y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer,           intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sub_si(r%m, x%m, int(y, kind=c_long), rnd_)
  end subroutine

  subroutine sub_int_mpfr(r, x, y, rnd)
    !! r = x - y
    type(mpfr),        intent(inout) :: r
    integer,           intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_si_sub(r%m, int(x, kind=c_long), y%m, rnd_)
  end subroutine

  subroutine mul_mpfr_mpfr(r, x, y, rnd)
    !! r = x * y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_mul(r%m, x%m, y%m, rnd_)
  end subroutine

  subroutine mul_mpfr_real(r, x, y, rnd)
    !! r = x * y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    real,              intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_mul_d(r%m, x%m, y, rnd_)
  end subroutine

  subroutine mul_mpfr_int(r, x, y, rnd)
    !! r = x * y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer,           intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_mul_si(r%m, x%m, int(y, kind=c_long), rnd_)
  end subroutine

  subroutine sqr_mpfr(r, x, rnd)
    !! r = x**2
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sqr(r%m, x%m, rnd_)
  end subroutine

  subroutine div_mpfr_mpfr(r, x, y, rnd)
    !! r = x / y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_div(r%m, x%m, y%m, rnd_)
  end subroutine

  subroutine div_mpfr_real(r, x, y, rnd)
    !! r = x / y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    real,              intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_div_d(r%m, x%m, y, rnd_)
  end subroutine

  subroutine div_real_mpfr(r, x, y, rnd)
    !! r = x / y
    type(mpfr),        intent(inout) :: r
    real,              intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_d_div(r%m, x, y%m, rnd_)
  end subroutine

  subroutine div_mpfr_int(r, x, y, rnd)
    !! r = x / y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer,           intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_div_si(r%m, x%m, int(y, kind=c_long), rnd_)
  end subroutine

  subroutine div_int_mpfr(r, x, y, rnd)
    !! r = x / y
    type(mpfr),        intent(inout) :: r
    integer,           intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_si_div(r%m, int(x, kind=c_long), y%m, rnd_)
  end subroutine

  subroutine pow_mpfr_mpfr(r, x, y, rnd)
    !! r = x ** y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_pow(r%m, x%m, y%m, rnd_)
  end subroutine

  subroutine pow_mpfr_int(r, x, y, rnd)
    !! r = x ** y
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer,           intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_pow_si(r%m, x%m, int(y, kind=c_long), rnd_)
  end subroutine

  subroutine pow_int_mpfr(r, x, y, rnd)
    !! r = x ** y
    type(mpfr),        intent(inout) :: r
    integer,           intent(in)    :: x
    type(mpfr),        intent(in)    :: y
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    if (x < 0) call program_error("x must be non-negative")
    ret = mpfr_ui_pow(r%m, int(x, kind=c_long), y%m, rnd_)
  end subroutine

  subroutine sqrt_mpfr(r, x, rnd)
    !! r = sqrt(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sqrt(r%m, x%m, rnd_)
  end subroutine

  subroutine neg_mpfr(r, x, rnd)
    !! r = -x
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_neg(r%m, x%m, rnd_)
  end subroutine

  subroutine abs_mpfr(r, x, rnd)
    !! r = abs(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)       :: ret
    integer(mpfr_rnd_t)  :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_abs(r%m, x%m, rnd_)
  end subroutine

  subroutine log_mpfr(r, x, rnd)
    !! r = log(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_log(r%m, x%m, rnd_)
  end subroutine

  subroutine log1p_mpfr(r, x, rnd)
    !! r = log1p(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_log1p(r%m, x%m, rnd_)
  end subroutine

  subroutine exp_mpfr(r, x, rnd)
    !! r = exp(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_exp(r%m, x%m, rnd_)
  end subroutine

  subroutine expm1_mpfr(r, x, rnd)
    !! r = expm1(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_expm1(r%m, x%m, rnd_)
  end subroutine

  subroutine sin_mpfr(r, x, rnd)
    !! r = sin(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sin(r%m, x%m, rnd_)
  end subroutine

  subroutine asin_mpfr(r, x, rnd)
    !! r = asin(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_asin(r%m, x%m, rnd_)
  end subroutine

  subroutine cos_mpfr(r, x, rnd)
    !! r = cos(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_cos(r%m, x%m, rnd_)
  end subroutine

  subroutine acos_mpfr(r, x, rnd)
    !! r = acos(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_acos(r%m, x%m, rnd_)
  end subroutine

  subroutine tan_mpfr(r, x, rnd)
    !! r = tan(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_tan(r%m, x%m, rnd_)
  end subroutine

  subroutine atan_mpfr(r, x, rnd)
    !! r = atan(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_atan(r%m, x%m, rnd_)
  end subroutine

  subroutine sinh_mpfr(r, x, rnd)
    !! r = sinh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_sinh(r%m, x%m, rnd_)
  end subroutine

  subroutine asinh_mpfr(r, x, rnd)
    !! r = asinh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_asinh(r%m, x%m, rnd_)
  end subroutine

  subroutine cosh_mpfr(r, x, rnd)
    !! r = cosh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_cosh(r%m, x%m, rnd_)
  end subroutine

  subroutine acosh_mpfr(r, x, rnd)
    !! r = acosh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_acosh(r%m, x%m, rnd_)
  end subroutine

  subroutine tanh_mpfr(r, x, rnd)
    !! r = tanh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_tanh(r%m, x%m, rnd_)
  end subroutine

  subroutine atanh_mpfr(r, x, rnd)
    !! r = atanh(x)
    type(mpfr),        intent(inout) :: r
    type(mpfr),        intent(in)    :: x
    integer, optional, intent(in)    :: rnd

    integer(c_int)      :: ret
    integer(mpfr_rnd_t) :: rnd_

    rnd_ = default_rnd
    if (present(rnd)) rnd_ = int(rnd, kind = mpfr_rnd_t)

    ret = mpfr_atanh(r%m, x%m, rnd_)
  end subroutine

end module

m4_divert(0)
