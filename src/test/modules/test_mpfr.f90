m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module test_mpfr_m

  use mpfr_c_m
  use mpfr_m
  use string_m,    only: string
  use test_case_m, only: test_case

  implicit none

  private
  public test_mpfr

contains

  subroutine test_mpfr()
    type(test_case) :: tc

    call mpfr_startup()

    call tc%init("mpfr")

    call tc%assert(mpfr_buildopt_tls_p() /= 0, "MPFR compiled with --enable-thread-safe")
    call tc%assert(mpfr_buildopt_float128_p() /= 0, "MPFR compiled with --enable-float128")
    call tc%assert(mpfr_buildopt_sharedcache_p() == 0, "MPFR not compiled with --enable-shared-cache")

    ! get precision
    block
      type(mpfr) :: m
      integer    :: prec

      call m%init(100)

      prec = m%get_prec()
      call tc%assert_eq(100, prec, "get_prec")

      call m%destruct()
    end block

    ! init_set
    block
      type(mpfr)   :: mr, mm, mi, ms
      type(string) :: s, s_exp

      call mr%init_set(1.0)
      call tc%assert_eq(1.0, mr%to_real(), 0.0, "init_set real")

      call mm%init_set(mr)
      call tc%assert_eq(1.0, mm%to_real(), 0.0, "init_set mpfr")

      call mi%init_set(1)
      call tc%assert_eq(1, mi%to_int(), "init_set int")

      call ms%init_set("16.0")
      s%s = ms%to_string()
      s_exp%s = "1.600000000000000000000000000000000000000E+1"
      call tc%assert_eq(s_exp, s, "init_set string")

      call mr%destruct()
      call mm%destruct()
      call mi%destruct()
      call ms%destruct()
    end block

    ! add
    block
      type(mpfr)   :: a, b, c
      type(string) :: s, s_exp

      call a%init()
      call b%init()
      call c%init()

      call a%set("3.14159")
      call b%set("1.00000")
      call add(c, a, b)
      s%s = c%to_string()
      s_exp%s = "4.141589999999999999999999999999999999988E+0"
      call tc%assert_eq(s_exp, s, "add mpfr mpfr")

      call a%set("3.1415926535897932384626433832795")
      call add(c, a, 1.0)
      s%s = c%to_string()
      s_exp%s = "4.141592653589793238462643383279499999991E+0"
      call tc%assert_eq(s_exp, s, "add mpfr real")

      call a%set("3.1415926535897932384626433832795")
      call add(c, a, 1)
      s%s = c%to_string()
      s_exp%s = "4.141592653589793238462643383279499999991E+0"
      call tc%assert_eq(s_exp, s, "add mpfr int")

      call a%destruct()
      call b%destruct()
      call c%destruct()
    end block

    ! sub
    block
      type(mpfr)   :: a, b, c
      type(string) :: s, s_exp

      call a%init()
      call b%init()
      call c%init()

      call a%set("3.14159")
      call b%set("1.00000")
      call sub(c, a, b)
      s%s = c%to_string()
      s_exp%s = "2.141590000000000000000000000000000000000E+0"
      call tc%assert_eq(s_exp, s, "sub mpfr mpfr")

      call a%set("3.1415926535897932384626433832795")
      call sub(c, a, 1.0)
      s%s = c%to_string()
      s_exp%s = "2.141592653589793238462643383279500000002E+0"
      call tc%assert_eq(s_exp, s, "sub mpfr real")

      call a%set("3.1415926535897932384626433832795")
      call sub(c, 1.0, a)
      s%s = c%to_string()
      s_exp%s = "-2.141592653589793238462643383279500000002E+0"
      call tc%assert_eq(s_exp, s, "sub real mpfr")

      call a%set("3.1415926535897932384626433832795")
      call sub(c, a, 1)
      s%s = c%to_string()
      s_exp%s = "2.141592653589793238462643383279500000002E+0"
      call tc%assert_eq(s_exp, s, "sub mpfr int")

      call a%set("3.1415926535897932384626433832795")
      call sub(c, 1, a)
      s%s = c%to_string()
      s_exp%s = "-2.141592653589793238462643383279500000002E+0"
      call tc%assert_eq(s_exp, s, "sub int mpfr")

      call a%destruct()
      call b%destruct()
      call c%destruct()
    end block

    ! mul
    block
      type(mpfr)   :: a, b, c
      type(string) :: s, s_exp

      call a%init()
      call b%init()
      call c%init()

      call a%set("3.1415926535897932384626433832795")
      call b%set("2.7182818284590452353602874713527")
      call mul(c, a, b)
      s%s = c%to_string()
      s_exp%s = "8.539734222673567065463550869546684471768E+0"
      call tc%assert_eq(s_exp, s, "mul mpfr mpfr")

      call a%set("3.1415926535897932384626433832795")
      call mul(c, a, 2.718281828459045)
      s%s = c%to_string()
      s_exp%s = "8.539734222673566611300185395393421040520E+0"
      call tc%assert_eq(s_exp, s, "mul mpfr real")

      call a%set("3.1415926535897932384626433832795")
      call mul(c, a, 2)
      s%s = c%to_string()
      s_exp%s = "6.283185307179586476925286766559000000005E+0"
      call tc%assert_eq(s_exp, s, "mul mpfr int")

      call a%destruct()
      call b%destruct()
      call c%destruct()
    end block

    ! div
    block
      type(mpfr)   :: a, b, c
      type(string) :: s, s_exp

      call a%init()
      call b%init()
      call c%init()

      call a%set("3.1415926535897932384626433832795")
      call b%set("2.7182818284590452353602874713527")
      call div(c, a, b)
      s%s = c%to_string()
      s_exp%s = "1.155727349790921717910093183312679293316E+0"
      call tc%assert_eq(s_exp, s, "div mpfr mpfr")

      call a%set("3.1415926535897932384626433832795")
      call div(c, a, 2.718281828459045)
      s%s = c%to_string()
      s_exp%s = "1.155727349790921779374420885450460118450E+0"
      call tc%assert_eq(s_exp, s, "div mpfr real")

      call a%set("3.1415926535897932384626433832795")
      call div(c, 2.718281828459045, a)
      s%s = c%to_string()
      s_exp%s = "8.652559794322650412014050328172537248468E-1"
      call tc%assert_eq(s_exp, s, "div real mpfr")

      call a%set("3.1415926535897932384626433832795")
      call div(c, a, 2)
      s%s = c%to_string()
      s_exp%s = "1.570796326794896619231321691639750000001E+0"
      call tc%assert_eq(s_exp, s, "div mpfr int")

      call a%set("3.1415926535897932384626433832795")
      call div(c, 2, a)
      s%s = c%to_string()
      s_exp%s = "6.366197723675813430755350534900580325970E-1"
      call tc%assert_eq(s_exp, s, "div int mpfr")

      call a%destruct()
      call b%destruct()
      call c%destruct()
    end block

    ! FIXME: add missing tests

    call tc%finish()
  end subroutine

end module

m4_divert(0)
