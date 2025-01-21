m4_include(util/macro.f90.inc)

module degen_m

  use gauss_m,       only: gauss_legendre, gauss_laguerre
  use gauss_table_m, only: gauss_table
  use mathm_m,       only: roots

  use fermi_m, only: fermi_dirac_integral_1h_reg, inv_fermi_dirac_integral_1h_reg

  implicit none

  private
  public :: DEGEN_DEBUG, DEGEN_TANH_SINH, degen_init, degen_get, fd12, ifd12, get_current, weierstrass

  logical :: DEGEN_DEBUG     = .false.
  logical :: DEGEN_TANH_SINH = .true.

  interface

    subroutine dist_(eta, F, dF1, dF2, dF3, dF4, dF5, dF6)
      !! distribution function
      real,           intent(in)  :: eta
        !! normalized Fermi level
      real, optional, intent(out) :: F
        !! output distribution function
      real, optional, intent(out) :: dF1
        !! output first derivative of F wrt eta
      real, optional, intent(out) :: dF2
        !! output second derivative of F wrt eta
      real, optional, intent(out) :: dF3
        !! output third derivative of F wrt eta
      real, optional, intent(out) :: dF4
        !! output fourth derivative of F wrt eta
      real, optional, intent(out) :: dF5
        !! output fifth derivative of F wrt eta
      real, optional, intent(out) :: dF6
        !! output sixth derivative of F wrt eta
    end subroutine

    subroutine idist_(F, eta, detadF)
      !! inverse distribution function
      real,           intent(in)  :: F
        !! distribution function
      real, optional, intent(out) :: eta
        !! output normalized Fermi level
      real, optional, intent(out) :: detadF
        !! output derivative of eta wrt F
    end subroutine

    subroutine integrand_(x, p, f, dfdx, dfdp)
      real, intent(in)  :: x
      real, intent(in)  :: p(:)
      real, intent(out) :: f
      real, intent(out) :: dfdx
      real, intent(out) :: dfdp(:)
    end subroutine

  end interface

  integer           :: NG = -1
  real, allocatable :: XLEG(:), WLEG(:), XLAG(:), WLAG(:)
  type(gauss_table) :: gtab

  integer,      parameter :: CASE1A  = 1
  integer,      parameter :: CASE1B  = 2
  integer,      parameter :: CASE1C  = 3
  integer,      parameter :: CASE1D  = 4
  integer,      parameter :: CASE1E  = 5
  integer,      parameter :: CASE2A  = 6
  integer,      parameter :: CASE2B  = 7
  character(2), parameter :: CASENAME(7) = ["1a", "1b", "1c", "1d", "1e", "2a", "2b"]
  real,         parameter :: CASETOL = 1e-3

  integer, parameter :: MAX_IT = 50
  real,    parameter :: RTOL   = 2e-14
  real,    parameter :: ATOL   = 1e-16

  real,    parameter :: F_GAMMA  = sqrt(0.125)
  real,    parameter :: DELTA    = 3.0

contains

  subroutine degen_init(ng_)
    integer, intent(in) :: ng_
      !! number of gauss nodes

    NG = ng_

    allocate (XLEG(NG), WLEG(NG), XLAG(NG), WLAG(NG), source = 0.0)

    call gauss_legendre(XLEG, WLEG)
    call gauss_laguerre(XLAG, WLAG)
    call gtab%init(NG, .false.)
  end subroutine

  subroutine degen_get(n, dpot, j, djdn, djddpot)
    real,               intent(in)  :: n(2)
      !! left/right normalized density
    real,               intent(in)  :: dpot
      !! normalized potential drop
    real,               intent(out) :: j
      !! output normalized edge current
    real,               intent(out) :: djdn(2)
      !! output derivatives of j wrt n
    real,               intent(out) :: djddpot
      !! output derivatives of j wrt dpot

    call get_current(fd12, ifd12, n(1), n(2), dpot, j, djdn(1), djdn(2), djddpot)
  end subroutine

  subroutine fd12(eta, F, dF1, dF2, dF3, dF4, dF5, dF6)
    !! fermi-dirac integral for j = 1/2 and derivatives
    real,           intent(in)  :: eta
      !! argument
    real, optional, intent(out) :: F
      !! output distribution function
    real, optional, intent(out) :: dF1
      !! output first derivative of F wrt eta
    real, optional, intent(out) :: dF2
      !! output second derivative of F wrt eta
    real, optional, intent(out) :: dF3
      !! output third derivative of F wrt eta
    real, optional, intent(out) :: dF4
      !! output fourth derivative of F wrt eta
    real, optional, intent(out) :: dF5
      !! output fifth derivative of F wrt eta
    real, optional, intent(out) :: dF6
      !! output sixth derivative of F wrt eta

    real, parameter :: G0 = 1.0 / gamma(1.5)
    real, parameter :: G1 = 1.0 / gamma(0.5)
    real, parameter :: G2 = 1.0 / gamma(-0.5)
    real, parameter :: G3 = 1.0 / gamma(-1.5)
    real, parameter :: G4 = 1.0 / gamma(-2.5)
    real, parameter :: G5 = 1.0 / gamma(-3.5)

    if (present(f  )) F   = fd1h( eta) * G0
    if (present(dF1)) dF1 = fdm1h(eta) * G1
    if (present(dF2)) dF2 = fdm3h(eta) * G2
    if (present(dF3)) dF3 = fdm5h(eta) * G3
    if (present(dF4)) dF4 = fdm7h(eta) * G4
    if (present(dF5)) dF5 = fdm9h(eta) * G5
    if (present(dF6)) dF6 = dfdm9h(eta) * G5
  end subroutine

  subroutine ifd12(F, eta, deta)
    !! inverse of fermi-dirac integral for j = 1/2 and derivative
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 698-707)
    real,           intent(in)  :: F
      !! fermi-dirac integral value (> 0)
    real, optional, intent(out) :: eta
      !! output argument of fermi-dirac integral
    real, optional, intent(out) :: deta
      !! output derivative of eta wrt F

    real, parameter :: C1(10) = [ &
      8.4973821066601840E-01, 1.5637783330562941E+05, 4.8177570589828698E+04, 5.8470721838381196E+03, 3.3539780796721942E+02, &
      7.8441186802991201E+00, 1.1776202905535090E+05,-1.9007269383703679E+04, 1.3762936928453139E+03,-5.4113726984817170E+01 ]
    real, parameter :: C2(17) = [ &
      3.7691787449019803E-01,-4.4356940732931460E-01, 4.8914044731041020E+02, 5.3350726931726194E+03, 2.0169073614044250E+04, &
      3.5247811559551090E+04, 3.0462366861471477E+04, 1.2567903242612896E+04, 2.1318678935739867E+03, 9.3652017208541949E+01, &
      6.5682620764306057E+02, 4.2748283105194159E+03, 1.0555758131015149E+04, 1.2341874209461188E+04, 6.9491885441319710E+03, &
      1.6921965063419400E+03, 1.2922177299158975E+02 ]
    real, parameter :: C3(17) = [ &
      1.0465156933592495E-01,-4.0080827720541695E-01, 1.0198488640664235E+03, 9.4401825500392206E+03, 3.3947661636376244E+04, &
      6.0256728098054278E+04, 5.5243004506305580E+04, 2.4769835480221085E+04, 4.5117728861766827E+03, 2.1143280633615015E+02, &
      3.5050207035358642E+02, 2.5310629620123404E+03, 6.9390985065943923E+03, 9.0054019797239653E+03, 5.6067361299413406E+03, &
      1.4887663456400508E+03, 1.2153702888941258E+02 ]
    real, parameter :: C4(17) = [ &
      2.5090716445082574E-02,-3.3585051328246379E-01, 1.1885877939839949E+04, 1.1322025082517879E+05, 4.0852437388119783E+05, &
      6.9567435748347593E+05, 5.6938991708850558E+05, 2.0643308201368144E+05, 2.7307253567197411E+04, 8.2443082679473071E+02, &
      1.6344049122086119E+03, 1.2218115855188402E+04, 3.2911786995779323E+04, 3.8934696303939934E+04, 2.0038835843822581E+04, &
      3.9494838089779696E+03, 2.1560740489099570E+02 ]
    real, parameter :: C5(17) = [ &
      7.3980341563880635E-03,-3.9387746247592931E-01, 1.1730701119043564E+04, 9.9421745579663359E+04, 3.2770696891070693E+05, &
      5.3042566801656317E+05, 4.3863190051655506E+05, 1.7532285566231585E+05, 2.8701960598881389E+04, 1.2582091446428640E+03, &
      6.3408047038302618E+02, 4.2956315986026584E+03, 1.0868526066891194E+04, 1.2781687199797707E+04, 7.0938073210076054E+03, &
      1.6750641705630003E+03, 1.2575090181775967E+02 ]
    real, parameter :: C6( 7) = [ &
      1.0801341205098402E+03, 1.1281349514482193E+07, 4.2036891115716088E+05, 1.6896947571453611E+03, 6.0880835083129587E+03, &
      2.2144523675946675E+02, 7.1821670869539778E-01 ]
    real, parameter :: FBND( 5) = [ &
      1.3279138832783541E+00, 4.3216142181900556E+00, 1.5103862150597896E+01, 6.0075840912753129E+01, 2.1260003088573106E+02 ]
    real, parameter :: G = gamma(1.5)

    real :: b, b1, n, n1, d, d1, s, s1, t

    if (F < FBND(1)) then
      t = C1(1) * G * F
      if (present(eta) .or. present(deta)) then
        n =       t*(C1(2)+t*(C1(3)+t*(C1( 4)+t*(C1(5)+t*C1(6)))))
        d = C1(7)+t*(C1(8)+t*(C1(9)+t*(C1(10)+t                )))
        b = n / d
        if (present(eta)) eta = log(b)
      end if
      if (present(deta)) then
        n1 = (C1(2)+t*(2*C1(3)+t*(3*C1( 4)+t*(4*C1(5)+t*5*C1(6))))) * G * C1(1)
        d1 = (C1(8)+t*(2*C1(9)+t*(3*C1(10)+t* 4                 ))) * G * C1(1)
        b1 = (n1 - n * d1 / d) / d
        deta = b1 / b
      end if
    elseif (F < FBND(2)) then
      t = C2(1) * G * F + C2(2)
      if (present(eta) .or. present(deta)) then
        n = (C2( 3)+t*(C2( 4)+t*(C2( 5)+t*(C2( 6)+t*(C2( 7)+t*(C2( 8)+t*(C2( 9)+t*C2(10))))))))
        d = (C2(11)+t*(C2(12)+t*(C2(13)+t*(C2(14)+t*(C2(15)+t*(C2(16)+t*(C2(17)+t       )))))))
        if (present(eta)) eta = n / d
      end if
      if (present(deta)) then
        n1 = (C2( 4)+t*(2*C2( 5)+t*(3*C2( 6)+t*(4*C2( 7)+t*(5*C2( 8)+t*(6*C2( 9)+t*7*C2(10))))))) * G * C2(1)
        d1 = (C2(12)+t*(2*C2(13)+t*(3*C2(14)+t*(4*C2(15)+t*(5*C2(16)+t*(6*C2(17)+t*7       )))))) * G * C2(1)
        deta = (n1 - n * d1 / d) / d
      end if
    elseif (F < FBND(3)) then
      t = C3(1) * G * F + C3(2)
      if (present(eta) .or. present(deta)) then
        n = C3( 3)+t*(C3( 4)+t*(C3( 5)+t*(C3( 6)+t*(C3( 7)+t*(C3( 8)+t*(C3( 9)+t*C3(10)))))))
        d = C3(11)+t*(C3(12)+t*(C3(13)+t*(C3(14)+t*(C3(15)+t*(C3(16)+t*(C3(17)+t       ))))))
        if (present(eta)) eta = n / d
      end if
      if (present(deta)) then
        n1 = (C3( 4)+t*(2*C3( 5)+t*(3*C3( 6)+t*(4*C3( 7)+t*(5*C3( 8)+t*(6*C3( 9)+t*7*C3(10))))))) * G * C3(1)
        d1 = (C3(12)+t*(2*C3(13)+t*(3*C3(14)+t*(4*C3(15)+t*(5*C3(16)+t*(6*C3(17)+t*7       )))))) * G * C3(1)
        deta = (n1 - n * d1 / d) / d
      end if
    elseif (F < FBND(4)) then
      t = C4(1) * G * F + C4(2)
      if (present(eta) .or. present(deta)) then
        n = C4( 3)+t*(C4( 4)+t*(C4( 5)+t*(C4( 6)+t*(C4( 7)+t*(C4( 8)+t*(C4( 9)+t*C4(10)))))))
        d = C4(11)+t*(C4(12)+t*(C4(13)+t*(C4(14)+t*(C4(15)+t*(C4(16)+t*(C4(17)+t       ))))))
        if (present(eta)) eta = n / d
      end if
      if (present(deta)) then
        n1 = (C4( 4)+t*(2*C4( 5)+t*(3*C4( 6)+t*(4*C4( 7)+t*(5*C4( 8)+t*(6*C4( 9)+t*7*C4(10))))))) * G * C4(1)
        d1 = (C4(12)+t*(2*C4(13)+t*(3*C4(14)+t*(4*C4(15)+t*(5*C4(16)+t*(6*C4(17)+t*7       )))))) * G * C4(1)
        deta = (n1 - n * d1 / d) / d
      end if
    elseif (F < FBND(5)) then
      t = C5(1) * G * F + C5(2)
      if (present(eta) .or. present(deta)) then
        n = C5( 3)+t*(C5( 4)+t*(C5( 5)+t*(C5( 6)+t*(C5( 7)+t*(C5( 8)+t*(C5( 9)+t*C5(10)))))))
        d = C5(11)+t*(C5(12)+t*(C5(13)+t*(C5(14)+t*(C5(15)+t*(C5(16)+t*(C5(17)+t       ))))))
        if (present(eta)) eta = n / d
      end if
      if (present(deta)) then
        n1 = (C5( 4)+t*(2*C5( 5)+t*(3*C5( 6)+t*(4*C5( 7)+t*(5*C5( 8)+t*(6*C5( 9)+t*7*C5(10))))))) * G * C5(1)
        d1 = (C5(12)+t*(2*C5(13)+t*(3*C5(14)+t*(4*C5(15)+t*(5*C5(16)+t*(6*C5(17)+t*7       )))))) * G * C5(1)
        deta = (n1 - n * d1 / d) / d
      end if
    else
      s = C6(1) * (G * F)**(- 4.0/3.0)
      t = 1 - s
      if (present(eta) .or. present(deta)) then
        n = C6(2)+t*(C6(3)+t*(C6(4)+t))
        d = s*(C6(5)+t*(C6(6)+t*C6(7)))
        b = n / d
        if (present(eta)) eta = sqrt(b)
      end if
      if (present(deta)) then
        s1 = C6(1) * (-4.0/3.0) * (G * F)**(-7.0/3.0) * G
        n1 = - (C6(3)+t*(2*C6(4)+t*3)) * s1
        d1 = s1*(C6(5)+t*(C6(6)+t*C6(7)) - s*(C6(6)+t*2*C6(7)))
        b1 = (n1 - n * d1 / d) / d
        deta = 0.5 / sqrt(b) * b1
      end if
    end if
  end subroutine

  subroutine get_current(dist, idist, n1, n2, dpot, j, djdn1, djdn2, djddpot)
    procedure(dist_)  :: dist
      !! distribution function
    procedure(idist_) :: idist
      !! inverse distribution function
    real, intent(in)  :: n1
      !! normalized density at the left point
    real, intent(in)  :: n2
      !! normalized density at the right point
    real, intent(in)  :: dpot
      !! normalized potential drop
    real, intent(out) :: j
      !! output normalized current density
    real, intent(out) :: djdn1
      !! output derivative of j wrt n1
    real, intent(out) :: djdn2
      !! output derivative of j wrt n2
    real, intent(out) :: djddpot
      !! output derivative of j wrt dpot

    integer :: cs, it
    logical :: flip
    real    :: n1_, n2_, dpot_, eta1, eta2, deta1dn1, deta2dn2, deta, djdeta, tmp
    real    :: jj, jjmin, jjmax, res, dresdjj, djj, err
    real    :: F0, F1, F2, F3, F4, F5, F6

    ! flip edge direction if potential drop is negative
    flip = (dpot < 0)
    if (flip) then
      n1_   = n2
      n2_   = n1
    else
      n1_   = n1
      n2_   = n2
    end if
    dpot_ = abs(dpot)

    ! get normalized Fermi levels
    call idist(n1_, eta = eta1, detadF = deta1dn1)
    call idist(n2_, eta = eta2, detadF = deta2dn2)
    deta = eta2 - eta1

    ! get case
    if (abs(dpot_) > CASETOL) then
      if (deta - dpot_ > dpot_ * CASETOL) then
        cs = CASE1A
      elseif (deta - dpot_ >= - dpot_ * CASETOL) then
        cs = CASE1B
      elseif (deta > CASETOL) then
        cs = CASE1C
      elseif (deta >= - CASETOL) then
        cs = CASE1D
      else
        cs = CASE1E
      end if
    else
      if (abs(deta) > CASETOL) then
        cs = CASE2A
      else
        cs = CASE2B
      end if
    end if

    ! debug info
    if (DEGEN_DEBUG) then
      print *
      print "(A,ES25.16E3)", "dpot = ", dpot_
      print "(A,ES25.16E3)", "eta1 = ", eta1
      print "(A,ES25.16E3)", "eta2 = ", eta2
      print "(A,ES25.16E3)", "deta = ", deta
      print "(2A)",          "case ", CASENAME(cs)
    end if

    ! calculate current based on case
    select case (cs)
    case (CASE1B)
      call case_1b()
    case (CASE1D)
      call case_1d()
    case (CASE2A)
      call case_2a()
    case (CASE2B)
      call case_2b()
    case default
      ! range for jj by mean value theorem
      jjmin = abs(1 - deta / dpot_) * min(n1_, n2_)
      jjmax = abs(1 - deta / dpot_) * max(n1_, n2_)
      if (deta > dpot_) then
        ! flip sign and order
        jj    = jjmin
        jjmin = - jjmax
        jjmax = - jj
      end if

      ! further reduce range if possible
      if (cs == CASE1A) then
        call integrate_dist(dist, etac, eta2 - dpot_, eta2, -1, II, dII, ncalls)
        nquad = nquad + ncalls
        jjmin = max(jjmin, (dpot_ - deta) / II)
      elseif (cs == CASE1C) then
        jjmax = min(jjmax, n1_)
      elseif (cs == CASE1E) then
        jjmin = max(jjmin, n1_)
      end if

      ! initial value: enhanced diffusion approximation
      if (n2_ == n1_) then
        tmp = 1.0 / (n1_ * deta1dn1)
      else
        tmp = dpot_ * log(n2_ / n1_) / deta
      end if
      jj = - n1_ / expm1(-tmp) - n2_ / expm1(tmp)
      if ((jj < jjmin) .or. (jj > jjmax)) jj = 0.5 * (jjmin + jjmax)

      if (DEGEN_DEBUG) then
        print "(A,ES25.16E3)", "jjmin = ", jjmin
        print "(A,ES25.16E3)", "jjmax = ", jjmax
        print "(A,ES25.16E3)", "jj    = ", jj
      end if

      ! Newton iteration
      err = huge(1.0)
      it  = 0
      do while ((it < 6) .or. (err > RTOL * abs(jj) + ATOL / abs(dpot_)))
        it = it + 1
        if (it > MAX_IT) call fatal_error("No convergence after " // int2str(MAX_IT) // " iterations")

        ! evaluate residual and get Newton update
        call residual()

        ! treat singularity of res at interval end-points
        do while (.not. ieee_is_finite(res))
          if (jj - jjmin < 1e-13 * (jjmax - jjmin)) then
            jj    = ieee_next_after(jj, huge(1.0))
            jjmin = jj
          elseif (jjmax - jj < 1e-13 * (jjmax - jjmin)) then
            jj    = ieee_next_after(jj, -huge(1.0))
            jjmax = jj
          else
            call fatal_error("residual not finite")
          end if
          call residual()
        end do

        ! Newton update
        djj = - res / dresdjj
        err = abs(djj)

        ! update bounds (assume monotonic behaviour)
        if (djj > 0) then
          jjmin = jj
        else
          jjmax = jj
        end if

        ! update solution
        jjold = jj
        jj    = jj + djj

        ! bisection
        if ((jj < jjmin) .or. (jj > jjmax) .or. ((jjold == jjmin) .and. (jj == jjmax))) then
          if (DEGEN_DEBUG) print "(A)", "bisection"
          jj = 0.5 * (jjmin + jjmax)
          err = min(err, jjmax - jjmin)
        end if

        if (DEGEN_DEBUG) then
          print "(I6,A,ES25.16E3,A,ES25.16E3,A,ES25.16E3)", it, ": jj = ", jj, "  +/-", err, "  ,", err / jj
        end if
      end do

      ! get current and derivatives with implicit differentiation
      call residual()
      j       = dpot_ * jj
      djdeta  =    - dpot_ * dresdeta  / dresdjj
      djddpot = jj - dpot_ * dresddpot / dresdjj

    end select

    ! derivatives wrt densities
    djdn1 = djdeta(1) * deta1dn1
    djdn2 = djdeta(2) * deta2dn2

    ! reverse edge direction flip (keep dj/ddpot = d(-j)/d(-dpot) unchanged)
    if (flip) then
      j     = - j
      tmp   = djdn1
      djdn1 = - djdn2
      djdn2 = - tmp
    end if

  contains

    subroutine case_1b()
      real          :: I(-5:0), dI(2,-5:0), jc(1:5), djc(2,1:5)
      real(kind=16) :: I16(-5:0), jc16(1:5), djc16(2,1:5), tmp16, dtmp16(2), j16

      do k = -5, -1
        call integrate_dist(dist, eta1, eta2, k, I(k), dI(:,k))
      end do
      I(0) = deta

      I16 = I

      dI(:,0) = [-1.0, 1.0]

      jc16(1)  = I16(0) / I16(-1)
      djc16(:,1) = (dI(:,0) - dI(:,-1) * I16(0) / I16(-1)) / I16(-1)

      tmp16    = I16(-1)**2 - I16(-2)*I16(0)
      dtmp16   = 2 * I16(-1) * dI(:,-1) - dI(:,-2)*I16(0) - I16(-2)*dI(:,0)
      jc16(2)  = tmp16 / I16(-1)**3
      djc16(:,2) = (dtmp16 - 3 * tmp16 * dI(:,-1) / I16(-1)) / I16(-1)**3

      tmp16    = 2*I16(-2)**2*I16(0) - I16(-1)**2*I16(-2) - I16(-1)*I16(-3)*I16(0)
      dtmp16   = 4*I16(-2)*dI(:,-2)*I16(0) + 2*I16(-2)**2*dI(:,0) - 2*I16(-1)*dI(:,-1)*I16(-2) - I16(-1)**2*dI(:,-2) &
        &      - dI(:,-1)*I16(-3)*I16(0) - I16(-1)*dI(:,-3)*I16(0) - I16(-1)*I16(-3)*dI(:,0)
      jc16(3)  = tmp16 / I16(-1)**5
      djc16(:,3) = (dtmp16 - 5 * tmp16 * dI(:,-1) / I16(-1)) / I16(-1)**5

      tmp16    = 2*I16(-1)**2*I16(-2)**2 - I16(-1)**3*I16(-3) - 5*I16(-2)**3*I16(0) - I16(-1)**2*I16(-4)*I16(0) + 5*I16(-1)*I16(-2)*I16(-3)*I16(0)
      dtmp16   = 4*I16(-1)*dI(:,-1)*I16(-2)**2 + 4*I16(-1)**2*I16(-2)*dI(:,-2) - 3*I16(-1)**2*dI(:,-1)*I16(-3) - I16(-1)**3*dI(:,-3) &
        &      - 15*I16(-2)**2*dI(:,-2)*I16(0) - 5*I16(-2)**3*dI(:,0) - 2*I16(-1)*dI(:,-1)*I16(-4)*I16(0) - I16(-1)**2*dI(:,-4)*I16(0) &
        &      - I16(-1)**2*I16(-4)*dI(:,0) + 5*dI(:,-1)*I16(-2)*I16(-3)*I16(0) + 5*I16(-1)*dI(:,-2)*I16(-3)*I16(0) &
        &      + 5*I16(-1)*I16(-2)*dI(:,-3)*I16(0) + 5*I16(-1)*I16(-2)*I16(-3)*dI(:,0)
      jc16(4)  = tmp16 / I16(-1)**7
      djc16(:,4) = (dtmp16 - 7 * tmp16 * dI(:,-1) / I16(-1)) / I16(-1)**7

      tmp16    = 14*I16(0)*I16(-2)**4 - I16(-1)**4*I16(-4) - 5*I16(-1)**2*I16(-2)**3 - I16(0)*I16(-1)**3*I16(-5) &
        &      + 5*I16(-1)**3*I16(-2)*I16(-3) + 3*I16(0)*I16(-1)**2*I16(-3)**2 - 21*I16(0)*I16(-1)*I16(-2)**2*I16(-3) &
        &      + 6*I16(0)*I16(-1)**2*I16(-2)*I16(-4)
      dtmp16   = 14*dI(:,0)*I16(-2)**4 + 56*I16(0)*I16(-2)**3*dI(:,-2) - 4*I16(-1)**3*dI(:,-1)*I16(-4) - I16(-1)**4*dI(:,-4) &
        &      - 10*I16(-1)*dI(:,-1)*I16(-2)**3 - 15*I16(-1)**2*I16(-2)**2*dI(:,-2) - dI(:,0)*I16(-1)**3*I16(-5) &
        &      - 3*I16(0)*I16(-1)**2*dI(:,-1)*I16(-5) - I16(0)*I16(-1)**3*dI(:,-5) + 15*I16(-1)**2*dI(:,-1)*I16(-2)*I16(-3) &
        &      + 5*I16(-1)**3*dI(:,-2)*I16(-3) + 5*I16(-1)**3*I16(-2)*dI(:,-3) + 3*dI(:,0)*I16(-1)**2*I16(-3)**2 &
        &      + 6*I16(0)*I16(-1)*dI(:,-1)*I16(-3)**2 + 6*I16(0)*I16(-1)**2*I16(-3)*dI(:,-3) - 21*dI(:,0)*I16(-1)*I16(-2)**2*I16(-3) &
        &      - 21*I16(0)*dI(:,-1)*I16(-2)**2*I16(-3) - 42*I16(0)*I16(-1)*I16(-2)*dI(:,-2)*I16(-3) - 21*I16(0)*I16(-1)*I16(-2)**2*dI(:,-3) &
        &      + 6*dI(:,0)*I16(-1)**2*I16(-2)*I16(-4) + 12*I16(0)*I16(-1)*dI(:,-1)*I16(-2)*I16(-4) + 6*I16(0)*I16(-1)**2*dI(:,-2)*I16(-4) &
        &      + 6*I16(0)*I16(-1)**2*I16(-2)*dI(:,-4)
      jc16(5)    = tmp16 / I16(-1)**9
      djc16(:,5) = (dtmp16 - 9 * tmp16 * dI(:,-1) / I16(-1)) / I16(-1)**9

      j16     = 0
      djdeta  = 0
      djddpot = 0
      jc  = real(jc16)
      djc = real(djc16)
      do k = 1, 5
        j16     = j16 + jc16(k) * (dpot_ - I16(0))**k
        djdeta  = djdeta + (djc(:,k) * (dpot_ - I(0)) - k * jc(k) * dI(:,0)) * (dpot_ - I(0))**(k-1)
        djddpot = djddpot + jc(k) * k * (dpot_ - I(0))**(k-1)
      end do
      j = real(j16)
    end subroutine

    subroutine case_1d()
      use math_m, only: expm1

      real :: etam, A, B, C, dAdetam, dBdetam, dCdetam, dCddpot, u, dudetam, duddpot, v, dvdetam, dvddpot, w
      real :: jc(0:3), djcdetam(0:3), djcddpot(0:3)

      etam = 0.5 * (eta1 + eta2)
      call dist(etam, F = F0, dF1 = F1, dF2 = F2, dF3 = F3, dF4 = F4)

      A       = F1**2 - F0 * F2
      dAdetam = F1 * F2 - F0 * F3
      B       = F2**2 - F1 * F3
      dBdetam = F2 * F3 - F1 * F4
      C       = dpot_ * F1
      dCdetam = dpot_ * F2
      dCddpot = F1

      jc(0)       = dpot_ * F0
      djcdetam(0) = C
      djcddpot(0) = F0

      u       = C / F0
      dudetam = dCdetam / F0 - C * F1 / F0**2
      duddpot = dCddpot / F0
      v       = exp(u)
      dvdetam = v * dudetam
      dvddpot = v * duddpot
      w       = expm1(u)
      jc(1)       = - 0.5 * C * (v + 1) / w
      djcdetam(1) = - 0.5 * (dCdetam * (v + 1) + C * dvdetam - C * (v + 1) * dvdetam / w) / w
      djcddpot(1) = - 0.5 * (dCddpot * (v + 1) + C * dvddpot - C * (v + 1) * dvddpot / w) / w

      u       = C**2 * F0 * (A + F1**2) + 2 * A * jc(1) * (C**2 - 4 * jc(1) * (F0 + jc(1)))
      dudetam = (2 * C * dCdetam * F0 + C**2 * F1) * (A + F1**2) + C**2 * F0 * (dAdetam + 2 * F1 * F2) &
        &     + 2 * (dAdetam * jc(1) + A * djcdetam(1)) * (C**2 - 4 * jc(1) * (F0 + jc(1))) &
        &     + 2 * A * jc(1) * (2 * C * dCdetam - 4 * djcdetam(1) * (F0 + jc(1)) - 4 * jc(1) * (F1 + djcdetam(1)))
      duddpot = 2 * C * dCddpot * F0 * (A + F1**2) + 2 * A * djcddpot(1) * (C**2 - 4 * jc(1) * (F0 + jc(1))) &
        &     + 2 * A * jc(1) * (2 * C * dCddpot - 4 * djcddpot(1) * F0 - 8 * jc(1) * djcddpot(1))
      v       = 8 * C * F0**2 * F1
      dvdetam = 8 * dCdetam * F0**2 * F1 + 16 * C * F0 * F1**2 + 8 * C * F0**2 * F2
      dvddpot = 8 * dCddpot * F0**2 * F1
      jc(2)       = u / v
      djcdetam(2) = (dudetam - u * dvdetam / v) / v
      djcddpot(2) = (duddpot - u * dvddpot / v) / v

      u       = 192*A**2*jc(1)**5 &
              + 48*F0*jc(1)**4*(8*A**2 - B*F0**2 + 2*A*F0*F2) &
              - 24*jc(1)**3*(2*B*F0**4 + A**2*(3*C**2 - 8*F0**2) - 4*A*F0**3*F2) &
              - 12*C**2*F0*jc(1)**2*(10*A**2 - B*F0**2 + 3*A*F0*F2) &
              - 4*C**2*F0**2*jc(1)*(12*A**2 - F0**2*(3*B + F1*F3) + 6*F0*F2*A) &
              + 3*C**4*(A**2*(F0 + 2*jc(1)) + A*F0*F1**2)
      dudetam = 384*A*dAdetam*jc(1)**5 + 960*A**2*jc(1)**4*djcdetam(1) &
        &     + 48*(F1*jc(1)**4 + 4*F0*jc(1)**3*djcdetam(1)) * (8*A**2 - B*F0**2 + 2*A*F0*F2) &
        &     + 48*F0*jc(1)**4*(16*A*dAdetam - dBdetam*F0**2 - 2*B*F0*F1 + 2*dAdetam*F0*F2 + 2*A*F1*F2 + 2*A*F0*F3) &
        &     - 72*jc(1)**2*djcdetam(1)*(2*B*F0**4 + A**2*(3*C**2 - 8*F0**2) - 4*A*F0**3*F2) &
        &     - 24*jc(1)**3*(2*dBdetam*F0**4 + 8*B*F0**3*F1 + 2*A*dAdetam*(3*C**2 - 8*F0**2) + A**2*(6*C*dCdetam - 16*F0*F1) - 4*dAdetam*F0**3*F2 - 12*A*F0**2*F1*F2 - 4*A*F0**3*F3) &
        &     - 12*(2*C*dCdetam*F0*jc(1)**2 + C**2*F1*jc(1)**2 + 2*C**2*F0*jc(1)*djcdetam(1))*(10*A**2 - B*F0**2 + 3*A*F0*F2) &
        &     - 12*C**2*F0*jc(1)**2*(20*A*dAdetam - dBdetam*F0**2 - 2*B*F0*F1 + 3*(dAdetam*F0*F2 + A*F1*F2 + A*F0*F3)) &
        &     - 4*C*F0*(2*dCdetam*F0*jc(1) + 2*C*F1*jc(1) + C*F0*djcdetam(1))*(12*A**2 - F0**2*(3*B + F1*F3) + 6*F0*F2*A) &
        &     - 4*C**2*F0**2*jc(1)*(24*A*dAdetam - 2*F0*F1*(3*B + F1*F3) - F0**2*(3*dBdetam + F2*F3 + F1*F4) + 6*F1*F2*A + 6*F0*F3*A + 6*F0*F2*dAdetam) &
        &     + 12*C**3*dCdetam*(A**2*(F0 + 2*jc(1)) + A*F0*F1**2) &
        &     + 3*C**4*(2*A*dAdetam*(F0 + 2*jc(1)) + A**2*(F1 + 2*djcdetam(1)) + dAdetam*F0*F1**2 + A*F1**3 + 2*A*F0*F1*F2)
      duddpot = 960*A**2*jc(1)**4*djcddpot(1) &
        &     + 192*F0*jc(1)**3*djcddpot(1)*(8*A**2 - B*F0**2 + 2*A*F0*F2) &
        &     - 72*jc(1)**2*djcddpot(1)*(2*B*F0**4 + A**2*(3*C**2 - 8*F0**2) - 4*A*F0**3*F2) &
        &     - 144*jc(1)**3*A**2*C*dCddpot &
        &     - 24*F0*C*jc(1)*(dCddpot*jc(1) + C*djcddpot(1))*(10*A**2 - B*F0**2 + 3*A*F0*F2) &
        &     - 4*F0**2*(2*C*dCddpot*jc(1) + C**2*djcddpot(1))*(12*A**2 - F0**2*(3*B + F1*F3) + 6*F0*F2*A) &
        &     + 12*C**3*dCddpot*(A**2*(F0 + 2*jc(1)) + A*F0*F1**2) &
        &     + 6*C**4*A**2*djcddpot(1)
      v       = 96 * C**2 * F0**4 * F1**2
      dvdetam = 192 * C * dCdetam * F0**4 * F1**2 + 384 * C**2 * F0**3 * F1**3 + 192 * C**2 * F0**4 * F1 * F2
      dvddpot = 192 * C * dCddpot * F0**4 * F1**2
      jc(3)       = u / v
      djcdetam(3) = (dudetam - u * dvdetam / v) / v
      djcddpot(3) = (duddpot - u * dvddpot / v) / v

      j = jc(0) + jc(1) * deta + jc(2) * deta**2 + jc(3) * deta**3
      djdeta(1) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) - jc(1) - 2 * jc(2) * deta - 3 * jc(3) * deta**2
      djdeta(2) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) + jc(1) + 2 * jc(2) * deta + 3 * jc(3) * deta**2
      djddpot   = djcddpot(0) + djcddpot(1) * deta + djcddpot(2) * deta**2 + djcddpot(3) * deta**3
    end subroutine

    subroutine case_2a()
      integer :: ncalls
      real    :: I(5), dI(2,5), jc(0:4), djc(2,0:4), dtmp(2)

      do k = 1, 5
        call integrate_dist(dist, etac, eta1, eta2, k, I(k), dI(:,k), ncalls)
        nquad = nquad + ncalls
      end do
      jc(0)    = - I(1)
      djc(:,0) = - dI(:,1)

      jc(1) = I(2) / I(1)
      djc(:,1) = (dI(:,2) - I(2) * dI(:,1) / I(1)) / I(1)

      tmp      = I(2)**2 - I(1)*I(3)
      dtmp     = 2*I(2)*dI(:,2) - dI(:,1)*I(3) - I(1)*dI(:,3)
      jc(2)    = tmp / I(1)**3
      djc(:,2) = (dtmp - 3 * tmp * dI(:,1) / I(1)) / I(1)**3

      tmp      = I(1)**2*I(4) + 2*I(2)**3 - 3*I(1)*I(2)*I(3)
      dtmp     = 2*I(1)*dI(:,1)*I(4) + I(1)**2*dI(:,4) + 6*I(2)**2*dI(:,2) - 3*dI(:,1)*I(2)*I(3) - 3*I(1)*dI(:,2)*I(3) - 3*I(1)*I(2)*dI(:,3)
      jc(3)    = tmp / I(1)**5
      djc(:,3) = (dtmp - 5 * tmp * dI(:,1) / I(1)) / I(1)**5

      tmp      = 5*I(2)**4 - I(1)**3*I(5) + 2*I(1)**2*I(3)**2 - 10*I(1)*I(2)**2*I(3) + 4*I(1)**2*I(2)*I(4)
      dtmp     = 20*I(2)**3*dI(:,2) - 3*I(1)**2*dI(:,1)*I(5) - I(1)**3*dI(:,5) + 4*I(1)*dI(:,1)*I(3)**2 &
        &      + 4*I(1)**2*I(3)*dI(:,3) - 10*dI(:,1)*I(2)**2*I(3) - 20*I(1)*I(2)*dI(:,2)*I(3) - 10*I(1)*I(2)**2*dI(:,3) &
        &      + 8*I(1)*dI(:,1)*I(2)*I(4) + 4*I(1)**2*dI(:,2)*I(4) + 4*I(1)**2*I(2)*dI(:,4)
      jc(4)    = tmp / I(1)**7
      djc(:,4) = (dtmp - 7 * tmp * dI(:,1) / I(1)) / I(1)**7

      j       = jc(0)
      djdeta  = djc(:,0)
      djddpot = 0
      do k = 1, 4
        j       = j + jc(k) * dpot_**k
        djdeta  = djdeta + djc(:,k) * dpot_**k
        djddpot = djddpot + jc(k) * k * dpot_**(k-1)
      end do
    end subroutine

    subroutine case_2b()
      real :: etam, dtmpdetam, dtmpddeta, dtmpddpot, djdetam, djddeta

      etam = 0.5 * (eta1 + eta2)
      call dist(etam, F = F0, dF1 = F1, dF2 = F2, dF3 = F3)

      tmp = F0 - F1**2 / (12 * F0) * deta * dpot_ + F2 / 24 * deta**2

      dtmpdetam = F1 - (2 * F2 - F1**2 / F0) * F1 / (12 * F0) * deta * dpot_ + F3 / 24 * deta**2
      dtmpddeta = - F1**2 / (12 * F0) * dpot_ + F2 / 12 * deta
      dtmpddpot = - F1**2 / (12 * F0) * deta

      j       = (dpot_ - deta) * tmp
      djdetam = (dpot_ - deta) * dtmpdetam
      djddeta = - tmp + (dpot_ - deta) * dtmpddeta
      djdeta  = djdetam * [0.5, 0.5] + djddeta * [-1.0, 1.0]
      djddpot = tmp + (dpot_ - deta) * dtmpddpot
    end subroutine

    subroutine residual()
      integer :: k
      real    :: dresdjj1(1)

      if (DEGEN_TANH_SINH) then
        call quad(integrand_u, eta1, eta2, [jj], res, dresdeta(1), dresdeta(2), dresdjj1)
        res       = res - dpot_
        dresdjj   = dresdjj1(1)
        dresddpot = -1
      else
        if (cs == CASE1A) then
          kmax = 3
          eta_min(1) = - huge(1.0)
          eta_max(1) = 0.0
          eta_min(2) = 0.0
          eta_max(2) = 2.0
          eta_min(3) = 2.0
          eta_max(3) = huge(1.0)
          d_eta_min  = 0.0
          d_eta_max  = 0.0
        else
          ! pole position
          call idist(jj, eta = eta0, detadF = deta0djj)
          F0 = jj
          F1 = 1.0 / deta0djj

          kmax = 7
          d_eta_min = 0.0
          d_eta_max = 0.0

          eta_min(1) = - huge(1.0)
          eta_max(1) = eta0 - DELTA
          if (eta_max(1) > 0) then
            eta_max(1) = 0
          else
            d_eta_max(1) = 1
          end if

          eta_min(2) = eta0 + DELTA
          if (eta_min(2) < 0) then
            eta_min(2) = 0
          else
            d_eta_min(2) = 1
          end if
          eta_max(2) = eta0 + DELTA
          if (eta_max(2) < 2) then
            eta_max(2) = 2
          else
            d_eta_max(2) = 1
          end if

          eta_min(3) = eta0 + DELTA
          if (eta_min(3) < 2) then
            eta_min(3) = 2
          else
            d_eta_min(3) = 1
          end if
          eta_max(3) = huge(1.0)

          eta_min(4)   = eta0 - DELTA
          d_eta_min(4) = 1
          eta_max(4)   = eta0 + DELTA
          d_eta_max(4) = 1

          eta_min(5)   = eta0 + DELTA
          d_eta_min(5) = 1
          eta_max(5) = eta0 + DELTA
          if (eta_max(5) < 0) then
            eta_max(5) = 0
          else
            d_eta_max(5) = 1
          end if

          eta_min(6) = eta0 - DELTA
          if (eta_min(6) > 0) then
            eta_min(6) = 0
          else
            d_eta_min(6) = 1
          end if
          eta_max(6) = eta0 - DELTA
          if (eta_max(6) > 2) then
            eta_max(6) = 2
          else
            d_eta_max(6) = 1
          end if

          eta_min(7) = eta0 - DELTA
          if (eta_min(7) > 2) then
            eta_min(7) = 2
          else
            d_eta_min(7) = 1
          end if
          eta_max(7) = eta0 - DELTA
          d_eta_max(7) = 1
        end if

        res = 0
        do k = 1, kmax
          if (eta_min(k) == eta_max(k)) cycle

          ! merge bounds with eta1, eta2
          mask1 = ([eta1, eta2] <= eta_min(k))
          mask2 = ([eta1, eta2] <= eta_max(k))
          bnd = merge([eta_min(k), eta_min(k)], merge([eta1, eta2], [eta_max(k), eta_max(k)], mask2), mask1)
          if (bnd(1) == bnd(2)) cycle
          dbnddeta = merge([0.0, 0.0], merge([1.0, 1.0], [0.0, 0.0], mask2), mask1)

          ! principal part (analytic integration)
          select case (k)
          case (1) ! blake for eta <= 0
            tmp = jj * F_GAMMA - 1
            res = res + log((exp(bnd(1)) * tmp + jj) / (exp(bnd(2)) * tmp + jj)) / tmp
            ....
          case (2) ! blake for 0 <= eta <= 2
          case (3) ! blake for eta >= 2
            ! find roots of polynomial
            omega = roots([0.125*PI**2, - jj, 0, 0])


            ! function roots(p) result(rts)
            !   !! Find complex roots of polynomial
            !   real, intent(in) :: p(:)
            !     !! coefficients of polynomial in reduced form: f(x) = p(1) + p(2) * x + p(3) * x**2 + ... + x**n
            !   complex          :: rts(size(p))
            !     !! return roots of polynomial
          case (4)
          case (5)
          case (6)
          case (7)
          end select

          ! rest (numerical integration)

        end do
      end if
    end subroutine

  end subroutine

  subroutine integrate_dist(dist, eta1, eta2, k, I, dIdeta)
    !! int_eta1^eta2 dist(eta)^k deta using Gauss-Laguerre, Gauss-Legendre
    procedure(dist_)     :: dist
      !! distribution function (e.g. fermi-dirac integral)
    real,    intent(in)  :: eta1
      !! lower integration bound
    real,    intent(in)  :: eta2
      !! upper integration bound
    integer, intent(in)  :: k
      !! exponent
    real,    intent(out) :: I
      !! output integral over dist^k
    real,    intent(out) :: dIdeta(2)
      !! output derivative of I wrt eta1, eta2

    logical :: mask(2)
    real    :: bnd(2), dbnd(2), I1, dI1dbnd(2), dum(0), dum2(0), dum3

    ! special case: k = 0
    if (k == 0) then
      I      = eta2 - eta1
      dIdeta = [-1.0, 1.0]
      return
    end if

    if (DEGEN_TANH_SINH) then
      call quad(dist_k, eta1, eta2, dum, I, dIdeta(1), dIdeta(2), dum2)
      return
    end if

    ! reset output
    I      = 0
    dIdeta = 0
    ncalls = 0

    ! eta < 0: Gauss with Exponential weight
    mask = [eta1, eta2] < 0.0
    bnd  = merge([eta1, eta2], [0.0, 0.0], mask)
    if (bnd(1) /= bnd(2)) then
      dbnd = merge([1.0, 1.0], [0.0, 0.0], mask)
      call integrate_gauss_exp(dist_k, bnd(1), bnd(2), dum, real(k), I1, dI1dbnd(1), dI1dbnd(2), dum2, dum3)
      I      = I + I1
      dIdeta = dIdeta + dI1dbnd * dbnd
    end if

    ! eta > 0: Gauss-Legendre
    bnd = merge([0.0, 0.0], [eta1, eta2], mask)
    if (bnd(1) /= bnd(2)) then
      dbnd = merge([0.0, 0.0], [1.0, 1.0], mask)
      call integrate_gauss_legendre(dist_k, bnd(1), bnd(2), dum, I1, dI1dbnd(1), dI1dbnd(2), dum2)
      I      = I + I1
      dIdeta = dIdeta + dI1dbnd * dbnd
    end if

  contains

    subroutine dist_k(x, p, f, dfdx, dfdp)
      real, intent(in)  :: x
      real, intent(in)  :: p(:)
      real, intent(out) :: f
      real, intent(out) :: dfdx
      real, intent(out) :: dfdp(:)

      real :: FF, dFF1

      m4_ignore(p)
      m4_ignore(dfdp)

      call dist(x, F = FF, dF1 = dFF1)
      f    = FF**k
      dfdx = k * FF**(k-1) * dFF1
    end subroutine

  end subroutine

  subroutine integrate_gauss_legendre(integrand, a, b, p, I, dIda, dIdb, dIdp)
    !! integrate function using Gauss-Legendre quadrature
    procedure(integrand_) :: integrand
      !! integrand
    real,    intent(in)   :: a
      !! lower bound
    real,    intent(in)   :: b
      !! upper bound
    real,    intent(in)   :: p(:)
      !! parameters
    real,    intent(out)  :: I
      !! output value of integral
    real,    intent(out)  :: dIda
      !! output derivative of I wrt lower bound
    real,    intent(out)  :: dIdb
      !! output derivative of I wrt upper bound
    real,    intent(out)  :: dIdp(:)
      !! output derivative of I wrt parameters

    integer :: l
    real    :: x, w, f, dfdx, dfdp(size(p))

    ! reset output
    I    = 0
    dIda = 0
    dIdb = 0
    dIdp = 0

    ! quadrature
    do l = 1, NG
      x = 0.5 * (a + b) + 0.5 * (b - a) * XLEG(l)
      w = 0.5 * (b - a) * WLEG(l)

      call integrand(x, p, f, dfdx, dfdp)
      I = I + w * f
      dIda = dIda + 0.5 * (w * dfdx * (1 - XLEG(l)) - WLEG(l) * f)
      dIdb = dIdb + 0.5 * (w * dfdx * (1 + XLEG(l)) + WLEG(l) * f)
      dIdp = dIdp + w * dfdp
    end do
  end subroutine

  subroutine integrate_gauss_laguerre(integrand, a, b, p, k, I, dIda, dIdb, dIdp)
    !! integrate function using two Gauss-Laguerre quadratures
    procedure(integrand_) :: integrand
      !! integrand
    real,    intent(in)   :: a
      !! lower bound
    real,    intent(in)   :: b
      !! upper bound
    real,    intent(in)   :: p(:)
      !! parameter
    real,    intent(in)   :: k
      !! exponent (can be positive or negative)
    real,    intent(out)  :: I
      !! output value of integral
    real,    intent(out)  :: dIda
      !! output derivative of I wrt lower bound
    real,    intent(out)  :: dIdb
      !! output derivative of I wrt upper bound
    real,    intent(out)  :: dIdp(:)
      !! output derivative of I wrt parameter

    integer :: l
    real    :: x, w, f, dfdx, dfdp(size(p))

    m4_assert(k /= 0.0)

    ! reset output
    I    = 0
    dIda = 0
    dIdb = 0
    dIdp = 0

    ! quadrature
    do l = 1, NG
      w = WLAG(l) * exp(XLAG(l)) / k

      x    = a - XLAG(l) / k
      call integrand(x, p, f, dfdx, dfdp)
      I    = I    - w * f
      dIda = dIda - w * dfdx
      dIdp = dIdp - w * dfdp

      x = b - XLAG(l) / k
      call integrand(x, p, f, dfdx, dfdp)
      I    = I    + w * f
      dIdb = dIdb + w * dfdx
      dIdp = dIdp + w * dfdp
    end do
  end subroutine

  subroutine integrate_gauss_exp(integrand, a, b, p, k, I, dIda, dIdb, dIdp, dIdk)
    !! integrate function using custom exponentially weighted quadrature rule
    procedure(integrand_) :: integrand
      !! integrand
    real,    intent(in)   :: a
      !! lower bound
    real,    intent(in)   :: b
      !! upper bound
    real,    intent(in)   :: p(:)
      !! parameters
    real,    intent(in)   :: k
      !! exponent (can be positive or negative)
    real,    intent(out)  :: I
      !! output value of integral
    real,    intent(out)  :: dIda
      !! output derivative of I wrt lower bound
    real,    intent(out)  :: dIdb
      !! output derivative of I wrt upper bound
    real,    intent(out)  :: dIdp(:)
      !! output derivative of I wrt parameter
    real,    intent(out)  :: dIdk
      !! output derivative of I wrt exponent

    integer :: l
    real :: alpha, t(NGEXP), dtdalpha(NGEXP), w(NGEXP), dwdalpha(NGEXP)
    real :: x, dxda, dxdb, dxdk, v, dvda, dvdb, dvdk, f, dfdx, dfdp(size(p))

    ! reset output
    I    = 0
    dIda = 0
    dIdb = 0
    dIdp = 0
    dIdk = 0

    alpha = k * (b - a)
    call gtab%get(alpha, t, dtdalpha, w, dwdalpha)

    do l = 1, NGEXP
      x    = a + (b - a) * t(l)
      dxda = 1 - t(l) - alpha * dtdalpha(l)
      dxdb =     t(l) + alpha * dtdalpha(l)
      dxdk = dtdalpha(l) * (b - a)**2
      v    =   w(l) * (b - a)
      dvda = - w(l) - alpha * dwdalpha(l)
      dvdb =   w(l) + alpha * dwdalpha(l)
      dvdk = dwdalpha(l) * (b - a)**2
      call integrand(x, p, f, dfdx, dfdp)

      I = I + v * f
      dIda = dIda + dvda * f + v * dfdx * dxda
      dIdb = dIdb + dvdb * f + v * dfdx * dxdb
      dIdp = dIdp + v * dfdp
      dIdk = dIdk + dvdk * f + v * dfdx * dxdk
    end do
  end subroutine

end module
