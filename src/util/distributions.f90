m4_include(macro.f90.inc)

module distributions_m

  use, intrinsic :: ieee_arithmetic
  use math_m, only: expm1, PI
  m4_ifdef({m4_quadpack},{use quadpack_m, only: quadpack_int})

  implicit none

  private
  public bose_einstein,        d_bose_einstein
  public fermi_dirac,          d_fermi_dirac
  public maxwell_boltzmann,    d_maxwell_boltzmann
  public fermi_dirac_integral_approx
  m4_ifdef({m4_quadpack},{public fermi_dirac_integral})
  public fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h, fermi_dirac_integral_m1h

contains

  real function maxwell_boltzmann(E) result(f)
    !! maxwell-boltzmann distribution
    real, intent(in) :: E
      !! energy

    f = exp(-E)
  end function

  real function d_maxwell_boltzmann(E) result(dfdE)
    !! derivative of maxwell-boltzmann distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(-E)
  end function

  real function bose_einstein(E) result(f)
    !! bose-einstein distribution
    real, intent(in) :: E
      !! energy

    f = 1 / expm1(E)
  end function

  real function d_bose_einstein(E) result(dfdE)
    !! derivative of bose-einstein distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(E) / expm1(E)**2
  end function

  real function fermi_dirac(E) result(f)
    !! fermi-dirac distribution
    real, intent(in) :: E
      !! energy

    f = 1 / (exp(E) + 1)
  end function

  real function d_fermi_dirac(E) result(dfdE)
    !! derivative of fermi-dirac distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(E) / (exp(E) + 1)**2
  end function

  subroutine fermi_dirac_integral_approx(x, y, dydx)
    !! Approximation of Fermi-Dirac integral F_/2
    !!  wiki: https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
    real,           intent(in)  :: x
    real,           intent(out) :: y
    real, optional, intent(out) :: dydx
      !! derivative

    if (x < 1.3) then
      y = 1/(exp(-x) + 0.27)
      if (present(dydx)) dydx = exp(-x)/(exp(-x) + 0.27)**2
    else
      y = 4/(3*sqrt(PI)) * (x*x + PI*PI/6)**0.75
      if (present(dydx)) dydx =  2*x/sqrt(PI) / (x*x + PI*PI/6)**0.25
    end if
  end subroutine

  m4_ifdef({m4_quadpack},{recursive subroutine fermi_dirac_integral(x, y, j, atol, rtol, dydx)
    !! Fermi-Dirac integral
    !!  wiki: https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
    !! y = y(x) = int_0^\infty \frac{t}{1+\exp(t-x)} \dd{t}

    use quadpack_m, only: quadpack_int

    real,           intent(in)  :: x
    real,           intent(out) :: y
    real, optional, intent(in)  :: j
      !! index j. default: 1/2
    real, optional, intent(in)  :: atol
      !! optional absolute tolerance. default: 1e-10
    real, optional, intent(in)  :: rtol
      !! optional relative tolerance. default: 1e-6
    real, optional, intent(out) :: dydx
      !! derivative

    real :: atol_, rtol_, j_

    j_ = 0.5
    if (present(j)) j_ = j
    atol_ = 1e-10
    if (present(atol)) atol_ = atol
    rtol_ = 1e-6
    if (present(rtol)) rtol_ = rtol

    call quadpack_int(f, 0.0, ieee_value(0.0, ieee_positive_inf), atol_, rtol_, y)
    y = y / gamma(j_+1)

    ! derivative: d(F_j(x))/dx = F_{j-1}(x)
    if (present(dydx)) call fermi_dirac_integral(x, dydx, j=j_-1, atol=atol, rtol=rtol)

  contains

    real function f(t) result(y)
      !! integrand: f(t) = f(t; x) = \frac{t^j}{1+\exp(t-x)}
      real, intent(in) :: t
        !! integration variable

      y = t**j_ / (1 + exp(t-x))
    end function

  end subroutine})

  subroutine fermi_dirac_integral_1h(x, y, dydx)
    !! fermi-dirac integral for j = 1/2
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 708-729)
    real, intent(in)  :: x
      !! argument
    real, intent(out) :: y
      !! value of fermi-dirac integral
    real, intent(out) :: dydx
      !! output derivative of y wrt x (computed with complex step derivative)

    real, parameter :: A(11) = [ &
      7.3890560989306504E+00, 8.8622692545275805E-01, 1.9894455338695167E+04, 4.5096432995594860E+03, 3.0346178903514237E+02, &
      5.7574879114754740E+00, 2.7508898684976261E-03, 6.3493915041308050E+04, 1.9070117824360394E+04, 1.9621936214123509E+03, &
      7.9250704958640156E+01 ]
    real, parameter :: B(16) = [ &
      1.4946258776886523E+02, 2.2812588988505016E+01,-6.2925639553428547E-01, 9.0812044151599522E+00, 3.3535747840183530E+00, &
     -4.7367769691555578E-01,-4.6719091355618597E-01,-8.8061031727233077E-02, 2.6220808049157267E-03, 2.6994660938022645E+02, &
      3.4364199263362468E+02, 3.2390494709019413E+02, 2.1889170769294023E+02, 1.0231331350098316E+02, 3.6319337289702666E+01, &
      8.3317401231389461E+00 ]
    real, parameter :: C(17) = [ &
      7.1652717119215551E+04, 1.3495473407022376E+05, 1.5369383335031563E+05, 1.2324728074570341E+05, 7.2886293647930725E+04, &
      3.2081249942236296E+04, 1.0210996733776292E+04, 2.1527111038132080E+03, 2.3290658816520505E+02, 1.0566783985429880E+05, &
      3.1946075298931446E+04, 7.1158788776422211E+04, 1.5650899013818742E+04, 1.3521803365778344E+04, 1.6469825828352789E+03, &
      6.1890691969249406E+02,-3.3631959175539472E+00 ]
    real, parameter :: D(15) = [ &
      2.3744870699331430E+04, 6.8257858985562299E+04, 8.9327446768333466E+04, 6.2766341560044253E+04, 2.0093662260990201E+04, &
     -2.2138908411977795E+03,-3.9016605726757739E+03, 9.4864289594485888E+02, 9.4886197291956578E+03, 1.2514812552695308E+04, &
      9.9034408820745102E+03, 2.1381542091033430E+03,-5.2839486373083821E+02,-6.6103363399544969E+02,-5.1448147025096233E+01 ]
    real, parameter :: E(15) = [ &
      3.1133745266158256E+05, 1.1126707441664820E+06, 1.7563862889567174E+06, 1.5963085580377246E+06, 9.1081893545618374E+05, &
      3.2649273355070123E+05, 6.5507262497285294E+04, 4.8094564952728688E+03, 3.9721664162508969E+04, 8.6424752910766241E+04, &
      8.8163725525215181E+04, 5.0615736351115738E+04, 1.7334977480500820E+04, 2.7121317080904255E+03, 8.2220582835462906E+01 ]
    real, parameter :: F(14) = [ &
      7.2687006300305976E+06, 2.7904973485477604E+07, 4.4279176775974236E+07, 3.6373501751236334E+07, 1.5576634246367980E+07, &
      2.9746935708529949E+06, 1.5451644703159839E+05, 3.4054254436020972E+05, 8.0502146864762006E+05, 7.5908823545500264E+05, &
      3.0468667137164035E+05, 3.9289406140054234E+04, 5.8242613812639831E+02, 1.1272819458158603E+01 ]
    real, parameter :: G(13) = [ &
      4.8144979754196312E+06, 1.8516285071312759E+07, 2.7763096752257444E+07, 2.0327593768807061E+07, 7.4157887158936933E+06, &
      1.2119311359618905E+06, 6.3211954514464487E+04, 8.0492776597523742E+04, 1.8932867815265484E+05, 1.5115589065148256E+05, &
      4.8146324225383723E+04, 5.4070887839418056E+03, 1.1219504441077558E+02 ]
    real, parameter :: H( 5) = [ &
      8.1097939074447795E+03, 3.4206986745470408E+02, 1.0714170229350459E+00, 6.5699847253282906E+03, 2.8070646585168379E+02 ]

    real, parameter :: K = 0.999999999999999877
      !! 1.0 - epsilon(1.0)

    real, parameter :: dx = 1e-32
      !! complex step size

    complex :: xx, ex, yy, s, t, w

    xx = cmplx(x, dx)

    if (x < -2) then
      ex = exp(xx)
      t  = ex * A(1)
      yy = ex * (A(2)-ex*(A(3)+t*(A(4)+t*(A(5)+t*(A(6)+t*A(7)))))/(A(8)+t*(A(9)+t*(A(10)+t*(A(11)+t)))))
    elseif (x < 0) then
      s  = - 0.5 * xx
      t  = 1 - s
      yy = (B(1)+t*(B(2)+t*(B(3)+t*(B(4)+t*(B(5)+t*(B(6)+t*(B(7)+t*(B(8)-t*B(9)))))))))/(B(10)+s*(B(11)+s*(B(12)+s*(B(13)+s*(B(14)+s*(B(15)+s*(B(16)+s)))))))
    elseif (x < 2) then
      t  = 0.5 * xx
      yy = (C(1)+t*(C(2)+t*(C(3)+t*(C(4)+t*(C(5)+t*(C(6)+t*(C(7)+t*(C(8)+t*C(9)))))))))/(C(10)+t*(C(11)+t*(C(12)+t*(C(13)+t*(C(14)+t*(C(15)+t*(C(16)+t*(C(17)+t))))))))
    elseif (x < 5) then
      t  = (xx - 2) / 3
      yy = (D(1)+t*(D(2)+t*(D(3)+t*(D(4)+t*(D(5)+t*(D(6)+t*(D(7)-t*D(8))))))))/(D(9)+t*(D(10)+t*(D(11)+t*(D(12)+t*(D(13)+t*(D(14)+t*(D(15)+t)))))))
    elseif (x < 10) then
      t  = 0.2 * xx - 1
      yy = (E(1)+t*(E(2)+t*(E(3)+t*(E(4)+t*(E(5)+t*(E(6)+t*(E(7)+t*E(8))))))))/(E(9)+t*(E(10)+t*(E(11)+t*(E(12)+t*(E(13)+t*(E(14)+t*(E(15)-t)))))))*K
    elseif (x < 20) then
      t  = 0.1 * xx - 1
      yy = (F(1)+t*(F(2)+t*(F(3)+t*(F(4)+t*(F(5)+t*(F(6)+t*F(7)))))))/(F(8)+t*(F(9)+t*(F(10)+t*(F(11)+t*(F(12)+t*(F(13)+t*(F(14)-t)))))))*K
    elseif (x < 40) then
      t  = 0.05 * xx - 1
      yy = (G(1)+t*(G(2)+t*(G(3)+t*(G(4)+t*(G(5)+t*(G(6)+t*G(7)))))))/(G(8)+t*(G(9)+t*(G(10)+t*(G(11)+t*(G(12)+t*(G(13)-t))))))*K
    else
      w  = 1 / (xx * xx)
      s  = 1 - 1600 * w
      yy = 2 * xx * sqrt(xx) / 3 * (1 + w*(H(1)+s*(H(2)+s*H(3)))/(H(4)+s*(H(5)+s)))
    end if

    yy = yy / gamma(1.5)

    y    = real(yy)
    dydx = imag(yy) / dx
  end subroutine

  subroutine inv_fermi_dirac_integral_1h(y, x, dxdy)
    !! inverse of fermi-dirac integral for j = 1/2
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 698-707)
    real, intent(in) :: y
      !! fermi-dirac integral value
    real, intent(out) :: x
      !! output argument of fermi-dirac integral
    real, intent(out) :: dxdy
      !! output derivative of x wrt y (computed with complex step derivative)

    real, parameter :: A(10) = [ &
      8.4973821066601840E-01, 1.5637783330562941E+05, 4.8177570589828698E+04, 5.8470721838381196E+03, 3.3539780796721942E+02, &
      7.8441186802991201E+00, 1.1776202905535090E+05,-1.9007269383703679E+04, 1.3762936928453139E+03,-5.4113726984817170E+01 ]
    real, parameter :: B(17) = [ &
      3.7691787449019803E-01,-4.4356940732931460E-01, 4.8914044731041020E+02, 5.3350726931726194E+03, 2.0169073614044250E+04, &
      3.5247811559551090E+04, 3.0462366861471477E+04, 1.2567903242612896E+04, 2.1318678935739867E+03, 9.3652017208541949E+01, &
      6.5682620764306057E+02, 4.2748283105194159E+03, 1.0555758131015149E+04, 1.2341874209461188E+04, 6.9491885441319710E+03, &
      1.6921965063419400E+03, 1.2922177299158975E+02 ]
    real, parameter :: C(17) = [ &
      1.0465156933592495E-01,-4.0080827720541695E-01, 1.0198488640664235E+03, 9.4401825500392206E+03, 3.3947661636376244E+04, &
      6.0256728098054278E+04, 5.5243004506305580E+04, 2.4769835480221085E+04, 4.5117728861766827E+03, 2.1143280633615015E+02, &
      3.5050207035358642E+02, 2.5310629620123404E+03, 6.9390985065943923E+03, 9.0054019797239653E+03, 5.6067361299413406E+03, &
      1.4887663456400508E+03, 1.2153702888941258E+02 ]
    real, parameter :: D(17) = [ &
      2.5090716445082574E-02,-3.3585051328246379E-01, 1.1885877939839949E+04, 1.1322025082517879E+05, 4.0852437388119783E+05, &
      6.9567435748347593E+05, 5.6938991708850558E+05, 2.0643308201368144E+05, 2.7307253567197411E+04, 8.2443082679473071E+02, &
      1.6344049122086119E+03, 1.2218115855188402E+04, 3.2911786995779323E+04, 3.8934696303939934E+04, 2.0038835843822581E+04, &
      3.9494838089779696E+03, 2.1560740489099570E+02 ]
    real, parameter :: E(17) = [ &
      7.3980341563880635E-03,-3.9387746247592931E-01, 1.1730701119043564E+04, 9.9421745579663359E+04, 3.2770696891070693E+05, &
      5.3042566801656317E+05, 4.3863190051655506E+05, 1.7532285566231585E+05, 2.8701960598881389E+04, 1.2582091446428640E+03, &
      6.3408047038302618E+02, 4.2956315986026584E+03, 1.0868526066891194E+04, 1.2781687199797707E+04, 7.0938073210076054E+03, &
      1.6750641705630003E+03, 1.2575090181775967E+02 ]
    real, parameter :: F( 7) = [ &
      1.0801341205098402E+03, 1.1281349514482193E+07, 4.2036891115716088E+05, 1.6896947571453611E+03, 6.0880835083129587E+03, &
      2.2144523675946675E+02, 7.1821670869539778E-01 ]
    real, parameter :: G( 5) = [ &
      1.3279138832783541E+00, 4.3216142181900556E+00, 1.5103862150597896E+01, 6.0075840912753129E+01, 2.1260003088573106E+02 ]

    real    :: dy
    complex :: xx, yy, s, t, v, w, z

    dy = y * 1e-32

    yy = cmplx(y, dy)
    yy = yy * gamma(1.5)

    if (y < G(1)) then
      t  = yy * A(1)
      z  = t*(A(2)+t*(A(3)+t*(A(4)+t*(A(5)+t*A(6)))))/(A(7)+t*(A(8)+t*(A(9)+t*(A(10)+t))))
      xx = log(z)
    elseif (y < G(2)) then
      t  =  B(1) * yy + B(2)
      xx = (B(3)+t*(B(4)+t*(B(5)+t*(B(6)+t*(B(7)+t*(B(8)+t*(B(9)+t*B(10))))))))/(B(11)+t*(B(12)+t*(B(13)+t*(B(14)+t*(B(15)+t*(B(16)+t*(B(17)+t)))))))
    elseif (y < G(3)) then
      t  = C(1) * yy + C(2)
      xx = (C(3)+t*(C(4)+t*(C(5)+t*(C(6)+t*(C(7)+t*(C(8)+t*(C(9)+t*C(10))))))))/(C(11)+t*(C(12)+t*(C(13)+t*(C(14)+t*(C(15)+t*(C(16)+t*(C(17)+t)))))))
    elseif (y < G(4)) then
      t  = D(1) * yy + D(2)
      xx = (D(3)+t*(D(4)+t*(D(5)+t*(D(6)+t*(D(7)+t*(D(8)+t*(D(9)+t*D(10))))))))/(D(11)+t*(D(12)+t*(D(13)+t*(D(14)+t*(D(15)+t*(D(16)+t*(D(17)+t)))))))
    elseif (y < G(5)) then
      t  = E(1) * yy + E(2)
      xx = (E(3)+t*(E(4)+t*(E(5)+t*(E(6)+t*(E(7)+t*(E(8)+t*(E(9)+t*E(10))))))))/(E(11)+t*(E(12)+t*(E(13)+t*(E(14)+t*(E(15)+t*(E(16)+t*(E(17)+t)))))))
    else
      v  = yy**(- 4.0 / 3.0)
      s  = F(1) * v
      t  = 1 - s
      w  = (F(2)+t*(F(3)+t*(F(4)+t)))/(s*(F(5)+t*(F(6)+t*F(7))))
      xx = sqrt(w)
    end if

    x    = real(xx)
    dxdy = imag(xx) / dy
  end subroutine

  subroutine fermi_dirac_integral_m1h(x, y, dydx)
    !! fermi-dirac integral for j = -1/2
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 708-729)
    real, intent(in)  :: x
      !! argument
    real, intent(out) :: y
      !! value of fermi-dirac integral
    real, intent(out) :: dydx
      !! output derivative of y wrt x (computed with complex step derivative)

    real, parameter :: A(11) = [ &
      7.3890560989306504E+00, 1.7724538509055161E+00, 4.0641453751028443E+04, 9.3957080940846445E+03, 6.4996168315267300E+02, &
      1.2797229580475896E+01, 1.5386435076758546E-03, 3.2427188476529293E+04, 1.1079920566127479E+04, 1.3229662700147885E+03, &
      6.3738361029333468E+01 ]
    real, parameter :: B(16) = [ &
      2.7277009213193270E+02, 3.0884565384468285E+01,-6.4353763238036610E+00, 1.4874747309821787E+01, 4.8692886284214261E+00, &
     -1.5326583455067366E+00,-1.0269889831559749E+00,-1.7768682092860594E-01, 3.7714132550924644E-03, 2.9307537818766787E+02, &
      3.0581816268627080E+02, 2.9996239544929762E+02, 2.0764083408749426E+02, 9.2038480318185179E+01, 3.7016491411279119E+01, &
      7.8850095027142055E+00 ]
    real, parameter :: C(15) = [ &
      3.5315036056824306E+03, 6.0775339658420035E+03, 6.1997700433981327E+03, 4.4127870191956763E+03, 2.2522734309281091E+03, &
      8.1184098649224086E+02, 1.9183640105363713E+02, 2.3288183895918380E+01, 3.2938370258479627E+03, 1.5289747402978910E+03, &
      2.5684856281498605E+03, 9.2564264653555824E+02, 5.7423248354035991E+02, 1.3280385932066727E+02, 2.9844716655210213E+01 ]
    real, parameter :: D(15) = [ &
      4.0607075340411825E+03, 1.0812729133305276E+04, 1.3897564948224259E+04, 1.0628474985274002E+04, 5.1077067019067899E+03, &
      1.5408433012600337E+03, 2.8445272011297033E+02, 2.9521441735848416E+01, 1.5645819561263354E+03, 2.8257517227785042E+03, &
      3.1891606616998156E+03, 1.9550397906903256E+03, 8.2800033369181472E+02, 1.8149811108951837E+02, 3.2035285779480375E+01 ]
    real, parameter :: E(14) = [ &
      1.1984171902955750E+03, 3.2635145455490865E+03, 3.8749758847137650E+03, 2.6231306031719982E+03, 1.1004135563712123E+03, &
      2.6746953249050358E+02, 2.5420767181271835E+01, 3.8988775423455579E-01, 2.7340795779255700E+02, 5.9591831895205860E+02, &
      6.0520245226166082E+02, 3.4318330273561998E+02, 1.2218762201569572E+02, 2.0901635907985593E+01 ]
    real, parameter :: F(15) = [ &
      9.4460016943523769E+03, 3.6843444847402861E+04, 6.3710111541992621E+04, 6.2985219736107480E+04, 3.7634523139570090E+04, &
      1.2810989862780776E+04, 1.9815689613892096E+03, 8.1493017189766761E+01, 1.5000469781013367E+03, 5.0869138105279408E+03, &
      7.7300159374762188E+03, 6.6408337623936059E+03, 3.3389959030082641E+03, 8.6049904388680295E+02, 7.8856582418692668E+01 ]
    real, parameter :: G(15) = [ &
      2.2977965785536722E+04, 1.2341661681388778E+05, 2.6115376517235511E+05, 2.7461889451409579E+05, 1.4971071838992485E+05, &
      4.0129337170018458E+04, 4.4704649588141510E+03, 1.3268434683100298E+02, 2.5716884252533569E+03, 1.2521498229077535E+04, &
      2.3268157432505534E+04, 2.0477232011975815E+04, 8.7265257796226815E+03, 1.6474289689676991E+03, 1.0647527514207663E+02 ]
    real, parameter :: H(6) = [ &
      4.1123351671200997E-01, 1.1098041003408895E-03, 1.1368929899017368E-05, 2.5693179067943678E-07, 9.9789778675544616E-09, &
      8.6766769879110857E-10 ]

    real, parameter :: dx = 1e-32
      !! complex step size

    complex :: xx, ex, yy, s, t, w

    xx = cmplx(x, dx)

    if(x < -2) then
      ex = exp(xx)
      t  = ex * A(1)
      yy = ex*(A(2)-ex*(A(3)+t*(A(4)+t*(A(5)+t*(A(6)+t*A(7)))))/(A(8)+t*(A(9)+t*(A(10)+t*(A(11)+t)))))
    elseif(x < 0) then
      s  = -0.5 * xx
      t  = 1.0 - s
      yy = (B(1)+t*(B(2)+t*(B(3)+t*(B(4)+t*(B(5)+t*(B(6)+t*(B(7)+t*(B(8)-t*B(9)))))))))/(B(10)+s*(B(11)+s*(B(12)+s*(B(13)+s*(B(14)+s*(B(15)+s*(B(16)+s)))))))
    elseif(x < 2) then
      t  = 0.5 * xx
      yy = (C(1)+t*(C(2)+t*(C(3)+t*(C(4)+t*(C(5)+t*(C(6)+t*(C(7)+t*C(8))))))))/(C(9)+t*(C(10)+t*(C(11)+t*(C(12)+t*(C(13)+t*(C(14)+t*(C(15)+t)))))))
    elseif(x < 5) then
      t  = (xx - 2.0) / 3.0
      yy = (D(1)+t*(D(2)+t*(D(3)+t*(D(4)+t*(D(5)+t*(D(6)+t*(D(7)+t*D(8))))))))/(D(9)+t*(D(10)+t*(D(11)+t*(D(12)+t*(D(13)+t*(D(14)+t*(D(15)+t)))))))
    elseif(x < 10) then
      t  = 0.2 * xx - 1.0
      yy = (E(1) +t*(E(2) +t*(E(3)+t*(E(4)+t*(E(5)+t*(E(6)+t*(E(7)+t*E(8))))))))/(E(9)+t*(E(10)+t*(E(11)+t*(E(12)+t*(E(13)+t*(E(14)+t))))))
    elseif(x < 20) then
      t  = 0.1 * xx - 1.0
      yy = (F(1)+t*(F(2)+t*(F(3)+t*(F(4)+t*(F(5)+t*(F(6)+t*(F(7)+t*F(8))))))))/(F(9)+t*(F(10)+t*(F(11)+t*(F(12)+t*(F(13)+t*(F(14)+t*(F(15)+t)))))))
    elseif(x < 40) then
      t  = 0.05 * xx - 1.0
      yy = (G(1)+t*(G(2)+t*(G(3)+t*(G(4)+t*(G(5)+t*(G(6)+t*(G(7)+t*G(8))))))))/(G(9)+t*(G(10)+t*(G(11)+t*(G(12)+t*(G(13)+t*(G(14)+t*(G(15)+t)))))))
    else
      w  = 1.0 / (xx * xx)
      t  = 1600.0 * w
      yy = sqrt(xx) * 2.0 * (1.0 - w*(H(1)+t*(H(2)+t*(H(3)+t*(H(4)+t*(H(5)+t*H(6)))))))
    end if

    yy = yy / gamma(0.5)

    y    = real(yy)
    dydx = imag(yy) / dx
  end subroutine

end module


