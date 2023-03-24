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
  public fermi_dirac_integral_one_half, inverse_fermi_dirac_integral_one_half

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

  subroutine fermi_dirac_integral_one_half(x, y, dydx)
    !! fermi-dirac integral for j = 1/2
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 708-729)
    real, intent(in)  :: x
      !! argument
    real, intent(out) :: y
      !! value of fermi-dirac integral
    real, intent(out) :: dydx
      !! output derivative of y wrt x (computed with complex step derivative)

    real, parameter :: A(11) = [7.38905609893065023, 0.886226925452758014, 19894.4553386951666, 4509.64329955948557, &
      303.461789035142376, 5.7574879114754736, 0.00275088986849762610, 63493.915041308052, 19070.1178243603945, &
      1962.19362141235102, 79.250704958640158]
    real, parameter :: B(16) = [149.462587768865243, 22.8125889885050154, -0.629256395534285422, 9.08120441515995244, &
      3.35357478401835299, -0.473677696915555805, -0.467190913556185953, -0.0880610317272330793, 0.00262208080491572673, &
      269.94660938022644, 343.6419926336247, 323.9049470901941, 218.89170769294024, 102.31331350098315, 36.319337289702664, &
      8.3317401231389461]
    real, parameter :: C(17) = [71652.717119215557, 134954.734070223743, 153693.833350315645, 123247.280745703400, &
      72886.293647930726, 32081.2499422362952, 10210.9967337762918, 2152.71110381320778, 232.906588165205042, &
      105667.839854298798, 31946.0752989314444, 71158.788776422211, 15650.8990138187414, 13521.8033657783433, &
      1646.98258283527892, 618.90691969249409, -3.36319591755394735]
    real, parameter :: D(15) = [23744.8706993314289, 68257.8589855623002, 89327.4467683334597, 62766.3415600442563, &
      20093.6622609901994, -2213.89084119777949, -3901.66057267577389, 948.642895944858861, 9488.61972919565851, &
      12514.8125526953073, 9903.44088207450946, 2138.15420910334305, -528.394863730838233, -661.033633995449691, &
      -51.4481470250962337]
    real, parameter :: E(15) = [311337.452661582536, 1.11267074416648198e6, 1.75638628895671735e6, 1.59630855803772449e6, &
      910818.935456183774, 326492.733550701245, 65507.2624972852908, 4809.45649527286889, 39721.6641625089685, &
      86424.7529107662431, 88163.7255252151780, 50615.7363511157353, 17334.9774805008209, 2712.13170809042550, &
      82.2205828354629102]
    real, parameter :: F(14) = [7.26870063003059784e6, 2.79049734854776025e7, 4.42791767759742390e7, 3.63735017512363365e7, &
      1.55766342463679795e7, 2.97469357085299505e6, 154516.447031598403, 340542.544360209743, 805021.468647620047, &
      759088.235455002605, 304686.671371640343, 39289.4061400542309, 582.426138126398363, 11.2728194581586028 ]
    real, parameter :: G(13) = [4.81449797541963104e6, 1.85162850713127602e7, 2.77630967522574435e7, 2.03275937688070624e7, &
      7.41578871589369361e6, 1.21193113596189034e6, 63211.9545144644852, 80492.7765975237449, 189328.678152654840, &
      151155.890651482570, 48146.3242253837259, 5407.08878394180588, 112.195044410775577]
    real, parameter :: H(5) = [8109.79390744477921, 342.069867454704106, 1.07141702293504595, 6569.98472532829094, 280.706465851683809 ]

    real, parameter :: K = 0.999999999999999877
      !! 1.0 - epsilon(1.0)

    real, parameter :: dx = 1e-32
      !! complex step size

    complex :: xx, ex, yy, s, t, w

    xx = cmplx(x, dx)

    if (x < -2) then
      ex = exp(xx)
      t  = ex * A(1)
      yy  = ex * (A(2)-ex*(A(3)+t*(A(4)+t*(A(5)+t*(A(6)+t*A(7)))))/(A(8)+t*(A(9)+t*(A(10)+t*(A(11)+t)))))
    elseif (x < 0) then
      s  = - 0.5 * xx
      t  = 1 - s
      yy  = (B(1)+t*(B(2)+t*(B(3)+t*(B(4)+t*(B(5)+t*(B(6)+t*(B(7)+t*(B(8)-t*B(9)))))))))/(B(10)+s*(B(11)+s*(B(12)+s*(B(13)+s*(B(14)+s*(B(15)+s*(B(16)+s)))))))
    elseif (x < 2) then
      t = 0.5 * xx
      yy = (C(1)+t*(C(2)+t*(C(3)+t*(C(4)+t*(C(5)+t*(C(6)+t*(C(7)+t*(C(8)+t*C(9)))))))))/(C(10)+t*(C(11)+t*(C(12)+t*(C(13)+t*(C(14)+t*(C(15)+t*(C(16)+t*(C(17)+t))))))))
    elseif (x < 5) then
      t = (xx - 2) / 3
      yy = (D(1)+t*(D(2)+t*(D(3)+t*(D(4)+t*(D(5)+t*(D(6)+t*(D(7)-t*D(8))))))))/(D(9)+t*(D(10)+t*(D(11)+t*(D(12)+t*(D(13)+t*(D(14)+t*(D(15)+t)))))))
    elseif (x < 10) then
      t = 0.2 * xx - 1
      yy = (E(1)+t*(E(2)+t*(E(3)+t*(E(4)+t*(E(5)+t*(E(6)+t*(E(7)+t*E(8))))))))/(E(9)+t*(E(10)+t*(E(11)+t*(E(12)+t*(E(13)+t*(E(14)+t*(E(15)-t)))))))*K
    elseif (x < 20) then
      t = 0.1 * xx - 1
      yy = (F(1)+t*(F(2)+t*(F(3)+t*(F(4)+t*(F(5)+t*(F(6)+t*F(7)))))))/(F(8)+t*(F(9)+t*(F(10)+t*(F(11)+t*(F(12)+t*(F(13)+t*(F(14)-t)))))))*K
    elseif (x < 40) then
      t = 0.05 * xx - 1
      yy = (G(1)+t*(G(2)+t*(G(3)+t*(G(4)+t*(G(5)+t*(G(6)+t*G(7)))))))/(G(8)+t*(G(9)+t*(G(10)+t*(G(11)+t*(G(12)+t*(G(13)-t))))))*K
    else
      w = 1 / (xx * xx)
      s = 1 - 1600 * w
      yy = 2 * xx * sqrt(xx) / 3 * (1 + w*(H(1)+s*(H(2)+s*H(3)))/(H(4)+s*(H(5)+s)))
    end if

    yy = yy / gamma(1.5)

    y    = real(yy)
    dydx = imag(yy) / dx
  end subroutine

  subroutine inverse_fermi_dirac_integral_one_half(y, x, dxdy)
    !! inverse of fermi-dirac integral for j = 1/2
    !! Reference: Fukushima, T. (2015, App. Math. Comp., 259, 698-707)
    real, intent(in) :: y
      !! fermi-dirac integral value
    real, intent(out) :: x
      !! output argument of fermi-dirac integral
    real, intent(out) :: dxdy
      !! output derivative of x wrt y (computed with complex step derivative)

    real, parameter :: A(10) = [0.849738210666018375, 156377.8333056294, 48177.5705898287, 5847.07218383812, &
      335.3978079672194, 7.84411868029912, 117762.02905535089, -19007.26938370368, 1376.2936928453140, -54.11372698481717 ]
    real, parameter :: B(17) = [0.376917874490198033, -0.443569407329314587, 489.140447310410217, 5335.07269317261966, &
      20169.0736140442509, 35247.8115595510907, 30462.3668614714761, 12567.9032426128967, 2131.86789357398657, &
      93.6520172085419439, 656.826207643060606, 4274.82831051941605, 10555.7581310151498, 12341.8742094611883, &
      6949.18854413197094, 1692.19650634194002, 129.221772991589751]
    real, parameter :: C(17) = [0.104651569335924949, -0.400808277205416960, 1019.84886406642351, 9440.18255003922075, &
      33947.6616363762463, 60256.7280980542786, 55243.0045063055787, 24769.8354802210838, 4511.77288617668292, &
      211.432806336150141, 350.502070353586442, 2531.06296201234050, 6939.09850659439245, 9005.40197972396592, &
      5606.73612994134056, 1488.76634564005075, 121.537028889412581]
    real, parameter :: D(17) = [0.0250907164450825724, -0.335850513282463787, 11885.8779398399498, 113220.250825178799, &
      408524.373881197840, 695674.357483475952, 569389.917088505552, 206433.082013681440, 27307.2535671974100, &
      824.430826794730740, 1634.40491220861182, 12218.1158551884025, 32911.7869957793233, 38934.6963039399331, &
      20038.8358438225823, 3949.48380897796954, 215.607404890995706]
    real, parameter :: E(17) = [0.00739803415638806339, -0.393877462475929313, 11730.7011190435638, 99421.7455796633651, &
      327706.968910706902, 530425.668016563224, 438631.900516555072, 175322.855662315845, 28701.9605988813884, &
      1258.20914464286403, 634.080470383026173, 4295.63159860265838, 10868.5260668911946, 12781.6871997977069, &
      7093.80732100760563, 1675.06417056300026, 125.750901817759662]
    real, parameter :: F(7) = [1080.13412050984017, 1.12813495144821933e7, 420368.911157160874, 1689.69475714536117, &
      6088.08350831295857, 221.445236759466761, 0.718216708695397737]
    real, parameter :: G(5) = [1.17683303804380831, 3.82993088157949761, 13.3854493161866553, 53.2408277860982205, 188.411871723022843] / gamma(1.5)

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

end module


