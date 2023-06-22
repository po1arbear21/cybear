m4_include(util/macro.f90.inc)

module current_density_m

  use density_m,        only: density
  use device_params_m,  only: device_params
  use distributions_m,  only: fermi_dirac_integral_1h, fermi_dirac_integral_m1h, inv_fermi_dirac_integral_1h
  use dual_m
  use equation_m,       only: equation
  use error_m,          only: assert_failed, program_error
  use grid_m,           only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m, only: DIR_NAME
  use high_precision_m
  use ieee_arithmetic,  only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_value
  use jacobian_m,       only: jacobian
  use math_m,           only: ber, dberdx, log1p
  use mobility_m,       only: mobility
  use potential_m,      only: potential
  use radau5_m,         only: ode_options, ode_result, radau5
  use semiconductor_m,  only: CR_NAME, CR_CHARGE
  use stencil_m,        only: dirichlet_stencil, near_neighb_stencil
  use variable_m,       only: variable_real

  implicit none

  private
  public current_density, calc_current_density

  type, extends(variable_real) :: current_density
    !! electron/hole current density

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_current_density
    !! calculate current density by drift-diffusion model

    type(device_params),   pointer :: par   => null()
    type(potential),       pointer :: pot   => null()
    type(density),         pointer :: dens  => null()
    type(current_density), pointer :: cdens => null()
    type(mobility),        pointer :: mob   => null()

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer :: jaco_pot  => null()
    type(jacobian), pointer :: jaco_dens => null()
    type(jacobian), pointer :: jaco_mob  => null()
  contains
    procedure :: init => calc_current_density_init
    procedure :: eval => calc_current_density_eval

    procedure, private :: get_curr       => calc_current_density_get_curr
    procedure, private :: get_curr_sg    => calc_current_density_get_curr_sg
    procedure, private :: get_curr_degen => calc_current_density_get_curr_degen
  end type

contains

  subroutine current_density_init(this, par, ci, idx_dir)
    !! initialize current density
    class(current_density), intent(out) :: this
    type(device_params),    intent(in)  :: par
      !! device parameters
    integer,                intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    integer,                intent(in)  :: idx_dir
      !! edge direction

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init(CR_NAME(ci)//"cdens"//DIR_NAME(idx_dir), "1/cm^2/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
    this%ci = ci

    ! get pointer to data
    select case (par%g%idx_dim)
    case (1)
      p1 => this%data%get_ptr1()
      this%x1 => p1%data
    case (2)
      p2 => this%data%get_ptr2()
      this%x2 => p2%data
    case (3)
      p3 => this%data%get_ptr3()
      this%x3 => p3%data
    case default
      call program_error("Maximal 3 dimensions allowed")
    end select
  end subroutine

  subroutine calc_current_density_init(this, par, pot, dens, cdens, mob)
    !! initialize drift-diffusion equation
    class(calc_current_density),   intent(out) :: this
    type(device_params),   target, intent(in)  :: par
      !! device parameters
    type(potential),       target, intent(in)  :: pot
      !! potential variable
    type(density),         target, intent(in)  :: dens
      !! density variable
    type(current_density), target, intent(in)  :: cdens
      !! current density variable
    type(mobility),        target, intent(in)  :: mob
      !! mobility variable

    integer :: idx_dir, ci, iprov

    idx_dir = cdens%idx_dir
    ci      = cdens%ci

    ! init base
    call this%equation_init("calc_"//cdens%name)
    this%par   => par
    this%pot   => pot
    this%dens  => dens
    this%cdens => cdens
    this%mob   => mob

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_nn%init(par%g, IDX_EDGE, idx_dir, IDX_VERTEX, 0)

    ! provide current density
    iprov = this%provide(cdens, par%transport(IDX_EDGE, idx_dir))

    ! depend on potential
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on density
    this%jaco_dens => this%init_jaco(iprov, this%depend(dens, par%transport(IDX_VERTEX, 0)), st = [this%st_nn%get_ptr()])

    ! depend on mobility
    this%jaco_mob => this%init_jaco(iprov, this%depend(mob, par%transport(IDX_EDGE, idx_dir)), st = [this%st_dir%get_ptr()])

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    !! evaluate drift-diffusion equation
    class(calc_current_density), intent(inout) :: this

    integer               :: i, idx_dir, idx_dim
    real                  :: pot(2), dens(2), j, djdpot(2), djddens(2), djdmob, len, mob
    integer, allocatable  :: idx(:), idx1(:), idx2(:)
    logical               :: status

    idx_dim = this%par%g%idx_dim
    idx_dir = this%cdens%idx_dir
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim))

    ! loop over transport edges in parallel
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(i,pot,dens,j,djdpot,djddens,djdmob,len,mob,idx,idx1,idx2,status) &
    !$omp shared(this,idx_dir,idx_dim)
    do i = 1, this%par%transport(IDX_EDGE, idx_dir)%n
      idx = this%par%transport(IDX_EDGE, idx_dir)%get_idx(i)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
      call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

      ! parameters
      len     = this%par%g%get_len(idx1, idx_dir)
      pot( 1) = this%pot%get( idx1)
      pot( 2) = this%pot%get( idx2)
      dens(1) = this%dens%get(idx1)
      dens(2) = this%dens%get(idx2)
      mob     = this%mob%get( idx1)

      ! get current along edge
      call this%get_curr(this%par%smc%degen, len, pot, dens, mob, j, djdpot, djddens, djdmob)

      ! set current density + derivatives
      call this%cdens%set(idx, j)
      call this%jaco_pot%set( idx, idx1, djdpot(1))
      call this%jaco_pot%set( idx, idx2, djdpot(2))
      call this%jaco_dens%set(idx, idx1, djddens(1))
      call this%jaco_dens%set(idx, idx2, djddens(2))
      call this%jaco_mob%set( idx, idx,  djdmob)
    end do
    !$omp end parallel do
  end subroutine

  subroutine calc_current_density_get_curr(this, degen, len, pot, dens, mob, j, djdpot, djddens, djdmob)
    !! get current along edge
    class(calc_current_density), intent(in)  :: this
    logical,                     intent(in)  :: degen
      !! degeneracy flag (true: FD statistics, false: MB statistics)
    real,                        intent(in)  :: len
      !! edge length
    real,                        intent(in)  :: pot(2)
      !! potential at edge endpoints
    real,                        intent(in)  :: dens(2)
      !! density at edge endpoints
    real,                        intent(in)  :: mob
      !! mobility
    real,                        intent(out) :: j
      !! output current density
    real,                        intent(out) :: djdpot(2)
      !! output derivatives of j wrt pot
    real,                        intent(out) :: djddens(2)
      !! output derivatives of j wrt dens
    real,                        intent(out) :: djdmob
      !! output derivatives of j wrt mob

    integer :: ci
    real    :: ch, dpot, edos, n(2), djdn(2), djddpot

    ! abbreviations
    ci   = this%cdens%ci
    ch   = CR_CHARGE(ci)
    edos = this%par%smc%edos(ci)

    ! normalize
    dpot = - ch * (pot(2) - pot(1))
    n    = dens / edos

    if (degen) then
      ! generalized Scharfetter-Gummel
      call this%get_curr_degen(n, dpot, j, djdn, djddpot)
    else
      ! Scharfetter-Gummel
      call this%get_curr_sg(n, dpot, j, djdn, djddpot)
    end if

    ! denormalize
    j       = j * mob * edos / len
    djdpot  = - ch * [-1.0, 1.0] * djddpot * mob * edos / len
    djddens = djdn * mob / len
    djdmob  = j / mob
  end subroutine

  subroutine calc_current_density_get_curr_sg(this, n, dpot, j, djdn, djddpot)
    use math_m, only: ber

    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real :: B1, B2, dB1, dB2

    m4_ignore(this)

    B1  = ber(    dpot)
    B2  = ber(   -dpot)
    dB1 = dberdx( dpot)
    dB2 = dberdx(-dpot)

    j       = B2 * n(1) - B1 * n(2)
    djdn    = [B2, -B1]
    djddpot = - dB2 * n(1) - dB1 * n(2)
  end subroutine

  recursive subroutine calc_current_density_get_curr_degen(this, n, dpot, j, djdn, djddpot)
    class(calc_current_density), intent(in)  :: this
    real,                        intent(in)  :: n(2)
    real,                        intent(in)  :: dpot
    real,                        intent(out) :: j
    real,                        intent(out) :: djdn(2)
    real,                        intent(out) :: djddpot

    real, parameter :: etaF = -16, alpha = sqrt(0.125)

    integer           :: dir
    logical           :: small_eta
    real              :: eta(2), deta, detadn(2), jmin, jmax, jsgn, djdeta(2), nc, t
    type(ode_options) :: opt

    ! use symmetry to confine dpot to <= 0
    if (dpot > 0) then
      call this%get_curr_degen([n(2), n(1)], -dpot, j, djdn, djddpot)
      j    = - j
      djdn = - [djdn(2), djdn(1)]
      return
    end if

    ! get eta with F12(eta) = n
    call inv_fermi_dirac_integral_1h(n(1), eta(1), detadn(1))
    call inv_fermi_dirac_integral_1h(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

    ! small eta flag
    small_eta = all(eta <= etaF)

    if (small_eta) then
      ! special case: 1/F12(eta) ~ exp(-eta) + sqrt(1/8) -> Bernoulli iteration
      call this%get_curr_sg(n, dpot, j, djdn, djddpot)
      call newton_iteration(j, djdeta, djddpot)
      djdn = djdeta * detadn
    else
      ! get jmin, jmax by slope (detadx must be equal to deta for some x in [0, 1])
      jmin = min(n(1), n(2)) * abs(dpot - deta)
      jmax = max(n(1), n(2)) * abs(dpot - deta)
      jsgn = sign(1.0, dpot - deta)
      if (jsgn < 0) then
        ! correct sign and swap jmin, jmax
        t    = jmin
        jmin = - jmax
        jmax = - t
      end if

      ! determine shooting direction and reduce j range further if possible
      if (eta(2) > eta(1)) then
        dir  = 1
        jmax = min(jmax, dpot * n(2))
      elseif (eta(2) >= eta(1) + dpot) then
        dir  = 0
        jmin = max(jmin, dpot * n(2))
        jmax = min(jmax, 0.0)
      else
        dir  = -1
        call fermi_dirac_integral_1h(eta(1) + 0.5 * dpot, nc, t)
        jmin = max(jmin, 0.0)
        jmax = min(jmax, abs(dpot - deta) * nc)
      end if

      if (abs(deta) < 1e-6) then
        ! eta approximately constant
        block
          use math_m, only: ber, dberdx
          real :: Fm1h, dFm1h, B, dB, dBdn
          call fermi_dirac_integral_m1h(eta(1), Fm1h, dFm1h)

          B    = ber(   dpot * Fm1h / n(1))
          dB   = dberdx(dpot * Fm1h / n(1))
          dBdn = dB * dpot / n(1) * dFm1h * detadn(1) - dB * dpot * Fm1h / n(1)**2

          j  = n(1) * (dpot - deta * B)

          ! j  = n(1) * (dpot + (eta(1) - eta(2)) * B)
          djddpot  = n(1) - deta * dB * Fm1h
          djdn(1) = (dpot - deta * B) + n(1) * (detadn(1) * B + deta * dBdn)
          djdn(2) = - n(1) * detadn(2) * B
        end block
      elseif (abs(deta - dpot) < 1e-6) then
        ! eta changes approximately linear
        block
          real, parameter :: xk(21) = [-1.0, &
            &                          -0.9825722966045480282345, &
            &                          -0.9419762969597455342961, &
            &                          -0.8792947553235904644512, &
            &                          -0.7960019260777124047443, &
            &                          -0.6940510260622232326273, &
            &                          -0.575831960261830686927, &
            &                          -0.4441157832790021011945, &
            &                          -0.3019898565087648872754, &
            &                          -0.1527855158021854660064, &
            &                          0.0, &
            &                          0.1527855158021854660064, &
            &                          0.3019898565087648872754, &
            &                          0.444115783279002101195, &
            &                          0.575831960261830686927, &
            &                          0.6940510260622232326273, &
            &                          0.7960019260777124047443, &
            &                          0.8792947553235904644512, &
            &                          0.9419762969597455342961, &
            &                          0.9825722966045480282345, &
            &                          1.0]
          real, parameter :: wk(21) = [0.004761904761904761904762, &
            &                          0.0291848400985054586095, &
            &                          0.0518431690008496250727, &
            &                          0.0732739181850741442525, &
            &                          0.0929854679578860653011, &
            &                          0.110517083219123335267, &
            &                          0.1254581211908689480152, &
            &                          0.1374584628600413435809, &
            &                          0.1462368624479774592673, &
            &                          0.151587575111681384453, &
            &                          0.1533851903321749485516, &
            &                          0.1515875751116813844533, &
            &                          0.1462368624479774592673, &
            &                          0.1374584628600413435809, &
            &                          0.1254581211908689480152, &
            &                          0.110517083219123335267, &
            &                          0.092985467957886065301, &
            &                          0.0732739181850741442525, &
            &                          0.0518431690008496250727, &
            &                          0.0291848400985054586095, &
            &                          0.004761904761904761904762]

          integer :: ii
          real    :: navg, dnavg(2), ee, dee(2), ff, dff

          navg = 0
          dnavg = 0
          do ii = 1, size(xk)
            ee = 0.5 * (eta(1) + eta(2)) + 0.5 * (eta(2) - eta(1)) * xk(ii)
            dee = 0.5 + [-0.5, 0.5] * xk(ii)
            call fermi_dirac_integral_1h(ee, ff, dff)
            navg  =  navg + 0.5 * wk(ii) / ff
            dnavg = dnavg - 0.5 * wk(ii) / ff**2 * dff * dee
          end do
          dnavg = - dnavg / navg**2
          navg  = 1.0 / navg

          j       = navg * (dpot - deta)
          djddpot = navg
          djdn(1) = (dnavg(1) * (dpot - deta) + navg) * detadn(1)
          djdn(2) = (dnavg(2) * (dpot - deta) - navg) * detadn(2)
        end block
      else
        if (abs(deta) < 5) then
          call this%get_curr_sg(n, dpot, j, djdn, djddpot)
          if ((j < jmin) .or. (j > jmax)) then
            j = 0.5 * (jmin + jmax)
          end if
        else
          j = jmax
          if (deta < 0) j = jmin
        end if

        call newton_iteration(j, djdeta, djddpot)
        djdn = djdeta * detadn
      end if
    end if

  contains

    subroutine newton_iteration(j, djdeta, djddpot)
      real,    intent(inout) :: j
      real,    intent(out)   :: djdeta(2)
      real,    intent(out)   :: djddpot

      integer, parameter :: max_it = 1000

      integer :: it, it2
      logical :: status
      real    :: atol, rtol, dj, err, err0, f, dfdj, dfdeta(2), dfddpot, fmin, fmax, jmin0, jmax0, smin, smax, s

      ! clear
      djdeta  = 0
      djddpot = 0

      if (.not. small_eta) then
        ! init ode options
        call opt%init(1, atol = [1e-16], rtol = [1e-14], max_rejected = 200)

        ! lower bound
        jmin0 = - huge(1.0)
        it = 0
        do while (.true.)
          call residual(jmin, status, f = fmin)
          if (status) exit
          it = it + 1
          if (it > max_it) goto 200
          if (j - jmin < jmax - j) then
            jmin = jmin + max(abs(jmin) * 1e-15, 1e-100)
          else
            jmin = 0.5 * (jmin + jmax)
          end if
          jmin0 = jmin
        end do
        if (fmin == 0) then
          j = jmin
          goto 100
        end if

        ! upper bound
        jmax0 = huge(1.0)
        it = 0
        do while (.true.)
          call residual(jmax, status, f = fmax)
          if (status) exit
          it = it + 1
          if (it > max_it) goto 200
          if (j - jmin < jmax - j) then
            jmax = 0.5 * (jmin + jmax)
          else
            jmax  = jmax - max(abs(jmax) * 1e-15, 1e-100)
          end if
          jmax0 = jmax
        end do
        if (fmax == 0) then
          j = jmax
          goto 100
        end if

        ! widen bounds if necessary
        smin = sign(1.0, fmin)
        smax = sign(1.0, fmax)
        it2 = 0
        do while (smin == smax)
          it2 = it2 + 1
          if (it2 > max_it) goto 200
          if (abs(fmin) < abs(fmax)) then
            if (jmin == jmin0) then
              j = jmin
              goto 100
            end if

            j = jmin - max(abs(jmin) * 1e-12, 1e-100)
            it = 0
            do while (.true.)
              call residual(j, status, f = f)
              if (status) exit

              if (it > max_it) goto 200
              j     = j + max(abs(j) * 1e-15, 1e-100)
              jmin0 = j
            end do
            if (f == 0) goto 100

            jmax = jmin
            fmax = fmin
            jmin = j
            fmin = f
            smin = sign(1.0, fmin)
          else
            if (jmax == jmax0) then
              j = jmax
              goto 100
            end if

            j = jmax + max(abs(jmax) * 1e-10, 1e-100)
            it = 0
            do while (.true.)
              call residual(j, status, f = f)
              if (status) exit
              it = it + 1
              if (it > max_it) goto 200
              j     = j - max(abs(j) * 1e-15, 1e-100)
              jmax0 = j
            end do
            if (f == 0) goto 100

            jmin = jmax
            fmin = fmax
            jmax = j
            fmax = f
            smax = sign(1.0, fmax)
          end if
          j = 0.5 * (jmin + jmax)
        end do
      end if

      ! tolerances
      atol = max(2e-16 * abs(j), 1e-100)
      rtol = 1e-13
      err0 = huge(1.0)
      if (.not. small_eta) then
        err  = 0.5 * (jmax - jmin)
      else
        err = 0.5 * err0
      end if

      ! newton iteration
      it = 0
      do while ((err > atol) .and. (err > abs(j) * rtol))
        it = it + 1

        ! bisection
        if ((.not. small_eta) .and. ((j < jmin) .or. (j > jmax) .or. (err0 <= err))) j = 0.5 * (jmin + jmax)

        ! evaluate residual
        call residual(j, status, f = f, dfdj = dfdj)
        if (.not. status) call program_error("could not evaluate residual, even though jmin <= j <= jmax")
        if (f == 0) goto 100

        ! calculate newton update and new error
        dj   = f / dfdj
        err0 = err
        err  = abs(dj)

        ! update bounds
        if (.not. small_eta) then
          s = sign(1.0, f)
          if (s * smax > 0) then
            jmax = j
            fmax = f
          else
            jmin = j
            fmin = f
          end if
        end if

        ! update solution
        j = j - dj

        ! error if stuck
        if (it > max_it) goto 200

        ! exit if close to solution
        if (.not. small_eta) then
          if (((jmax - jmin) < 0.5 * abs(jmin + jmax) * rtol) .or. (0.5 * (jmax - jmin) < atol)) then
            j = 0.5 * (jmax + jmin)
            goto 100
          end if
        end if
      end do

      ! calculate derivatives and return
      100 call residual(j, status, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      if (.not. status) call program_error("could not evaluate residual at solution")
      djdeta  = - dfdeta / dfdj
      djddpot = - dfddpot / dfdj
      return

      ! error
      200 print "(A,ES25.16E3)", "eta1 = ", eta(1)
      print "(A,ES25.16E3)", "eta2 = ", eta(2)
      print "(A,ES25.16E3)", "dpot = ", dpot
      call program_error("stuck in iteration")
    end subroutine

    subroutine residual(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      if (small_eta) then
        call residual_small_eta(j, status, f = f, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      else
        call residual_shooting(j, status, f = f, dfdj = dfdj, dfdeta = dfdeta, dfddpot = dfddpot)
      end if
    end subroutine

    subroutine residual_small_eta(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      type(hp_real) :: harg, hB1, hB2, he(2), hf
      real          :: arg, B1, B2, dB1, dB2, e(2)

      status = ieee_is_finite(j)

      he%x = eta
      he%y = 0
      he = exp(he)
      e = hp_to_real(he)

      ! shifted bernoulli argument
      harg = dpot - TwoProduct(alpha, j)
      arg  = hp_to_real(harg)
      hB1  = ber(-harg)
      hB2  = ber( harg)
      B1   = hp_to_real(hB1)
      B2   = hp_to_real(hB2)
      dB1 = dberdx(-arg)
      dB2 = dberdx( arg)

      if (present(f)) then
        hf = j - hB1 * he(1) + hB2 * he(2)
        f  = hp_to_real(hf)
      end if
      if (present(dfdj)) then
        dfdj = 1 + alpha * (dB1 * e(1) + dB2 * e(2))
      end if
      if (present(dfdeta)) then
        dfdeta(1) = - B1 * e(1)
        dfdeta(2) =   B2 * e(2)
      end if
      if (present(dfddpot)) then
        dfddpot = - dB1 * e(1) + dB2 * e(2)
      end if
    end subroutine

    subroutine residual_shooting(j, status, f, dfdj, dfdeta, dfddpot)
      real,           intent(in)  :: j
      logical,        intent(out) :: status
      real, optional, intent(out) :: f
      real, optional, intent(out) :: dfdj
      real, optional, intent(out) :: dfdeta(2)
      real, optional, intent(out) :: dfddpot

      real :: e, e2, dedeta, de2de, dedj, de2dj, deddpot, de2ddpot, xknee

      status = .true.

      if (dir > 0) then
        call solve_ode(0.0, 1.0, eta(1), j, status, e, dedeta, deddpot, dedj)
        if (present(f)) then
          f = e - eta(2)
        end if
        if (present(dfdj)) then
          dfdj = dedj
        end if
        if (present(dfdeta)) then
          dfdeta(1) = dedeta
          dfdeta(2) = - 1.0
        end if
        if (present(dfddpot)) then
          dfddpot = deddpot
        end if
      elseif (dir < 0) then
        call solve_ode(1.0, 0.0, eta(2), j, status, e, dedeta, deddpot, dedj)

        if (present(f)) then
          f = e - eta(1)
        end if
        if (present(dfdj)) then
          dfdj = dedj
        end if
        if (present(dfdeta)) then
          dfdeta(1) = - 1.0
          dfdeta(2) = dedeta
        end if
        if (present(dfddpot)) then
          dfddpot = deddpot
        end if
      else
        if (abs(dpot) < 1e-3) then
          xknee = 0.5
        else
          xknee = (eta(2) - eta(1)) / dpot
        end if
        call solve_ode(0.0, xknee, eta(1), j, status, e, dedeta, deddpot, dedj)
        if (.not. status) return
        call solve_ode(xknee, 1.0, e, j, status, e2, de2de, de2ddpot, de2dj)
        if (.not. status) return

        if (present(f)) then
          f = e2 - eta(2)
        end if
        if (present(dfdj)) then
          dfdj = de2dj + de2de * dedj
        end if
        if (present(dfdeta)) then
          dfdeta(1) = de2de * dedeta
          dfdeta(2) = - 1.0
        end if
        if (present(dfddpot)) then
          dfddpot = de2ddpot + de2de * deddpot
        end if
      end if
    end subroutine

    recursive subroutine solve_ode(x0, x1, eta0, j, status, eta1, deta1deta0, deta1ddpot, deta1dj)
      use math_m, only: ber

      real,    intent(in)  :: x0
      real,    intent(in)  :: x1
      real,    intent(in)  :: eta0
      real,    intent(in)  :: j
      logical, intent(out) :: status
      real,    intent(out) :: eta1
      real,    intent(out) :: deta1deta0
      real,    intent(out) :: deta1ddpot
      real,    intent(out) :: deta1dj

      real             :: deta1dx(1)
      type(dual_3)     :: dl_B, dl_dpot, dl_eta0, dl_eta1, dl_j, dl_t, dl_t1, dl_t2, dl_xF
      type(ode_result) :: result

      status = .true.

      if (eta0 < etaF) then
        ! use dual numbers to avoid most of the manual derivative calculations
        call dl_eta0%init(eta0, i = 1)
        call dl_dpot%init(dpot, i = 2)
        call dl_j%init(      j, i = 3)
        dl_t = dl_dpot - alpha * dl_j

        ! find x where eta crosses etaF
        dl_t1 = dl_t * exp(   etaF) - dl_j
        dl_t2 = dl_t * exp(dl_eta0) - dl_j
        if (((dl_t1%x > 0) .and. (dl_t2%x > 0)) .or. ((dl_t1%x < 0) .and. (dl_t2%x < 0))) then
          if (abs(dl_t%x) < 1e-6) then
            dl_xF = x0 +      (exp(    dl_eta0) - exp(    etaF)) /      dl_j     &
              &   + dl_t    * (exp(2 * dl_eta0) - exp(2 * etaF)) / (2 * dl_j**2) &
              &   + dl_t**2 * (exp(3 * dl_eta0) - exp(3 * etaF)) / (3 * dl_j**3)
          else
            dl_xF = x0 + log(dl_t1 / dl_t2) / dl_t
          end if

          if (((x1 > x0) .and. (dl_xF%x >= x0) .and. (dl_xF%x < x1)) .or. &
            & ((x1 < x0) .and. (dl_xF%x <= x0) .and. (dl_xF%x > x1))) then
            ! go from (xF,etaF) to (x1,eta1) using ode solver
            call solve_ode(dl_xF%x, x1, etaF, j, status, eta1, deta1deta0, deta1ddpot, deta1dj)
            if (.not. status) return
            call detadx(x1, [eta1], [dpot, j], status, f = deta1dx)
            if (.not. status) return
            deta1deta0 =            - deta1dx(1) * dl_xF%dx(1)
            deta1ddpot = deta1ddpot - deta1dx(1) * dl_xF%dx(2)
            deta1dj    = deta1dj    - deta1dx(1) * dl_xF%dx(3)
          else
            ! stay below -16
            dl_B%x  = ber(   dl_t%x * (x1 - x0))
            dl_B%dx = dberdx(dl_t%x * (x1 - x0)) * dl_t%dx
            dl_eta1 = dl_eta0 + log1p((dl_t - dl_j * exp(-dl_eta0)) * (x1 - x0) / dl_B)

            eta1       = dl_eta1%x
            deta1deta0 = dl_eta1%dx(1)
            deta1ddpot = dl_eta1%dx(2)
            deta1dj    = dl_eta1%dx(3)
          end if
        else
          ! stay constant (avoid numerical issues)
          eta1       = eta0
          deta1deta0 = 1
          deta1ddpot = 0
          deta1dj    = 0
        end if
      else
        ! go from (x0,eta0) to (x1,eta1) using ode solver
        call radau5(detadx, x0, x1, [x1], [eta0], [dpot, j], opt, status, result)
        if (.not. status) return
        eta1       = result%Usmp(    1,  1)
        deta1deta0 = result%dUsmpdU0(1,1,1)
        deta1ddpot = result%dUsmpdP( 1,1,1)
        deta1dj    = result%dUsmpdP( 1,2,1)
      end if
    end subroutine

    subroutine detadx(x, U, p, status, f, dfdU, dfdp)
      !! ode right-hand side
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (eta)
      real,           intent(in)  :: p(:)
        !! parameters (dpot, j)
      logical,        intent(out) :: status
        !! success/fail
      real, optional, intent(out) :: f(:)
        !! output deta/dx
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivative of f wrt eta
      real, optional, intent(out) :: dfdp(:,:)
        !! output derivative of f wrt P

      real :: F12, dF12

      m4_ignore(x)

      call fermi_dirac_integral_1h(U(1), F12, dF12)
      status = ieee_is_finite(F12)

      if (present(f)) then
        f(1) = p(1) - p(2) / F12
      end if
      if (present(dfdU)) then
        dfdU(1,1) = p(2) / F12**2 * dF12
      end if
      if (present(dfdp)) then
        dfdp(1,1) = 1.0
        dfdp(1,2) = - 1.0 / F12
      end if
    end subroutine

  end subroutine

end module
