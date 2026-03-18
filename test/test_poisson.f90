module test_poisson_m
  use charge_density_m, only: charge_density
  use device_params_m,  only: device_params
  use electric_field_m, only: calc_electric_field, calc_electric_field_component, electric_field, electric_field_component
  use esystem_m,        only: esystem
  use grid_m,           only: IDX_EDGE, IDX_VERTEX
  use grid_generator_m, only: DIR_NAME
  use input_m,          only: input_file, input_section
  use input_src_m,      only: const_src
  use math_m,           only: PI
  use normalization_m,  only: denorm, init_normconst, norm
  use poisson_m,        only: poisson
  use potential_m,      only: potential
  use steady_state_m,   only: steady_state
  use test_case_m,      only: test_case
  use voltage_m,        only: voltage

  implicit none

  public test_poisson

  ! temperature, used for normalization
  real, parameter :: T = 300.0

  ! output simulation errors in console
  logical, parameter :: OUT = .false.

contains

  subroutine test_poisson
    type(test_case) :: tc

    call init_normconst(T) ! may be more useful inside test
    call tc%init("poisson")

    call test_poisson_oneDim_cap(tc)
    call test_poisson_oneDim_cap_permShift(tc)
    call test_poisson_oneDim_rho(tc)
    call test_poisson_twoDim_cap(tc)
    call test_poisson_twoDim_rho(tc)
    call test_poisson_threeDim_cap(tc)

    call tc%finish()
  end subroutine

  subroutine test_poisson_oneDim_cap(tc)
    ! One-dimensional parallel plate capacitor
    ! Expected solution: linear potential dependant on x-coordinate and constant field
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(2)

    type(calc_electric_field_component) :: calc_ef_dir(1)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(1)
    type(electric_field)                :: ef

    integer           :: i, idx_bnd(2,1)
    real              :: atol, c(2), p(3,1), rtol, v(2)
    real, allocatable :: sol(:), tmp(:)

    if(OUT) print "(A)", "test poisson: 1D capacitor"

    call file%init("cap1D.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    do i = 1, 2
      call volt(i)%init("V_"//par%contacts(i)%name)
    end do
    call ef%init(par)

    ! init electric field
    call ef_dir(1)%init(par, 1)
    call ef%init(par)
    call calc_ef_dir(1)%init(par, ef_dir(1), pot)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    call sys%add_equation(calc_ef_dir(1))
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    do i = 1, 2
      call sys%provide(volt(i), input = .true.)
    end do
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! set input
    tmp = 0.0
    call rho%set(tmp)
    v(1) = norm(1.0, "V")
    v(2) = norm(2.0, "V")
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    call calc_ef_dir(1)%eval()
    call calc_ef%eval()

    ! test potential
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    c(1) = (v(1) - v(2)) / (p(1,1) - p(2,1))
    c(2) = v(1) - c(1) * p(1,1)
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
      sol(i) = p(3,1)
    end do
    sol = c(1) * sol + c(2)
    atol = norm(1e-14, "V")
    rtol = 1e-12

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 1D capacitor - potential")

    ! test directional electric field
    sol = -c(1)
    atol = norm(1e-6, "V/m")
    rtol = 1e-10

    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. relative error: ", max(maxval(abs((sol - ef_dir(1)%get()) / sol), mask = abs(sol - ef_dir(1)%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. absolute error: ", denorm(maxval(abs((sol - ef_dir(1)%get()))), "V/m")

    call tc%assert_eq(sol, ef_dir(1)%get(), rtol, atol, "poisson: 1D capacitor - electric field in "//DIR_NAME(1)//" direction")

    ! test absolute electric field
    sol = abs(c(1))
    atol = norm(1e-6, "V/m")
    rtol = 1e-10

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, atol, "poisson: 1D capacitor - electric field")

  end subroutine

  subroutine test_poisson_oneDim_cap_permShift(tc)
    ! One-dimensional parallel plate capacitor with 2 permittivity regions
    ! Expected potential solution: linear potential for each region of constant permittivity, dependant on x-coordinate
    ! Expected field solution: constant for each region and a middle point with an average of both, due to permittivity being defined on cells
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(2)

    type(calc_electric_field_component) :: calc_ef_dir(1)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(1)
    type(electric_field)                :: ef

    integer           :: i, idx1(1), idx2(1), idx_bnd(2,1)
    logical           :: status(2)
    real              :: atol, c(3), eps(2), p(4,1), rtol, v(2)
    real, allocatable :: sol(:), tmp(:)

    if(OUT) print "(A)", "test poisson: 1D capacitor with a shift in permittivity"

    call file%init("cap1DpermShift.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    do i = 1, 2
      call volt(i)%init("V_"//par%contacts(i)%name)
    end do

    ! init electric field
    call ef_dir(1)%init(par, 1)
    call calc_ef_dir(1)%init(par, ef_dir(1), pot)
    call ef%init(par)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    call sys%add_equation(calc_ef_dir(1))
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    do i = 1, 2
      call sys%provide(volt(i), input = .true.)
    end do
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! set input
    tmp = 0.0
    call rho%set(tmp)
    v(1) = norm(1.0, "V")
    v(2) = norm(2.0, "V")
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    call calc_ef_dir(1)%eval()
    call calc_ef%eval()

    ! get point of permittivity change
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    do i = 2, par%poisson(IDX_VERTEX, 0)%n-1
      call par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, par%poisson(IDX_VERTEX, 0)%get_idx(i), 1, idx1, status(1))
      call par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, 1, par%poisson(IDX_VERTEX, 0)%get_idx(i), 2, idx2, status(2))
      if (all(status)) then
        if(par%eps(IDX_EDGE, 1)%get(idx1) /= par%eps(IDX_EDGE, 1)%get(idx2)) then
          call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
          exit
        end if
      end if
    end do

    ! test potential
    call par%g%get_idx_bnd(IDX_EDGE, 1, idx_bnd)
    eps(1) = par%eps(IDX_EDGE, 1)%get(idx_bnd(1,:))
    eps(2) = par%eps(IDX_EDGE, 1)%get(idx_bnd(2,:))
    c(3) = (v(2) - v(1)) / ((p(3,1) - p(1,1)) / (eps(1)) - (p(3,1) - p(2,1)) / (eps(2)))
    c(1) = v(1) - c(3) * p(1,1) / eps(1)
    c(2) = v(2) - c(3) * p(2,1) / eps(2)
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(4,:))
      if (p(4,1) < p(3,1)) then
        sol(i) = p(4,1) * c(3) / eps(1) + c(1)
      else
        sol(i) = p(4,1) * c(3) / eps(2) + c(2)
      end if
    end do
    atol = norm(1e-14, "V")
    rtol = 1e-12

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 1D capacitor with a shift in permittivity - potential")

    ! test directional electric field
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(4,:))
      if (p(4,1) < p(3,1)) then
        sol(i) = - c(3) / eps(1)
      else if(p(4,1) == p(3,1)) then
        sol(i) = - ((c(3) / eps(1))+(c(3) / eps(2)))/2
      else
        sol(i) = - c(3) / eps(2)
      end if
    end do
    atol = norm(1e-6, "V/m")
    rtol = 1e-11

    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. relative error: ", max(maxval(abs((sol - ef_dir(1)%get()) / sol), mask = abs(sol - ef_dir(1)%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. absolute error: ", denorm(maxval(abs((sol - ef_dir(1)%get()))), "V/m")

    call tc%assert_eq(sol, ef_dir(1)%get(), rtol, atol, "poisson: 1D capacitor with a shift in permittivity - electric field in "//DIR_NAME(1)//" direction")

    ! test absolute electric field
    sol = abs(sol)
    atol = norm(1e-6, "V/m")
    rtol = 1e-11

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, atol, "poisson: 1D capacitor with a shift in permittivity - electric field")

  end subroutine

  subroutine test_poisson_oneDim_rho(tc)
    ! One-dimensional parallel plate capacitor with sinusoidal charge, assume device to be symmetric around x=0
    ! Expected potential solution: linear and cosine dependance on x
    ! Expected field solution: linear and sine dependence on x
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(2)

    type(calc_electric_field_component) :: calc_ef_dir(1)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(1)
    type(electric_field)                :: ef

    integer           :: i, idx(1), idx_bnd(2,1)
    real              :: atol, eps, len, p(3,1), rho0, rtol, v(2), vol, x(2)
    real, allocatable :: sol(:), tmp(:)

    if(OUT) print "(A)", "test poisson: 1D capacitor with sinusoidal charge"

    call file%init("1Dtransport.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    do i = 1, 2
      call volt(i)%init("V_"//par%contacts(i)%name)
    end do

    ! init electric field
    call ef_dir(1)%init(par, 1)
    call calc_ef_dir(1)%init(par, ef_dir(1), pot)
    call ef%init(par)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    call sys%add_equation(calc_ef_dir(1))
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    do i = 1, 2
      call sys%provide(volt(i), input = .true.)
    end do
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! calculate and set charge density
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    rho0 = norm(1e-2, "C/cm^3")
    len  = abs(p(2,1) - p(1,1))
    vol  = len / (idx_bnd(2,1) - 1)
    eps  = par%eps(IDX_VERTEX, 0)%get(idx_bnd(1,:))
    v(1) = norm(1.0, "V")
    v(2) = norm(2.0, "V")
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      idx = par%poisson(IDX_VERTEX, 0)%get_idx(i)
      call par%g%get_vertex(idx, p(3,:))
      x = p(3,1)
      if(idx(1) /= idx_bnd(1,1)) x(1) = x(1) - vol / 2
      if(idx(1) /= idx_bnd(2,1)) x(2) = x(2) + vol / 2
      tmp(i) = rho0 / 2 * (1 + cos(2 * PI * p(3,1) / len))
    end do

    call rho%set(tmp)
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    call calc_ef_dir(1)%eval()
    call calc_ef%eval()

    ! test potential
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
      tmp(i) = p(3,1)
    end do
    sol =  rho0 / (4 * eps) * (len**2 / 4 + len**2 / (2 * PI**2) - tmp**2 + len**2 / (2 * PI**2) * cos(2 * PI * tmp / len)) + tmp * (v(2) - v(1)) / len + (v(1) + v(2)) / 2
    atol = norm(1e-7, "V")
    rtol = 1e-5

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 1D capacitor with sinusoidal charge - potential")

    ! test directional electric field
    sol = rho0 / eps * (tmp / 2 + len / (4 * PI) * sin(2 * PI * tmp / len)) + (v(1) - v(2)) / len
    atol = norm(1e0, "V/m")
    rtol = 1e-2

    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. relative error: ", max(maxval(abs((sol - ef_dir(1)%get()) / sol), mask = abs(sol - ef_dir(1)%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(1)//" max. absolute error: ", denorm(maxval(abs((sol - ef_dir(1)%get()))), "V/m")

    call tc%assert_eq(sol, ef_dir(1)%get(), rtol, atol, "poisson: 1D capacitor with sinusoidal charge - electric field in "//DIR_NAME(1)//" direction")

    ! test absolute electric field
    sol = abs(sol)
    atol = norm(1e0, "V/m")
    rtol = 1e-2

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, norm(2e3, "V"), "poisson: 1D capacitor with sinusoidal charge - electric field")

  end subroutine

  subroutine test_poisson_twoDim_cap(tc)
    ! Two-dimensional parallel plate capacitor
    ! Expected solution: linear potential dependant on x-coordinate and constant field
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(2)

    type(calc_electric_field_component) :: calc_ef_dir(2)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(2)
    type(electric_field)                :: ef

    integer           :: i, idx_bnd(2,2)
    real              :: atol, c(2), p(3,2), rtol, v(2)
    real, allocatable :: sol(:), sol2(:,:), tmp(:)

    if(OUT) print "(A)", "test poisson: 2D capacitor"

    call file%init("cap2D.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    do i = 1, 2
      call volt(i)%init("V_"//par%contacts(i)%name)
    end do

    ! init electric field
    do i = 1, 2
      call ef_dir(i)%init(par, i)
      call calc_ef_dir(i)%init(par, ef_dir(i), pot)
    end do
    call ef%init(par)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    do i = 1, 2
      call sys%add_equation(calc_ef_dir(i))
    end do
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    do i = 1, 2
      call sys%provide(volt(i), input = .true.)
    end do
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), sol2(2, par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! set input
    tmp = 0.0
    call rho%set(tmp)
    v(1) = norm(1.0, "V")
    v(2) = norm(2.0, "V")
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    do i = 1, 2
      call calc_ef_dir(i)%eval()
    end do
    call calc_ef%eval()

    ! test potential
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    c(1) = (v(1) - v(2)) / (p(1,1) - p(2,1))
    c(2) = v(1) - c(1) * p(1,1)
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
      tmp(i) = p(3,1)
    end do
    sol = c(1) * tmp + c(2)
    atol = norm(1e-14, "V")
    rtol = 1e-12

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 2D capacitor - potential")

    ! test directional electric field
    sol2(1,:) = - c(1)
    sol2(2,:) = 0
    atol = norm(1e-4, "V/m")
    rtol = 1e-11

    do i = 1, 2
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. relative error: ", max(maxval(abs((sol2(i,:) - ef_dir(i)%get()) / sol2(i,:)), mask = abs(sol2(i,:) - ef_dir(i)%get()) > atol), 0.0)
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. absolute error: ", denorm(maxval(abs((sol2(i,:) - ef_dir(i)%get()))), "V/m")

      call tc%assert_eq(sol2(i,:), ef_dir(i)%get(), rtol, atol, "poisson: 2D capacitor - electric field in "//DIR_NAME(i)//" direction")
    end do

    ! test absolute electric field
    sol = abs(sol2(1,:))
    atol = norm(1e-4, "V/m")
    rtol = 1e-11

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, atol, "poisson: 2D capacitor - electric field")

  end subroutine

  subroutine test_poisson_twoDim_rho(tc)
    ! Two-dimensional square transport region with charge density that satisfies Neumann BC; Contact at corner for potential reference
    ! Expected potential solution: phi = C * (x^2 * (L-x)^2 + y^2 * (L-y)^2)) + V_c with C = rho_0 / (L^2 * eps)
    ! Expected field solution:     E_x = C * (6 * L * x^2 - 4 * x^3 - 2 * L^2 * x) similar for y
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(1)

    type(calc_electric_field_component) :: calc_ef_dir(2)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(2)
    type(electric_field)                :: ef

    integer           :: i, idx(2), idx_bnd(2,2)
    real              :: atol, eps, len, p(3,2), rho0, rtol, v(1), x(2), y(2)
    real, allocatable :: sol(:), sol2(:,:), tmp(:)

    if(OUT) print "(A)", "test poisson: 2D device with polynomial charge"

    call file%init("2Dtransport.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    call volt(1)%init("V_"//par%contacts(1)%name)

    ! init electric field
    do i = 1, 2
      call ef_dir(i)%init(par, i)
      call calc_ef_dir(i)%init(par, ef_dir(i), pot)
    end do
    call ef%init(par)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    do i = 1, 2
      call sys%add_equation(calc_ef_dir(i))
    end do
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    call sys%provide(volt(1), input = .true.)
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), sol2(2, par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! calculate and set charge density
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    len  = abs(p(2,1) - p(1,1))
    rho0 = norm(1e-1, "C/cm^3")
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      idx = par%poisson(IDX_VERTEX, 0)%get_idx(i)
      call par%g%get_vertex(idx, p(3,:))
      x = p(3,1)
      if(idx(1) /= idx_bnd(1,1)) x(1) = x(1) - len / (idx_bnd(2,1) - 1) / 2
      if(idx(1) /= idx_bnd(2,1)) x(2) = x(2) + len / (idx_bnd(2,1) - 1) / 2
      y = p(3,2)
      if(idx(2) /= idx_bnd(1,2)) y(1) = y(1) - len / (idx_bnd(2,2) - 1) / 2
      if(idx(2) /= idx_bnd(2,2)) y(2) = y(2) + len / (idx_bnd(2,2) - 1) / 2
      tmp(i) = - rho0 / len**2 * (4 * (x(1)**2 + x(1) * x(2) + x(2) **2) - 6 * len * (x(1) + x(2)) &
                                + 4 * (y(1)**2 + y(1) * y(2) + y(2) **2) - 6 * len * (y(1) + y(2)) + 4 * len**2)
    end do
    call rho%set(tmp)

    v(1) = norm(1.0, "V")
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    do i = 1, 2
      call calc_ef_dir(i)%eval()
    end do
    call calc_ef%eval()

    ! test potential
    eps  = par%eps(IDX_EDGE, 1)%get(idx_bnd(1,:))
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
      x(1) = p(3,1)
      y(1) = p(3,2)
      sol(i) = rho0 / (eps * len**2) * (y(1)**2 * (len - y(1))**2 + x(1)**2 * (len - x(1))**2) + v(1)
    end do
    atol = norm(1e-6, "V")
    rtol = 1e1

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 2D device with polynomial charge - potential")

    ! test directional electric field
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      idx = par%poisson(IDX_VERTEX, 0)%get_idx(i)
      call par%g%get_vertex(idx, p(3,:))
      x = p(3,1)
      if(idx(1) /= idx_bnd(1,1)) x(1) = x(1) - len / (idx_bnd(2,1) - 1) / 2
      if(idx(1) /= idx_bnd(2,1)) x(2) = x(2) + len / (idx_bnd(2,1) - 1) / 2
      y = p(3,2)
      if(idx(2) /= idx_bnd(1,2)) y(1) = y(1) - len / (idx_bnd(2,2) - 1) / 2
      if(idx(2) /= idx_bnd(2,2)) y(2) = y(2) + len / (idx_bnd(2,2) - 1) / 2
      sol2(1,i) = rho0 / (len**2 * eps * abs(x(2) - x(1))) * (len**2 * x(1)**2 - 2 * len * x(1)**3 + x(1)**4 &
                                                             -len**2 * x(2)**2 + 2 * len * x(2)**3 - x(2)**4)
      sol2(2,i) = rho0 / (len**2 * eps * abs(y(2) - y(1))) * (len**2 * y(1)**2 - 2 * len * y(1)**3 + y(1)**4 &
                                                             -len**2 * y(2)**2 + 2 * len * y(2)**3 - y(2)**4)
    end do
    atol = norm(1e2, "V/m")
    rtol = 1e-2

    ! current implementation of electric field yields constant relative errors on boundaries, set solution on boundaries to simulated value
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      idx = par%poisson(IDX_VERTEX, 0)%get_idx(i)
      if((idx(1) == idx_bnd(1,1)) .or. (idx(1) == idx_bnd(2,1))) sol2(1,i) = ef_dir(1)%get(idx)
      if((idx(2) == idx_bnd(1,2)) .or. (idx(2) == idx_bnd(2,2))) sol2(2,i) = ef_dir(2)%get(idx)
    end do

    do i = 1, 2
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. relative error: ", max(maxval(abs((sol2(i,:) - ef_dir(i)%get()) / sol2(i,:)), mask = abs(sol2(i,:) - ef_dir(i)%get()) > atol), 0.0)
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. absolute error: ", denorm(maxval(abs((sol2(i,:) - ef_dir(i)%get()))), "V/m")

      call tc%assert_eq(sol2(i,:), ef_dir(i)%get(), rtol, atol, "poisson: 2D device with polynomial charge - electric field in "//DIR_NAME(i)//" direction")
    end do

    ! test absolute electric field
    sol = sqrt(sol2(1,:)**2 + sol2(2,:)**2)
    atol = norm(1e2, "V/m")
    rtol = 1e-2

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, atol, "poisson: 2D device with polynomial charge - electric field")

  end subroutine

  subroutine test_poisson_threeDim_cap(tc)
    ! Three-dimensional parallel plate capacitor
    ! Expected solution: linear potential dependant on x-coordinate and constant field
    type(test_case), intent(inout) :: tc

    type(const_src)     :: ss_input
    type(device_params) :: par
    type(esystem)       :: sys
    type(input_file)    :: file
    type(input_section) :: ss_params
    type(poisson)       :: poiss
    type(steady_state)  :: ss

    type(charge_density) :: rho
    type(potential)      :: pot
    type(voltage)        :: volt(2)

    type(calc_electric_field_component) :: calc_ef_dir(3)
    type(calc_electric_field)           :: calc_ef
    type(electric_field_component)      :: ef_dir(3)
    type(electric_field)                :: ef

    integer           :: i, idx_bnd(2,3)
    real              :: atol, c(2), p(3,3), rtol, v(2)
    real, allocatable :: sol(:), sol2(:,:), tmp(:)

    if(OUT) print "(A)", "test poisson: 3D capacitor"

    call file%init("cap3D.ini")
    call par%init(file, T)

    ! init variables
    call pot%init(par)
    call rho%init(par)
    do i = 1, 2
      call volt(i)%init("V_"//par%contacts(i)%name)
    end do

    ! init electric field
    do i = 1, 3
      call ef_dir(i)%init(par, i)
      call calc_ef_dir(i)%init(par, ef_dir(i), pot)
    end do
    call ef%init(par)
    call calc_ef%init(par, ef_dir, ef)

    ! init system
    call poiss%init(par, pot, rho, volt)
    call sys%init("poisson")
    call sys%add_equation(poiss)
    do i = 1, 3
      call sys%add_equation(calc_ef_dir(i))
    end do
    call sys%add_equation(calc_ef)
    call sys%provide(rho)
    do i = 1, 2
      call sys%provide(volt(i), input = .true.)
    end do
    call sys%init_final()

    allocate(sol(par%poisson(IDX_VERTEX, 0)%n), sol2(3, par%poisson(IDX_VERTEX, 0)%n), tmp(par%poisson(IDX_VERTEX, 0)%n))

    ! set input
    tmp = 0.0
    call rho%set(tmp)
    v(1) = norm(1.0, "V")
    v(2) = norm(2.0, "V")
    call ss_input%init(v)

    ! run simulation
    call ss_params%init("")
    call ss_params%set("log", .true.)
    call ss%init(sys)
    call ss%set_params(ss_params)
    call ss%run(input = ss_input)

    ! evaluate electric field
    do i = 1, 3
      call calc_ef_dir(i)%eval()
    end do
    call calc_ef%eval()

    ! test potential
    call par%g%get_idx_bnd(IDX_VERTEX, 0, idx_bnd)
    call par%g%get_vertex(idx_bnd(1,:), p(1,:))
    call par%g%get_vertex(idx_bnd(2,:), p(2,:))
    c(1) = (v(1) - v(2)) / (p(1,1) - p(2,1))
    c(2) = v(1) - c(1) * p(1,1)
    do i = 1, par%poisson(IDX_VERTEX, 0)%n
      call par%g%get_vertex(par%poisson(IDX_VERTEX, 0)%get_idx(i), p(3,:))
      tmp(i) = p(3,1)
    end do
    sol = c(1) * tmp + c(2)
    atol = norm(1e-14, "V")
    rtol = 1e-12

    if(OUT) print "(A,ES25.16E3)", "potential max. relative error: ", max(maxval(abs((sol - pot%get()) / sol), mask = abs(sol - pot%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "potential max. absolute error: ", denorm(maxval(abs((sol - pot%get()))), "V")

    call tc%assert_eq(sol, pot%get(), rtol, atol, "poisson: 3D capacitor - potential")

    ! test directional electric field
    sol2(1,:) = - c(1)
    sol2(2,:) = 0
    sol2(3,:) = 0
    atol = norm(1e-3, "V/m")
    rtol = 1e-10

    do i = 1, 3
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. relative error: ", max(maxval(abs((sol2(i,:) - ef_dir(i)%get()) / sol2(i,:)), mask = abs(sol - sol2(i,:)) > atol), 0.0)
      if(OUT) print "(A,ES25.16E3)", "electric field "//DIR_NAME(i)//" max. absolute error: ", denorm(maxval(abs((sol2(i,:) - ef_dir(i)%get()))), "V/m")
      call tc%assert_eq(sol2(i,:), ef_dir(i)%get(), rtol, atol, "poisson: 3D capacitor - electric field in "//DIR_NAME(i)//" direction")
    end do

    ! test absolute electric field
    sol = abs(sol2(1,:))
    atol = norm(1e-3, "V/m")
    rtol = 1e-10

    if(OUT) print "(A,ES25.16E3)", "total electric field max. relative error: ", max(maxval(abs((sol - ef%get()) / sol), mask = abs(sol - ef%get()) > atol), 0.0)
    if(OUT) print "(A,ES25.16E3)", "total electric field max. absolute error: ", denorm(maxval(abs((sol - ef%get()))), "V/m")

    call tc%assert_eq(sol, ef%get(), rtol, atol, "poisson: 3D capacitor - electric field")

  end subroutine

end module
