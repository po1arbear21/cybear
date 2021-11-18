module test_circuit_m

  use test_case_m,    only: test_case
  use capacitor_m,    only: capacitor, new_capacitor
  use circuit_m,      only: circuit, component
  use csource_m,      only: csource, new_csource
  use ground_m,       only: ground, new_ground
  use inductor_m,     only: inductor, new_inductor
  use input_src_m,    only: const_src
  use newton_m,       only: newton_opt
  use resistor_m,     only: resistor, new_resistor
  use steady_state_m, only: steady_state
  use vsource_m,      only: vsource, new_vsource

  implicit none

  private
  public test_circuit

contains

  subroutine test_circuit()
    type(test_case) :: tc

    type(circuit)             :: crt
    type(ground),     pointer :: GND
    type(vsource),    pointer :: V1
    type(resistor),   pointer :: R1, R2, R3
    type(capacitor),  pointer :: C1, C2
    type(inductor),   pointer :: L1
    class(component), pointer :: comp
    real                      :: I2

    type(const_src)    :: input
    type(steady_state) :: ss
    type(newton_opt)   :: nopt

    call tc%init("circuit")

    ! create circuit components
    call crt%init("circuit")
    GND => new_ground(   crt, "GND")
    V1  => new_vsource(  crt, "V1")
    R1  => new_resistor( crt, "R1", 0.2)
    R2  => new_resistor( crt, "R2", 0.5)
    R3  => new_resistor( crt, "R3", 1.0)
    C1  => new_capacitor(crt, "C1", 0.1)
    C2  => new_capacitor(crt, "C2", 0.3)
    L1  => new_inductor( crt, "L1", 2.5)
    call tc%assert_eq(8, crt%comps%n, "circuit add components")

    ! connect
    call crt%connect(GND%terms(1), V1%terms(1))
    call crt%connect(V1%terms(2), R1%terms(1))
    call crt%connect(C1%terms(1), R1%terms(2))
    call crt%connect(C2%terms(1), L1%terms(1))
    call crt%connect(C1%terms(2), GND%terms(1))
    call crt%connect(C2%terms(2), GND%terms(1))
    call crt%connect(L1%terms(1), R1%terms(2))
    call crt%connect(R2%terms(1), R3%terms(1))
    call crt%connect(R2%terms(2), GND%terms(1))
    call crt%connect(R3%terms(2), GND%terms(1))
    call crt%connect(L1%terms(2), R2%terms(1))
    call tc%assert_eq(4, crt%nodes%n-crt%free_nodes%n, "circuit connect components")

    ! get
    comp => crt%get("R2")
    select type (comp)
      type is (resistor)
        call tc%assert_eq(0.5, comp%R, 0.0, "circuit get component")
      class default
        call tc%assert(.false., "circuit get component")
    end select

    call crt%init_final()

    ! solve steady state; linear => maximal 2 iterations
    call input%init([2.0]) ! V1
    call nopt%init(crt%sys%n, atol = 1e-15, rtol = 1e-12, log = .false., max_it = 2)
    call ss%run(crt%sys, nopt = nopt, input = input)

    I2 = V1%V%x / (R1%R + R2%R + R1%R * R2%R / R3%R)
    call tc%assert_eq( I2, R2%terms(1)%curr%x, 1e-14, "circuit DC R2.1.I")
    call tc%assert_eq(-I2, R2%terms(2)%curr%x, 1e-14, "circuit DC R2.2.I")

    call tc%finish()
  end subroutine

end module
