module test_circuit_m

  use test_case_m,    only: test_case
  use capacitor_m,    only: capacitor
  use circuit_m,      only: circuit
  use csource_m,      only: csource
  use ground_m,       only: ground
  use inductor_m,     only: inductor
  use input_src_m,    only: const_src
  use newton_m,       only: newton_opt
  use resistor_m,     only: resistor
  use steady_state_m, only: steady_state
  use vsource_m,      only: vsource

  implicit none

  private
  public test_circuit

contains

  subroutine test_circuit()
    type(test_case) :: tc

    type(circuit)   :: crt
    type(ground)    :: GND
    type(vsource)   :: V1
    type(resistor)  :: R1, R2, R3
    type(capacitor) :: C1, C2
    type(inductor)  :: L1
    real            :: I2

    type(const_src)    :: input
    type(steady_state) :: ss
    type(newton_opt)   :: nopt

    call tc%init("circuit")

    ! create circuit
    call crt%init("circuit")
    call GND%init(crt, "GND", "N0")
    call V1%init(crt, "V1", "N0", "N1")
    call R1%init(crt, "R1", "N1", "N2", 0.2)
    call R2%init(crt, "R2", "N3", "N0", 0.5)
    call R3%init(crt, "R3", "N3", "N0", 1.2)
    call C1%init(crt, "C1", "N2", "N0", 0.1)
    call C2%init(crt, "C2", "N2", "N0", 0.3)
    call L1%init(crt, "L1", "N2", "N3", 2.5)
    call tc%assert_eq(8, crt%comps%n, "circuit add components")
    call crt%init_final()

    ! solve steady state; linear => maximal 2 iterations
    call input%init([2.0]) ! V1
    call nopt%init(crt%sys%n, atol = 1e-15, rtol = 1e-12, log = .false., max_it = 2)
    call ss%run(crt%sys, nopt = nopt, input = input)

    ! I2 = V1 / (R1 + R2 + R1 * R2 / R3)
    I2 = V1%V%x / (0.2 + 0.5 + 0.2 * 0.5 / 1.2)
    call tc%assert_eq( I2, R2%comp%terms(1)%curr%x, 1e-14, "circuit DC R2.I1")
    call tc%assert_eq(-I2, R2%comp%terms(2)%curr%x, 1e-14, "circuit DC R2.I2")

    call tc%finish()
  end subroutine

end module
