m4_include(util/macro.f90.inc)

module electric_field_m

  use device_params_m, only: device_params
  use equation_m,      only: equation
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use grid_generator_m, only: DIR_NAME
  use jacobian_m,      only: jacobian
  use potential_m,     only: potential
  use stencil_m,       only: near_neighb_stencil
  use variable_m,      only: variable_real
  use vselector_m,     only: vselector
  use normalization_m, only: denorm

  implicit none

  private
  public electric_field, calc_efield

  type, extends(variable_real) :: electric_field
    !! electric field component at vertices

    integer :: dir
      !! field component direction (1=x, 2=y, 3=z)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => electric_field_init
  end type

  type, extends(equation) :: calc_efield
    !! calculate electric field component from potential gradient using finite volume method

    type(device_params),  pointer :: par    => null()
      !! pointer to device parameters

    type(vselector) :: pot
      !! variable selector for potential in equation order
    type(vselector) :: efield
      !! variable selector for electric field in equation order

    type(near_neighb_stencil) :: st_nn
      !! near neighbor stencil for vertex-to-vertex connectivity

    type(jacobian), pointer :: jaco_pot => null()
      !! jacobian for potential dependency
  contains
    procedure :: init => calc_efield_init
    procedure :: eval => calc_efield_eval
  end type

contains

  subroutine electric_field_init(this, par, dir)
    !! initialize electric field component
    class(electric_field), intent(out) :: this
    type(device_params),   intent(in)  :: par
      !! device parameters
    integer,               intent(in)  :: dir
      !! field component direction (1=x, 2=y, 3=z)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("E"//DIR_NAME(dir), "V/cm", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%dir = dir

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

  subroutine calc_efield_init(this, par, pot, efield)
    !! initialize electric field calculation equation
    class(calc_efield),    intent(out) :: this
    type(device_params),   target, intent(in)  :: par
      !! device parameters
    type(potential),       target, intent(in)  :: pot
      !! potential variable
    type(electric_field),  target, intent(in)  :: efield
      !! electric field component variable

    integer :: iprov, ict, i, j, idx_dir, dir, max_neighb
    real    :: geometric_factor, A_kl, L_kl, Omega_k
    real    :: p_k(par%g%dim), p_l(par%g%dim)
    integer, allocatable :: idx_k(:), idx_l(:), edge_idx(:)
    logical :: status

    print "(A)", "calc_efield_init for component "//DIR_NAME(efield%dir)

    ! init base
    call this%equation_init("calc_E"//DIR_NAME(efield%dir))
    this%par => par

    ! init stencil
    call this%st_nn%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)

    ! create variable selectors for equation order
    call this%pot%init(pot, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    call this%efield%init(efield, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])

    ! provide electric field component at all vertices (interior + contacts)
    iprov = this%provide(efield, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])

    ! depend on potential at all vertices
    this%jaco_pot => this%init_jaco(iprov, this%depend(pot, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)]), &
      & st = [(this%st_nn%get_ptr(), ict = 0, par%nct)], const = .true.)

    ! Set up geometric factors in Jacobian
    dir = efield%dir
    allocate(idx_k(par%g%idx_dim), idx_l(par%g%idx_dim), edge_idx(par%g%idx_dim))

    ! loop over all transport vertices
    do ict = 0, par%nct
      do i = 1, par%transport_vct(ict)%n
        idx_k = par%transport_vct(ict)%get_idx(i)

        ! get control volume at vertex k (same for all vertices)
        Omega_k = par%tr_vol%get(idx_k)
        if (Omega_k <= 0.0) cycle  ! skip if zero volume

        ! get coordinates of vertex k
        call par%g%get_vertex(idx_k, p_k)

        ! loop over all edge directions
        do idx_dir = 1, par%g%idx_dim
          max_neighb = par%g%get_max_neighb(IDX_VERTEX, 0, IDX_EDGE, idx_dir)

          ! loop over neighboring edges in this direction
          do j = 1, max_neighb
            ! get edge index
            call par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, idx_dir, idx_k, j, edge_idx, status)
            if (.not. status) cycle

            ! check if edge is in transport region
            ! if (.not. par%transport(IDX_EDGE, idx_dir)%contains(idx_edge)) cycle
            if (.not. par%transport(IDX_EDGE, idx_dir)%flags%get(edge_idx)) cycle

            ! get the other vertex of this edge
            call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, edge_idx, 1, idx_l, status)
            if (all(idx_l == idx_k)) then
              call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, edge_idx, 2, idx_l, status)
            end if

            ! get edge properties (use transport surface)
            A_kl = par%tr_surf(idx_dir)%get(edge_idx)
            L_kl = par%g%get_len(edge_idx, idx_dir)

            ! get coordinates of vertex l
            call par%g%get_vertex(idx_l, p_l)

            ! calculate geometric factor
            ! E_k,dir = -sum_l [A_kl * (coord_l - coord_k)_dir / (Omega_k * L_kl)] * (phi_l - phi_k)
            geometric_factor = -0.5 * A_kl * (p_l(dir) - p_k(dir)) / (Omega_k * L_kl)

            ! set Jacobian entries (accumulate contributions)
            call this%jaco_pot%add(idx_k, idx_l, geometric_factor)
            call this%jaco_pot%add(idx_k, idx_k, -geometric_factor)
          end do
        end do
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_efield_eval(this)
      !! evaluate electric field calculation equation
      class(calc_efield), intent(inout) :: this

      real, allocatable :: tmp(:)

      allocate(tmp(this%efield%n))

      call this%jaco_pot%matr%mul_vec(this%pot%get(), tmp)
      call this%efield%set(tmp)

      deallocate(tmp)
  end subroutine

end module

