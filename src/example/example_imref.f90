module example_imref_m

  use example_contact_m,   only: contacts, uncontacted
  use example_density_m,   only: dens
  use example_device_m,    only: dop, grd, n_intrin
  use grid_m,              only: grid_data1_real, IDX_VERTEX
  use variable_m,          only: variable

  implicit none

  private
  public iref

  type, extends(variable) :: imref
  !! quasi-fermi-potential
    real, pointer :: x(:) => null()
  contains
    procedure :: init => imref_init
    procedure :: calc => imref_calc
  end type

  type(imref)     :: iref

contains

  subroutine imref_init(this)
    class(imref), intent(out) :: this

    type(grid_data1_real), pointer :: p

    call this%variable_init("imref", "V", g = grd, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine imref_calc(this)
    class(imref), intent(inout) :: this

    real, allocatable :: dens_arr(:)

    allocate(dens_arr(size(grd%x)))
    dens_arr = dens%get()
    call this%set(-log(dens_arr/n_intrin)+pot%get())
  end subroutine

end module
