module density_m

  use device_params_m, only: device_params, CR_NAME
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data2_real
  use variable_m,      only: variable_real

  implicit none

  private
  public density

  type, extends(variable_real) :: density
    !! electron/hole density
    integer       :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => density_init
  end type

contains

  subroutine density_init(this, par, ci)
    !! initialize density
    class(density),      intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data2_real), pointer :: p => null()

    call this%variable_init(CR_NAME(ci)//"dens", "1/cm^3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
  end subroutine

end module
