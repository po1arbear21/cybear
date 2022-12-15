module density_m

  use device_params_m, only: device_params, CR_NAME
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use variable_m,      only: variable_real
  use error_m,         only: assert_failed, program_error

  implicit none

  private
  public density

  type, extends(variable_real) :: density
    !! electron/hole density
    integer       :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    real, pointer :: x1(:)     => null()
    real, pointer :: x2(:,:)   => null()
    real, pointer :: x3(:,:,:) => null()
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
    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    call this%variable_init(CR_NAME(ci)//"dens", "1/cm^3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
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

end module
