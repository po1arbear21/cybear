module mobility_m

  use device_params_m, only: device_params, CR_NAME, DIR_NAME
  use dual_m
  use equation_m,      only: equation
  use grid_m,          only: IDX_VERTEX, IDX_EDGE
  use grid_data_m,     only: grid_data2_real
  use imref_m,         only: imref
  use jacobian_m,      only: jacobian
  use stencil_m,       only: near_neighb_stencil
  use variable_m,      only: variable_real

  implicit none

  private
  public :: mobility, calc_mobility

  type, extends(variable_real) :: mobility
    !! electron/hole mobility
    integer       :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => mobility_init
  end type

  type, extends(equation) :: calc_mobility
    !! calculate electron/hole mobility using Caughey-Thomas model

    type(device_params), pointer :: par  => null()
    type(imref),         pointer :: iref => null()
    type(mobility),      pointer :: mob  => null()

    type(near_neighb_stencil) :: st

    type(jacobian), pointer :: jaco_iref => null()
  contains
    procedure :: init => calc_mobility_init
    procedure :: eval => calc_mobility_eval
  end type

contains

  subroutine mobility_init(this, par, ci, idx_dir)
    !! initialize mobility variable
    class(mobility),     intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters
    integer,             intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)
    integer,             intent(in)  :: idx_dir
      !! edge direction

    type(grid_data2_real), pointer :: p

    ! init base
    call this%variable_init(CR_NAME(ci)//"mob"//DIR_NAME(idx_dir), "cm^2/V/s", g = par%g, idx_type = IDX_EDGE, idx_dir = idx_dir)
    this%ci = ci

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
  end subroutine

  subroutine calc_mobility_init(this, par, iref, mob)
    !! initialize mobility calculation equation
    class(calc_mobility),        intent(out) :: this
    type(device_params), target, intent(in)  :: par
      !! device parameters
    type(imref),         target, intent(in)  :: iref
      !! imref variable
    type(mobility),      target, intent(in)  :: mob
      !! mobility variable

    ! init base
    call this%equation_init("calc_"//mob%name)
    this%par  => par
    this%iref => iref
    this%mob  => mob

    ! init stencil
    call this%st%init(par%g, IDX_EDGE, mob%idx_dir, IDX_VERTEX, 0)

    ! provide mobility, depend on imref
    this%jaco_iref => this%init_jaco(this%provide(mob,  par%transport(IDX_EDGE,mob%idx_dir)), &
      &                              this%depend( iref, par%transport(IDX_VERTEX,0)), &
      &                              st = [this%st%get_ptr()], const = .false.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_mobility_eval(this)
    !! evaluate mobility calculation equation
    class(calc_mobility), intent(inout) :: this

    integer, parameter :: eye(2,2) = reshape([1, 0, 0, 1], [2, 2])

    integer      :: i, idx1(2), idx2(2)
    real         :: mob0, v_sat, beta, len
    type(dual_1) :: mob, diref

    call diref%init(1.0, 1)
    do i = 1, this%par%transport(IDX_EDGE,this%mob%idx_dir)%n
      idx1 = this%par%transport(IDX_EDGE,this%mob%idx_dir)%get_idx(i)
      idx2 = idx1 + eye(:,this%mob%idx_dir)

      ! parameters
      mob0  = this%par%mob0(IDX_EDGE,this%mob%idx_dir,this%mob%ci)%get(idx1)
      v_sat = this%par%v_sat(this%mob%ci)
      beta  = this%par%beta( this%mob%ci)
      len   = this%par%g%get_len(idx1, this%mob%idx_dir)

      ! delta imref
      diref%x = this%iref%get(idx2) - this%iref%get(idx1)

      ! calculate mobility including derivative wrt delta imref
      mob = mob0 / (1.0 + ((mob0/(v_sat*len)) * abs(diref))**beta)**(1.0/beta)

      ! extract result and derivative from dual number
      call this%mob%set(idx1, mob%x)
      call this%jaco_iref%set(idx1, idx1, - mob%dx(1))
      call this%jaco_iref%set(idx1, idx2,   mob%dx(1))
    end do
  end subroutine

end module
