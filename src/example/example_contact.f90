module example_contact_m

  use bin_search_m,      only: bin_search
  use example_current_m, only: current
  use example_device_m,  only: dop, dop_v, grd, n_intrin
  use example_voltage_m, only: voltage
  use grid_m,            only: grid_data1_int, grid_table, IDX_CELL, IDX_VERTEX
  use input_m,           only: input_file

  implicit none

  private
  public init_contacts
  public contacts, uncontacted, grd_contacts

  type contact
    character(:), allocatable :: name
    real                      :: phi_ms
    type(voltage)             :: volt
    type(current)             :: curr
    type(grid_table)          :: conts
    integer                   :: idx(1)
  contains
    procedure :: init => contact_init
  end type

  type(contact), allocatable :: contacts(:)
  type(grid_table)           :: uncontacted
  type(grid_data1_int)       :: grd_contacts

contains

subroutine init_contacts(f)
  type(input_file), intent(in) :: f
    !! input file

  integer              :: i, num_ct
  integer, allocatable :: sid(:)

  ! get different sections of contact
  call f%get_sections("contact", sid)
  num_ct = size(sid)

  ! initialise uncontacted with default True
  call uncontacted%init("uncontacted", grd, IDX_VERTEX, 0, initial_flags = .true.)
  call grd_contacts%init(grd, IDX_VERTEX, 0)

  ! getting input for the contact at each section
  allocate(contacts(num_ct))
  do i = 1, num_ct
    call contacts(i)%init(f, sid(i))
    call grd_contacts%set(contacts(i)%idx, i)
  end do

  ! finalize uncontacted
  call uncontacted%init_final()
end subroutine

subroutine contact_init(this, f, sid)
  class(contact),   intent(out) :: this
  type(input_file), intent(in)  :: f
  integer,          intent(in)  :: sid

  real :: x, v

  ! set name
  call f%get(sid, "name", this%name)

  ! initialise voltage and current
  call this%volt%init(this%name)
  call this%curr%init(this%name)

  ! set voltage
  call f%get(sid, "V", v)
  call this%volt%set([v])

  ! init grid_table
  call this%conts%init(this%name, grd, IDX_VERTEX, 0)

  ! set the flag on x on true
  call f%get(sid, "x", x)
  this%idx(1) = bin_search(grd%x, x)
  call this%conts%flags%set( this%idx, .true.)
  call uncontacted%flags%set(this%idx, .false.)

  ! finalize grid_table init
  call this%conts%init_final()

  ! set phi_ms
  this%phi_ms = log(dop_v%get(this%idx) / n_intrin)
end subroutine

end module
