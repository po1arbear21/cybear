module example_contact_m

  use bin_search_m,     only: bin_search
  use example_device_m, only: grd, dop, dop_v, n_intrin
  use grid_m,           only: grid_table, IDX_CELL, IDX_VERTEX
  use input_m,          only: input_file

  implicit none

  private
  public contact, contacts, uncontacted, init_contacts

  type contact
    character(:), allocatable :: name
    real                      :: phi_ms
    type(grid_table)          :: conts
  contains
    procedure :: init => contact_init
  end type

  type(contact), allocatable :: contacts(:)
  type(grid_table)           :: uncontacted

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

  ! getting input for the contact at each section
  allocate(contacts(num_ct))
  do i = 1, num_ct
    call contacts(i)%init(f, sid(i))
  end do

  ! finalize uncontacted
  call uncontacted%init_final()

end subroutine

subroutine contact_init(this, f, sid)
  class(contact),   intent(out) :: this
  type(input_file), intent(in)  :: f
  integer,          intent(in)  :: sid

  integer :: indx
  real    :: x

  ! set name
  call f%get(sid, "name", this%name)

  ! init grid_table
  call this%conts%init(this%name, grd, IDX_VERTEX, 0)

  ! set the flag on x on true
  call f%get(sid, "x", x)
  indx = bin_search(grd%x, x)
  call this%conts%flags%set( [indx], .true.)
  call uncontacted%flags%set([indx], .false.)

  ! finalize grid_table init
  call this%conts%init_final()

  ! set phi_ms
  this%phi_ms = log(dop_v%get([indx])/n_intrin)

end subroutine

end module
