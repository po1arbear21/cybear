module contact_m

  implicit none

  private
  public :: CT_OHMIC, CT_GATE, contact

  ! contact types
  integer, parameter :: CT_OHMIC = 1
  integer, parameter :: CT_GATE  = 2

  type contact
    !! device contact

    character(:), allocatable :: name
      !! contact name
    integer                   :: type
      !! type of contact (CT_OHMIC, CT_GATE)
    real                      :: phims
      !! metal-semiconductor workfunction difference
  end type

end module
