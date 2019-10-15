module string_m
  implicit none

  type string
    character(len=:), allocatable :: s
  end type

end module