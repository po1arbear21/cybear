module carrier_m

  implicit none

  integer,      parameter :: CR_ELEC      = 1
  integer,      parameter :: CR_HOLE      = 2
  character(*), parameter :: CR_NAME(2)   = [ "n", "p" ]
  real,         parameter :: CR_CHARGE(2) = [ -1.0, +1.0 ]

end module
