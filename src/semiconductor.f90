module semiconductor_m

  implicit none

  integer,      parameter :: CR_ELEC      = 1
  integer,      parameter :: CR_HOLE      = 2
  character(*), parameter :: CR_NAME(2)   = [ "n", "p" ]
  real,         parameter :: CR_CHARGE(2) = [ -1.0, +1.0 ]

  type semiconductor
    real              :: n_intrin
      !! intrinsic carrier density
    real              :: edos(2)
      !! effective density of states Nc, Nv
    real              :: band_gap
      !! band gap
    real              :: band_edge(2)
      !! conduction/valence band edge
    real, allocatable :: mass(:)
      !! effective mass
    logical           :: degen
      !! degenerate case: Fermi-Dirac statistics
    real, allocatable :: alpha(:)
      !! Caughey-Thomas alpha parameter
    real, allocatable :: beta(:)
      !! Caughey-Thomas beta parameter
    real, allocatable :: mob_min(:)
      !! Caughey-Thomas minimal mobility
    real, allocatable :: mob_max(:)
      !! Caughey-Thomas maximal mobility
    real, allocatable :: N_ref(:)
      !! Caughey-Thomas reference density
    real, allocatable :: v_sat(:)
      !! Caughey-Thomas saturation velocity
  end type

end module
