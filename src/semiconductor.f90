module semiconductor_m

  implicit none

  integer,      parameter :: CR_ELEC      = 1
  integer,      parameter :: CR_HOLE      = 2
  character(*), parameter :: CR_NAME(2)   = [ "n", "p" ]
  real,         parameter :: CR_CHARGE(2) = [ -1.0, 1.0 ]

  integer,      parameter :: DOP_DCON      = 1
  integer,      parameter :: DOP_ACON      = 2
  character(*), parameter :: DOP_NAME(2)   = [ "D", "A" ]
  real,         parameter :: DOP_CHARGE(2) = [ 1.0, -1.0]

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

    real :: rec_tau(2,2)
      !! generation/recombination time constants (generation/recombination, carrier index)

    logical           :: incomp_ion
      !! enable/disable incomplete ionization (Altermatt-Schenk model)
    real, allocatable :: ii_E_dop0(:)
      !! dopant energy relative to the carrier band for a single dopant
    real, allocatable :: ii_N_ref(:)
      !! dopant energy reference density for dopant energy
    real, allocatable :: ii_c(:)
      !! dopant energy exponent
    real, allocatable :: ii_N_b(:)
      !! reference density for fraction of bound states in dopant clusters
    real, allocatable :: ii_d(:)
      !! exponent for fitting the fraction of bound states
    real, allocatable :: ii_g(:)
      !! degeneracy factor of dopant states
  end type

end module
