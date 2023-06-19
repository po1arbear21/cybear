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
      !! enable/disable incomplete ionization
    real, allocatable :: edop(:)
      !! dopant energy relative to the carrier band for a single dopant
    real, allocatable :: N_asr(:)
      !! Altermatt-Schenk reference density for dopant energy
    real, allocatable :: asc(:)
      !! Altermatt-Schenk exponent for fitting of the dopant energy
    real, allocatable :: N_asb(:)
      !! Altermatt-Schenk reference density for fraction of bound states in dopant clusters
    real, allocatable :: asd(:)
      !! Altermatt-Schenk exponent for fitting the fraction of bound states
    real, allocatable :: g_dop(:)
      !! degeneracy of dopant states

  end type

end module
