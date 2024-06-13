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

    real, allocatable :: genrec_tau(:)
      !! generation/recombination time constants (carrier index)

    logical           :: incomp_ion
      !! enable/disable incomplete ionization (Altermatt-Schenk model)
    real, allocatable :: ii_E_dop0(:)
      !! dopant energy relative to the carrier band for a single dopant
    real, allocatable :: ii_g(:)
      !! degeneracy factor of dopant states
    real, allocatable :: ii_N_crit(:)
      !! cricital concentration (alternative to altermatt-schenk)
    real, allocatable :: ii_dop_th(:)
      !! full ionization if doping > threshold
  end type

end module
