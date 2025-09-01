module schottky_m

  use device_params_m,  only: device_params
  use normalization_m,  only: norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE

  implicit none

  private
  public :: schottky_injection_mb, schottky_velocity

contains

  subroutine schottky_injection_mb(par, ci, ict, ninj)
    !! Calculate equilibrium density n0 for Schottky contact Robin BC
    !! This is the thermionic emission density at zero bias
    !! Electrons: n0 = Nc * exp(-phi_Bn)
    !! Holes:     p0 = Nv * exp(-phi_Bp) where phi_Bp = E_g - phi_Bn

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(out) :: ninj    ! Equilibrium density n0 (normalized)

    real :: phi_Bn, phi_Bp

    ! Get normalized barrier height (already converted in device_params)
    phi_Bn = par%contacts(ict)%phi_b

    if (ci == CR_ELEC) then
      ! Electrons: n0 = Nc * exp(-phi_Bn)
      ninj = par%smc%edos(CR_ELEC) * exp(-phi_Bn)
    else  ! CR_HOLE
      ! Holes: barrier from valence band
      phi_Bp = par%smc%band_gap - phi_Bn
      ninj = par%smc%edos(CR_HOLE) * exp(-phi_Bp)
    end if
  end subroutine

  function schottky_velocity(par, ci, ict) result(s)
    !! Calculate thermionic emission velocity at Schottky contact
    !! Using Richardson constant: v_surf = A*T^2/(q*Nc)
    !! The normalization system handles the elementary charge q

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! Carrier index
    integer,             intent(in) :: ict  ! Contact index
    real                            :: s

    ! Check if Richardson constant is provided and > 0
    if (par%contacts(ict)%A_richardson > 0.0) then
      ! Calculate normalized surface velocity
      ! v_surf = A*T^2/(q*Nc) where q is handled by normalization
      ! T must be normalized, Nc is already normalized
      s = par%contacts(ict)%A_richardson * norm(par%T, "K") * norm(par%T, "K") / par%smc%edos(ci)

      ! Debug output
      print *, "DEBUG schottky_velocity for contact ", ict, ":"
      print *, "  A_richardson = ", par%contacts(ict)%A_richardson, " A/cm^2/K^2"
      print *, "  Temperature = ", par%T, " K"
      print *, "  Nc (normalized) = ", par%smc%edos(ci)
      print *, "  v_surf (normalized) = ", s
      print *, "  v_surf (physical) = ", denorm(s, "cm/s"), " cm/s"
    else
      ! Default thermal velocity estimate (v_th/4)
      s = 0.25  ! v_th/4 in normalized units
    end if
  end function

end module schottky_m
