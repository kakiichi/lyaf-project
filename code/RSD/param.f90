module param
  implicit none

  real*8,parameter :: om_m0=0.26
  real*8,parameter :: om_v0=0.74

  real*8,parameter :: redshift=3.0
  real*8,parameter :: Hubble=100.d0*dsqrt(om_v0+om_m0*(1.d0+redshift)**3) ! km/s/(Mpc/h)

  real*8,parameter :: pi=3.14159265359d0

  ! model parameters
  real*8 :: r0, slope_xi                   ! [cMpc/h], [-] 
  real*8 :: v_inflow, r_inflow, slope_v12  ! [km/s], [cMpc/h], [-] 
  real*8 :: v_outflow, r_outflow           ! [km/s], [cMpc/h]
  real*8 :: s12_0                          ! [km/s]
  real*8 :: r_eq, betaN                    ! [cMpc/h], [-]


end module param
