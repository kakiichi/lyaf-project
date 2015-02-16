module param
  implicit none

  real*8,parameter :: om_m0=0.3
  real*8,parameter :: om_v0=0.7

  real*8,parameter :: redshift=3.0
  real*8,parameter :: Hubble=100.d0*dsqrt(om_v0+om_m0*(1.d0+redshift)**3) ! km/s/(Mpc/h)

  real*8,parameter :: pi=3.14159265359d0

  ! astrophysical model parameters
  real*8 :: r0
  real*8 :: slope
end module param
