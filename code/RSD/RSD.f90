program RSD
  implicit none
  integer, parameter :: N=500
  real*8, parameter :: s_min=-8.d0, s_max=+8.d0 ! cMpc/h
  real*8, parameter :: r_max=+8.d0              ! cMpc/h
  integer :: i,j
  real*8 :: xi_s, s_perp, s_para, ds
  real*8 :: xi, r, dr
  real*8 :: v12,s12

   ! write the real-space correlation function and velocity statistics
  dr=r_max/N
  open(9,file='output/xi_v12_s12.dat')
  do i=1,N
     r=i*dr
     write(9,*) r,xi(r),v12(r),s12(r)
  end do
  close(9)

  ! write the redshift-space correlation function
  open(10,file='output/RSD.dat')
  ds=(s_max-s_min)/N
  do i=1,N
     do j=1,N

        s_perp=s_min+i*ds
        s_para=s_min+j*ds
        r=sqrt(s_perp**2+s_para**2)
        write(10,*) s_perp, s_para, xi_s(s_perp,s_para)

     enddo
     write(10,*) ''
  enddo
  close(10)


end program RSD

function xi_s(s_perp,s_para)
  implicit none
  real*8 :: xi_s,s_perp,s_para
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-3
  real*8 :: result,abserr,y_min,y_max
  integer :: neval,ier
  real*8 :: s1,s2
  common /svar/ s1,s2
  external intg

  s1=s_perp
  s2=s_para
  y_min=-500.d0
  y_max=+500.d0
  call qags(intg, y_min, y_max, epsabs, epsrel, result, abserr, neval, ier)
  xi_s=result-1.d0

  return
end function xi_s

function intg(y)
  use param
  implicit none
  real*8 :: intg,y
  real*8 :: r,mu,f_v12
  real*8 :: xi,v12,s12
  real*8 :: s1,s2
  common /svar/ s1,s2

  r=sqrt(s1**2+y**2)
  mu=y/r
  f_v12=1.d0/(dsqrt(2.d0*pi)*s12(r)/Hubble) * &
        exp( -(s2-y-mu*v12(r)/Hubble)**2/(2.d0*(s12(r)/Hubble)**2) )       
  !f_v12=1.d0/(dsqrt(2.d0)*s12(r)/Hubble) * &
  !      exp( -sqrt(2.0)*abs(s2-y-mu*v12(r)/Hubble)/(s12(r)/Hubble) )       
  intg=(1.d0+xi(r))*f_v12

  return
end function intg

! 2-point galaxy-absorber correlation function 
function xi(r) 
  implicit none
  real*8 :: xi,r
  real*8, parameter :: rc=0.1d0 
  real*8, parameter :: r0=3.d0 ! cMpc/h
  real*8, parameter :: slope=1.74d0

  xi=(r/r0)**(-slope) * ( 1.d0-exp(-(r/rc)**2) )

  return
end function xi

! pairwise velocity field
function v12(r)
  implicit none
  real*8 :: v12,r
  real*8, parameter :: v0=400.d0 ! km/s
  real*8, parameter :: r0=3.d0 ! cMpc/h
  real*8, parameter :: slope=1.74d0
  v12=-v0/(1.d0+(r/r0)**slope)
  return
end function v12

! pairwise velocity dispersion
function s12(r)
  implicit none
  real*8 :: s12,r
  real*8, parameter :: s12_0=300.d0 ! km/s
  s12=s12_0
  return
end function s12
  
