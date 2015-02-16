from pylab import *
import scipy.integrate
import scipy.special

# load parameter file
param=genfromtxt('param.input',comments='#')
h=param[0]
om_v0=param[1]
om_m0=param[2]
om_b0=param[3]
z=param[4]
xi_on=int(param[5])
r0=param[6]
slope=param[7]
phr_on=int(param[8])
req=param[9]
betaN=param[10]
v0=param[11]
s0=param[12]
outflow_on=int(param[13])
v_out=param[14]
r_out=param[15]
form=param[16]

def g_v12(v):
    def Hubble(z):
        Hubble=100.0*sqrt(om_v0+om_m0*(1.+z)**3) # (km/s)/(Mpc/h)
        return Hubble
    def xi(r):
        if xi_on==1:
            xi=(r/r0)**(-slope)
        if xi_on==0:
            xi=0.
        return xi
    def PDF_r(r):
        if phr_on==1:
            factor=( (r/req)**(-2)+1. )**(-betaN+1)
        if phr_on==0:
            factor=1.0
        PDF_r=factor*(1.+xi(r))
        return PDF_r
    def v12(r):
        if outflow_on==1:
            v12=v_out*2.0*exp(-r/r_out)-v0/(1.0+(r/r0)**slope)
        if outflow_on==0:
            v12=-v0/(1.0+(r/r0)**slope)
        return v12
    def s12(r):
        s12=s0
        return s12
    def PDF_v(v,r):
        PDF_v=( 1./(sqrt(2.*pi)*s12(r)) )*exp( -(v-v12(r))**2/(2.*s12(r)**2) )
        return PDF_v
    def f_v12(v):
        if form == 1:    # Gaussian
            r_min=0.001  # comoving radius min. and max. limit of integral
            r_max=100.0  # cMpc/h
            integrand=lambda r: PDF_r(r)*PDF_v(v-Hubble(z)*r/(1.+z),r)
            result,err=scipy.integrate.quad(integrand,r_min,r_max,limit=1000)
            f_v12=result*Hubble(z)/(1.+z)
        if form == 0:    # Hubble only
            f_v12=PDF_r( v*(1.+z)/Hubble(z) )
        return f_v12

    return f_v12(v)

# main
N=1000
v=linspace(-2000,2000,N)
f=zeros([N])
for i in range(N):
    f[i]=g_v12(v[i])

# output
file=open('f_v12.output','w')
for i in range(N):
    file.write( str(v[i])+'\t'+str(f[i])+'\n' )
file.close()

