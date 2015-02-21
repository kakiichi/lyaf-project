# computation of effective optical depth
import sys
from pylab import *
from scipy.interpolate import interp1d
from scipy.integrate import quad,romberg

# parameters
G=6.674e-8                # cgs
c=2.998e5                 # km/s
Mpc2cm=3.0857e24          # 1Mpc=3.0857cm
cross_sec=0.011           # [cm2*Hz]
freq_lya=2.466e15         # Hz lya resonance frequency
wavelength_lya=(1.0e5*c)/freq_lya # cm
damp_coef=6.25e8          # 1/s

z_s=3.0

# min. and max. HI column density input
if sys.argv[3] == 'LAF':
    Nmin=10.**14.0
    Nmax=10.**17.0
if sys.argv[3] == 'LLS':
    Nmin=10.**17.0
    Nmax=10.**20.3
if sys.argv[3] == 'DLA':
    Nmin=10.**20.3
    Nmax=10.**22.0

# load velocity PDF
infile=sys.argv[1]
print 'input velocity PDF from ... ', infile
data=genfromtxt(infile,comments='#')
vmin=min(data[:,0])
vmax=max(data[:,0])
g_v12_data=interp1d(data[:,0],data[:,1])

def g_v12(v):
    if v<vmin:
       g_v12=0.0
    elif v>vmax:
        g_v12=1.0
    else:
        g_v12=g_v12_data(v)
    return g_v12

# functions
def line_profile(freq,form):
    def voigt(T,x):      # dimensionless voigt profile
        a=0.0472/sqrt(T) # T in [K]
        xx=x * x
        y=(xx-0.855)/(xx+3.42)
        if y <= 0.0:
            q=0.0
        elif y > 0.0:
            q=y*(1.+21./xx)*(a/(pi*(xx+1.0)))*(0.1117+y*(4.421+y*(-9.207+5.674*y)))
        voigt=q+exp(-xx)/sqrt(pi)
        return voigt
    if form == 'Voigt':
        temperature=1.e4 # K
        doppler_width=1.057e9*sqrt(temperature) # Hz
        x=(freq-freq_lya)/doppler_width
        line_profile=1./doppler_width*voigt(temperature,x) # [1/Hz]
    return line_profile

# Lya optical depth for a single absorbing cloud
def optdpt(freq,NHI,form):
    optdpt=cross_sec*NHI*line_profile(freq,form)
    return optdpt

#column density distribution function from Haardt&Madau 2014
#def dNdNHIdz(NHI,z): # [cm2]
#    A=8.7e18
#    betaN=2.0
#    betaZ=1.27
#    dNdNHIdz=A*NHI**(-betaN)*(1.+z)**betaZ
#    return dNdNHIdz

def dNdNHIdz(NHI,z):
    if   ( (log10(NHI) >= 11.0) & (log10(NHI) < 15.0) ):
        A=1.2e7
        beta=1.5
        gamma=3.0
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 15.0) & (log10(NHI) < 17.5) ):
        A=3.8e14
        beta=2.0
        gamma=3.0
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 19.0) & (log10(NHI) < 20.3) ):
        A=0.45
        beta=1.05
        gamma=1.27
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 20.3) & (log10(NHI) < 21.55) ):
        A=8.7e18
        beta=2.0
        gamma=1.27
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 17.5) & (log10(NHI) < 19.0) ):
        A1=3.8e14
        beta1=2.0
        gamma1=3.0
        C1=A1*(1.+z)**gamma1

        A2=0.45
        beta2=1.05
        gamma2=1.27
        C2=A2*(1.+z)**gamma2
            
        logy1=log10(C1)-17.5*beta1
        logy2=log10(C2)-19.0*beta2
        
        beta=(logy1-logy2)/(19.0-17.5)
        C=(10.**logy1)*(10.**17.5)**beta

        dNdNHIdz=C*NHI**(-beta)
    else:
        dNdNHIdz=0.0
        
    return dNdNHIdz


def EW_eff(NHI,freq,beta_max,form):
    #integrand=lambda beta: g_v12(beta*c)*( 1.-exp(-optdpt(freq*(1.0-beta),NHI,form)) )
    #EW_eff1=quad(integrand,  -0.01, -0.001,limit=100)[0]
    #EW_eff2=quad(integrand, -0.001,   0.0, limit=100)[0]
    #EW_eff3=quad(integrand,    0.0, +0.001,limit=100)[0]
    #EW_eff4=quad(integrand, +0.001,  +0.01,limit=100)[0]
    #EW_eff5=quad(integrand, +0.01,    +0.1,limit=100)[0]
    #EW_eff6=quad(integrand, +0.1,    +0.5,limit=100)[0]
    #EW_eff=EW_eff1+EW_eff2+EW_eff3+EW_eff4+EW_eff5+EW_eff6

    integrand=lambda x: g_v12( c*(1.-freq_lya/freq*(1.+x)) ) * ( 1.-exp(-optdpt(freq_lya*(1.0+x),NHI,form)) )
    xmin=-1.0
    xmax=2.0*freq/freq_lya-1.0

    EW_eff=freq_lya/freq*( quad(integrand, -1.0  , -1.e-1,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, -1.e-1, -1.e-2,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, -1.e-2, -1.e-3,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, -1.e-4, -1.e-5,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, -1.e-5, +1.e-5,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, +1.e-5, +1.e-4,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, +1.e-4, +1.e-3,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, +1.e-3, +1.e-2,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] + \
                           quad(integrand, +1.e-1,   xmax,limit=100,epsabs=1.e-7,epsrel=1.e-7)[0] )   

    return EW_eff

# Lya effective optical depth
def optdpt_eff(freq,beta_max,Nmin,Nmax,z_s,form):
    W_eff=interp1d(NHI,W_eff_data)

    # linear space integration
    #integrand=lambda NHI : dNdNHIdz(NHI,z_s)*(1.+z_s)*result1(NHI,freq,beta_max,form)
    #optdpt_eff=quad(integrand,Nmin,Nmax,limit=100)[0]
    # logspace integration
    zz=freq_lya/freq*(1.+z_s)-1.
    intgrd=lambda lnNHI : exp(lnNHI)*dNdNHIdz(exp(lnNHI),zz)*(1.+z_s)*W_eff(exp(lnNHI))
    optdpt_eff=quad(intgrd,log(Nmin),log(Nmax),limit=100)[0]

    return optdpt_eff

# main: plot
wavelength_LL1=912.0  * 1.e-8  # cm
wavelength_LL2=1200.0 * 1.e-8  # cm
wavelength_min=1201.0 * 1.e-8  # cm
wavelength_max=1220.0 * 1.e-8  # cm
#npoints=100 # 100
#wavelength=linspace(wavelength_min,wavelength_max,npoints)
npoints1=100 # 100
npoints2=100 # 100
npoints=npoints1+npoints2
wavelength1=linspace(wavelength_LL1,wavelength_LL2,npoints1)
wavelength2=linspace(wavelength_min,wavelength_max,npoints2)
wavelength=zeros(npoints1+npoints2)
for i in range(npoints1):
    wavelength[i]=wavelength1[i]
for i in range(npoints2):
    wavelength[i+npoints1]=wavelength2[i]

freq=1.e5*c/wavelength
vel=c*(1.0-wavelength_lya/wavelength)  # km/s

N_NHI=100
#Nmin=10.**(14.0)
#Nmax=10.**(21.55)
NHI=logspace(log10(Nmin),log10(Nmax),N_NHI)

#beta_max=vmax/c #0.01
beta_max=5000.0/c #0.01
#beta_max=2000.0/c #0.01

tau_eff=zeros(npoints)
print 'computing the optical depth from Lya-absorbing clouds ...'
for i in range(npoints):

    # precompute the effective equivalent width
    W_eff_data=zeros(N_NHI)
    for j in range(N_NHI):
        W_eff_data[j]=EW_eff(NHI[j],freq[i],beta_max,'Voigt')

    # compute the effective optical depth
    tau_eff[i]=optdpt_eff(freq[i],beta_max,Nmin,Nmax,z_s,'Voigt')
    print i,freq[i],vel[i],tau_eff[i]



# output
outfile=sys.argv[2]
print 'output effective Lya opacity from ... ', outfile
file=open(outfile,'w')
for i in range(npoints):
    output=str(wavelength[i]/1.e-8)+'\t'+str(freq[i])+'\t'+str(vel[i])+'\t'+ \
           str(tau_eff[i])+'\t'+str(exp(-tau_eff[i]))+'\n'
    file.write(output)
file.close()
