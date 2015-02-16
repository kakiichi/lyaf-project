# fast computation of effective optical depth
from pylab import *
from scipy.interpolate import interp1d
from scipy.integrate import quad

# parameters
G=6.674e-8                # cgs
c=2.998e10                # cm/s
Mpc2cm=3.0857e24          # 1Mpc=3.0857cm
cross_sec=0.011           # [cm2*Hz]
freq_lya=2.466e15         # Hz lya resonance frequency
wavelength_lya=c/freq_lya # cm

z_s=3.0

# read the oscillator strength and Lyn resonant frequancy
freq_lyn,osci_lyn=genfromtxt('Lyman_parameters.input',usecols=(1,3),unpack=True)

# functions
def line_profile(freq,form,lyn):
    def voigt(T,x,lyn):      # dimensionless voigt profile
        dampingwidth=6.67e-5*(osci_lyn[lyn]/3.)*(0.01*c/freq_lyn[lyn])**(-2) # Hz
        dopplerwidth=4.28e-7*freq_lyn[lyn]*sqrt(T)                               # Hz,  T in [K]
        a=dampingwidth/(4.*pi*dopplerwidth)
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
        dopplerwidth=4.28e-7*freq_lyn[lyn]*sqrt(temperature)                               # Hz,  T in [K]
        x=(freq-freq_lyn[lyn])/dopplerwidth
        line_profile=1./dopplerwidth*voigt(temperature,x,lyn) # [1/Hz]
    return line_profile

# Lya optical depth for a single absorbing cloud
def optdpt(freq,NHI,form,lyn):
    cross_sec=0.026502*osci_lyn[lyn] # cm^2
    optdpt=cross_sec*NHI*line_profile(freq,form,lyn)
    return optdpt

#column density distribution function from Haardt&Madau 2014
def dNdNHIdz(NHI,z): # [cm2]
    A=8.7e18
    betaN=2.0
    betaZ=1.27
    dNdNHIdz=A*NHI**(-betaN)*(1.+z)**betaZ
    return dNdNHIdz

def EW(NHI,form,lyn):
    def voigt(T,x,lyn):      # dimensionless voigt profile
        dampingwidth=6.67e-5*(osci_lyn[lyn]/3.)*(0.01*c/freq_lyn[lyn])**(-2) # Hz
        dopplerwidth=4.28e-7*freq_lyn[lyn]*sqrt(T)                               # Hz,  T in [K]
        a=dampingwidth/(4.*pi*dopplerwidth)
        xx=x * x
        y=(xx-0.855)/(xx+3.42)
        if y <= 0.0:
            q=0.0
        elif y > 0.0:
            q=y*(1.+21./xx)*(a/(pi*(xx+1.0)))*(0.1117+y*(4.421+y*(-9.207+5.674*y)))
        voigt=q+exp(-xx)/sqrt(pi)
        return voigt
    T=1.e4                                      # K
    cross_sec=0.026502*osci_lyn[lyn]            # cm^2.Hz
    dopplerwidth=4.28e-7*freq_lyn[lyn]*sqrt(T)  # Hz,  T in [K]
    integrand=lambda x: ( 1.-exp(-cross_sec*NHI/dopplerwidth*voigt(T,x,lyn) ) )
    EW=quad(integrand, -inf, inf, limit=100)[0]
    return EW

n=100
lyn_max=10
logNHI_min=12.0
logNHI_max=21.0
NHI=logspace(logNHI_min,logNHI_max,n)
F=zeros([n,lyn_max])
for i in range(n):
    for lyn in range(lyn_max):
        F[i,lyn]=EW(NHI[i],'Voigt',lyn)

# Lya effective optical depth
def optdpt_eff(freq,Nmin,Nmax,z_s,form,lyn):
    T=1.e4                                      # K
    dopplerwidth=4.28e-7*freq_lyn[lyn]*sqrt(T)  # Hz,  T in [K]
    F_EW=interp1d(NHI,F[:,lyn])
    integrand=lambda NHI : dNdNHIdz(NHI,freq_lyn[lyn]/freq*(1.+z_s)-1.)*F_EW(NHI)*(1.+z_s)*dopplerwidth/freq
    optdpt_eff=quad(integrand,Nmin,Nmax,limit=100)[0]
    return optdpt_eff

# main: plot
wavelength_lyn=c/freq_lyn*1.e8

#wavelength_min=1213.5 * 1.e-8  # cm
wavelength_min=912.0 * 1.e-8  # cm
wavelength_max=1220.0 * 1.e-8  # cm
npoints=100
wavelength=linspace(wavelength_min,wavelength_max,npoints)
freq=c/wavelength
vel=1.e-5*c*(1.0-wavelength_lya/wavelength)  # km/s

Nmin=10.**(18.0) # 10**18
Nmax=10.**(21.0)

tau_eff=zeros([npoints,lyn_max])
print 'computing the optical depth from Lya-absorbing clouds ...'
for lyn in range(lyn_max):
    for i in range(npoints):
        tau_eff[i,lyn]=optdpt_eff(freq[i],Nmin,Nmax,z_s,'Voigt',lyn)
        print lyn,i,freq[i],vel[i],tau_eff[i,lyn]

tau_eff_total=zeros(npoints)
for i in range(npoints):
    for lyn in range(lyn_max):
        if freq[i]>=freq_lyn[lyn]:
            tau_eff_total[i]=tau_eff_total[i]+tau_eff[i,lyn]

# output
file=open('tau_eff_Lyn.output','w')
for i in range(npoints):
    output=str(wavelength[i]/1.e-8)+'\t'+str(freq[i])+'\t'+str(vel[i])+'\t'+str(exp(-tau_eff_total[i]))+'\n'
    file.write(output)
file.close()
