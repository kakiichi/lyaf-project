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
damp_coef=6.25e8          # 1/s

z_s=3.0

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
    if form == 'Lorentz':
        line_profile=(damp_coef/(4.*pi**2))/((freq-freq_lya)**2+(damp_coef/(4.*pi))**2) # [1/Hz]
    elif form == 'Doppler':
        temperature=1.e4 # K
        doppler_width=1.057e9*sqrt(temperature) # Hz
        x=(freq-freq_lya)/doppler_width
        line_profile=1./(sqrt(pi)*doppler_width)*exp(-x*x) # 1/Hz
    elif form == 'Rayleigh':
        f=freq/freq_lya
        line_profile=(damp_coef/(4.*pi**2))*f**4/((freq-freq_lya)**2+(damp_coef/(4.*pi))**2*f**6) # [1/Hz]
    elif form == 'Voigt':
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
def dNdNHIdz(NHI,z): # [cm2]
    A=8.7e18
    betaN=2.0
    betaZ=1.27
    dNdNHIdz=A*NHI**(-betaN)*(1.+z)**betaZ
    return dNdNHIdz

# Lya effective optical depth
def optdpt_eff(freq,Nmin,Nmax,z_s,form):
    def result(NHI,freq,z_s,form):
        integrand=lambda nu: (1.-exp(-optdpt(nu,NHI,form)))*(1.+z_s)/freq
        result=quad(integrand, 0.0, freq_lya+freq_lya*1000,limit=100)[0]
        return result
    integrand=lambda NHI : dNdNHIdz(NHI,freq_lya/freq*(1.+z_s)-1.)*result(NHI,freq,z_s,form)
    optdpt_eff=quad(integrand,Nmin,Nmax,limit=100)[0]
    return optdpt_eff

# main: plot
#wavelength_min=1213.5 * 1.e-8  # cm
wavelength_min=1000.0 * 1.e-8  # cm
wavelength_max=1220.0 * 1.e-8  # cm
npoints=100
wavelength=linspace(wavelength_min,wavelength_max,npoints)
freq=c/wavelength
vel=1.e-5*c*(1.0-wavelength_lya/wavelength)  # km/s

Nmin=10.**(18.0)
Nmax=10.**(21.0)

tau_eff=zeros(npoints)
print 'computing the optical depth from Lya-absorbing clouds ...'
for i in range(npoints):
    tau_eff[i]=optdpt_eff(freq[i],Nmin,Nmax,z_s,'Voigt')
    print i,freq[i],vel[i],tau_eff[i]

# output
file=open('tau_eff_global.output','w')
for i in range(npoints):
    output=str(wavelength[i]/1.e-8)+'\t'+str(freq[i])+'\t'+str(vel[i])+'\t'+ \
           str(tau_eff[i])+'\t'+str(exp(-tau_eff[i]))+'\n'
    file.write(output)
file.close()
