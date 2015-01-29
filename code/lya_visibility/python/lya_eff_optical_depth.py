# compute the effective optical depth with pairwise velocity statistics
import scipy.integrate

# parameters
G=6.674e-8                # cgs
c=2.998e10                # cm/s

cross_sec=0.011           # [cm2*Hz]
freq_lya=2.466e15         # Hz lya resonance frequency
wavelength_lya=c/freq_lya # cm
damp_coef=6.25e8          # 1/s

Mpc2cm=3.0857e24          # 1Mpc=3.0857cm

h=0.72
om_v0=0.73
om_m0=0.27
om_b0=0.044
Y=0.25                                  # He abundance in mass fraction
H0=3.240e-18*h                          # 1/s

def g_v12(v,form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out):
###  r0=3.32    #cMpc/h correlation length
###  slope=1.74 # slope of correlation function
###  
###  req=1.0
###  betaN=2.
###
###  v0=200.
###  r0=3.32
###  slope=1.74
###  
###  s0=200.
###
    def Hubble(z):
        Hubble=100.0*sqrt(om_v0+om_m0*(1.+z)**3) # (km/s)/(Mpc/h)
        return Hubble
    def xi(r):
        if xi_on==True:
            xi=(r/r0)**(-slope)
        if xi_on==False:
            xi=0.
        return xi
    def PDF_r(r):
        if phr_on==True:
            factor=( (r/req)**(-2)+1. )**(-betaN+1)
        if phr_on==False:
            factor=1.0
        PDF_r=factor*(1.+xi(r))
        return PDF_r
    def v12(r):
        if outflow_on==True:
            v12=v_out*2.0*exp(-r/r_out)-v0/(1.0+(r/r0)**slope)
        if outflow_on==False:
            v12=-v0/(1.0+(r/r0)**slope)
        return v12
    def s12(r):
        s12=s0
        return s12
    def PDF_v(v,r):
        PDF_v=( 1./(sqrt(2.*pi)*s12(r)) )*exp( -(v-v12(r))**2/(2.*s12(r)**2) )
        return PDF_v
    def f_v12(v,form):
        if form == 'Gaussian':
            r_min=0.001  # comoving radius min. and max. limit of integral
            r_max=100.0 # cMpc/h
            integrand=lambda r: PDF_r(r)*PDF_v(v-Hubble(z_s)*r/(1.+z_s),r)
            result,err=scipy.integrate.quad(integrand,r_min,r_max)
            f_v12=result*Hubble(z_s)/(1.+z_s)
        if form == 'Hubble':
            f_v12=PDF_r( v*(1.+z_s)/Hubble(z_s) )
        return f_v12

    return f_v12(v,form)

# line profile
def line_profile(freq,form):
    def voigt(T,x):     # dimensionless voigt profile
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
def optdpt_eff(freq,beta_max,Nmin,Nmax,z_s,form):
    
    def result(NHI,freq,beta_max,form):
        integrand=lambda beta: 1.-exp(-optdpt(freq*(1.0-beta),NHI,form))
        result=scipy.integrate.quad(integrand, 0.0, beta_max,limit=100)[0]
        return result

    integrand=lambda NHI : dNdNHIdz(NHI,z_s)*(1.+z_s)*result(NHI,freq,beta_max,form)
    optdpt_eff=scipy.integrate.quad(integrand,Nmin,Nmax,limit=100)[0]
    return optdpt_eff

# Lya effective optical depth with velocity space clustering
def optdpt_eff_clustering(freq,beta_max,Nmin,Nmax,z_s,form,v_form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out):
    
    def result(NHI,freq,beta_max,form,v_form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out):
        integrand=lambda beta: g_v12(beta*c*1e-5,v_form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out) \
                   *( 1.-exp(-optdpt(freq*(1.0-beta),NHI,form)) )
        result=scipy.integrate.quad(integrand, -beta_max, beta_max,limit=100)[0]
        return result

    integrand=lambda NHI : dNdNHIdz(NHI,z_s)*(1.+z_s)*result(NHI,freq,beta_max,form, \
                                                             v_form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out)
    optdpt_eff_clustering=scipy.integrate.quad(integrand,Nmin,Nmax,limit=100)[0]
    return optdpt_eff_clustering

# main: plot
# parameters for your interest
z_s=3.0

wavelength_min=1213.5 * 1.e-8  # cm
wavelength_max=1220.0 * 1.e-8  # cm
npoints=100
wavelength=linspace(wavelength_min,wavelength_max,npoints)
freq=c/wavelength
vel=1.e-5*c*(1.0-wavelength_lya/wavelength)  # km/s

Nmin=10.**(18.0)
Nmax=10.**(21.0)
beta_max=0.01

optdpt_total=zeros(npoints)
print 'computing the optical depth from Lya-absorbing clouds ...'
for i in range(npoints):
    optdpt_total[i]=optdpt_eff(freq[i],beta_max,Nmin,Nmax,z_s,'Voigt')
    print i,freq[i],vel[i],optdpt_total[i]

xi_on=True
r0=3.32     # cMpc/h, correlation length
slope=1.74  # power-law slope of correlation function
phr_on=True
req=0.32     # cMpc/h, local-global UV background equality radius
betaN=2.0   # CDDF slope
v0=200.0    # km/s
s0=200.0    # km/s
outflow_on=False
v_out=0.0
r_out=0.0

optdpt_total_clustering=zeros(npoints)
print 'computing the optical depth from Lya-absorbing clouds ...'
for i in range(npoints):
    optdpt_total_clustering[i]=optdpt_eff_clustering(freq[i],beta_max,Nmin,Nmax,z_s,'Voigt', \
                                                     'Gaussian',xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out)
    print i,freq[i],vel[i],optdpt_total_clustering[i]


figure()
plt.rc('font',**{'family':'serif', 'size':13})
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=12)

ylim(-0.1,1.05)
xlim(-500,1100)
ylabel('transmission $e^{-\\tau}$',fontsize=15)
xlabel('$\Delta v$ [km/s]',fontsize=15)

plot(vel,exp(-optdpt_total),'b-',linewidth=4,label='$\\tau_{\\rm{eff}}$')
plot(vel,exp(-optdpt_total_clustering),'r-',linewidth=4,label='$\\tau_{\\rm{eff}}$')

lg=plt.legend(loc='lower right',fontsize=13)
lg.draw_frame(False)

figure()
plot(vel,1.+exp(-(vel-400)**2/(2.*200**2)),'k-')
plot(vel,1.+exp(-(vel-400)**2/(2.*50**2)),'k--')
plot(vel,exp(-optdpt_total_clustering)*(1.+exp(-(vel-400)**2/(2.*200**2))),'r-')
plot(vel,exp(-optdpt_total_clustering)*(1.+exp(-(vel-400)**2/(2.*50**2))),'r--')
plot(vel,exp(-optdpt_total)*(1.+exp(-(vel-400)**2/(2.*200**2))),'b-')
plot(vel,exp(-optdpt_total)*(1.+exp(-(vel-400)**2/(2.*50**2))),'b--')

