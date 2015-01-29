# python script to plot Lya red damping wing due to the collection of many lya-absorbing clouds

import scipy.integrate

# parameters
G=6.674e-8                # cgs
c=2.998e10                # cm/s
mH=1.673e-24              # g

cross_sec=0.011           # [cm2*Hz]
freq_lya=2.466e15         # Hz lya resonance frequency
wavelength_lya=c/freq_lya # cm
damp_coef=6.25e8          # 1/s

Mpc2cm=3.0857e24          # 1Mpc=3.0857cm

h=0.72
om_m=0.27
om_b=0.044
Y=0.25                                  # He abundance in mass fraction
H0=3.240e-18*h                          # 1/s
nH0=(3.*H0**2)/(8.*pi*G)*(1.-Y)*om_b/mH # comoving hydrogen number density 1/cm3

def Hubble(z):
    Hubble=H0*sqrt(om_m*(1.+z)**3)
    return Hubble

def corrfunc(beta,z):
    corr_slope=1.74
    r0=3.32        # cMpc/h
    l0=r0/h/(1.+z) # pMpc
    beta0=Hubble(z)*(l0*Mpc2cm)/c
    lmin=0.1       # pMpc (=300pkpc) CGM scale
    beta_min=Hubble(z)*(lmin*Mpc2cm)/c
    if beta > beta_min:
        corrfunc=(beta/beta0)**(-corr_slope)
    else:
        corrfunc=0.0   
    #print r0,l0,beta0
    return corrfunc

# Prochaska et al 2005 HI CDDF various best fit function
def CDDF(N,form):
    # single-power law
    if form == 'single':
        logk1=23.16
        a1=-2.19
        k1=10.**logk1
        CDDF=k1*N**a1
    # Gamma function
    elif form == 'gamma':
        logk2=-23.52
        logN2=21.48
        a2=-1.80
        k2=10.**logk2
        N2=10.**logN2
        CDDF=k2*(N/N2)**a2*exp(-N/N2)
    # double power-law
    elif form == 'double':
        CDDF=zeros(len(N))
        logk3=-23.83
        logN3=21.50
        a3=-2.00
        a4=-6.00
        k3=10.**logk3
        N3=10.**logN3
        for i in range(len(N)):
            if N[i] < N3:
                CDDF[i]=k3*(N[i]/N3)**a3
            else:
                CDDF[i]=k3*(N[i]/N3)**a4
    elif form =='4powlaw':
        k12=10**(-9.52)
        b12=-1.67
        b17p5=-1.07
        b20p3=-1.71
        b21p5=-11.10
        N12p0=10.**12.0
        N17p5=10.**17.5
        k17p5=k12*(N17p5/N12p0)**b12
        N20p3=10.**20.3
        k20p3=k17p5*(N20p3/N17p5)**b17p5
        N21p5=10.**21.5
        k21p5=k20p3*(N21p5/N20p3)**b20p3
        if log10(N) < 17.5 :
            CDDF=k12*(N/N12p0)**b12
        elif (log10(N) >= 17.5) & (log10(N) < 20.3) :
            CDDF=k17p5*(N/N17p5)**b17p5
        elif (log10(N) >= 20.3) & (log10(N) < 21.5) :
            CDDF=k20p3*(N/N20p3)**b20p3
        elif log10(N) >= 21.5 :
            CDDF=k21p5*(N/N21p5)**b21p5
    else:
        print 'The funtional form must be specified.'
    return CDDF


# dimensionless voigt profile
def voigt(T,x):
    a=0.0472/sqrt(T) # T in [K]
    xx=x * x
    y=(xx-0.855)/(xx+3.42)
    if y <= 0.0:
        q=0.0
    elif y > 0.0:
        q=y*(1.+21./xx)*(a/(pi*(xx+1.0)))*(0.1117+y*(4.421+y*(-9.207+5.674*y)))
    
    voigt=q+exp(-xx)/sqrt(pi)
    return voigt

# line profile
def line_profile(freq,form):
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

# lya optical depth for single absorber at Rabs
def optdpt_selfsheild(freq,NHI,form,Rabs,z_s):
    freq_obs=freq*(1.0-Hubble(z_s)*Rabs*Mpc2cm/c)
    optdpt_selfsheild=optdpt(freq_obs,NHI,form)
    return optdpt_selfsheild

#column density distribution function from Becker & Bolton 2013
def dNdNHIdz(NHI,z,CDDFform):
    if CDDFform=='Becker13':
        A=0.93
        betaN=1.33
        betaZ=1.92
        N_LL=10.**(17.2)
        dNdNHIdz=(A/N_LL)*(NHI/N_LL)**(-betaN)*((1.+z)/4.5)**betaZ
    elif CDDFform=='OMeara13':
        dNdNHIdz=CDDF(NHI,'4powlaw')*H0*(1.+z)**2/Hubble(z)
    return dNdNHIdz

# Lya effective optical depth with local approximation
def optdpt_eff_local(freq,beta_max,Nmin,Nmax,z_s,form,CDDFform):
    
    def result(NHI,freq,beta_max,form):
        lmin=0.0       # pMpc (=300pkpc) CGM scale
        beta_min=Hubble(z_s)*(lmin*Mpc2cm)/c

        integrand=lambda beta: 1.-exp(-optdpt(freq*(1.0-beta),NHI,form))
        result=scipy.integrate.quad(integrand, beta_min, beta_max,limit=1000)[0]
        return result
    
    integrand=lambda NHI : dNdNHIdz(NHI,z_s,CDDFform)*(1.+z_s)*result(NHI,freq,beta_max,form)
    optdpt_eff_local=scipy.integrate.quad(integrand,Nmin,Nmax,limit=1000)[0]
    return optdpt_eff_local

# Lya effective optical depth with local approximation
def optdpt_eff_local_clustering(freq,beta_max,Nmin,Nmax,z_s,form,CDDFform):
    
    def result(NHI,freq,beta_max,form):
        lmin=0.0       # pMpc (=300pkpc) CGM scale
        beta_min=Hubble(z_s)*(lmin*Mpc2cm)/c
 
        integrand=lambda beta: (1.-exp(-optdpt(freq*(1.0-beta),NHI,form)))*(1.+corrfunc(beta,z_s))
        result=scipy.integrate.quad(integrand, beta_min, beta_max,limit=1000)[0]
        return result
    
    integrand=lambda NHI : dNdNHIdz(NHI,z_s,CDDFform)*(1.+z_s)*result(NHI,freq,beta_max,form)
    optdpt_eff_local=scipy.integrate.quad(integrand,Nmin,Nmax,limit=1000)[0]
    return optdpt_eff_local

# main: plot
# parameters for your interest
z_s=6.6

wavelength_min=1213.5 * 1.e-8  # cm
wavelength_max=1220.0 * 1.e-8  # cm
npoints=10
vel=linspace(100,700,npoints)  # km/s
wavelength=wavelength_lya/(1.-vel*.1e5/c)
freq=c/wavelength

Nmin=10.**(19.0) #10.**(15.0)
Nmax=logspace(20,22.0,3)
beta_max=0.1

optdpt_totalBB=zeros([3,npoints])
optdpt_totalOM=zeros([3,npoints])
optdpt_totalBB_clustering=zeros([3,npoints])
optdpt_totalOM_clustering=zeros([3,npoints])
print 'computing the optical depth from Lya-absorbing clouds ...'
for n in range(npoints):
    for i in range(3):
        optdpt_totalBB[i,n]=optdpt_eff_local(freq[n],beta_max,Nmin,Nmax[i],z_s,'Voigt','Becker13')
        optdpt_totalOM[i,n]=optdpt_eff_local(freq[n],beta_max,Nmin,Nmax[i],z_s,'Voigt','OMeara13')
        optdpt_totalBB_clustering[i,n]=optdpt_eff_local_clustering(freq[n],beta_max,Nmin,Nmax[i],z_s,'Voigt','Becker13')
        optdpt_totalOM_clustering[i,n]=optdpt_eff_local_clustering(freq[n],beta_max,Nmin,Nmax[i],z_s,'Voigt','OMeara13')
        #print optdpt_eff_local_clustering(freq[n],beta_max,Nmax[0],Nmax[i],z_s,'Voigt','Becker13')

        print n,i,vel[n],Nmax[i],optdpt_totalOM[i,n],optdpt_totalOM_clustering[i,n],optdpt_totalOM_clustering[i,n]/optdpt_totalOM[i,n]




figure()
plt.rc('font',**{'family':'serif', 'size':13})
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=12)
#yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])2014MNRAS.438..529R

#ylim(-0.05,1.15)
#xlim(-500,1100)
ylabel('$\\tau_{\\rm{eff}}(N_{HI}<N_{HI}^{max})/\\tau_{\\rm{eff}}(N_{HI}<N_{HI}^{max}=10^{21}cm^{-2})$',fontsize=15)
xlabel('$\\Delta v$ $[km/s]$',fontsize=15)

#plot(vel,optdpt_totalBB[0,:],'b--',linewidth=4,label='BB13 CDDF, $\Delta v=-300km/s$')
#plot(vel,optdpt_totalBB[1,:],'k--',linewidth=4,label='BB13 CDDF, $\Delta v=0km/s$')
#plot(vel,optdpt_totalBB[2,:],'r--',linewidth=4,label='BB13 CDDF, $\Delta v=+300km/s$')
#plot(vel,optdpt_totalBB_clustering[0,:],'b--',linewidth=4,label='BB13 CDDF, $\Delta v=-300km/s$')
#plot(vel,optdpt_totalBB_clustering[1,:],'k--',linewidth=4,label='BB13 CDDF, $\Delta v=0km/s$')
#plot(vel,optdpt_totalBB_clustering[2,:],'r--',linewidth=4,label='BB13 CDDF, $\Delta v=+300km/s$')

#plot(vel,optdpt_totalOM[0,:],'bo',linewidth=4,label='BB13 CDDF, $\Delta v=-300km/s$')
#plot(vel,optdpt_totalOM[1,:],'ko',linewidth=4,label='BB13 CDDF, $\Delta v=0km/s$')
#plot(vel,optdpt_totalOM[2,:],'ro',linewidth=4,label='BB13 CDDF, $\Delta v=+300km/s$')
#plot(vel,optdpt_totalOM_clustering[0,:],'bs',linewidth=4,label='BB13 CDDF, $\Delta v=-300km/s$')
#plot(vel,optdpt_totalOM_clustering[1,:],'ks',linewidth=4,label='BB13 CDDF, $\Delta v=0km/s$')
#plot(vel,optdpt_totalOM_clustering[2,:],'rs',linewidth=4,label='BB13 CDDF, $\Delta v=+300km/s$')

plot(vel,optdpt_totalOM_clustering[0,:]/optdpt_totalOM[0,:],'bo',linewidth=4,label='BB13 CDDF, $\Delta v=-300km/s$')
plot(vel,optdpt_totalOM_clustering[1,:]/optdpt_totalOM[1,:],'ko',linewidth=4,label='BB13 CDDF, $\Delta v=0km/s$')
plot(vel,optdpt_totalOM_clustering[2,:]/optdpt_totalOM[2,:],'ro',linewidth=4,label='BB13 CDDF, $\Delta v=+300km/s$')



lg=plt.legend(loc='upper left',fontsize=10)
lg.draw_frame(False)

