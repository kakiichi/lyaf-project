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

# HI CDDF various best fit function (modified to be integrated by quad)
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
        logk3=-23.83
        logN3=21.50
        a3=-2.00
        a4=-6.00
        k3=10.**logk3
        N3=10.**logN3
        for i in range(len(N)):
            if N < N3:
                CDDF=k3*(N/N3)**a3
            else:
                CDDF=k3*(N/N3)**a4
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
    # Peroux et al 2003 best fit HI CDDF gamma function
    elif form == 'Peroux_z2.35':
        fc=3.25e-2
        beta=1.08
        logNc=21.27
        Nc=10.**(logNc)
        CDDF=(fc/Nc)*(N/Nc)**(-beta)*exp(-N/Nc)
    elif form == 'Peroux_z3.1':
        fc=4.06e-2
        beta=1.10
        logNc=21.18
        Nc=10.**(logNc)
        CDDF=(fc/Nc)*(N/Nc)**(-beta)*exp(-N/Nc)
    elif form == 'Peroux_z3.9':
        fc=25.1e-2
        beta=0.80
        logNc=20.46
        Nc=10.**(logNc)
        CDDF=(fc/Nc)*(N/Nc)**(-beta)*exp(-N/Nc)

    # Kim et al 2002 best fit HI CDDF
    # single-power law
    elif form == 'Kim_z2.1':
        logA=11.00
        beta=1.71
        CDDF=10.**(logA)*N**(-beta)
    elif form == 'Kim_z3.3':
        logA=9.18
        beta=1.55
        CDDF=10.**(logA)*N**(-beta)
    elif form == 'Kim_z3.8':
        logA=7.77
        beta=1.44
        CDDF=10.**(logA)*N**(-beta)
    else:
        print 'The funtional form must be specified.'
    return CDDF


# column density distribution function in dNdNHIdz
def dNdNHIdz(NHI,z,form):
    H0=1.02275e-10 # [h/yr]
    dNdNHIdz=CDDF(NHI,form)*H0*(1.+z)**2/Hubble(z)
    return dNdNHIdz

# absorber number counts
def Number_counts_absorber(logNHI_min,logNHI_max,z1,z2,form):
    from scipy.integrate import quad
    z_integral=quad(lambda z: H0*(1.+z)**2/Hubble(z), z1, z2)[0]
    CDDF_integral=quad(lambda logNHI: log(10.)*(10.**logNHI)*CDDF((10.**logNHI),form),logNHI_min,logNHI_max)[0]
    Number_counts_absorber=z_integral*CDDF_integral
    return Number_counts_absorber


# check
#N=logspace(12.3,22.0,100)
#DLAindex=find(log10(N)>=20.3)
#LLSindex=find((log10(N)>=17.0)&(log10(N)<20.3))
#LLSDLAindex=find(log10(N)>=17.2)
#Lyafindex=find(log10(N)<17.0)

                   
# plot
# plot in velocity space
#figure()
#plt.rc('font',**{'family':'serif', 'size':13})
#plt.rc('xtick', labelsize=13)
#plt.rc('ytick', labelsize=13)
#plt.rc('legend', fontsize=12)
#
##ylim(-26,-10)
##xlim(-500,1100)
#ylabel('$\log_{10} d^2N/dN_{HI}dz [cm^2]$',fontsize=15)
#xlabel('$\log_{10} N_{HI} [cm^{-2}]$',fontsize=15)
#
## O'Meara fits at z~2.4
#z=2.4
#plot(log10(N[:]),log10(N[:]*dNdNHIdz(N[:],2.4,'4powlaw')),'k-',linewidth=2,label='$z\sim2.4$ O\'Meara+ (2013)')
##plot(log10(N[DLAindex]),log10(dNdNHIdz(N[DLAindex],3.06,'double')),'g-.',linewidth=1,label='$z\sim3.06$ Prochaska+ (2005)')
##plot(log10(N[LLSDLAindex]),log10(dNdNHIdz(N[LLSDLAindex],2.35,'Peroux_z2.35')), 'm:',linewidth=1,label='$z\sim2.35$ Peroux+ (2005)')
##plot(log10(N[LLSDLAindex]),log10(dNdNHIdz(N[LLSDLAindex],3.1,'Peroux_z3.1')), 'm-.',linewidth=1,label='$z\sim3.1$ Peroux+ (2005)')
##plot(log10(N[LLSDLAindex]),log10(dNdNHIdz(N[LLSDLAindex],3.9,'Peroux_z3.9')), 'm--',linewidth=1,label='$z\sim3.9$ Peroux+ (2005)')
#
#
#
##text(700, 0.93,'$\\tau(<R_b)$ only')
##text(700, 0.65,'$R_b=1.0$pMpc')
#lg=plt.legend(loc='lower left',fontsize=8)
#lg.draw_frame(False)

