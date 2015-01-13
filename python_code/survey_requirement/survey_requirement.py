execfile('cosmo_toolkits.py')

# cosmological parameter for distance measure
h=0.70
om_m0=0.3
om_v0=0.7

# Ouchi+ 2008 LAE luminosity function
def dndL(L,z):

    if z==3.1:
        phi_c=9.2e-4 # 1/cMpc3
        Lc=5.8e42    #  erg/s
        slope=-1.5

    dndL=phi_c/Lc*(L/Lc)**slope*exp(-L/Lc)

    return dndL

# number density of LAE down to L_lim luminosity [(cMpc/h)^-3]
def n_LAE(L_lim,z,fit):
    import mpmath
    if z==3.1:
        if fit=='1p5':
            phi_c=9.2e-4/h**3 # (h/cMpc)^3
            Lc=5.8e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=14.9e-4/h**3 # (h/cMpc)^3
            Lc=4.1e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=3.9e-4/h**3 # (h/cMpc)^3
            Lc=9.1e42         #  erg/s
            slope=-2.0
    if z==3.7:
        if fit=='1p5':
            phi_c=3.4e-4/h**3 # (h/cMpc)^3
            Lc=10.2e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=5.7e-4/h**3 # (h/cMpc)^3
            Lc=7.2e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=1.3e-4/h**3 # (h/cMpc)^3
            Lc=16.2e42         #  erg/s
            slope=-2.0
    if z==5.7:
        if fit=='1p5':
            phi_c=7.7e-4/h**3 # (h/cMpc)^3
            Lc=6.8e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=10.4e-4/h**3 # (h/cMpc)^3
            Lc=5.4e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=3.6e-4/h**3 # (h/cMpc)^3
            Lc=9.5e42         #  erg/s
            slope=-2.0


    n_LAE=phi_c*mpmath.fp.gammainc(slope+1.0, L_lim/Lc)
    return n_LAE

# limiting luminosity [erg/s] corresponding to flux limit [erg/s/cm2]
def L_limit(F_lim,z):
    Mpc2cm=3.08567758e24                # 1Mpc=3.08567758e24cm
    DL=luminosity_distance(z)/h*Mpc2cm  # Mpc/h to cm
    L_limit=4.*pi*DL**2*F_lim
    return L_limit

def Number_counts(F_lim,deg,z1,z2,z_fixed,fit):
    from scipy.integrate import quad
    FoV=(deg*(pi/180.))**2 # FoV in radian squared
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)*angular_distance(z)**2*n_LAE(L_limit(F_lim,z),z_fixed,fit)
    Number_counts=FoV*quad(integrand,z1,z2,limit=100)[0]
    return Number_counts

wavelength_lya=1216.0 # [Angstrom]
wavelength_lyb=1026.0 # [Angstrom]

deg=0.5 # 30arcmin FoV # hypothetical QSO field 

Np=50
F_lim=logspace(-19,-15,Np)
N_LAE=zeros([Np,3,3])

z_fixed=3.1 # galaxy redshift
z_QSO=3.5  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LAE[i,0,0]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p0')
    N_LAE[i,1,0]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p5')
    N_LAE[i,2,0]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'2p0')

z_fixed=3.7
z_QSO=4.5  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LAE[i,0,1]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p0')
    N_LAE[i,1,1]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p5')
    N_LAE[i,2,1]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'2p0')

z_fixed=5.7
z_QSO=6.1  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LAE[i,0,2]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p0')
    N_LAE[i,1,2]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'1p5')
    N_LAE[i,2,2]=Number_counts(F_lim[i],deg,z1,z2,z_fixed,'2p0')

F_limit_Subaru  =1.6e-17 # erg/s/cm2 Subaru SDF flux limit
F_limit_Magellan=3.0e-18 # erg/s/cm2 Magellan/IMACS flux limit

# plot
from matplotlib import gridspec

fig = plt.figure(figsize=(6.5, 9)) 
plt.rc('font',**{'family':'serif', 'size':15})
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
#plt.rc('legend', fontseize=15)

gs = gridspec.GridSpec(3, 1)

subplots_adjust(hspace=0.00)
ax1=subplot(gs[0])
ax1.loglog(F_lim,N_LAE[:,0,0],'k--',linewidth=2)
ax1.loglog(F_lim,N_LAE[:,1,0],'k-',linewidth=2)
ax1.loglog(F_lim,N_LAE[:,2,0],'k:',linewidth=2)
ax1.xaxis.set_ticklabels([]) #xticks([])

ylabel(' $N_{LAE}(>F_{lim})\\times$FoV$_{900}$')
vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(2.5e-17,3e3,'$<z>=3.1$ galaxies',rotation='horizontal',fontsize=15)
text(2.5e-17,1e3,'$z_Q=3.5$ QSO field',rotation='horizontal',fontsize=15)
ylim(2,2e4)

ax2=subplot(gs[1])
ax2.loglog(F_lim,N_LAE[:,0,1],'k--',linewidth=2)
ax2.loglog(F_lim,N_LAE[:,1,1],'k-',linewidth=2)
ax2.loglog(F_lim,N_LAE[:,2,1],'k:',linewidth=2)
ax2.xaxis.set_ticklabels([]) #xticks([])

ylabel(' $N_{LAE}(>F_{lim})\\times$FoV$_{900}$')
vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(2.5e-17,3e3,'$<z>=3.7$ galaxies',rotation='horizontal',fontsize=15)
text(2.5e-17,1e3,'$z_Q=4.5$ QSO field',rotation='horizontal',fontsize=15)

ylim(2,2e4)

ax3=subplot(gs[2])
ax3.loglog(F_lim,N_LAE[:,0,2],'k--',linewidth=2)
ax3.loglog(F_lim,N_LAE[:,1,2],'k-',linewidth=2)
ax3.loglog(F_lim,N_LAE[:,2,2],'k:',linewidth=2)

ylabel(' $N_{LAE}(>F_{lim})\\times$FoV$_{900}$')
xlabel('Flux limit $F_{lim}$ $[erg/s/cm^2]$')
vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(1.1e-17,1e4,'Subaru',rotation='vertical',fontsize=12)
text(2.1e-18,1e4,'Magellan',rotation='vertical',alpha=1.0,fontsize=12)
text(2.5e-17,3e3,'$<z>=5.7$ galaxies',rotation='horizontal',fontsize=15)
text(2.5e-17,1e3,'$z_Q=6.1$ QSO field',rotation='horizontal',fontsize=15)
ylim(2,2e4)

