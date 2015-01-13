execfile('cosmo_toolkits.py')

# cosmological parameter for distance measure
h=0.70
om_m0=0.3
om_v0=0.7


# Number density of LBG down to M_lim [AB mag]
# using Bouwen+2014 fit
# n_LBG [(cMpc/h)^-3]
def n_LBG(M_lim,z,fit):
    import mpmath
    if fit == 'Bouwens':
        slope=-1.85-0.09*(z-6.)
        M_c=-20.89+0.12*(z-6.)                        # AB mag of UV
        phi_c=0.48*10.**(-0.19*(z-6.))*10.**(-3)/h**3 # (h/cMpc)^3
    n_LBG=phi_c*mpmath.fp.gammainc( slope+1.0, 10.**(-0.4*(M_lim-M_c)) )
    return n_LBG

## limiting luminosity [erg/s] corresponding to flux limit [erg/s/cm2]
#def L_limit(F_lim,z):
#    Mpc2cm=3.08567758e24                # 1Mpc=3.08567758e24cm
#    DL=luminosity_distance(z)/h*Mpc2cm  # Mpc/h to cm
#    L_limit=4.*pi*DL**2*F_lim
#    return L_limit

def Number_counts(M_lim,deg,z1,z2,fit):
    from scipy.integrate import quad
    FoV=(deg*(pi/180.))**2 # FoV in radian squared
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)*angular_distance(z)**2*n_LBG(M_lim,z,fit)
    Number_counts=FoV*quad(integrand,z1,z2,limit=100)[0]
    return Number_counts

wavelength_lya=1216.0 # [Angstrom]
wavelength_lyb=1026.0 # [Angstrom]

deg=0.5 # 30arcmin FoV # hypothetical QSO field 

Np=50
M_lim=linspace(-23,-15,Np)
N_LBG=zeros([Np,3,3])

z_QSO=3.5  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LBG[i,0,0]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,1,0]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,2,0]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')

z_QSO=4.5  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LBG[i,0,1]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,1,1]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,2,1]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')

z_QSO=6.1  # redshift of background QSO
dz_QSO=(1.+z_QSO)*(wavelength_lya-wavelength_lyb)/wavelength_lya
z1=z_QSO-dz_QSO
z2=z_QSO
for i in range(Np):
    N_LBG[i,0,2]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,1,2]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')
    N_LBG[i,2,2]=Number_counts(M_lim[i],deg,z1,z2,'Bouwens')

#F_limit_Subaru  =1.6e-17 # erg/s/cm2 Subaru SDF flux limit
#F_limit_Magellan=3.0e-18 # erg/s/cm2 Magellan/IMACS flux limit

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
ax1.semilogy(M_lim,N_LBG[:,0,0],'k--',linewidth=2)
ax1.semilogy(M_lim,N_LBG[:,1,0],'k-',linewidth=2)
ax1.semilogy(M_lim,N_LBG[:,2,0],'k:',linewidth=2)
ax1.xaxis.set_ticklabels([]) #xticks([])

ylabel(' $N_{LBG}(>M_{lim})\\times$FoV$_{900}$')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(-22.5,3e3,'$<z>=3.1$ galaxies',rotation='horizontal',fontsize=15)
text(-22.5,1e3,'$z_Q=3.5$ QSO field',rotation='horizontal',fontsize=15)
ylim(2,2e4)

ax2=subplot(gs[1])
ax2.semilogy(M_lim,N_LBG[:,0,1],'k--',linewidth=2)
ax2.semilogy(M_lim,N_LBG[:,1,1],'k-',linewidth=2)
ax2.semilogy(M_lim,N_LBG[:,2,1],'k:',linewidth=2)
ax2.xaxis.set_ticklabels([]) #xticks([])

ylabel(' $N_{LBG}(>M_{lim})\\times$FoV$_{900}$')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(-22.5,3e3,'$<z>=3.7$ galaxies',rotation='horizontal',fontsize=15)
text(-22.5,1e3,'$z_Q=4.5$ QSO field',rotation='horizontal',fontsize=15)

ylim(2,2e4)

ax3=subplot(gs[2])
ax3.semilogy(M_lim,N_LBG[:,0,2],'k--',linewidth=2)
ax3.semilogy(M_lim,N_LBG[:,1,2],'k-',linewidth=2)
ax3.semilogy(M_lim,N_LBG[:,2,2],'k:',linewidth=2)

ylabel(' $N_{LBG}(>M_{lim})\\times$FoV$_{900}$')
xlabel('Flux limit $M_{UV,lim}$ [AB mag]')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
#text(1.1e-17,1e4,'Subaru',rotation='vertical',fontsize=12)
#text(2.1e-18,1e4,'Magellan',rotation='vertical',alpha=1.0,fontsize=12)
text(-22.5,3e3,'$<z>=5.7$ galaxies',rotation='horizontal',fontsize=15)
text(-22.5,1e3,'$z_Q=6.1$ QSO field',rotation='horizontal',fontsize=15)
ylim(2,2e4)

