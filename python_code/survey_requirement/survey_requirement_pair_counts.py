from scipy.integrate import quad
execfile('cosmo_toolkits.py')

# cosmological parameter for distance measure
h=0.70
om_m0=0.3
om_v0=0.7

# survey volume [(cMpc/h)^3]
def survey_volume(z1,z2,deg):
    from scipy.integrate import quad
    FoV=(deg*(pi/180.))**2 # FoV in radian squared
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)*angular_distance(z)**2
    survey_volume=FoV*quad(integrand,z1,z2,limit=100)[0]
    return survey_volume

execfile('galaxy_counts.py')
execfile('CDDF.py')


# galaxy-absorber correlation function
def corrfunc(r):
    slope=-1.74 # Cooke+2003 best fit DLA-LBG correlation function z~3
    r0=3.32 # cMpc/h
    corrfunc=0.0 #(r/r0)**slope
    return corrfunc


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

Nr=20
r,dr=linspace(0.1,10.0,Nr,retstep=True)
N_pairs=zeros([Nr,Np,3,3])
for n in range(Nr):
    for i in range(Np):
        J0=quad(lambda r: 4.*pi*r**2*(1.+0.0),         r[n]-0.5*dr, r[n]+0.5*dr)[0]
        J1=quad(lambda r: 4.*pi*r**2*(1.+corrfunc(r)), r[n]-0.5*dr, r[n]+0.5*dr)[0]
        N_absorber=Number_counts_absorber(17.0,23.0,z1,z2,'4powlaw')

        N_pairs[n,i,0,0]=N_absorber*N_LAE[i,0,0]/survey_volume(z1,z2,deg)*J0
        N_pairs[n,i,1,0]=N_absorber*N_LAE[i,0,0]/survey_volume(z1,z2,deg)*J1
        #N_pairs[n,i,2,0]=N_LAE[i,0,0]/survey_volume(z1,z2,deg)*J

# plot
from matplotlib import gridspec

fig = plt.figure(figsize=(7, 9)) 
plt.rc('font',**{'family':'serif', 'size':15})
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
#plt.rc('legend', fontseize=15)

gs = gridspec.GridSpec(3, 2)

subplots_adjust(hspace=0.00)
ax1=subplot(gs[0,0])

ax1.semilogy(r,N_pairs[:,10,0,0],'k--',drawstyle='steps-mid',linewidth=1)
ax1.semilogy(r,N_pairs[:,10,1,0],'k-',drawstyle='steps-mid',linewidth=1)
ax1.semilogy(r,N_pairs[:,5,1,0],'k-',drawstyle='steps-mid',linewidth=1)
ax1.semilogy(r,N_pairs[:,15,1,0],'k-',drawstyle='steps-mid',linewidth=1)
ax1.semilogy(r,N_pairs[:,20,1,0],'k-',drawstyle='steps-mid',linewidth=1)

#ax1.loglog(F_lim,N_pairs[:,10,2,0],'k:',linewidth=2)
ax1.xaxis.set_ticklabels([]) #xticks([])

ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
text(2.5e-17,3e3,'$<z>=3.1$ galaxies',rotation='horizontal',fontsize=15)
text(2.5e-17,1e3,'$z_Q=3.5$ QSO field',rotation='horizontal',fontsize=15)

ax2=subplot(gs[1,0])
ax2.semilogy(r,N_pairs[:,10,0,0],'k--',drawstyle='steps-mid',linewidth=2)
ax2.semilogy(r,N_pairs[:,10,1,0],'k-',drawstyle='steps-mid',linewidth=2)
#ax2.loglog(F_lim,N_LAE[:,2,1],'k:',linewidth=2)
ax2.xaxis.set_ticklabels([]) #xticks([])
#
ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
#text(2.5e-17,3e3,'$<z>=3.7$ galaxies',rotation='horizontal',fontsize=15)
#text(2.5e-17,1e3,'$z_Q=4.5$ QSO field',rotation='horizontal',fontsize=15)
#
ax3=subplot(gs[2,0])
ax3.semilogy(r,N_pairs[:,10,0,0],'k--',drawstyle='steps-mid',linewidth=2)
ax3.semilogy(r,N_pairs[:,10,1,0],'k-',drawstyle='steps-mid',linewidth=2)
#ax3.loglog(F_lim,N_LAE[:,2,2],'k:',linewidth=2)
#
ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
xlabel('Pair seperation $r$ [$h^{-1}$cMpc]')
#vlines(F_limit_Subaru,1,2e4,linestyles='dotted')
#vlines(F_limit_Magellan,1,2e4,linestyles='dashed')
#text(1.1e-17,1e4,'Subaru',rotation='vertical',fontsize=12)
#text(2.1e-18,1e4,'Magellan',rotation='vertical',alpha=1.0,fontsize=12)
#text(2.5e-17,3e3,'$<z>=5.7$ galaxies',rotation='horizontal',fontsize=15)
#text(2.5e-17,1e3,'$z_Q=6.1$ QSO field',rotation='horizontal',fontsize=15)

subplots_adjust(wspace=0.50)
ax4=subplot(gs[0,1])
ax4.loglog(F_lim,N_pairs[5,:,0,0],'k--',drawstyle='steps-mid',linewidth=2)
ax4.loglog(F_lim,N_pairs[5,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax4.loglog(F_lim,N_pairs[1,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax4.loglog(F_lim,N_pairs[10,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax4.xaxis.set_ticklabels([]) #xticks([])
ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
ylim(1e-2,1e2)

ax5=subplot(gs[1,1])
ax5.loglog(F_lim,N_pairs[5,:,0,0],'k--',drawstyle='steps-mid',linewidth=2)
ax5.loglog(F_lim,N_pairs[5,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax5.loglog(F_lim,N_pairs[1,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax5.loglog(F_lim,N_pairs[10,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax5.xaxis.set_ticklabels([]) #xticks([])
ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
ylim(1e-2,1e2)

ax6=subplot(gs[2,1])
ax6.loglog(F_lim,N_pairs[5,:,0,0],'k--',drawstyle='steps-mid',linewidth=2)
ax6.loglog(F_lim,N_pairs[5,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax6.loglog(F_lim,N_pairs[1,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ax6.loglog(F_lim,N_pairs[10,:,1,0],'k-',drawstyle='steps-mid',linewidth=2)
ylabel(' $N_{pairs}(r|>F_{lim})\\times$FoV$_{900}$')
xlabel('Flux limit $F_{lim}$ $[erg/s/cm^2]$')
ylim(1e-2,1e2)
