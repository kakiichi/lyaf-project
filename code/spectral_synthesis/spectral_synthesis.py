from scipy.interpolate import interp1d
from scipy.integrate import quad


# parameter 
c=2.998e10            # cm/s
wavelength_lya=1216.0 # Angstrom, Lyman-alpha line wavelength
freq_lya=c/(wavelength_lya*1e-8)      # ~2.47e15Hz, Lyman-alpha line frequency
Lsun2cgs=3.827e33     # 1Lsun=3.827e33 erg/s

## model galaxy parameter ##
# stellar 
metal='z020'    # metal=z040(Z=0.04),z020(Z=0.02),z008(Z=0.008),z004(Z=0.004),z001(Z=0.001)
mass=1.0e11    # [Msun] stellar mass
tage=1.e9      # [yr] age
SFH='dec'      # SFH='const' or 'dec': constant and exponential decline profile exp(-t/tSF)
tSF=1.e9       # [yr] star formation time scale
EBV=0.05        # E(B-V) dust extinction
# HII region of ISM
EW=0.1        # HII region EW of lya line
# HI region of ISM
T=100.         # HI ISM temperature
a=0.0472/sqrt(T)
tau0=1.0e8
doppler_width=1.057e9*sqrt(T) # Hz

############################

def find_nearest(array,value):
    idx = (abs(array-value)).argmin()
    return array[idx],idx

# Lya RT (static spherical symmetry solution by Dijkstra+2006)
def Ja(x,a,tau0):
    Ja=sqrt(pi/24.)/(a*tau0)*x*x/(1.+cosh(sqrt(2.*pi**3/27.)*abs(x**3)/(a*tau0)))
    return Ja

# load SED library
path='/home/koki/work/MBII-project/BPASS/SEDS/'
type='sin'       # choose single stellar population or binary type='sin' or type='bin'

damp=loadtxt(path+'sed.bpass.instant.cloudy.'+type+'.'+metal)
wavelength=damp[:,0]            # Angstrom
frequency=c/(1.e-8*wavelength)  # Hz
SSP=damp[:,1:]/1.e6             # Lsun/Angstrom/stellar mass (BPASS is per 1.e6 stellar mass)
age=logspace(6.0,10.0,41)       # yr

def SFR(t):
    if SFH=='const':
        SFR=mass/tage
    if SFH=='dec':
        SFR=(mass/tSF)*exp(-t/tSF)/(1.-exp(-tage/tSF))
    return SFR

# spectral synthesis    
Nt=1000
dt=tage/Nt
SED=zeros(wavelength.size)
for i in range(Nt+1):
    tt=i*dt
    age_idx=find_nearest(age,tt)[1]
    SED[:]=SED[:]+SFR(tt)*SSP[:,age_idx]*dt
spectrum=SED*wavelength/frequency*Lsun2cgs # erg/s/Hz
flux=interp1d(wavelength,spectrum)     # erg/s/Hz

# subtract Lya line from spectral synthesis
UV1500_idx=find_nearest(wavelength,1500.)
lya_idx=find_nearest(wavelength,wavelength_lya)
spectrum_no_lya=zeros(wavelength.size)
spectrum_no_lya[:]=spectrum[:]
spectrum_no_lya[lya_idx[1]]=spectrum[lya_idx[1]+1]
flux_no_lya=interp1d(wavelength,spectrum_no_lya)     # erg/s/Hz

# add Lya RT
lya_idx=find_nearest(wavelength,wavelength_lya)

La=SED[lya_idx[1]]*EW*Lsun2cgs # Lya luminosity erg/s
norm=quad(lambda x: Ja(x,a,tau0),-1000,1000)[0]

def lya_spectrum(freq):   # erg/s/Hz
    x=(freq-freq_lya)/doppler_width
    lya_spectrum=La/doppler_width*Ja(x,a,tau0)/norm
    return lya_spectrum

N=100000
wavelength=linspace(800,1700,N)
wavelength_optical=linspace(3500,5000,N)
frequency=c/(1.e-8*wavelength)  # Hz
lya_flux=zeros(N)
for i in range(N):
    lya_flux[i]=lya_spectrum(frequency[i])

# include IGM transmission factor
path='/home/koki/work/lyaf-project/code/lya_visibility/python/'
tau_data=genfromtxt(path+'tau_eff_with_clustering.output')
TIGM_eff=interp1d(tau_data[:,0],tau_data[:,4],bounds_error=False,fill_value=1.0)

# include Calzetti dust extinction law
# Calzetti dust extinction law (Extinction = E(B-V), wavelength in Angstrom )
def Calzetti(Extinction,wavelength):
    w=wavelength/1.e4 # Anstrom to micron
    #if ((wavelength >= 1200) and (wavelength<6300)):
    if ((wavelength >= 500) and (wavelength<6300)): # simple extrapolation to 500A
        k=2.656*(-2.156+1.509/w-0.198/w**2+0.011/w**3)+4.88
    elif ((wavelength >= 6300) and (wavelength<16000)):
        k=2.656*(-2.310+1.315/w)+4.88
    else:
        k=0.0
    Calzetti=10.**(-0.4*Extinction*k)
    return Calzetti

dust=zeros(N)
for i in range(N):
    dust[i]=Calzetti(EBV,wavelength[i])

####################################
# plot nicely
plt.rc('font',**{'family':'serif', 'size':19})
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('legend', fontsize=12)

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def plotinplot():
    fig = plt.figure(figsize=(15,7))
    axes = []
    subpos = [0.54,0.45,0.4,0.4]
    x = np.linspace(-np.pi,np.pi)
    for i in range(1):
        axes.append(fig.add_subplot(1,1,i))
    for axis in axes:
        axis.set_xlim(850,1700)
        axis.set_ylim(0,7)      
        ylabel('Relative flux $f_\\nu/f_\\nu(1500\\AA)$',fontsize=22)
        xlabel('Wavelength [$\\AA$]',fontsize=22)
        axis.plot(wavelength,(lya_flux+flux_no_lya(wavelength))/(flux_no_lya(1500)),'k-')
        axis.plot(wavelength,dust*TIGM_eff(wavelength)*(lya_flux+flux_no_lya(wavelength))/(flux_no_lya(1500)),'r-')

        subax1 = add_subplot_axes(axis,subpos)
        subax1.plot(wavelength-wavelength_lya,dust*(lya_flux+flux_no_lya(wavelength))/(flux_no_lya(1500)),'k-')
        subax1.plot(wavelength-wavelength_lya,dust*TIGM_eff(wavelength)*(lya_flux+flux_no_lya(wavelength))/(flux_no_lya(1500)),'r-')

        xlabel('$\\lambda-\\lambda_{\\alpha}$ [$\\AA$]',fontsize=20)
        ylabel('$f_\\nu/f_\\nu(1500\\AA)$',fontsize=20)
        subax1.set_ylim(0,7)
        subax1.set_xlim(-2,10)


if __name__ == '__main__':
    plotinplot()
    plt.show()


figure(figsize=(15,7))
xlim(3500,5000)
ylim(0,7)      
ylabel('$f_\\nu/f_\\nu(1500\\AA)$',fontsize=22)
xlabel('Wavelength [$\\AA$]',fontsize=22)
plot(wavelength_optical,(lya_flux+flux_no_lya(wavelength_optical))/flux_no_lya(1500),'k-')
plot(wavelength_optical,TIGM_eff(wavelength_optical)*(lya_flux+flux_no_lya(wavelength_optical))/flux_no_lya(1500),'r-')
