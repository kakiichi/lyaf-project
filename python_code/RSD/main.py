# cosmological parameters
ns=0.96
sig8=0.8
h=0.72
om_v0=0.70
om_m0=0.30
om_b0=0.044
Y=0.25

from scipy.interpolate import interp1d
execfile('cosmo_toolkits.py')
execfile('lin_powspec.py')

### load power spectrum
input_path='/home/koki/work/lyaf-survey-project/CLPT_GSRSD-master/data/PL_z0.55.dat'
k,P=genfromtxt(input_path,unpack=True )
P_L=interp1d(k,P,kind='linear',bounds_error=False,fill_value=0.0)
def powspec(k,z):
    powspec=P_L(k)*(growth_rate(z)/growth_rate(0.55))**2
    return powspec

execfile('2PCF.py')
execfile('velocity_statistics.py')

# main
# parameters
Mpc2km=3.08567758e19          # 1Mpc = 3.0857e19km
yr2s=3.15569e7                # 1yr = 3.15569e7 s

# redshift of interest
z=4.0
a=1./(1.+z)

# grid setup
Nr=50
r_min=0.1   # cMpc/h
r_max=100.0 # cMpc/h
r=logspace(log10(r_min),log10(r_max),Nr)

# compute 2PCF and velocity statistics <v12> & sigma12^2
vH=Hubble(z)/yr2s*(r*Mpc2km)*a
v12=zeros([Nr])
sig12=zeros([Nr])
xi=zeros([Nr])
for i in range(Nr):
    v12[i]=v12_linear(r[i],z,'quad')
    sig12[i]=sig12_linear(r[i],z,'quad')
    xi[i]=corrfunc_linear(r[i],z)

# plot two point correlation function
figure()
plot(r,r**2*xi,'r.-')

# plot pairwise velocity statistics
figure()
subplot(2,1,1)
semilogx(r,v12,'k.-')
semilogx(r,sqrt(sig12),'b-')

subplot(2,1,2)
xlim(0,20)
ylim(-100,2000)
plot(r,vH+v12,'r-')
plot(r,vH,'b--')
hlines(500,0,20,colors='k',linestyles='dashed')
text(13,400,'Region of influence')
