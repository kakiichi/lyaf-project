# redshift-space 2PCF using the streaming model

# parameters
Mpc2km=3.08567758e19          # 1Mpc = 3.0857e19km
yr2s=3.15569e7                # 1yr = 3.15569e7 s

# cosmological parameters
ns=0.96
sig8=0.8
h=0.72
om_v0=0.70
om_m0=0.30
om_b0=0.044
Y=0.25

# redshift of interest
z=0.55

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

Nr=100
r_min=0.   # cMpc/h
r_max=10.0 # cMpc/h
r=linspace(r_min,r_max,Nr)

execfile('2PCF.py')
execfile('velocity_statistics.py')
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

mean_v12=interp1d(r,v12,kind='linear',bounds_error=False,fill_value=0.0)
sig12_v=interp1d(r,sig12,kind='linear',bounds_error=False,fill_value=0.0)
CF=interp1d(r,xi,kind='linear',bounds_error=False,fill_value=0.0)

# function
#def mean_v12(r):
    #v0=500.0 # km/s
    #r0=1.0   # cMpc/h
    #mean_v12=-v0/(1+(r/r0)**2)
    #return mean_v12

#def sig12_v(r):
    #sigma=10.0 # km/s
    #sig12_v=sigma**2
    #return sig12_v(r)

def velocity_PDF(v12,r):
    velocity_PDF=1./sqrt(2.*pi*sig12_v(r))*exp( -(v12-mean_v12(r))**2/(2.*sig12_v(r)) )
    return velocity_PDF

#def CF(r):
    #slope=1.8
    #r0=3.0   # cMpc/h
    #CF=(r/r0)**(-slope)
    #return CF

# main
Ns=100
s_para=linspace(0,10,Ns) # cMpc/h
s_perp=linspace(0,10,Ns)

Np=1000
v_min=-500.0
v_max=+500.0
v,dv=linspace(v_min,v_max,Np,retstep=True) # km/s

def Hubble(z):
    H0=1.02275e-10 # [h/yr]
    Hubble=H0*sqrt(om_m0*(1.+z)**3+om_v0)
    return Hubble

corrfunc=zeros([Ns,Ns])
for i in range(Ns):
    for j in range(Ns):

        intg_result=0.0
        for n in range(Np):
            r_para=sqrt( s_perp[j]**2 + (s_para[i]-v[n]/Mpc2km/(Hubble(z)/yr2s))**2 )
            r_perp=s_perp[j]
            r=sqrt(r_para**2+r_perp**2)
            intg_result=intg_result+(1.+CF(r))*velocity_PDF(v[n],r)*dv
        corrfunc[i,j]=intg_result-1.0

figure()
imshow(corrfunc[:,:],origin='lower',extent=[0,10,0,10])
colorbar()
contour(corrfunc[:,:],colors='black',linewidth=2.0,extent=[0,10,0,10])
