h=0.72
om_m0=0.26
om_v0=0.74

H0=100.0*h # km/s/Mpc

tau=1.0      # optical depth
NHI=10**20.0 # cm^2

def D_infl(z,tau,NHI):
    v_infl=507.3/sqrt(tau)*sqrt(NHI/10**20.0) # km/s, region of infl. in velocity
    D_infl=v_infl/H0*(1.+z)/sqrt(om_m0*(1.+z)**3+om_v0) # cMpc/h
    return D_infl

execfile('cosmo_toolkits.py')
N=100
z=linspace(2.0,8.0,N)
DA=zeros(N)
for i in range(N):
    DA[i]=angular_distance(z[i]) # cMpc/h


fig= plt.figure(figsize=(6.5,3.3)) 
plt.rc('font',**{'family':'serif', 'size':15})
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

subplot(1,2,1)
lw=0.5
plot(z,D_infl(z,1.0,10**19.0),'r-',linewidth=lw+0.0)
plot(z,D_infl(z,1.0,10**20.0),'r-',linewidth=lw+0.5)
plot(z,D_infl(z,1.0,10**21.0),'r-',linewidth=lw+1.0)
plot(z,D_infl(z,1.0,10**22.0),'r-',linewidth=lw+1.5)
minorticks_on()
ylabel(' $D_{infl}$ [$h^{-1}$cMpc]')
xlabel('z')

subplot(1,2,2)
rad2arcmin=60.*180./pi # radian to arcmin
semilogy(z,(D_infl(z,1.0,10**19.0)/DA)*rad2arcmin,'r-',linewidth=lw+0.0)
semilogy(z,(D_infl(z,1.0,10**20.0)/DA)*rad2arcmin,'r-',linewidth=lw+0.5)
semilogy(z,(D_infl(z,1.0,10**21.0)/DA)*rad2arcmin,'r-',linewidth=lw+1.0)
semilogy(z,(D_infl(z,1.0,10**22.0)/DA)*rad2arcmin,'r-',linewidth=lw+1.5)
minorticks_on()
hlines(2.0,7,8,linestyles='dashed')
hlines(30.0,7,8,linestyles='solid')
xlim(2,8)
ylabel(' $\\theta_{infl}$ [arcmin]')
xlabel('z')

plt.tight_layout()
