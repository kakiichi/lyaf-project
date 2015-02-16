# python script to calculate and plot the velocity factor
import scipy.integrate
import scipy.special

# simulation set up
z=3.0

# parameters
h=0.72
om_v0=0.73
om_m0=0.27
om_b0=0.044

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
def f_v12(v,r):
    f_v12=PDF_r(r)*PDF_v(v-Hubble(z)*r/(1.+z),r)
    #f_v12=PDF_r(r)*PDF_v(v,r)
    return f_v12

# main
N1=100
N2=200
v=linspace(-1000,1000,N1)
r=linspace(0.001,10,N2)
f=zeros([N1,N2])


# model class  (Gaussian velocity PDF + inflow 200km/s + local ionization model)
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
for i in range(N1):
    for j in range(N2):
        f[i,j]=f_v12(v[i],r[j])


# plot
figure()
plt.rc('font',**{'family':'serif', 'size':13})
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=12)

imshow(log10(f),vmin=-5,vmax=-2,origin='lower',extent=(0.01,10,-1000,1000),aspect=0.005)
colorbar()
contour(f,extent=(0,10,-1000,1000))
plot(r,Hubble(z)*r/(1+z),'k--')
ylim(-1000,1000)

