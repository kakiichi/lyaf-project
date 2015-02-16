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
def g_v12(v,form,r0,slope,req,betaN,v_in,r_in,slope_v,s0,v_out,r_out):
    def Hubble(z):
        Hubble=100.0*sqrt(om_v0+om_m0*(1.+z)**3) # (km/s)/(Mpc/h)
        return Hubble
    def xi(r):
        xi=(r/r0)**(-slope)
        return xi
    def PDF_r(r):
        factor=( (r/req)**(-2)+1. )**(-betaN+1)
        PDF_r=factor*(1.+xi(r))
        return PDF_r
    def v12(r):
        v12=v_out*2.0*exp(-r/r_out)-v_in/(1.0+(r/r_in)**slope_v)
        return v12
    def s12(r):
        s12=s0
        return s12
    def PDF_v(v,r):
        PDF_v=( 1./(sqrt(2.*pi)*s12(r)) )*exp( -(v-v12(r))**2/(2.*s12(r)**2) )
        return PDF_v
    def f_v12(v,form):
        if form == 'Gaussian':
            r_min=0.001  # comoving radius min. and max. limit of integral
            r_max=100.0 # cMpc/h
            integrand=lambda r: PDF_r(r)*PDF_v(v-Hubble(z)*r/(1.+z),r)
            result,err=scipy.integrate.quad(integrand,r_min,r_max)
            f_v12=result*Hubble(z)/(1.+z)
        if form == 'Hubble':
            f_v12=PDF_r( v*(1.+z)/Hubble(z) )
        return f_v12

    return f_v12(v,form)

# main
N=500
N_model=10
v=linspace(-800,2000,N)
f=zeros([N,N_model+1])

# fiducial model
r0=3.32     # cMpc/h, correlation length
slope=1.74  # power-law slope of correlation function
req=0.32     # cMpc/h, local-global UV background equality radius
betaN=2.0   # CDDF slope
v_in=143.0  # km/s
r_in=8.4    # cMpc/h
slope_v=2.26
s0=200.0    # km/s
v_out=0.0   # km/s
r_out=1.0   # cMpc/h

for i in range(N):
    # fiducial model
    f[i,1]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,s0,v_out,r_out)
    # local photoionization effect
    f[i,2]=g_v12(v[i],'Gaussian',r0,slope,0.10,betaN,v_in,r_in,slope_v,s0,v_out,r_out)
    f[i,3]=g_v12(v[i],'Gaussian',r0,slope,1.01,betaN,v_in,r_in,slope_v,s0,v_out,r_out)
    # inflow
    f[i,4]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,100.,r_in,slope_v,s0,v_out,r_out)
    f[i,5]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,200.,r_in,slope_v,s0,v_out,r_out)
    # velocity dispersion
    f[i,6]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,180.,v_out,r_out)
    f[i,7]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,220.,v_out,r_out)
    # outflow
    f[i,8]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,s0,200.,r_out)
    f[i,9]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,s0,400.,r_out)
    f[i,10]=g_v12(v[i],'Gaussian',r0,slope,req, betaN,v_in,r_in,slope_v,s0,600.,r_out)


# plot
figure(figsize=(10,10))
plt.rc('font',**{'family':'serif', 'size':13})
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=12)

subplot(2,2,1)
ylabel('$g(v_{12})$',fontsize=15)
xlabel('$v_{12}$ $[km/s]$',fontsize=15)
ylim(0.01,100)
xlim(-800,1800)
semilogy(v,f[:,1],'k-',  linewidth=2,label='$r_{eq}=0.32$')
semilogy(v,f[:,2],'k--', linewidth=2,label='$r_{eq}=0.10$')
semilogy(v,f[:,3],'k-.', linewidth=2,label='$r_{eq}=1.01$')
NHI=1.e19                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,1000,linestyles='dotted',colors='gray')   
NHI=1.e20                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,1000,linestyles='dotted',colors='gray')
NHI=1.e21                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,1000,linestyles='dotted')
text(120,700,'$v_{infl}(N_{HI}=10^{19}cm^2)$',rotation=90)
text(450,80,'$v_{infl}(N_{HI}=10^{20}cm^2)$',rotation=90)
lg=plt.legend(loc='upper right',fontsize=13)
#lg.draw_frame(False)

subplot(2,2,2)
ylabel('$g(v_{12})$',fontsize=15)
xlabel('$v_{12}$ $[km/s]$',fontsize=15)
ylim(0.01,100)
xlim(-800,1800)
semilogy(v,f[:,1],'r-',  linewidth=2,label='$v_{in}=143km/s$')
semilogy(v,f[:,4],'r--', linewidth=2,label='$v_{in}=100km/s$')
semilogy(v,f[:,5],'r-.', linewidth=2,label='$v_{in}=200km/s$')
NHI=1.e19                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')   
NHI=1.e20                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
NHI=1.e21                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
text(0.0,0.3,'$v_{infl}(N_{HI}=10^{19}cm^2)$',rotation=90)
text(300,0.3,'$v_{infl}(N_{HI}=10^{20}cm^2)$',rotation=90)
text(1400,0.3,'$v_{infl}(N_{HI}=10^{21}cm^2)$',rotation=90)
lg=plt.legend(loc='upper right',fontsize=13)

subplot(2,2,3)
ylabel('$g(v_{12})$',fontsize=15)
xlabel('$v_{12}$ $[km/s]$',fontsize=15)
ylim(0.01,100)
xlim(-800,1800)
semilogy(v,f[:,1],'b-',  linewidth=2,label='$\\sigma_{12}=200km/s$')
semilogy(v,f[:,6],'b--', linewidth=2,label='$\\sigma_{12}=180km/s$')
semilogy(v,f[:,7],'b-.', linewidth=2,label='$\\sigma_{12}=220km/s$')
NHI=1.e19                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')   
NHI=1.e20                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
NHI=1.e21                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
text(0.0,0.3,'$v_{infl}(N_{HI}=10^{19}cm^2)$',rotation=90)
text(300,0.3,'$v_{infl}(N_{HI}=10^{20}cm^2)$',rotation=90)
text(1400,0.3,'$v_{infl}(N_{HI}=10^{21}cm^2)$',rotation=90)
lg=plt.legend(loc='upper right',fontsize=13)

subplot(2,2,4)
ylabel('$g(v_{12})$',fontsize=15)
xlabel('$v_{12}$ $[km/s]$',fontsize=15)
ylim(0.01,100)
xlim(-800,1800)
semilogy(v,f[:,1],'g-',linewidth=2,label='$v_{out}=0km/s$')
semilogy(v,f[:,8],'g--', linewidth=2,label='$v_{out}=200km/s$')
semilogy(v,f[:,9],'g-.', linewidth=2,label='$v_{out}=400km/s$')
semilogy(v,f[:,10],'g:', linewidth=2,label='$v_{out}=600km/s$')
NHI=1.e19                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')   
NHI=1.e20                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
NHI=1.e21                    # cm2, HI column density of absorbers
v_crit=507.3*sqrt(NHI/1.e20) # km/s
vlines(v_crit,0.001,100,linestyles='dotted',colors='gray')
text(0.0,0.3,'$v_{infl}(N_{HI}=10^{19}cm^2)$',rotation=90)
text(300,0.3,'$v_{infl}(N_{HI}=10^{20}cm^2)$',rotation=90)
text(1400,0.3,'$v_{infl}(N_{HI}=10^{21}cm^2)$',rotation=90)


lg=plt.legend(loc='upper right',fontsize=13)

plt.tight_layout()

#semilogy(v,0.5*(1+scipy.special.erf( v/(sqrt(2.)*s0) )), 'k--')


