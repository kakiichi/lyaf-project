# functions to compute the linear power spectrum
from scipy.integrate import quad

# E(z) function H(z)=H0*E(z)
def E(z):
    E=sqrt(om_m0*(1.+z)**3+om_v0)
    return E

def Omega_v(z):
    Omega_v=om_v0/(E(z)**2)
    return Omega_v

def Omega_m(z):
    Omega_m=om_m0*(1.+z)**3/(E(z)**2)
    return Omega_m

# Window function
def W2(x,filter): # x=kR
    if filter == 'gauss':
        W=exp(-0.5*x*x)
    elif filter == 'sharp-k':
        if (1.0-x) > 0.0 :
            W=1.0
        elif (1.0-x) == 0.0 :
            W=0.5
        else:
            W=0.0
    elif filter == 'tophat':
        W=3./(x**3)*(sin(x)-x*cos(x))
    elif filter == 'tophat+cutoff':
        W_TH=3./(x**3)*(sin(x)-x*cos(x))
        W_EXP=(1.0+(0.4*x)**2)**(-2)
        W=W_TH*W_EXP
    elif filter == 'gauss2':
        W=exp(-0.5*x*x)**2
    elif filter == 'powerlaw':
        W=1./(1.+x*x)
    else:
        print 'Window function is not speficied.'
    W2=W*W
    return W2

# Carroll, Press, Turner 1992 fit to growth factor g=D/a
def growth_factor(z):
    om_mz=Omega_m(z)
    om_vz=Omega_v(z)
    growth_factor=(5.*om_mz/2.)/(om_mz**(4./7.)-om_vz+(1.+om_mz/2.)*(1.+om_vz/70.))
    return growth_factor

# growth rate D(z)
def growth_rate(z):
    growth_rate=growth_factor(z)/(1.+z)
    return growth_rate

# BBKS transfer function: k in [h/cMpc]
def tranfer_func(k):
    q=k/(om_m0*h)
    tranfer_func=log(1.+2.34*q)/(2.34*q)*(1.+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-1./4.)
    return tranfer_func

# dummy linear power spectrum at z=0 with Ap=1 unit normalisation
def P_dummy(k,z):
    P_dummy=(growth_rate(z)**2)*(tranfer_func(k)**2)*(k**ns)
    return P_dummy

# Normalisation to linear power spectrum
integrand=lambda k: ((k**3)*P_dummy(k,0.)/(2.*pi**2))*W2(k*8.,'tophat')/k
normfactor=quad(integrand,0.0,inf,limit=500)[0]
Ap=sig8/normfactor

# linear power spectrum
def P_linear(k,z):
    P_linear=Ap*P_dummy(k,z)
    return P_linear

# linear dimensionless power spectrum
def Delta2(k,z):
    Delta2=k**3*P_linear(k,z)/(2.*pi**2)
    return Delta2
