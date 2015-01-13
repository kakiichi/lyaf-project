#from pylab import *

### pairwise velocity statistics <v12> and sigma12^2
# spherical Bessel function
def j1(x):
    j1=(sin(x)-x*cos(x))/x/x
    return j1

# mean pairwise velocity field
def v12_linear(r,z,method):
    from scipy.integrate import quad,romberg
    Mpc2km=3.08567758e19          # 1Mpc = 3.0857e19km
    yr2s=3.15569e7                # 1yr = 3.15569e7 s
    kmin=1.0e-4                   # minimum k [h/Mpc]
    kmax=1.0e2                    # maximum k [h/Mpc]
    integrand=lambda lnk : powspec(exp(lnk),z)*j1(exp(lnk)*r)*exp(lnk)**2
    if method == 'quad':
        intg_result=quad(integrand,log(kmin),log(kmax),limit=5000)[0]
    if method == 'romberg':
        intg_result=romberg(integrand,log(kmin),log(kmax),divmax=10)

    dlnDdlna=om_m0**0.6
    Hz=Hubble(z)/yr2s
    v12_linear=-Hz*dlnDdlna/(pi**2)*intg_result * Mpc2km # km/s

    return v12_linear

# pairwise velocity dispersion (parallel component)
def sig12_linear(r,z,method):
    from scipy.integrate import quad,romberg

    Mpc2km=3.08567758e19          # 1Mpc = 3.0857e19km
    yr2s=3.15569e7                # 1yr = 3.15569e7 s
    kmin=1.0e-4                   # minimum k [h/Mpc]
    kmax=1.0e2                    # maximum k [h/Mpc]

    dlnDdlna=om_m0**0.6
    Hz=Hubble(z)/yr2s    

    if method == 'quad':
        integrand=lambda lnk: powspec(exp(lnk),z)*exp(lnk)
        intg_result=quad(integrand,log(kmin),log(kmax),limit=5000)[0]
        sig2_v=(Hz*dlnDdlna)**2/(6.*pi**2)*intg_result * Mpc2km**2

        integrand=lambda lnk: powspec(exp(lnk),z)*j1(exp(lnk)*r)/(exp(lnk)*r)*exp(lnk)
        intg_result=quad(integrand,log(kmin),log(kmax),limit=5000)[0]
        psi=(Hz*dlnDdlna)**2/(2.*pi**2)*intg_result * Mpc2km**2

    sig12_linear=2.*(sig2_v-psi)

    return sig12_linear
