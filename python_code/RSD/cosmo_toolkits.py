# cosmological parameter for distance measure
#h=0.70
#om_m0=0.3
#om_v0=0.7

def Hubble(z):
    H0=1.02275e-10 # [h/yr]
    Hubble=H0*sqrt(om_m0*(1.+z)**3+om_v0)
    return Hubble

# comoving interval distance between redshift z2 to z1 [cMpc/h]
def comoving_interval(z2,z1): 
    import scipy.integrate
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)
    comoving_interval=scipy.integrate.quad(integrand,z1,z2)[0]
    return comoving_interval

# comoving distance to redshift z [cMpc/h]
def comoving_distance(z):  
    import scipy.integrate
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)
    comoving_distance=scipy.integrate.quad(integrand,0.0,z)[0]
    return comoving_distance

# comoving size of object with angular size in degree at redshift z [cMpc/h]
def comoving_size(degree,z):
    radian=(pi/180.)*degree
    comoving_size=comoving_distance(z)*radian
    return comoving_size

# time interval between two redshift [yr]
def time_redshift_interval(z2,z1):
    import scipy.integrate
    integrand=lambda z: 1./((1.+z)*Hubble(z)*h)
    time_redshift_interval=scipy.integrate.quad(integrand,z1,z2)[0]
    return time_redshift_interval

# angular diameter distance to redshift z [cMpc/h]
def angular_distance(z):
    angular_distance=comoving_distance(z)/(1.+z)
    return angular_distance

# luminosity distance to redshift z [cMpc/h]
def luminosity_distance(z):
    luminosity_distance=(1.+z)**2*angular_distance(z)
    return luminosity_distance
