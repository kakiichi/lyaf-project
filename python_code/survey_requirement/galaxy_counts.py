
# number density of LAE down to L_lim luminosity [(cMpc/h)^-3]
def n_LAE(L_lim,z,fit):
    import mpmath
    if z==3.1:
        if fit=='1p5':
            phi_c=9.2e-4/h**3 # (h/cMpc)^3
            Lc=5.8e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=14.9e-4/h**3 # (h/cMpc)^3
            Lc=4.1e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=3.9e-4/h**3 # (h/cMpc)^3
            Lc=9.1e42         #  erg/s
            slope=-2.0
    if z==3.7:
        if fit=='1p5':
            phi_c=3.4e-4/h**3 # (h/cMpc)^3
            Lc=10.2e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=5.7e-4/h**3 # (h/cMpc)^3
            Lc=7.2e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=1.3e-4/h**3 # (h/cMpc)^3
            Lc=16.2e42         #  erg/s
            slope=-2.0
    if z==5.7:
        if fit=='1p5':
            phi_c=7.7e-4/h**3 # (h/cMpc)^3
            Lc=6.8e42         #  erg/s
            slope=-1.5
        if fit=='1p0':
            phi_c=10.4e-4/h**3 # (h/cMpc)^3
            Lc=5.4e42          #  erg/s
            slope=-1.0
        if fit=='2p0':
            phi_c=3.6e-4/h**3 # (h/cMpc)^3
            Lc=9.5e42         #  erg/s
            slope=-2.0


    n_LAE=phi_c*mpmath.fp.gammainc(slope+1.0, L_lim/Lc)
    return n_LAE

# limiting luminosity [erg/s] corresponding to flux limit [erg/s/cm2]
def L_limit(F_lim,z):
    Mpc2cm=3.08567758e24                # 1Mpc=3.08567758e24cm
    DL=luminosity_distance(z)/h*Mpc2cm  # Mpc/h to cm
    L_limit=4.*pi*DL**2*F_lim
    return L_limit

def Number_counts(F_lim,deg,z1,z2,z_fixed,fit):
    from scipy.integrate import quad
    FoV=(deg*(pi/180.))**2 # FoV in radian squared
    c=3.0659e-7 # Mpc/yr
    integrand=lambda z: c/Hubble(z)*angular_distance(z)**2*n_LAE(L_limit(F_lim,z),z_fixed,fit)
    Number_counts=FoV*quad(integrand,z1,z2,limit=100)[0]
    return Number_counts
