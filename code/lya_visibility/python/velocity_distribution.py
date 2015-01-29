def g_v12(v,form,xi_on,r0,slope,phr_on,req,betaN,v0,s0,outflow_on,v_out,r_out):
###  r0=3.32    #cMpc/h correlation length
###  slope=1.74 # slope of correlation function
###  
###  req=1.0
###  betaN=2.
###
###  v0=200.
###  r0=3.32
###  slope=1.74
###  
###  s0=200.
###
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
