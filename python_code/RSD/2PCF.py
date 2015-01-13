def corrfunc_linear(r,z):
    from scipy.integrate import quadrature,quad,romberg
    kmin=1.0e-4                    # minimum k [h/Mpc]
    kmax=1.0e2                    # maximum k [h/Mpc]

    integrand=lambda lnk : exp(lnk)**3*powspec(exp(lnk),z)/(2.*pi**2)*sin(exp(lnk)*r)/(exp(lnk)*r)
    corrfunc_linear=quad(integrand,log(kmin),log(kmax),limit=5000)[0]

    return corrfunc_linear
