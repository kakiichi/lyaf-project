# plot correlation function

N=100
r=logspace(-2,1.5,N)

def xi(r,r0,slope):
    xi=(r/r0)**(-slope)
    return xi

def S(r,req,betaN):
    S=( (r/req)**(-2)+1. )**(-betaN+1)
    return S

r_data,xi_data=genfromtxt('xi_linear_z3p0.output',unpack=True)

#def xi_reion(r,Rb):
#    xi_reion=-1/(1+(r/Rb)**2)+0.5
#    return xi_reion

figure()
plt.rc('font',**{'family':'serif', 'size':22})
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=23)
ylabel('$\\xi(r)$',fontsize=22)
xlabel('$r$ [$h^{-1}$cMpc]',fontsize=22)

xlim(0.01,20)
ylim(0.01,1e5)
b_LBG=1.0
b_DLA=1.0
loglog(r_data,b_LBG*b_DLA*xi_data,'k:',label='linear theory')

r0=3.32
slope=1.74
loglog(r,xi(r,r0,slope),'k',label='$r_0=3.32$, $\gamma=1.74$')

# variation over equality radius
req=0.320
betaN=2
loglog(r,S(r,req,betaN)*(1+xi(r,r0,slope))-1,'r-',linewidth=3,label='$r_{eq}=0.32$, $\\beta_N=2.0$')
req=1.012
betaN=2
loglog(r,S(r,req,betaN)*(1+xi(r,r0,slope))-1,'r--',linewidth=3,label='$r_{eq}=1.01$, $\\beta_N=2.0$')
req=0.101
betaN=2
loglog(r,S(r,req,betaN)*(1+xi(r,r0,slope))-1,'r:',linewidth=3,label='$r_{eq}=0.10$, $\\beta_N=2.0$')

# variation over betaN
req=0.320
betaN=1.3
loglog(r,S(r,req,betaN)*(1+xi(r,r0,slope))-1,'b-',linewidth=3,label='$r_{eq}=0.32$, $\\beta_N=1.3$')
req=0.320
betaN=1.7
loglog(r,S(r,req,betaN)*(1+xi(r,r0,slope))-1,'b--',linewidth=3,label='$r_{eq}=0.32$, $\\beta_N=1.7$')

lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)


plt.tight_layout()


#loglog(r,0.5*xi(r,r0,slope)+0.5*xi_reion(r,10.0)+1)


