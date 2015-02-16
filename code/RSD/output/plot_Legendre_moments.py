from pylab import *

# read from files
nfiles=3
s_perp=range(nfiles)
s_para=range(nfiles)
xi_s=range(nfiles)

s_perp[0],s_para[0],xi_s[0]=genfromtxt('RSD_01.out',unpack=True,usecols=(2,3,4))
#s_perp[1],s_para[1],xi_s[1]=genfromtxt('RSD_02.out',unpack=True,usecols=(2,3,4))
#s_perp[2],s_para[2],xi_s[2]=genfromtxt('RSD_03.out',unpack=True,usecols=(2,3,4))
#s_perp[1],s_para[1],xi_s[1]=genfromtxt('RSD_04.out',unpack=True,usecols=(2,3,4))
#s_perp[2],s_para[2],xi_s[2]=genfromtxt('RSD_05.out',unpack=True,usecols=(2,3,4))
#s_perp[1],s_para[1],xi_s[1]=genfromtxt('RSD_06.out',unpack=True,usecols=(2,3,4))
#s_perp[2],s_para[2],xi_s[2]=genfromtxt('RSD_07.out',unpack=True,usecols=(2,3,4))
#s_perp[3],s_para[3],xi_s[3]=genfromtxt('RSD_08.out',unpack=True,usecols=(2,3,4))
s_perp[1],s_para[1],xi_s[1]=genfromtxt('RSD_09.out',unpack=True,usecols=(2,3,4))
s_perp[2],s_para[2],xi_s[2]=genfromtxt('RSD_10.out',unpack=True,usecols=(2,3,4))


title=range(nfiles)
#title[0]='$r_{eq}=0.32h^{-1}$cMpc (fid.)'
#title[1]='$r_{eq}=0.10h^{-1}$cMpc'
#title[2]='$r_{eq}=1.01h^{-1}$cMpc'
#title[0]='$v_{in}=143$km/s (fid.)'
#title[1]='$v_{in}=100$km/s'
#title[2]='$v_{in}=200$km/s'
#title[0]='$v_{out}=0$km/s (fid.)'
#title[1]='$v_{out}=200$km/s'
#title[2]='$v_{out}=400$km/s'
#title[3]='$v_{out}=600$km/s'
title[0]='$\\sigma_{12}=200$km/s (fid.)'
title[1]='$\\sigma_{12}=180$km/s'
title[2]='$\\sigma_{12}=220$km/s'


# reshape to Fortran like order [si,sj] (perp,para)

s_perp_map=range(nfiles)
s_para_map=range(nfiles)
xi_s_map=range(nfiles)
s_sample=range(nfiles)

Ns=int(sqrt(xi_s[0].size))
for ID in range(nfiles):
    s_perp_map[ID]=reshape(s_perp[ID],(Ns,Ns),order='F')
    s_para_map[ID]=reshape(s_para[ID],(Ns,Ns),order='F')
    xi_s_map[ID]=reshape(xi_s[ID],(Ns,Ns),order='F') 
    s_sample[ID]=s_perp_map[ID][0,:]

# Legendre moment of redshift space correlation function
def xi_moment(s,l,ID):   

    def find_nearest(array,value):
        idx = (abs(array-value)).argmin()
        return array[idx],idx

    def Legendre(x,l):
        if l == 0:
            Legendre=1.
        elif l == 2:
            Legendre=1./2.*(3.*x*x-1.)
        elif l == 4:
            Legendre=1./8.*(35.*x**4-30.*x**2+3.)
        else:
            print 'no Legendre moment implemented yet'
        return Legendre

    Ng=200
    mu=zeros(Ng)
    dmu=1./Ng
    for i in range(Ng):
        mu[i]=(i+0.5)*dmu

    xi_moment=0.
    for i in range(Ng):
        s_perp=s*sqrt(1.-mu[i]**2)
        s_para=s*mu[i]    
        s1_value,s1_idx=find_nearest(s_sample[ID],s_perp)
        s2_value,s2_idx=find_nearest(s_sample[ID],s_para)
        idx=where((s_perp_map[ID]==s1_value)&(s_para_map[ID]==s2_value))
        xi_moment=xi_moment+xi_s_map[ID][idx]*Legendre(mu[i],l)*dmu
    xi_moment=(2.*l+1.)/2. * (xi_moment * 2.) # x2 because we only integrated 
                                              # firstquandrant mu=(0,1) 

    return xi_moment



# main
Nbins=50
s_min=0.1 # cMpc/h
s_max=7.0 # cMpc/h
s=linspace(s_min,s_max,Nbins)

xi_0=range(nfiles)
xi_2=range(nfiles)
#xi_4=range(nfiles)
for ID in range(nfiles):
    print 'computing the Ledengre moments for model ID= ',ID
    xi_0[ID]=zeros(Nbins)
    xi_2[ID]=zeros(Nbins)
    #xi_4[ID]=zeros(Nbins)
    for i in range(Nbins):
        xi_0[ID][i]=xi_moment(s[i],0,ID)
        xi_2[ID][i]=xi_moment(s[i],2,ID)
        #xi_4[ID][i]=xi_moment(s[i],4,ID)

# plot
figure(figsize=(15,6))
plt.rc('font',**{'family':'serif', 'size':22})
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=23)

subplot(1,3,1)
ylabel('$\\xi_0(s)$',fontsize=22)
xlabel('$s$ [$h^{-1}$cMpc]',fontsize=22)
semilogy(s,xi_0[0],'r-', linewidth=2,label=title[0])
semilogy(s,xi_0[1],'r--',linewidth=2,label=title[1])
semilogy(s,xi_0[2],'r-.',linewidth=2,label=title[2])
#semilogy(s,xi_0[3],'r:', linewidth=2,label=title[3])
lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)

subplot(1,3,2)
ylabel('$s^2\\xi_2(s)$',fontsize=22)
xlabel('$s$ [$h^{-1}$cMpc]',fontsize=22)
plot(s,s**2*xi_2[0],'b-', linewidth=2,label=title[0])
plot(s,s**2*xi_2[1],'b--',linewidth=2,label=title[1])
plot(s,s**2*xi_2[2],'b-.',linewidth=2,label=title[2])
#plot(s,s**2*xi_2[3],'b:', linewidth=2,label=title[3])
lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)

subplot(1,3,3)
ylabel('$\\xi_2/\\xi_0$',fontsize=22)
xlabel('$s$ [$h^{-1}$cMpc]',fontsize=22)
plot(s,xi_2[0]/xi_0[0],'g-', linewidth=2,label=title[0])
plot(s,xi_2[1]/xi_0[1],'g--',linewidth=2,label=title[1])
plot(s,xi_2[2]/xi_0[2],'g-.',linewidth=2,label=title[2])
#plot(s,xi_2[3]/xi_0[3],'g:', linewidth=2,label=title[3])
lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)

plt.tight_layout()

#figure()
#imshow(xi_s_map)
