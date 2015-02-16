def dNdNHIdz(NHI,z):
    #if   ( (log10(NHI) >= 11.0) & (log10(NHI) < 15.0) ):
    if   ( (log10(NHI) >= 8.0) & (log10(NHI) < 15.0) ):
        A=1.2e7
        beta=1.5
        gamma=3.0
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 15.0) & (log10(NHI) < 17.5) ):
        A=3.8e14
        beta=2.0
        gamma=3.0
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 19.0) & (log10(NHI) < 20.3) ):
        A=0.45
        beta=1.05
        gamma=1.27
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    #elif ( (log10(NHI) >= 20.3) & (log10(NHI) < 21.55) ):
    elif ( (log10(NHI) >= 20.3) & (log10(NHI) < 23.0) ):
        A=8.7e18
        beta=2.0
        gamma=1.27
        dNdNHIdz=A*NHI**(-beta)*(1.+z)**gamma
    elif ( (log10(NHI) >= 17.5) & (log10(NHI) < 19.0) ):
        A1=3.8e14
        beta1=2.0
        gamma1=3.0
        C1=A1*(1.+z)**gamma1

        A2=0.45
        beta2=1.05
        gamma2=1.27
        C2=A2*(1.+z)**gamma2
            
        logy1=log10(C1)-17.5*beta1
        logy2=log10(C2)-19.0*beta2
        
        beta=(logy1-logy2)/(19.0-17.5)
        C=(10.**logy1)*(10.**17.5)**beta

        dNdNHIdz=C*NHI**(-beta)
    else:
        dNdNHIdz=0.0
        
    return dNdNHIdz

npoints=100
NHI=logspace(11,22.0,npoints)
f=zeros([npoints,4])
f1=zeros([npoints,4])
for i in range(npoints):
    f[i,0]=dNdNHIdz(NHI[i],2.0)
    f[i,1]=dNdNHIdz(NHI[i],3.0)
    f[i,2]=dNdNHIdz(NHI[i],4.0)
    f[i,3]=dNdNHIdz(NHI[i],5.0)
    f1[i,0]=0.01*dNdNHIdz(0.01*NHI[i],2.0)
    f1[i,1]=0.01*dNdNHIdz(0.01*NHI[i],3.0)
    f1[i,2]=0.01*dNdNHIdz(0.01*NHI[i],4.0)
    f1[i,3]=0.01*dNdNHIdz(0.01*NHI[i],5.0)


figure()
xlim(1e9,1e23)
loglog(NHI,f[:,0],'r-')
loglog(NHI,f[:,1],'r-')
loglog(NHI,f[:,2],'r-')
loglog(NHI,f[:,3],'r-')

loglog(NHI,f1[:,0],'b-')
loglog(NHI,f1[:,1],'b-')
loglog(NHI,f1[:,2],'b-')
loglog(NHI,f1[:,3],'b-')


