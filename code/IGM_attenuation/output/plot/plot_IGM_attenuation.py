from scipy.interpolate import interp1d

# read the oscillator strength and Lyn resonant frequancy
c=2.998e8 * 1.e10 # A/s
freq_lyn=genfromtxt('../../input/Lyman_parameters.input',usecols=(1),unpack=True)
w_lyn=c/freq_lyn[0:49]

ID='01'
lynmax=31
tau_lyn_LAF_data=range(lynmax)
tau_lyn_LLS_data=range(lynmax)
tau_lyn_DLA_data=range(lynmax)
for n in range(lynmax):
    # Lya forest absorber 11<log10(NHI/[1/cm2])<17
    infile='../IGM_attenuation_'+ID+'_LAF_ly'+str(n)+'.output'
    tau_lyn_LAF_data[n]=genfromtxt(infile)
    # LLS absorber 17<log10(NHI/[1/cm2])<20.3
    infile='../IGM_attenuation_'+ID+'_LLS_ly'+str(n)+'.output'
    tau_lyn_LLS_data[n]=genfromtxt(infile)
    # DLA absorber 20.3<log10(NHI/[1/cm2])<22.0
    infile='../IGM_attenuation_'+ID+'_DLA_ly'+str(n)+'.output'
    tau_lyn_DLA_data[n]=genfromtxt(infile)


tau_lyn_LAF=range(lynmax)
tau_lyn_LLS=range(lynmax)
tau_lyn_DLA=range(lynmax)
for n in range(lynmax):
    tau_lyn_LAF[n]=interp1d(tau_lyn_LAF_data[n][:,0],tau_lyn_LAF_data[n][:,3],bounds_error=False,fill_value=0.0)
    tau_lyn_LLS[n]=interp1d(tau_lyn_LLS_data[n][:,0],tau_lyn_LLS_data[n][:,3],bounds_error=False,fill_value=0.0)
    tau_lyn_DLA[n]=interp1d(tau_lyn_DLA_data[n][:,0],tau_lyn_DLA_data[n][:,3],bounds_error=False,fill_value=0.0)

w0=range(lynmax)
for n in range(lynmax):
    idx=where((tau_lyn_LAF_data[n][:,0]<=w_lyn[n]+6)&(tau_lyn_LAF_data[n][:,0]>=w_lyn[n]-6))
    w0[n]=tau_lyn_LAF_data[n][idx,0]

def tau_total(wavelength,lynmax,type):
    if type=='LAF':
        tau_total=0.
        for n in range(lynmax):
            tau_total=tau_total+tau_lyn_LAF[n](wavelength)
    if type=='LLS':
        tau_total=0.
        for n in range(lynmax):
            tau_total=tau_total+tau_lyn_LLS[n](wavelength)
    if type=='DLA':
        tau_total=0.
        for n in range(lynmax):
            tau_total=tau_total+tau_lyn_DLA[n](wavelength)
    return tau_total
Ng=100000
w=linspace(912,1230,Ng)

fig = plt.figure(figsize=(13,4))                      
ax = fig.add_subplot(1,1,1) 
#figure()
plt.rc('font',**{'family':'serif','size':15})
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('legend',fontsize=15)
major_ticks = numpy.arange(910, 1240, 30)                                              
minor_ticks = numpy.arange(910, 1240, 1)                                               
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                      
xlim(909,1230)
ylim(0,1)

xlabel('Rest frame wavelength $\\lambda$ $[\\AA]$',fontsize=18)
ylabel('$e^{-\\tau_{eff}}$',fontsize=18)

plot(w,exp(-tau_total(w,lynmax,'LAF')-tau_total(w,lynmax,'LLS')-tau_total(w,lynmax,'DLA')),'r-',linewidth=1.5)
#plot(w,exp(-tau_total(w,lynmax,'LLS')),'g-')
#plot(w,exp(-40*tau_total(w,lynmax,'DLA')),'r-')

vlines(912,0,1,colors='k',linestyles='solid')
#text(905,0.05,'LL')

for n in range(lynmax):
    vlines(w_lyn[n],0,1,colors='k',linestyles='dotted')
text(w_lyn[0],0.05,'Ly$\\alpha$')
text(w_lyn[1],0.05,'Ly$\\beta$')
text(w_lyn[2],0.05,'Ly$\\gamma$')
text(w_lyn[3],0.05,'Ly$\\delta$')
text(w_lyn[4],0.05,'Ly$\\epsilon$')

plt.tight_layout()

# plot-in-plot ####################################
# plot IGM attenuation

plt.rc('font',**{'family':'serif', 'size':11})
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
plt.rc('legend', fontsize=12)

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def plotinplot():
    fig = plt.figure(figsize=(14,4))
    axes = []
    subpos = [0.5,0.3,0.3,0.5]
    x = np.linspace(-np.pi,np.pi)
    for i in range(1):
        axes.append(fig.add_subplot(1,1,i))
    for axis in axes:
        major_ticks = numpy.arange(910, 1240, 30)                                              
        minor_ticks = numpy.arange(910, 1240, 1)                                               
        axis.set_xticks(major_ticks)                                                       
        axis.set_xticks(minor_ticks, minor=True)                      

        axis.set_xlim(909,1230)      
        axis.set_ylim(0,1)      
        xlabel('Rest frame wavelength $\\lambda$ $[\\AA]$',fontsize=18)
        ylabel('$e^{-\\tau_{eff}}$',fontsize=18)
        axis.plot(w,exp(-tau_total(w,lynmax,'LAF')-tau_total(w,lynmax,'LLS')-tau_total(w,lynmax,'DLA')),'r-',linewidth=1.5)
        text(1180,0.9,'$z_s=3$',fontsize=18)
        #plot(w,exp(-tau_total(w,lynmax,'LLS')),'g-')
        #plot(w,exp(-40*tau_total(w,lynmax,'DLA')),'r-')
        vlines(912,0,1,colors='k',linestyles='solid')
        for n in range(lynmax):
            vlines(w_lyn[n],0,1,colors='k',linestyles='dotted')
        text(w_lyn[0],0.05,'Ly$\\alpha$')
        text(w_lyn[1],0.05,'Ly$\\beta$')
        text(w_lyn[2],0.05,'Ly$\\gamma$')
        text(w_lyn[3],0.05,'Ly$\\delta$')
        text(w_lyn[4],0.05,'Ly$\\epsilon$')


        subax1 = add_subplot_axes(axis,subpos)
        subax1.plot(w0[0][0]-w_lyn[0],exp(-tau_total(w0[0][0],lynmax,'LAF')-tau_total(w0[0][0],lynmax,'LLS')-tau_total(w0[0][0],lynmax,'DLA')),'r-')
        subax1.plot(w0[1][0]-w_lyn[1],exp(-tau_total(w0[1][0],lynmax,'LAF')-tau_total(w0[1][0],lynmax,'LLS')-tau_total(w0[1][0],lynmax,'DLA')),'y-')
        subax1.plot(w0[2][0]-w_lyn[2],exp(-tau_total(w0[2][0],lynmax,'LAF')-tau_total(w0[2][0],lynmax,'LLS')-tau_total(w0[2][0],lynmax,'DLA')),'g-')
        subax1.plot(w0[3][0]-w_lyn[3],exp(-tau_total(w0[3][0],lynmax,'LAF')-tau_total(w0[3][0],lynmax,'LLS')-tau_total(w0[3][0],lynmax,'DLA')),'c-')
        subax1.plot(w0[4][0]-w_lyn[4],exp(-tau_total(w0[4][0],lynmax,'LAF')-tau_total(w0[4][0],lynmax,'LLS')-tau_total(w0[4][0],lynmax,'DLA')),'b-')
        subax1.vlines(0,0,1,colors='k',linestyles='dotted')
        xlabel('$\\lambda-\\lambda_{Lyn}$ $[\\AA]$',fontsize=13)
        ylabel('$e^{-\\tau_{eff}}$',fontsize=13)
        subax1.set_ylim(0,1)
        subax1.set_xlim(-5.5,5.5)


if __name__ == '__main__':
    plotinplot()
    plt.show()

plt.tight_layout()
