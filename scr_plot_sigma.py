# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:06:26 2017

@author: dori
"""

#from scattDB import shape
#from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from scattDB import scattering
from scattDB import psd
from scipy import integrate

plt.close('all')

c = 299792458000. # mm/s
lamx = c/(9.6*1e9)
lamu = c/(13.4*1e9)
lama = c/(35.6*1e9)
lamw = c/(94*1e9)
coeffx = lamx**4./(0.93*np.pi**5.)
coeffu = lamu**4./(0.95*np.pi**5.)
coeffa = lama**4./(0.95*np.pi**5.)
coeffw = lamw**4./(0.75*np.pi**5.)

tablesfolder = '/work/DBs/scattDB/tables/'
authors = ['DO','BJ','JL']
author = 'DO'
melt = '00'

dataDO = pd.read_csv(tablesfolder+author+'_'+melt+'.csv')
dataDO.drop('Unnamed: 0',axis=1,inplace=True)

author = 'BJ'
dataBJ2 = pd.read_csv(tablesfolder+author+'_agg2_0'+'.csv')
dataBJ2.drop('Unnamed: 0',axis=1,inplace=True)

dataBJ3 = pd.read_csv(tablesfolder+author+'_agg3_0'+'.csv')
dataBJ3.drop('Unnamed: 0',axis=1,inplace=True)

author = 'dataJL_A0.0'
dataJL = pd.read_csv(tablesfolder+author+'.csv')

#plt.figure()
#ax = plt.gca()
#ax.scatter(dataBJ2.Dmax,dataBJ2.Ku,label='BJ2')
#ax.scatter(dataBJ3.Dmax,dataBJ3.Ku,label='BJ3')
#ax.scatter(dataJL.Dmax, 1e6*dataJL.Ub, label='JL')
#ax.scatter(dataDO.Dmax, dataDO.Ku, label='DO')
#ax.legend()
#ax.grid()

def Z(D,s,psd):
    conc = psd(D)
    intZ = np.multiply(conc,s)
    return integrate.trapz(intZ,D)
    
def f3plot(data,title='title'):
    D0s = np.linspace(0.1,20.,20)
    mus = [-1.,0.,1.,2.,3.]
    marks=[',','+','h','v','.']
    plt.figure()
    ax = plt.gca()
    for mu,m in zip(mus,marks):
        XKa = []
        KaW = []
        for D0 in D0s:
            conc  = psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=22)
            Zx = 10.*np.log10(coeffx*Z(data.Dmax,data.X,conc))
            Zu = 10.*np.log10(coeffu*Z(data.Dmax,data.Ku,conc))
            if np.isnan(Zx):
                Zx = Zu
            Za = 10.*np.log10(coeffa*Z(data.Dmax,data.Ka,conc))
            Zw = 10.*np.log10(coeffw*Z(data.Dmax,data.W,conc))
            #LDR= 10.*np.log10(Z(data.Dmax,data.ldr,conc))
            XKa.append(Zx-Za)
            KaW.append(Za-Zw)
        s = ax.scatter(KaW,XKa,c=D0s,marker=m,label='$\mu=$ '+str(mu))
    ax.legend()
    ax.grid()
    ax.set_xlabel('DWR$_{Ka,W}$')
    ax.set_ylabel('DWR$_{X,Ka}$')
    ax.set_xlim([-1,15])
    ax.set_ylim([0,22])
    colorbar = plt.colorbar(mappable=s,ax=ax)
    colorbar.set_label('$D_0$')
    ax.set_title(title)

f3plot(dataDO,'Davide dry')
f3plot(dataBJ2,'BJ2 dry')
f3plot(dataBJ3,'BJ3 dry')
dataJL.columns = ['Dmax','model','ELWP','mkg','Dmax.1','Rgyr','ar','riming', 'Xa','Xs', 'X', 'Xe', 'Ua', 'Us', 'Ku', 'Ue', 'Aa', 'As', 'Ka', 'Ae', 'Wa','Ws', 'W', 'We']
f3plot(dataJL,'Jussi unrimed')

#%%


cols = ['Ku_Ka','Ka_W','LDRka','melt','Dm']

plt.figure()
ax = plt.gca()
plt.figure()
ax2 = plt.gca()
datatot = pd.DataFrame(columns=cols)
for melt_frac in melt_fracs:
    data = pd.DataFrame(columns=cols)
    i=0
    dists = {}
    for freq in freqs:

        
        scatt_folders = sorted(glob(scattfolder+'*aggregate2*'+'_f'+melt_frac+'_*'+freq))
        #print(scatt_folders)
        dist = scattering.ScattDist()
        
        for fld in scatt_folders:
            shpstr, aeffstr = fld.split('/')[-1].split('AEFF')
            shape_folders = glob(shapefolder+shpstr+'*')
            shapefile = shape_folders[0] + '/shape.dat'
            shp = shape.shapeDDA(shapefile)
            
            avgfiles = sorted(glob(fld+'/*.avg'))
            if 'to' in aeffstr:
                pass
#                for avg in avgfiles:
#                    scatt = scattering.ScattDDSCAT(avg) 
#                    scatt.D = shp.find_dmax(aeff=scatt.aeff)
#                    scatt.mass = shp.find_mass(aeff=scatt.aeff)
#                    scatt.melt = shp.get_melted_fraction()
#                    dist.add_scatterer(scatt)
                    
            else:
#                pass
                scatt = scattering.ScattDDSCAT(avgfiles[0]) 
                scatt.D = shp.find_dmax(aeff=scatt.aeff)
                scatt.mass = shp.find_mass(aeff=scatt.aeff)
                scatt.melt = shp.get_melted_fraction()
                #scatt.melt = float(shpstr[-7:-1])/10000.
                dist.add_scatterer(scatt)
        dists[freq] = dist        
        #array = dist.get_distro(['sig_bk','mass','melt','ldr'])

    sizes = np.linspace(0.25,30.7,1024)
    #plt.figure()
    #ax2 = plt.subplot(3,1,1)
    #ax3 = plt.subplot(3,1,2)
    #ax4 = plt.subplot(3,1,3)
    for lam in lambdas:
#        conc  = psd.ExponentialPSD(Lambda=lam,D_max=max(sizes)+0.1)
        conc  = psd.GammaPSD(D0=(3.67+mu)/lam,Nw=1.,mu=mu,D_max=max(sizes)+0.1)
        concD = conc(sizes) 
        propsu = dists['13.4'](sizes,['radar_xsect','ldr'])
        propsa = dists['35.6'](sizes,['radar_xsect','ldr'])
        propsW = dists[  '94'](sizes,['radar_xsect','ldr'])
        intZu  = np.multiply(concD,propsu[:,0])
        intLu  = np.multiply(concD,propsu[:,1])
        intZa  = np.multiply(concD,propsa[:,0])
        intLa  = np.multiply(concD,propsa[:,1])
        intZW  = np.multiply(concD,propsW[:,0])
        intLW  = np.multiply(concD,propsW[:,1])
        Zu = 10*np.log10(coeffu*integrate.trapz(intZu,sizes))
        LDRu = 10*np.log10(integrate.trapz(intLu,sizes))
        Za = 10*np.log10(coeffa*integrate.trapz(intZa,sizes))
        LDRa = 10*np.log10(integrate.trapz(intLa,sizes))
        ZW = 10*np.log10(coeffW*integrate.trapz(intZW,sizes))
        LDRW = 10*np.log10(integrate.trapz(intLW,sizes))
        datapart = pd.DataFrame(index=[i],columns=cols)
        datapart.loc[i] = [Zu-Za,Za-ZW,LDRa,melt_frac,1.0/lam]
        data = data.append(datapart)
        print(Zu-Za,Za-ZW,LDRa)
        i = i+1
        #ax2.plot(sizes,intZu,label=str(1.0/lam))
        #ax3.plot(sizes,intZa,label=str(1.0/lam))
        #ax4.plot(sizes,intZW,label=str(1.0/lam))
        
    #ax2.legend()
    #ax3.legend()
    #ax4.legend()
    #ax2.set_title('Ku'+melt_frac)
    #ax3.set_title('Ka')
    #ax4.set_title('W')
    s=ax.scatter(data['Ka_W'],data['Ku_Ka'],label=melt2perc(melt_frac))#,c=data['LDRka'])
    #data.plot(x='Ka_W',y='Ku_Ka',kind='scatter',ax=ax,label=melt_frac)
    datatot = datatot.append(data)

########ax.legend()  
#plt.colorbar(mappable=s,ax=ax)

###############################################################################

freqs = {'000':'C','001':'X','002':'Ku','003':'Ka','004':'W','005':'89','006':'157'}
freqs = {'001':'X','002':'Ku','003':'Ka','004':'W'}
scattfolder = '/work/DBs/0/'

cols = ['Dmax','X','Ku','Ka','W','ldr']#,'XKa','KaW']
subfolders = glob(scattfolder+'*')
data00 = pd.DataFrame(index=[x[len(scattfolder):] for x in subfolders],columns=cols)
for subfolder in subfolders:
    Dstr = subfolder[len(scattfolder):]
    D = int(Dstr)*1e-3
    print(subfolder, D)
    data00.loc[Dstr,'Dmax'] = D
    for freqidx in freqs.keys():
        sfld = glob(subfolder+'/run'+freqidx+'*')[0]
        scatt = scattering.ScattADDA(logfile=sfld+'/log',muellerfile=sfld+'/mueller',
                             csfile=sfld+'/CrossSec', D=D)
        data00.loc[Dstr,freqs[freqidx]] = scatt.sig_bk
        if freqidx == '003':
            data00.loc[Dstr,'ldr'] = scatt.ldr

scattfolder = '/work/DBs/10/'

subfolders = glob(scattfolder+'*')
data10 = pd.DataFrame(index=[x[len(scattfolder):] for x in subfolders],columns=cols)
for subfolder in subfolders:
    Dstr = subfolder[len(scattfolder):]
    D = int(Dstr)*1e-3
    data10.loc[Dstr,'Dmax'] = D
    print(subfolder, D)
    for freqidx in freqs.keys():
        sfld = glob(subfolder+'/run'+freqidx+'*')[0]
        scatt = scattering.ScattADDA(logfile=sfld+'/log',muellerfile=sfld+'/mueller',
                             csfile=sfld+'/CrossSec', D=D)
        data10.loc[Dstr,freqs[freqidx]] = scatt.sig_bk
        if freqidx == '003':
            data10.loc[Dstr,'ldr'] = scatt.ldr

data00 = data00.astype(float)
data10 = data10.astype(float)
data00.sort(columns='Dmax',inplace=True)
data10.sort(columns='Dmax',inplace=True)
data00['XKa' ] = coeffx*data00.X/(coeffa*data00.Ka)
data00['KuKa'] = coeffu*data00.Ku/(coeffa*data00.Ka)
data00['KaW' ] = coeffa*data00.Ka/(coeffW*data00.W)
data10['XKa' ] = coeffx*data10.X/(coeffa*data10.Ka)
data10['KuKa'] = coeffu*data10.Ku/(coeffa*data10.Ka)
data10['KaW' ] = coeffa*data10.Ka/(coeffW*data10.W)

plt.figure()
plt.plot(data00.Dmax,data00.ldr)
plt.plot(data10.Dmax,data10.ldr)
axx = plt.gca()
axx.set_yscale('log')

plt.figure()
plt.plot(data00.Dmax,(data00.XKa),label='0 XKa')
plt.plot(data10.Dmax,(data10.XKa),label='10 XKa')
plt.plot(data00.Dmax,(data00.KuKa),label='0 KuKa')
plt.plot(data10.Dmax,(data10.KuKa),label='10 KuKa')
plt.plot(data00.Dmax,(data00.KaW),label='0 KaW')
plt.plot(data10.Dmax,(data10.KaW),label='10 KaW')
plt.legend()
axx = plt.gca()
axx.set_yscale('log')

plt.figure()
plt.plot(data00.Dmax,(data10.XKa)-(data00.XKa),label='10-0 XKa')
plt.plot(data00.Dmax,(data10.KuKa)-(data00.KuKa),label='10-0 KuKa')
plt.plot(data00.Dmax,(data10.KaW)-(data00.KaW),label='10-0 KaW')
plt.legend()
axx = plt.gca()
axx.set_yscale('log')

#lambdas = 1.0/np.linspace(0.05,4,5) #13
Nexp = lambda l,x: np.exp(-l*x)

Z00 = pd.DataFrame(index=lambdas,columns=['Dm','X','Ku','Ka','W','XKa','KuKa','KaW','ldr'])
Z10 = pd.DataFrame(index=lambdas,columns=['Dm','X','Ku','Ka','W','XKa','KuKa','KaW','ldr'])




for lam in lambdas:
#     conc  = psd.ExponentialPSD(Lambda=lam,D_max=max(sizes)+0.1)
     D0 = (3.67+mu)/lam
     print(D0)
     conc  = psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=max(sizes)+0.1)
#     conc00 = Nexp(lam,data00.Dmax)
#     conc10 = Nexp(lam,data10.Dmax)
     conc00 = conc(data00.Dmax)
     conc10 = conc(data10.Dmax)
     #plt.figure()
     #axx = plt.gca()
     #axx.plot(data00.Dmax,(data00.X*conc00*coeffx ), label='X ')
     #axx.plot(data00.Dmax,(data00.Ku*conc00*coeffu), label='Ku')
     #axx.plot(data00.Dmax,(data00.Ka*conc00*coeffa), label='Ka')
     #axx.plot(data00.Dmax,(data00.W*conc00*coeffW ), label='W ')
     #axx.set_yscale('log')
     #axx.legend()
     Z00.loc[lam,'X' ] = 10.0*np.log10((data00.X*conc00 ).sum()*coeffx)
     Z00.loc[lam,'Ku'] = 10.0*np.log10((data00.Ku*conc00).sum()*coeffu)
     Z00.loc[lam,'Ka'] = 10.0*np.log10((data00.Ka*conc00).sum()*coeffa)
     Z00.loc[lam,'W' ] = 10.0*np.log10((data00.W*conc00 ).sum()*coeffW)
     Z00.loc[lam,'ldr' ] = 10.0*np.log10((data00.ldr*conc00 ).sum())
     Z10.loc[lam,'X' ] = 10.0*np.log10((data10.X*conc10 ).sum()*coeffx)
     Z10.loc[lam,'Ku'] = 10.0*np.log10((data10.Ku*conc10).sum()*coeffu)
     Z10.loc[lam,'Ka'] = 10.0*np.log10((data10.Ka*conc10).sum()*coeffa)
     Z10.loc[lam,'W' ] = 10.0*np.log10((data10.W*conc10 ).sum()*coeffW)
     Z10.loc[lam,'ldr' ] = 10.0*np.log10((data10.ldr*conc10 ).sum())



Z00['XKa' ] = Z00.X-Z00.Ka
Z00['KuKa'] = Z00.Ku-Z00.Ka
Z00['KaW' ] = Z00.Ka-Z00.W
Z10['XKa' ] = Z10.X-Z10.Ka
Z10['KuKa'] = Z10.Ku-Z10.Ka
Z10['KaW' ] = Z10.Ka-Z10.W

plt.figure()
axx = plt.gca()
plt.plot(Z00.KaW,Z00.XKa,label='dry ')
plt.plot(Z10.KaW,Z10.XKa,label='10 %')
plt.legend()

###############################################################################

plt.figure()
plt.scatter(datatot['Ka_W'],datatot['Ku_Ka'])#,c=datatot['LDRka'],cmap=plt.cm.jet)
plt.scatter(Z00.KaW,Z00.XKa)#,c=Z00.ldr,cmap=plt.cm.jet)
plt.scatter(Z10.KaW,Z10.XKa)#,c=Z10.ldr,cmap=plt.cm.jet)
plt.scatter(Z00.KaW,Z00.KuKa)#,c=Z00.ldr,cmap=plt.cm.jet)
plt.scatter(Z10.KaW,Z10.KuKa)#,c=Z10.ldr,cmap=plt.cm.jet)
plt.title('3freq')# + LDR')
#plt.colorbar(label='LDR')
plt.xlabel('DWR$_{Ka,W}$')
plt.ylabel('DWR$_{RAY,Ka}$')
plt.savefig('3freq+LDR.png',dpi=600)

ax.plot(Z00.KaW,Z00.XKa,label='X dry')#,c=Z00.ldr,cmap=plt.cm.jet)
ax.plot(Z10.KaW,Z10.XKa,label='X 10%')#,c=Z10.ldr,cmap=plt.cm.jet)
ax.plot(Z00.KaW,Z00.KuKa,label='Ku dry')#,c=Z00.ldr,cmap=plt.cm.jet)
ax.plot(Z10.KaW,Z10.KuKa,label='Ku 10%')#,c=Z10.ldr,cmap=plt.cm.jet)
ax.legend()
ax.grid()
ax.set_xlim([0,20])
ax.set_ylim([0,22])
ax.set_title('mu =  '+str(mu))


plt.figure()
ag = plt.gca()
D0s = np.linspace(0.05,8,10)
Z00 = pd.DataFrame(index=D0s,columns=['Dm','X','Ku','Ka','W','XKa','KuKa','KaW','ldr'])
Z10 = pd.DataFrame(index=D0s,columns=['Dm','X','Ku','Ka','W','XKa','KuKa','KaW','ldr'])
for mu in [-1.,0.,1.,2.,3.,4.]:
    for D0 in D0s:
        print(D0)
        conc  = psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=max(sizes)+0.1)
        conc00 = conc(data00.Dmax)
        conc10 = conc(data10.Dmax)
        Z00.loc[D0,'X' ] = 10.0*np.log10((data00.X*conc00 ).sum()*coeffx)
        Z00.loc[D0,'Ku'] = 10.0*np.log10((data00.Ku*conc00).sum()*coeffu)
        Z00.loc[D0,'Ka'] = 10.0*np.log10((data00.Ka*conc00).sum()*coeffa)
        Z00.loc[D0,'W' ] = 10.0*np.log10((data00.W*conc00 ).sum()*coeffW)
        Z00.loc[D0,'ldr' ] = 10.0*np.log10((data00.ldr*conc00 ).sum())
#        Z10.loc[lam,'X' ] = 10.0*np.log10((data10.X*conc10 ).sum()*coeffx)
#        Z10.loc[lam,'Ku'] = 10.0*np.log10((data10.Ku*conc10).sum()*coeffu)
#        Z10.loc[lam,'Ka'] = 10.0*np.log10((data10.Ka*conc10).sum()*coeffa)
#        Z10.loc[lam,'W' ] = 10.0*np.log10((data10.W*conc10 ).sum()*coeffW)
#        Z10.loc[lam,'ldr' ] = 10.0*np.log10((data10.ldr*conc10 ).sum())
    Z00['XKa' ] = Z00.X-Z00.Ka
    Z00['KuKa'] = Z00.Ku-Z00.Ka
    Z00['KaW' ] = Z00.Ka-Z00.W
    Z10['XKa' ] = Z10.X-Z10.Ka
    Z10['KuKa'] = Z10.Ku-Z10.Ka
    Z10['KaW' ] = Z10.Ka-Z10.W
    ag.plot(Z00.KaW,Z00.XKa,label='mu= '+str(mu))#,c=Z00.ldr,cmap=plt.cm.jet)
    #ax.plot(Z10.KaW,Z10.XKa,label='X 10%')#,c=Z10.ldr,cmap=plt.cm.jet)
    #ax.plot(Z00.KaW,Z00.KuKa,label='Ku dry')#,c=Z00.ldr,cmap=plt.cm.jet)
    #ax.plot(Z10.KaW,Z10.KuKa,label='Ku 10%')#,c=Z10.ldr,cmap=plt.cm.jet)
ag.legend()
ag.grid()
ag.set_xlabel('DWR$_{Ka W}$')
ag.set_ylabel('DWR$_{X Ka}$')





        
#        plt.figure()
#        plt.scatter(array[:,0],array[:,1],c=np.log10(array[:,2]))
#        plt.colorbar(label='log$_{10}$   mass   [mg]')
#        ax = plt.gca()
#        ax.set_ylabel('backscattering xsect [mm$^2$]')
#        ax.set_xlabel('maximum dimension [mm]')
#        ax.set_xlim([1e-1,2e1])
#        ax.set_ylim([1e-11,1e1])
#        ax.set_xscale('log')
#        ax.set_yscale('log')
#        ax.set_title(freq+'   GHz')
#        
#        plt.figure()
#        plt.scatter(array[:,0],array[:,1],c=array[:,3])
#        plt.colorbar(label='melted fraction')
#        ax = plt.gca()
#        ax.set_ylabel('backscattering xsect [mm$^2$]')
#        ax.set_xlabel('maximum dimension [mm]')
#        ax.set_xlim([1e-1,2e1])
#        ax.set_ylim([1e-11,1e1])
#        ax.set_xscale('log')
#        ax.set_yscale('log')
#        
#        plt.figure()
#        plt.scatter(array[:,3],array[:,3])
#        DF = pd.DataFrame(array,columns=['D','s','m','f','l'])
#        print(DF.shape,len(DF.f.drop_duplicates()))
#        
#        plt.figure()
#        plt.scatter(array[:,0],array[:,4],c=array[:,3])
#        plt.colorbar(label='melted fraction')
#        ax = plt.gca()
#        ax.set_ylabel('ldr')
#        ax.set_xlabel('maximum dimension [mm]')
#        ax.set_xlim([1e-1,2e1])
#        ax.set_ylim([1e-3,1e0])
#        ax.set_xscale('log')
#        ax.set_yscale('log')


#mm = [melt2perc(x) for x in melt_fracs]
#
#plt.figure()        
#for melt in melt_fracs:
#    subdata = datatot[datatot['melt']==melt]
#    plt.plot(subdata['Dm'],subdata['LDRka'],label=melt2perc(melt))
#plt.xlabel('$\Lambda^{-1}   [mm]$')
#plt.ylabel('LDR   Ka   [dBZ]')
#plt.legend(title='melted fraction')
#plt.grid()
#plt.savefig('ldr_melt.pdf')
#plt.savefig('ldr_melt.png',dpi=600)
#
#plt.figure()        
#for lam in datatot.Dm.drop_duplicates()[0:-1:3]:
#    subdata = datatot[datatot['Dm']==lam]
#    plt.plot(mm,subdata['LDRka'],label=str(lam)[0:3])
#plt.xlabel('melted fraction    [%]')
#plt.ylabel('LDR   Ka   [dBZ]')
#plt.legend(title='Dm')
#plt.grid()
#plt.savefig('ldr_dm.pdf')
#plt.savefig('ldr_dm.png',dpi=600)