# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:06:26 2017

@author: dori
"""

from scattDB import shape
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scattDB import scattering
from scattDB import psd

from scipy import integrate

shapefolder = '/work/DBs/melted_aggregate_shape_files/'
scattfolder = '/work/DBs/melted_aggregate_scaled_reff_Ku_Ka_W_89_165_183/'

#melt_fracs = ['000001','001010','002004','004988','007029']
#melt_fracs = ['000001','001017','001315']
melt_fracs = ['000001','002004','004988','007029']
freqs = ['13.4','35.6','94']

Zlab = {'13.4':'ZKu','35.6':'ZKa','94':'ZW'}
Llab = {'13.4':'LKu','35.6':'LKa','94':'LW'}

c = 299792458000. # mm/s
lamu = c/(13.4*1e9)
lama = c/(35.6*1e9)
lamW = c/(94*1e9)
coeffu = lamu**4./(0.95*np.pi**5.)
coeffa = lama**4./(0.95*np.pi**5.)
coeffW = lamW**4./(0.75*np.pi**5.)

lambdas = 1.0/np.linspace(0.5,3,10) #13

cols = ['Ku_Ka','Ka_W','LDRka','melt','Dm']

plt.figure()
ax = plt.gca()
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
    plt.figure()
    ax2 = plt.subplot(3,1,1)
    #plt.figure()
    ax3 = plt.subplot(3,1,2)
    #plt.figure()
    ax4 = plt.subplot(3,1,3)
    for lam in lambdas:
        conc  = psd.ExponentialPSD(Lambda=lam,D_max=max(sizes)+0.1)
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
        ax2.plot(sizes,intZu,label=str(1.0/lam))
        ax3.plot(sizes,intZa,label=str(1.0/lam))
        ax4.plot(sizes,intZW,label=str(1.0/lam))
        
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax2.set_title('Ku'+melt_frac)
    ax3.set_title('Ka')
    ax4.set_title('W')
    s=ax.scatter(data['Ka_W'],data['Ku_Ka'],c=data['LDRka'],label=melt_frac)
    #data.plot(x='Ka_W',y='Ku_Ka',kind='scatter',ax=ax,label=melt_frac)
    datatot = datatot.append(data)

ax.legend()  
plt.colorbar(mappable=s,ax=ax)

plt.figure()
plt.scatter(datatot['Ka_W'],datatot['Ku_Ka'],c=datatot['LDRka'],cmap=plt.cm.jet)
plt.title(i + "   " for i in melt_fracs)
plt.colorbar(label='LDR')
plt.xlabel('DWR$_{Ka,W}$')
plt.xlabel('DWR$_{Ku,Ka}$')
        
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