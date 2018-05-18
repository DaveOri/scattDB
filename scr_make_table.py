# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 11:50:26 2017

@author: dori
"""

from scattDB import shape
from glob import glob
import pandas as pd
import numpy as np
from scattDB import scattering

shapefolder = '/work/DBs/melted_aggregate_shape_files/'
scattfolder = '/work/DBs/melted_aggregate_scaled_reff_Ku_Ka_W_89_165_183/'

#melt_fracs = ['000001','001010','002004','004988','007029']
#melt_fracs = ['000001','001017','001315']
melt_fracs = ['000001','001010','002004','003030','004073','004988','007029']
freqs = ['13.4','35.6','94']

melt2perc = lambda x: int(int(x)*0.01)

Zlab = {'13.4':'Ku','35.6':'Ka','94':'W'}
cols = ['Dmax','X','Ku','Ka','W','ldr','xpolKa','mkg','CextC','CscaC','CabsC','CextX','CscaX','CabsX','CextKu','CscaKu','CabsKu','CextKa','CscaKa','CabsKa','CextW','CscaW','CabsW']
indexes = np.arange(10000)
for melt_frac in melt_fracs:
    data = pd.DataFrame(index=indexes, columns=cols)
    for freq in freqs:
        scatt_folders = sorted(glob(scattfolder+'*aggregate3*'+'_f'+melt_frac+'_*'+freq))
        dist = scattering.ScattDist()
        freqidx = Zlab[freq]
        i=0
        for fld in scatt_folders:
            shpstr, aeffstr = fld.split('/')[-1].split('AEFF')
            shape_folders = glob(shapefolder+shpstr+'*')
            shapefile = shape_folders[0] + '/shape.dat'
            shp = shape.shapeDDA(shapefile)
            
            avgfiles = sorted(glob(fld+'/*.avg'))
            for avg in avgfiles:
                scatt = scattering.ScattDDSCAT(avg) 
                scatt.D = shp.find_dmax(aeff=scatt.aeff)
                scatt.mass = shp.find_mass(aeff=scatt.aeff)
                scatt.melt = shp.get_melted_fraction()
                data.iloc[i][Zlab[freq]] = scatt.sig_bk
                data.iloc[i]['Dmax'] = scatt.D
                data.iloc[i]['mkg'] = scatt.mass
                if freq == '35.6':
                    data.iloc[i]['ldr'] = scatt.ldr
                    data.iloc[i]['xpolKa'] = scatt.xpol_xsect
                i = i + 1
    data.sort('Dmax',inplace=True)
    data.dropna(how='all',inplace=True)
    data.to_csv('tables/BJ_agg3_'+str(melt2perc(melt_frac))+'.csv')

freqs = {'9.6':'X','13.6':'Ku','35.6':'Ka','94':'W'}
scattfolder = '/data/optimice/scattering_databases/DavideOri_2014/melted/'
subfolders = glob(scattfolder+'*')
subfolders = [x for x in subfolders if '11938.9' not in x]
subfolders = [x for x in subfolders if '12973' not in x]
subfolders = [x for x in subfolders if '13885.6' not in x]
melt_fracs = [10,20,30,40,50,60,70]
for melt_frac in melt_fracs:
  data = pd.DataFrame(index=[x[len(scattfolder):] for x in subfolders],columns=cols)
  for subfolder in subfolders:
    Dstr = subfolder[len(scattfolder):]
    D = float(Dstr)*1e-3
    print(subfolder, D)
    data.loc[Dstr,'Dmax'] = D
    for freqidx in freqs.keys():
      datafolder = subfolder + '/1/' + freqidx +'/'+ str(melt_frac)+'/'
      scatt = scattering.ScattADDA(logfile=datafolder+'log',
                                   muellerfile=datafolder+'mueller',
                                   csfile=datafolder+'CrossSec', D=D)
      data.loc[Dstr,freqs[freqidx]] = scatt.sig_bk*1e12
      data.loc[Dstr,'mkg'] = scatt.mass*1e18
      data.loc[Dstr,'Cext'+freqs[freqidx]] = scatt.sig_ext*1e12
      data.loc[Dstr,'Csca'+freqs[freqidx]] = scatt.sig_sca*1e12
      data.loc[Dstr,'Cabs'+freqs[freqidx]] = scatt.sig_abs*1e12
      if freqidx == '35.6':
        data.loc[Dstr,'ldr'] = scatt.ldr
        data.loc[Dstr,'xpolKa'] = scatt.xpol_xsect*1e12
  data = data.astype(float)
  data.sort(columns='Dmax',inplace=True)
  data.set_index('Dmax',inplace=True)
  data.to_csv('tables/meltedDO'+str(melt_frac)+'.csv')


freqs = {'000':'C','001':'X','002':'Ku','003':'Ka','004':'W','005':'89','006':'157'}
freqs = {'001':'X','002':'Ku','003':'Ka','004':'W'}
scattfolder = '/data/optimice/scattering_databases/DavideOri_2014/0/'

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
        data00.loc[Dstr,'mkg'] = scatt.mass
        data00.loc[Dstr,'Cext'+freqs[freqidx]] = scatt.sig_ext
        data00.loc[Dstr,'Csca'+freqs[freqidx]] = scatt.sig_sca
        data00.loc[Dstr,'Cabs'+freqs[freqidx]] = scatt.sig_abs
        if freqidx == '003':
            data00.loc[Dstr,'ldr'] = scatt.ldr
            data00.loc[Dstr,'xpolKa'] = scatt.xpol_xsect

scattfolder = '/data/optimice/scattering_databases/DavideOri_2014/10/'
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
        data10.loc[Dstr,'mkg'] = scatt.mass
        data10.loc[Dstr,'Cext'+freqs[freqidx]] = scatt.sig_ext
        data10.loc[Dstr,'Csca'+freqs[freqidx]] = scatt.sig_sca
        data10.loc[Dstr,'Cabs'+freqs[freqidx]] = scatt.sig_abs
        if freqidx == '003':
            data10.loc[Dstr,'ldr'] = scatt.ldr
            data10.loc[Dstr,'xpolKa'] = scatt.xpol_xsect

data00 = data00.astype(float)
data10 = data10.astype(float)
data00.sort(columns='Dmax',inplace=True)
data10.sort(columns='Dmax',inplace=True)
data00.set_index('Dmax',inplace=True)
data10.set_index('Dmax',inplace=True)
data00.to_csv('tables/DO_00.csv')
data10.to_csv('tables/DO_10.csv')

from sys import path
path.append('/home/dori/develop/pyPamtra2/singleScattering/singleScattering/')
path.append('/home/dori/develop/pyPamtra2/refractiveIndex/refractiveIndex/')
from self_similar_rayleigh_gans import backscattering
from refractive import ice
diameters = np.linspace(0.4,40,200)
dataSSRG = pd.DataFrame(index=diameters,columns=cols)
dataSSRG['Dmax'] = diameters
frequencies = {'X':9.6,'Ku':13.6,'Ka':35.6,'W':94.0}
for freqidx in frequencies.keys():
    frGHz = frequencies[freqidx]
    frHz = frGHz*1.0e9
    for diam in diameters:
        dataSSRG.loc[diam,freqidx],dataSSRG.loc[diam,'mkg'] = backscattering(frHz,diam*0.001,ice.n(270.,frHz))

dataSSRG.to_csv('tables/DO_SSRG.csv')

dataSSRG.index.name = 'Dmax'
import matplotlib.pyplot as plt
plt.close('all')
plt.figure()
ax = plt.gca()
data00[['X','Ku','Ka','W']].plot(ax=ax,linestyle='--')
dataSSRG[['X','Ku','Ka','W']].plot(ax=ax)
ax.set_yscale('log')
ax.set_ylabel('C$_{bk}$    [mm$^2$]')

#brandes = lambda D: 7.9e-5*D**2.1
#smalles = lambda D: 4.1e-5*D**2.5
plt.close('all')
plt.figure()
ax = plt.gca()
data00['mkg'].plot(ax=ax,linestyle='-.')
dataSSRG['mkg'].plot(ax=ax)
#ax.plot(data00.index.values,brandes(data00.index.values))
#ax.plot(data00.index.values,smalles(data00.index.values))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('mass     [g]')
ax.grid()
