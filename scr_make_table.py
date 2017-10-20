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
#from scattDB import psd

#from scipy import integrate
shapefolder = '/work/DBs/melted_aggregate_shape_files/'
scattfolder = '/work/DBs/melted_aggregate_scaled_reff_Ku_Ka_W_89_165_183/'

#melt_fracs = ['000001','001010','002004','004988','007029']
#melt_fracs = ['000001','001017','001315']
melt_fracs = ['000001','001010','002004','003030','004073','004988','007029']
freqs = ['13.4','35.6','94']

melt2perc = lambda x: int(int(x)*0.01)

Zlab = {'13.4':'Ku','35.6':'Ka','94':'W'}
Llab = {'13.4':'LKu','35.6':'LKa','94':'LW'}

c = 299792458000. # mm/s
lamx = c/(9.6*1e9)
lamu = c/(13.4*1e9)
lama = c/(35.6*1e9)
lamW = c/(94*1e9)
coeffx = lamx**4./(0.93*np.pi**5.)
coeffu = lamu**4./(0.95*np.pi**5.)
coeffa = lama**4./(0.95*np.pi**5.)
coeffW = lamW**4./(0.75*np.pi**5.)

cols = ['Dmax','X','Ku','Ka','W','ldr']
indexes = np.arange(10000)
for melt_frac in melt_fracs:
    data = pd.DataFrame(index=indexes, columns=cols)
    for freq in freqs:
        scatt_folders = sorted(glob(scattfolder+'*aggregate2*'+'_f'+melt_frac+'_*'+freq))
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
                if freq == '35.6':
                    data.iloc[i]['ldr'] = scatt.ldr
                i = i + 1
    data.sort('Dmax',inplace=True)
    data.dropna(how='all',inplace=True)
    data.to_csv('tables/BJ_agg2_'+str(melt2perc(melt_frac))+'.csv')


freqs = {'000':'C','001':'X','002':'Ku','003':'Ka','004':'W','005':'89','006':'157'}
freqs = {'001':'X','002':'Ku','003':'Ka','004':'W'}
scattfolder = '/work/DBs/0/'

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
data00.to_csv('tables/DO_00.csv')
data10.to_csv('tables/DO_10.csv')