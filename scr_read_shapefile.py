# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:28:05 2017

@author: dori
"""

from scattDB import shape
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#from importlib import reload

#filename = '/work/DBs/melted_aggregate_shape_files/melt3a_aggregate2_010_20091210_222748_halfscale_f000001_AEFF_1000_ROT535_13.4/shape.dat'
#shp = shape.shapeDDA(filename)
#shp.draw()

shapefolder = '/work/DBs/melted_aggregate_shape_files/'
folders = glob(shapefolder+'*')
foldersto = glob(shapefolder+'*to*')
foldersNOTto = [x for x in folders if x not in foldersto]

DF = pd.DataFrame(index=range(len(foldersNOTto)),columns=['dmax','mass','melt'])
i = 0
for fld in foldersNOTto:
    shp = shape.shapeDDA(fld+'/shape.dat')
    fldsplit = fld.split('_')
    mfrac = float(fldsplit[-5][1:])*0.0001
    aeffstr = fldsplit[-3] # conversion to mm
    if 'to' in aeffstr:
        amin, rest = aeffstr.split('to')
        amax, aN = rest.split('x')
        aeff = np.linspace(int(amin),int(amax),int(aN))
    else:
        aeff=int(aeffstr)
    mass = shp.find_mass(aeff)
    DF.loc[i,'dmax'] = shp.find_dmax()*shp.d*0.001
    DF.loc[i,'melt'] = shp.get_melted_fraction()
    DF.loc[i,'mass'] = mass
    print(i,DF.loc[i])
    i = i + 1

lgm = np.log10(DF.mass.values.astype(float))
lgd = np.log10(DF.dmax.values.astype(float))
b,a = np.polyfit(lgd,lgm,1)
a = 10**a
print(a,b)

plt.figure()
dx=np.linspace(0.1,18,100)
plt.scatter(DF.dmax,DF.mass,c=DF.melt)
ax = plt.gca()
ax.plot(dx,a*dx**b)
plt.colorbar(label='melted fraction')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('maximum dimension   [mm]')
ax.set_ylabel('mass     [mg]')


DF = pd.DataFrame(index=range(len(foldersto)*30),columns=['dmax','mass','melt'])
i = 0
for fld in foldersto:
    shp = shape.shapeDDA(fld+'/shape.dat')
    fldsplit = fld.split('_')
    mfrac = float(fldsplit[-5][1:])*0.0001
    aeffstr = fldsplit[-3] # conversion to mm
    if 'to' in aeffstr:
        amin, rest = aeffstr.split('to')
        amax, aN = rest.split('x')
        aeff = np.linspace(int(amin),int(amax),int(aN))
    else:
        aeff=int(aeffstr)
    mass = shp.find_mass(aeff)
    DF.loc[30*i:30*i+29,'dmax'] = shp.find_dmax()*shp.d*0.001
    DF.loc[30*i:30*i+29,'melt'] = shp.get_melted_fraction()
    DF.loc[30*i:30*i+29,'mass'] = mass
    print(i,DF.loc[i])
    i = i + 1


lgm = np.log10(DF.mass.values.astype(float))
lgd = np.log10(DF.dmax.values.astype(float))
b,a = np.polyfit(lgd,lgm,1)
a = 10**a
print(a,b)

plt.figure()
plt.scatter(DF.dmax,DF.mass,c=DF.melt)
ax = plt.gca()
ax.plot(dx,a*dx**b)
plt.colorbar(label='melted fraction')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('maximum dimension   [mm]')
ax.set_ylabel('mass     [mg]')
