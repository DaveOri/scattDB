# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:41:29 2017

@author: dori
"""

import pandas as pd
from glob import glob
from scattDB import scattering
import matplotlib.pyplot as plt
import numpy as np

freqs = {'000':'C','001':'X','002':'Ku','003':'Ka','004':'W','005':'89','006':'157'}
freqs = {'001':'X','003':'Ka','004':'W'}
scattfolder = '/work/DBs/0/'

cols = ['Dmax','X','Ka','W','ldr']#,'XKa','KaW']
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
data00['XKa'] = data00.X/data00.Ka
data00['KaW'] = data00.Ka/data00.W
data10['XKa'] = data10.X/data10.Ka
data10['KaW'] = data10.Ka/data10.W

plt.figure()
plt.plot(data00.Dmax,data00.ldr)
plt.plot(data10.Dmax,data10.ldr)
ax = plt.gca()
ax.set_yscale('log')

plt.figure()
plt.plot(data00.Dmax,data00.XKa,label='0 XKa')
plt.plot(data10.Dmax,data10.XKa,label='10 XKa')
plt.plot(data00.Dmax,data00.KaW,label='0 KaW')
plt.plot(data10.Dmax,data10.KaW,label='10 KaW')
plt.legend()
ax = plt.gca()
ax.set_yscale('log')

c = 299792458000. # mm/s
lamx = c/(9.6*1e9)
lama = c/(35.6*1e9)
lamW = c/(94*1e9)
coeffx = lamx**4./(0.93*np.pi**5.)
coeffa = lama**4./(0.93*np.pi**5.)
coeffW = lamW**4./(0.75*np.pi**5.)

lambdas = 1.0/np.linspace(0.05,4,50) #13
Nexp = lambda l,x: np.exp(-l*x)

Z00 = pd.DataFrame(index=lambdas,columns=['Dm','X','Ka','W','KuKa','KaW'])
Z10 = pd.DataFrame(index=lambdas,columns=['Dm','X','Ka','W','KuKa','KaW'])

for lam in lambdas:
     conc00 = Nexp(lam,data00.Dmax)
     conc10 = Nexp(lam,data10.Dmax)
     #Dm = (flakes.mkg*flakes.Dmax*conc).sum()/((flakes.mkg*conc).sum())
     #Z.loc[lam,'Dm'] = Dm
     Z00.loc[lam,'X' ] = 10.0*np.log10((data00.X*conc00 ).sum()*coeffx)
     Z00.loc[lam,'Ka'] = 10.0*np.log10((data00.Ka*conc00).sum()*coeffa)
     Z00.loc[lam,'W' ] = 10.0*np.log10((data00.W*conc00 ).sum()*coeffW)
     Z10.loc[lam,'X' ] = 10.0*np.log10((data10.X*conc10 ).sum()*coeffx)
     Z10.loc[lam,'Ka'] = 10.0*np.log10((data10.Ka*conc10).sum()*coeffa)
     Z10.loc[lam,'W' ] = 10.0*np.log10((data10.W*conc10 ).sum()*coeffW)

Z00['XKa'] = Z00.X-Z00.Ka
Z00['KaW'] = Z00.Ka-Z00.W
Z10['XKa'] = Z10.X-Z10.Ka
Z10['KaW'] = Z10.Ka-Z10.W

plt.figure()
ax = plt.gca()
plt.plot(Z00.KaW,Z00.XKa,label='dry ')
plt.plot(Z10.KaW,Z10.XKa,label='10 %')
plt.legend()