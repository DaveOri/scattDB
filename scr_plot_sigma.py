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

def IntPsd(D,s,psd):
    conc = psd(D)
    intZ = np.multiply(conc,s)
    return integrate.trapz(intZ,D)

def myGammaPSD(D0,mu):
    return psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=22)



def f3plot(data,title='title',color=None):
    D0s = np.linspace(0.1,20.,20)
    mus = [-1.,0.,1.,2.,3.]
    marks=[',','+','h','v','.']
    plt.figure()
    ax = plt.gca()
    for mu,m in zip(mus,marks):
        Zx, Za, Zw, XKa, KaW, LDR = ([] for i in range(6))
        
        for D0 in D0s:
#            conc  = psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=22)
            conc  = myGammaPSD(D0,mu)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax,data.X,conc))
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax,data.Ku,conc))
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            Za.append(10.*np.log10(coeffa*IntPsd(data.Dmax,data.Ka,conc)))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax,data.W,conc)))
            LDR.append( 10.*np.log10(IntPsd(data.Dmax,data.ldr,conc)))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)
        cbarlabel = color
        if color is None:
            colorv = D0s
            vmin=0
            vmax=20
            cbarlabel = '$D_0$'
        elif color == 'ldr':
            colorv = LDR
            vmin = -30
            vmax = -20
            cbarlabel = 'LDR$_Ka$'
        elif color == 'Zx':
            colorv = Zx
            print(max(Zx))
            vmin = -20
            vmax = 25
        elif color == 'Zka':
            colorv = Za
            print(max(Za))
            vmin = -20
            vmax = 5
        elif color == 'Zw':
            colorv = Zw
            print(max(Zw))
            vmin = -30
            vmax = 0

        s = ax.scatter(KaW,XKa,c=colorv,marker=m,label='$\mu=$ '+str(mu),vmin=vmin,vmax=vmax)
    ax.legend()
    ax.grid()
    ax.set_xlabel('DWR$_{Ka,W}$')
    ax.set_ylabel('DWR$_{X,Ka}$')
    ax.set_xlim([-1,15])
    ax.set_ylim([0,22])
    colorbar = plt.colorbar(mappable=s,ax=ax)
    colorbar.set_label(cbarlabel)
    ax.set_title(title)

f3plot(dataDO,'Davide dry',color='Zw')
f3plot(dataBJ2,'BJ2 dry',color='Zw')
f3plot(dataBJ3,'BJ3 dry',color='Zw')
dataJL.columns = ['Dmax','model','ELWP','mkg','Dmax.1','Rgyr','ldr','riming', 'Xa','Xs', 'X', 'Xe', 'Ua', 'Us', 'Ku', 'Ue', 'Aa', 'As', 'Ka', 'Ae', 'Wa','Ws', 'W', 'We']
dataJL.X  = dataJL.X*1.0e6
dataJL.Ku = dataJL.Ku*1.0e6
dataJL.Ka = dataJL.Ka*1.0e6
dataJL.W  = dataJL.W*1.0e6
f3plot(dataJL,'Jussi unrimed',color='Zw')
