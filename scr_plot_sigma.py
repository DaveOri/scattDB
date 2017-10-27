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

tablesfolder = 'tables/'
authors = ['DO','BJ','JL']
author = 'DO'
melt = '00'

dataDO = pd.read_csv(tablesfolder+author+'_'+melt+'.csv')
dataDO.drop('Unnamed: 0',axis=1,inplace=True)
meltDO = pd.read_csv(tablesfolder+author+'_'+'10'+'.csv')
meltDO.drop('Unnamed: 0',axis=1,inplace=True)

author = 'BJ'
dataBJ2 = pd.read_csv(tablesfolder+author+'_agg2_0'+'.csv')
dataBJ2.drop('Unnamed: 0',axis=1,inplace=True)

dataBJ3 = pd.read_csv(tablesfolder+author+'_agg3_0'+'.csv')
dataBJ3.drop('Unnamed: 0',axis=1,inplace=True)

author = 'dataJL_B1.0'
dataJL = pd.read_csv(tablesfolder+author+'.csv')

plt.figure()
ax = plt.gca()
#ax.scatter(dataBJ2.Dmax,dataBJ2.mkg,label='BJ2')
#ax.scatter(dataBJ3.Dmax,dataBJ3.mkg,label='BJ3')
ax.scatter(dataJL.Dmax, dataJL.mkg, label='JL')
ax.scatter(dataDO.Dmax, dataDO.mkg, label='DO')
ax.plot(dataDO.Dmax,8.9e-5*dataDO.Dmax**2.1,label='brandes')
ax.plot(dataDO.Dmax,0.02414e-3*dataDO.Dmax**1.9,label='BF95')
ax.legend()
ax.grid()
ax.set_yscale('log')
ax.set_xscale('log')

def IntPsd(D,s,psd):
    conc = psd(D)
    intZ = np.multiply(conc,s)
    return integrate.trapz(intZ,D)

def myGammaPSD(D0,mu):
    return psd.GammaPSD(D0=D0,Nw=1.,mu=mu,D_max=22)

def Gauss(mu,sig):
    def F(D):
        x = D-mu 
        return np.exp(x*x/(2.0*sig*sig))/(sig*np.sqrt(2.0*np.pi))
    return F

D0s = np.linspace(0.1,20.,30)
mus = [-1.,0.,1.,2.,3.]
#mus = [0.1,0.3,0.5,0.8,1]
marks=[',','+','h','v','.']

def f3plot(data,title='title',color=None):
    plt.figure()
    ax = plt.gca()
    for mu,m in zip(mus,marks):
        Zx, Za, Zw, XKa, KaW, LDR, IWC = ([] for i in range(7))
        for D0 in D0s:
            conc  = myGammaPSD(D0,mu)
            #conc = Gauss(D0,mu)
            iwc = IntPsd(data.Dmax,data.mkg,conc)
            IWC.append(iwc)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax,data.X,conc)/iwc)
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax,data.Ku,conc)/iwc)
            
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            Za.append(10.*np.log10(coeffa*IntPsd(data.Dmax,data.Ka,conc)/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax,data.W,conc)/iwc))
            LDR.append( 10.*np.log10(IntPsd(data.Dmax,data.ldr,conc)))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)

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
        else:
            cbarlabel = color + '/iwc    [dB/(g/m$^3$)]'
            if color == 'Zx':
                colorv = Zx
                print(max(Zx))
                vmin = 20
                vmax = 40
            elif color == 'Zka':
                colorv = Za
                print(max(Za))
                vmin = 10
                vmax = 25
            elif color == 'Zw':
                colorv = Zw
                print(max(Zw))
                vmin = 0
                vmax = 15

        s = ax.scatter(KaW,XKa,c=colorv,marker=m,label='$\mu=$ '+str(mu),vmin=vmin,vmax=vmax,cmap='jet')
    ax.legend()
    ax.grid()
    ax.set_xlabel('DWR$_{Ka,W}$')
    ax.set_ylabel('DWR$_{X,Ka}$')
    ax.set_xlim([-1,15])
    ax.set_ylim([0,22])
    colorbar = plt.colorbar(mappable=s,ax=ax)
    colorbar.set_label(cbarlabel)
    ax.set_title(title)

f3plot(dataDO,'Davide dry')
f3plot(dataDO,'Davide dry',color='Zx')
f3plot(dataBJ2,'BJ2 dry',color='Zx')
#f3plot(dataBJ3,'BJ3 dry',color='Zw')
dataJL.columns = ['Dmax','model','ELWP','mkg','Dmax.1','Rgyr','ldr','riming', 'Xa','Xs', 'X', 'Xe', 'Ua', 'Us', 'Ku', 'Ue', 'Aa', 'As', 'Ka', 'Ae', 'Wa','Ws', 'W', 'We']
dataJL.X  = dataJL.X*1.0e6
dataJL.Ku = dataJL.Ku*1.0e6
dataJL.Ka = dataJL.Ka*1.0e6
dataJL.W  = dataJL.W*1.0e6
f3plot(dataJL,'Jussi unrimed',color='Zx')
f3plot(dataJL,'Jussi unrimed')

f3plot(dataDO,'Davide dry',color='Zka')
f3plot(dataJL,'Jussi unrimed',color='Zka')

f3plot(dataDO,'Davide dry',color='Zw')
f3plot(dataJL,'Jussi unrimed',color='Zw')

f3plot(dataDO,'Davide dry',color='ldr')
f3plot(meltDO,'Davide melt',color='ldr')