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

RHfile = '/home/dori/develop/scatdb/share/scatdb.csv'
dataRH = pd.read_csv(RHfile)
dataRH = dataRH[dataRH.flaketype == 20]
dataRH.set_index('max_dimension_mm',inplace=True)
dataRH.sort_index(inplace=True)

tol = 0.3
dataRHx = dataRH[abs(dataRH.frequencyghz-10.6) < tol]
dataRHu = dataRH[abs(dataRH.frequencyghz-13.6) < tol]
dataRHk = dataRH[abs(dataRH.frequencyghz-23.8) < tol]
dataRHa = dataRH[abs(dataRH.frequencyghz-35.6) < tol]
dataRHw = dataRH[abs(dataRH.frequencyghz-94)   < tol]

#dataRHx = dataRHx.loc[dataRHa.index]
#dataRHu = dataRHu.loc[dataRHa.index]
#dataRHk = dataRHk.loc[dataRHa.index]
#dataRHw = dataRHw.loc[dataRHa.index]

plt.close('all')

c = 299792458000. # mm/s
lamx = c/(9.6*1e9)
lamu = c/(13.4*1e9)
lamk = c/(24.*1e9)
lama = c/(35.6*1e9)
lamw = c/(94.*1e9)
coeffx = lamx**4./(0.93*np.pi**5.)
coeffu = lamu**4./(0.95*np.pi**5.)
coeffk = lamk**4./(0.95*np.pi**5.)
coeffa = lama**4./(0.95*np.pi**5.)
coeffw = lamw**4./(0.75*np.pi**5.)

tablesfolder = '/work/DBs/scattDB/tables/'
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

dataBJ2_10 = pd.read_csv(tablesfolder+author+'_agg2_10'+'.csv')
dataBJ2_10.drop('Unnamed: 0',axis=1,inplace=True)
dataBJ2_20 = pd.read_csv(tablesfolder+author+'_agg2_20'+'.csv')
dataBJ2_20.drop('Unnamed: 0',axis=1,inplace=True)
dataBJ2_30 = pd.read_csv(tablesfolder+author+'_agg2_30'+'.csv')
dataBJ2_30.drop('Unnamed: 0',axis=1,inplace=True)
dataBJ2_40 = pd.read_csv(tablesfolder+author+'_agg2_40'+'.csv')
dataBJ2_40.drop('Unnamed: 0',axis=1,inplace=True)
dataBJ2_49 = pd.read_csv(tablesfolder+author+'_agg2_49'+'.csv')
dataBJ2_49.drop('Unnamed: 0',axis=1,inplace=True)
dataBJ2_70 = pd.read_csv(tablesfolder+author+'_agg2_70'+'.csv')
dataBJ2_70.drop('Unnamed: 0',axis=1,inplace=True)

dataBJ3 = pd.read_csv(tablesfolder+author+'_agg3_0'+'.csv')
dataBJ3.drop('Unnamed: 0',axis=1,inplace=True)

author = 'dataJL_A0.5'
author = 'dataJL_A0.0'
dataJL = pd.read_csv(tablesfolder+author+'.csv')

Br07 = lambda x:0.08900e-3*x**2.1
BF95 = lambda x:0.02414e-3*x**1.9

plt.figure()
ax = plt.gca()
#ax.scatter(dataBJ2.Dmax,dataBJ2.mkg,label='BJ2')
#ax.scatter(dataBJ3.Dmax,dataBJ3.mkg,label='BJ3')
ax.scatter(dataJL.Dmax, dataJL.mkg, label='JL')
ax.scatter(dataDO.Dmax, dataDO.mkg, label='DO')
ax.plot(dataDO.Dmax,Br07(dataDO.Dmax),label='brandes')
ax.plot(dataDO.Dmax,BF95(dataDO.Dmax),label='BF95')
ax.legend()
ax.grid()
ax.set_yscale('log')
ax.set_xscale('log')

dataRH=pd.DataFrame(index=dataRHa.index,columns=dataDO.columns)
dataRH.Dmax = dataRH.index
dataRH.mkg = Br07(dataRH.Dmax)
dataRH.X  = dataRHx.cbk
dataRH.Ku = dataRHu.cbk
dataRH['K']  = dataRHk.cbk
dataRH.Ka = dataRHa.cbk
dataRH.W  = dataRHw.cbk
#%%

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
marks =[',','+','h','v','.']
colors=['C0','C1','C2','C3','C4']
def f3profile(data,title='title',what=1.0,color=None,ax=None):
    set_newplot = False
    if ax is None:
        set_newplot = True
        plt.figure()
        ax = plt.gca()
    for mu,color in zip(mus,colors):
        Zx, Za, Zw, XKa, KaW, LDR, IWC, MRR = ([] for i in range(8))
        for D0 in D0s:
            conc  = myGammaPSD(D0,mu)
            #conc = Gauss(D0,mu)
            iwc = IntPsd(data.Dmax,data.mkg,conc)
            if type(what) is float:
                iwc = iwc/what
            IWC.append(iwc)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax,data.X,conc)/iwc)
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax,data.Ku,conc)/iwc)
            if 'K' in data.columns:
                mrr = 10.*np.log10(coeffk*IntPsd(data.Dmax,data.K,conc)/iwc)
            else:
                mrr = np.nan
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            MRR.append(mrr)
            zalin = IntPsd(data.Dmax,data.Ka,conc)
            ldrka = IntPsd(data.Dmax,data.xpolKa,conc)/zalin
            Za.append(10.*np.log10(coeffa*zalin/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax,data.W,conc)/iwc))
            LDR.append(10.*np.log10(ldrka))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)
        if type(what) is float:
            ax.plot(Zx,D0s,c=color,linestyle='-.')
            ax.plot(Za,D0s,c=color,linestyle='-',label='$\mu$='+str(mu))
            ax.plot(Zw,D0s,c=color,linestyle=':')
            ax.plot(MRR,D0s,c=color,linestyle='--')
            ax.set_xlabel('Z    [dBZ]')
        elif what=='DWR':
            ax.plot(XKa,D0s,c=color,linestyle='-.')
            ax.plot(KaW,D0s,c=color,linestyle='-',label='$\mu$='+str(mu))
            ax.set_xlabel('DWR     [dBZ]')
        elif what=='LDR':
            ax.plot(LDR,D0s,c=color,linestyle='-.',label='$\mu$='+str(mu))
            ax.set_xlim([-40,-15])
            ax.set_xlabel('LDR     [dB]')
        ax.invert_yaxis()
        ax.set_ylabel('D$_0$')
        ax.grid()
        ax.legend()
        ax.set_title(title)
    return ax


def f3plot(data,title='title',color=None,ax=None):
    set_colorbar = False
    if ax is None:
        set_colorbar = True
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
            zalin = IntPsd(data.Dmax,data.Ka,conc)
            ldrka = IntPsd(data.Dmax,data.xpolKa,conc)/zalin
            Za.append(10.*np.log10(coeffa*zalin/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax,data.W,conc)/iwc))
            LDR.append(10.*np.log10(ldrka))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)

        if color is None:
            colorv = D0s
            vmin=0
            vmax=20
            cbarlabel = '$D_0$'
        elif color == 'ldr':
            colorv = LDR
            vmin = -35
            vmax = -15
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
                vmax = 30
            elif color == 'Zw':
                colorv = Zw
                print(max(Zw))
                vmin = 0
                vmax = 20

        s = ax.scatter(KaW,XKa,c=colorv,marker=m,label='$\mu=$ '+str(mu),vmin=vmin,vmax=vmax,cmap='jet')
    ax.grid()
    ax.set_xlabel('DWR$_{Ka,W}$')
    ax.set_ylabel('DWR$_{X,Ka}$')
    ax.set_xlim([0,20])
    ax.set_ylim([0,22])
    if set_colorbar:
        colorbar = plt.colorbar(mappable=s,ax=ax)
        colorbar.set_label(cbarlabel)
        ax.legend()
    ax.set_title(title)
    return ax

f3plot(dataDO,'Davide dry')
#f3plot(dataDO,'Davide dry',color='Zx')
#f3plot(dataBJ2,'BJ2 dry',color='Zx')

dataJL.columns = ['Dmax','model','ELWP','mkg','Dmax.1','Rgyr','ldr','riming', 'Xa','Xs', 'X', 'Xe', 'Ua', 'Us', 'Ku', 'Ue', 'Aa', 'As', 'Ka', 'Ae', 'Wa','Ws', 'W', 'We']
dataJL.X  = dataJL.X*1.0e6
dataJL.Ku = dataJL.Ku*1.0e6
dataJL.Ka = dataJL.Ka*1.0e6
dataJL.W  = dataJL.W*1.0e6
#f3plot(dataJL,'Jussi rimed A0.5',color='Zx')
#f3plot(dataJL,'Jussi rimed A0.5')

#f3plot(dataDO,'Davide dry',color='Zka')
#f3plot(dataJL,'Jussi rimed A0.5',color='Zka')

#f3plot(dataDO,'Davide dry',color='Zw')
#f3plot(dataJL,'Jussi rimed A0.5',color='Zw')

ax = f3plot(dataDO,'Davide dry',color='ldr')
f3plot(meltDO,'Davide 0 - 10 %',color='ldr',ax=ax)
ax.grid()

ax = f3plot(dataBJ2,'BJ2 dry',color='ldr')
f3plot(dataBJ2_10,'BJ2 10%',color='ldr',ax=ax)
f3plot(dataBJ2_20,'BJ2 20%',color='ldr',ax=ax)
f3plot(dataBJ2_30,'BJ2 30%',color='ldr',ax=ax)
f3plot(dataBJ2_40,'BJ2 40%',color='ldr',ax=ax)
f3plot(dataBJ2_49,'BJ2 49%',color='ldr',ax=ax)
f3plot(dataBJ2_70,'BJ2 0 - 70%',color='ldr',ax=ax)
ax.grid()


#f3profile(dataJL,title='Jussi unrimed')
#f3profile(dataJL,title='Jussi unrimed 0.01 kg/m$^2$',what=0.01)
#f3profile(dataJL,title='Jussi unrimed',what='DWR')

#f3profile(dataDO,title='Davide dry')
#f3profile(dataDO,title='Davide dry 0.01 kg/m$^2$',what=0.01)
#f3profile(dataDO,title='Davide dry',what='DWR')
f3profile(dataDO,title='Davide dry',what='LDR')
#f3profile(dataRH,title='RH spherical')
#f3profile(dataRH,title='RH spherical',what='DWR')
#f3plot(dataRH,'RH spherical')