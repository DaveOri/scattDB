# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:06:26 2017

@author: dori
"""

# from scattDB import shape
# from glob import glob
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from scattDB import scattering
from sys import path
path.append('/home/dori/develop/pyPamtra2/libs/singleScattering/singleScattering/')
from self_similar_rayleigh_gans import backscattering

try:
    from scattDB import psd
except:
    import scattDB.psd as psd
from scipy import integrate

Br07 = lambda x:0.08900e-3*x**2.1
BF95 = lambda x:0.02414e-3*x**1.9

RHfile = '/home/dori/develop/scatdb/share/scatdb.csv'
dataRH = pd.read_csv(RHfile)
dataRH = dataRH[dataRH.flaketype == 20]
dataRH.set_index('max_dimension_mm', inplace=True)
dataRH.sort_index(inplace=True)

tol = 0.3
dataRHx = dataRH[abs(dataRH.frequencyghz-10.6) < tol]
dataRHu = dataRH[abs(dataRH.frequencyghz-13.6) < tol]
dataRHk = dataRH[abs(dataRH.frequencyghz-23.8) < tol]
dataRHa = dataRH[abs(dataRH.frequencyghz-35.6) < tol]
dataRHw = dataRH[abs(dataRH.frequencyghz-94) < tol]

# dataRHx = dataRHx.loc[dataRHa.index]
# dataRHu = dataRHu.loc[dataRHa.index]
# dataRHk = dataRHk.loc[dataRHa.index]
# dataRHw = dataRHw.loc[dataRHa.index]

plt.close('all')

c = 299792458000.  # mm/s
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

tablesfolder = 'tables/'
authors = ['DO', 'BJ', 'JL', 'meltedDO']


author = 'meltedDO'
meltedDO10 = pd.read_csv(tablesfolder+author+'10'+'.csv')
meltedDO20 = pd.read_csv(tablesfolder+author+'20'+'.csv')
meltedDO30 = pd.read_csv(tablesfolder+author+'30'+'.csv')
meltedDO40 = pd.read_csv(tablesfolder+author+'40'+'.csv')
meltedDO50 = pd.read_csv(tablesfolder+author+'50'+'.csv')
meltedDO60 = pd.read_csv(tablesfolder+author+'60'+'.csv')
meltedDO70 = pd.read_csv(tablesfolder+author+'70'+'.csv')


author = 'DO'
melt = '00'
dataDO = pd.read_csv(tablesfolder+author+'_'+melt+'.csv')
# dataDO.drop('Unnamed: 0',axis=1,inplace=True)
meltDO = pd.read_csv(tablesfolder+author+'_'+'10'+'.csv')
# meltDO.drop('Unnamed: 0',axis=1,inplace=True)


dataDOssrg = pd.read_csv(tablesfolder+author+'_'+'SSRG.csv')
# dataDOssrg.drop('Unnamed: 0',axis=1,inplace=True)

dataJT = pd.read_csv(tablesfolder+'jani_fractals/JT_fractals.csv')
dataJT.drop('Unnamed: 0', axis=1, inplace=True)
dataJT.mkg = Br07(dataJT.Dmax)

author = 'BJ'
dataBJ2 = pd.read_csv(tablesfolder+author+'_agg2_0'+'.csv')
dataBJ2.drop('Unnamed: 0', axis=1, inplace=True)

dataBJ2_10 = pd.read_csv(tablesfolder+author+'_agg2_10'+'.csv')
dataBJ2_10.drop('Unnamed: 0', axis=1, inplace=True)
dataBJ2_20 = pd.read_csv(tablesfolder+author+'_agg2_20'+'.csv')
dataBJ2_20.drop('Unnamed: 0', axis=1, inplace=True)
dataBJ2_30 = pd.read_csv(tablesfolder+author+'_agg2_30'+'.csv')
dataBJ2_30.drop('Unnamed: 0', axis=1, inplace=True)
dataBJ2_40 = pd.read_csv(tablesfolder+author+'_agg2_40'+'.csv')
dataBJ2_40.drop('Unnamed: 0', axis=1, inplace=True)
dataBJ2_49 = pd.read_csv(tablesfolder+author+'_agg2_49'+'.csv')
dataBJ2_49.drop('Unnamed: 0', axis=1, inplace=True)
dataBJ2_70 = pd.read_csv(tablesfolder+author+'_agg2_70'+'.csv')
dataBJ2_70.drop('Unnamed: 0', axis=1, inplace=True)

dataBJ3 = pd.read_csv(tablesfolder+author+'_agg3_0'+'.csv')
dataBJ3.drop('Unnamed: 0', axis=1, inplace=True)

#author = 'dataJL_B1.0'
author = 'dataJL_A0.1'
author = 'dataJL_A0.0'
dataJL = pd.read_csv(tablesfolder+author+'.csv')
dataJL['xpolKa'] = 0.0

author = 'dataJL_B1.0'
dataRimed = pd.read_csv(tablesfolder+author+'.csv')
dataRimed['xpolKa'] = 0.0

#bm, am = np.polyfit(np.log10(dataJL.Dmax),np.log10(dataJL.mkg),1)
#am = 10.0**am
bm = 2.08
am = 0.015

dataRH = pd.DataFrame(index=dataRHa.index, columns=dataDO.columns)
dataRH.Dmax = dataRH.index
dataRH.mkg = Br07(dataRH.Dmax)
dataRH.X = dataRHx.cbk
dataRH.Ku = dataRHu.cbk
dataRH['K'] = dataRHk.cbk
dataRH.Ka = dataRHa.cbk
dataRH.W = dataRHw.cbk

plt.figure()
ax = plt.gca()
# ax.scatter(dataBJ2.Dmax,dataBJ2.mkg,label='BJ2')
# ax.scatter(dataBJ3.Dmax,dataBJ3.mkg,label='BJ3')
ax.scatter(dataJT.Dmax, dataJT.mkg, label='JT')
ax.scatter(dataRH.Dmax, dataRH.mkg, label='RH')
ax.scatter(dataJL.Dmax, dataJL.mkg, label='JL')
ax.scatter(dataDO.Dmax, dataDO.mkg, label='DO')
ax.plot(dataDO.Dmax, Br07(dataDO.Dmax), label='brandes')
ax.plot(dataDO.Dmax, am*(dataDO.Dmax)**bm, label='JLfit')
ax.plot(dataDO.Dmax, BF95(dataDO.Dmax), label='BF95')
ax.legend()
ax.grid()
ax.set_yscale('log')
ax.set_xscale('log')


# %%

def IntPsd(D, s, psd):
    conc = psd(D)
    intZ = np.multiply(conc, s)
    return integrate.trapz(intZ, D)


def myGammaPSD(D0, mu, D_max=100.0):
    return psd.GammaPSD(D0=D0, Nw=1., mu=mu, D_max=D_max)


def Gauss(mu, sig):
    def F(D):
        x = D-mu
        return np.exp(x*x/(2.0*sig*sig))/(sig*np.sqrt(2.0*np.pi))
    return F

D0s = np.linspace(0.1, 20., 30)
mus = [-1., 0., 1., 2., 3.]
mus = [0.0]
marks = [',', '+', 'h', 'v', '.']
colors = ['C0', 'C1', 'C2', 'C3', 'C4']


def f3profile(data, title='title', what=1.0, color=None, ax=None):
    # set_newplot = False
    if ax is None:
        # set_newplot = True
        plt.figure()
        ax = plt.gca()
    for mu, color in zip(mus, colors):
        Zx, Za, Zw, XKa, KaW, LDR, IWC, MRR = ([] for i in range(8))
        for D0 in D0s:
            conc = myGammaPSD(D0, mu)
            # conc = Gauss(D0,mu)
            iwc = IntPsd(data.Dmax, data.mkg, conc)
            if type(what) is float:
                iwc = iwc/what
            IWC.append(iwc)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax, data.X, conc)/iwc)
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax, data.Ku, conc)/iwc)
            if 'K' in data.columns:
                mrr = 10.*np.log10(coeffk*IntPsd(data.Dmax, data.K, conc)/iwc)
            else:
                mrr = np.nan
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            MRR.append(mrr)
            zalin = IntPsd(data.Dmax, data.Ka, conc)
            ldrka = IntPsd(data.Dmax, data.xpolKa, conc)/zalin
            Za.append(10.*np.log10(coeffa*zalin/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax, data.W, conc)/iwc))
            LDR.append(10.*np.log10(ldrka))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)
        if type(what) is float:
            ax.plot(Zx, D0s, c=color, linestyle='-.')
            ax.plot(Za, D0s, c=color, linestyle='-', label='$\mu$='+str(mu))
            ax.plot(Zw, D0s, c=color, linestyle=':')
            ax.plot(MRR, D0s, c=color, linestyle='--')
            ax.set_xlabel('Z    [dBZ]')
        elif what == 'DWR':
            ax.plot(XKa, D0s, c=color, linestyle='-.')
            ax.plot(KaW, D0s, c=color, linestyle='-', label='$\mu$='+str(mu))
            ax.set_xlabel('DWR     [dBZ]')
        elif what == 'LDR':
            ax.plot(LDR, D0s, c=color, linestyle='-.', label='$\mu$='+str(mu))
            # ax.set_xlim([-50,-15])
            ax.set_xlabel('LDR     [dB]')
        ax.invert_yaxis()
        ax.set_ylabel('D$_0$')
        ax.grid()
        ax.legend()
        ax.set_title(title)
    return ax


def f3plot(data, title='title', color=None, ax=None, fig=None):
    set_colorbar = False
    if ax is None:
        set_colorbar = True
        fig = plt.figure()
        ax = plt.gca()
    for mu, m in zip(mus, marks):
        Zx, Za, Zw, XKa, KaW, LDR, IWC = ([] for i in range(7))
        for D0 in D0s:
            conc = myGammaPSD(D0, mu)
            # conc = Gauss(D0,mu)
            iwc = IntPsd(data.Dmax, data.mkg, conc)
            IWC.append(iwc)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax, data.X, conc)/iwc)
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax, data.Ku, conc)/iwc)
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            zalin = IntPsd(data.Dmax, data.Ka, conc)
            ldrka = IntPsd(data.Dmax, data.xpolKa, conc)/zalin
            Za.append(10.*np.log10(coeffa*zalin/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax, data.W, conc)/iwc))
            LDR.append(10.*np.log10(ldrka))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)
        print(max(LDR))
        if color is None:
            colorv = D0s
            vmin = 0
            vmax = 20
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

        s = ax.scatter(KaW, XKa, c=colorv, marker=m, label='$\mu=$ '+str(mu),
                       vmin=vmin, vmax=vmax, cmap='jet')
    ax.grid()
    ax.set_xlabel('DWR$_{Ka,W}$')
    ax.set_ylabel('DWR$_{X,Ka}$')
    ax.set_xlim([0, 20])
    ax.set_ylim([0, 20])
    if set_colorbar:
        colorbar = plt.colorbar(mappable=s, ax=ax)
        colorbar.set_label(cbarlabel)
        ax.legend()
    ax.set_title(title)
    return fig, ax

# f3plot(dataDOssrg,'Dave SSRG')
# f3plot(dataDO,'Davide dry')
# f3plot(dataDO,'Davide dry',color='Zx')
# f3plot(dataBJ2,'BJ2 dry',color='Zx')

# dataJL.columns = ['Dmax','model','ELWP','mkg','Dmax.1','Rgyr','ldr','riming', 'Xa','Xs', 'X', 'Xe', 'Ua', 'Us', 'Ku', 'Ue', 'Aa', 'As', 'Ka', 'Ae', 'Wa','Ws', 'W', 'We','xpolKa']
# dataJL.X  = dataJL.X*1.0e6
# dataJL.Ku = dataJL.Ku*1.0e6
# dataJL.Ka = dataJL.Ka*1.0e6
# dataJL.W  = dataJL.W*1.0e6
# f3plot(dataJL,'Jussi rimed A0.5',color='Zx')
# f3plot(dataJL,'Jussi rimed A0.5')

# f3plot(dataDO,'Davide dry',color='Zka')
# f3plot(dataJL,'Jussi rimed A0.5',color='Zka')

# f3plot(dataDO,'Davide dry',color='Zw')
# f3plot(dataJL,'Jussi rimed A0.5',color='Zw')

fig, ax = f3plot(dataDO, color='ldr')
f3plot(meltedDO30, color='ldr', ax=ax, fig=fig)
f3plot(meltedDO50, color='ldr', ax=ax, fig=fig)
f3plot(meltedDO70, color='ldr', title='DO melt = 0 30 50 70', ax=ax, fig=fig)
plt.grid()
fig.savefig('3f_melted_Davide.png', dpi=300)

fig, ax = f3plot(dataBJ2, color='ldr')
f3plot(dataBJ2_30, color='ldr', ax=ax, fig=fig)
f3plot(dataBJ2_49, color='ldr', ax=ax, fig=fig)
f3plot(dataBJ2_70, color='ldr', title='Ben2 melt = 0 30 50 70', ax=ax, fig=fig)
plt.grid()
fig.savefig('3f_melted_Ben.png', dpi=300)

# f3plot(meltDO,'Davide 0 - 10 %',color='ldr',ax=ax)
# ax.grid()

# f3plot(dataJT,'JT fractals')
# f3plot(dataJT,'JT fractals',color='ldr')
# f3profile(dataJT,title='JT fractals',what='LDR')
# f3profile(dataJT,title='JT fractals')
# f3profile(dataJT,title='JT fractals',what='DWR')

# ax = f3plot(dataBJ2,'BJ2 dry',color='ldr')
# f3plot(dataBJ2_10,'BJ2 10%',color='ldr',ax=ax)
# f3plot(dataBJ2_20,'BJ2 20%',color='ldr',ax=ax)
# f3plot(dataBJ2_30,'BJ2 30%',color='ldr',ax=ax)
# f3plot(dataBJ2_40,'BJ2 40%',color='ldr',ax=ax)
# f3plot(dataBJ2_49,'BJ2 49%',color='ldr',ax=ax)
# f3plot(dataBJ2_70,'BJ2 0 - 70%',color='ldr',ax=ax)
# ax.grid()

# f3profile(dataJL,title='Jussi unrimed')
# f3profile(dataJL,title='Jussi unrimed 0.01 kg/m$^2$',what=0.01)
# f3profile(dataJL,title='Jussi unrimed',what='DWR')

ax = f3profile(dataBJ2, what='LDR')
f3profile(dataBJ2_30, ax=ax, what='LDR')
f3profile(dataBJ2_49, ax=ax, what='LDR')
f3profile(dataBJ2_70, title='Ben J 0 30 50 70', ax=ax, what='LDR')
ax.invert_yaxis()
plt.grid()
plt.savefig('ldr_profile_BJ.png', dpi=400)

ax = f3profile(dataDO, what='LDR')
f3profile(meltedDO30, ax=ax, what='LDR')
f3profile(meltedDO50, ax=ax, what='LDR')
f3profile(meltedDO70, title='Davide 0 30 50 70', ax=ax, what='LDR')
ax.invert_yaxis()
plt.grid()
plt.savefig('ldr_profile_Davide.png', dpi=400)
# f3profile(dataDO,title='Davide dry 0.01 kg/m$^2$',what=0.01)
# f3profile(dataDO,title='Davide dry',what='DWR')
# f3profile(dataDO,title='Davide dry',what='LDR')
# f3profile(dataDOssrg,title='Dave SSRG 0.01 kg/m$^2$',what=0.01)
# f3profile(dataDOssrg,title='Dave SSRG',what='DWR')

# f3profile(dataRH,title='RH spherical')
# f3profile(dataRH,title='RH spherical',what='DWR')
# f3plot(dataRH,'RH spherical')

D0s = [10., 15., 20.]
marks = ['--', '-', '.']
mus = [0.0, 4.0]
Dmax = np.linspace(10, 50, 100)
data = 1.0*dataDOssrg
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)

for D0, m in zip(D0s, marks):
    for mu in mus:
        Zx, Za, Zw, XKa, KaW, LDR, IWC, MRR = ([] for i in range(8))
        for Dx in Dmax:
            conc = myGammaPSD(D0, mu, D_max=Dx)
            iwc = IntPsd(data.Dmax, data.mkg, conc)
            IWC.append(iwc)
            zx = 10.*np.log10(coeffx*IntPsd(data.Dmax, data.X, conc)/iwc)
            zu = 10.*np.log10(coeffu*IntPsd(data.Dmax, data.Ku, conc)/iwc)
            if np.isnan(zx):
                zx = zu
            Zx.append(zx)
            zalin = IntPsd(data.Dmax, data.Ka, conc)
            ldrka = IntPsd(data.Dmax, data.xpolKa, conc)/zalin
            Za.append(10.*np.log10(coeffa*zalin/iwc))
            Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax, data.W, conc)/iwc))
            LDR.append(10.*np.log10(ldrka))
        XKa = np.array(Zx)-np.array(Za)
        KaW = np.array(Za)-np.array(Zw)
        ax1.plot(Dmax, Zx, m, label='D$_0$='+str(D0)+' $\mu$='+str(mu))
        ax2.plot(Dmax, Za, m, label='D$_0$='+str(D0)+' $\mu$='+str(mu))
        ax3.plot(Dmax, Zw, m, label='D$_0$='+str(D0)+' $\mu$='+str(mu))
        ax4.plot(Dmax, KaW, m, label='D$_0$='+str(D0)+' $\mu$='+str(mu))
ax1.legend()
ax1.grid()
ax1.set_ylabel('Zx/iwc')
ax2.grid()
ax2.set_ylabel('Za/iwc')
ax3.grid()
ax3.set_ylabel('Zw/iwc')
ax3.set_xlabel('Truncation size')
ax4.grid()
ax4.set_ylabel('DWR Ka-W')
ax4.set_xlabel('Truncation size')


# %%
D0s = np.linspace(0.1, 20., 30)
mus = [-1., 0., 1., 2., 3.]
mus = [0.0]
marks = [',', '+', 'h', 'v', '.']
colors = ['C0', 'C1', 'C2', 'C3', 'C4']
colorv = D0s
vmin = 0
vmax = 20
cbarlabel = '$D_0$'

mu = 0.0


def integratePSD(data):
    Zx, Za, Zw, XKa, KaW, LDR, IWC = ([] for i in range(7))
    for D0 in D0s:
        conc = myGammaPSD(D0, mu)
        iwc = IntPsd(data.Dmax, data.mkg, conc)
        IWC.append(iwc)
        zx = 10.*np.log10(coeffx*IntPsd(data.Dmax, data.X, conc)/iwc)
        zu = 10.*np.log10(coeffu*IntPsd(data.Dmax, data.Ku, conc)/iwc)
        if np.isnan(zx):
            zx = zu
        Zx.append(zx)
        zalin = IntPsd(data.Dmax, data.Ka, conc)
        ldrka = IntPsd(data.Dmax, data.xpolKa, conc)/zalin
        Za.append(10.*np.log10(coeffa*zalin/iwc))
        Zw.append(10.*np.log10(coeffw*IntPsd(data.Dmax, data.W, conc)/iwc))
        LDR.append(10.*np.log10(ldrka))
    XKa = np.array(Zx)-np.array(Za)
    KaW = np.array(Za)-np.array(Zw)
    return Zx, Za, Zw, XKa, KaW, LDR, IWC

fig, ax = plt.subplots(1, 2, sharey=True, figsize=(10, 4))
fig2, ax2 = plt.subplots(1, 1, figsize=(6, 4))
fig3, ax3 = plt.subplots(1, 1, figsize=(6, 4))

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataDO)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker=',', label='DO dry', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker=',', label='DO dry', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='D0 dry')
ax3.scatter(KaW, XKa, c=LDR, marker=',', label='D0 dry', cmap='jet',
            vmin=-35, vmax=-15)

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(meltedDO40)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker='+', label='DO 40%', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker='+', label='DO 40%', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='D0 40%')
ax3.scatter(KaW, XKa, c=LDR, marker='+', label='D0 40%', cmap='jet',
            vmin=-35, vmax=-15)

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(meltedDO70)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker='h', label='DO 70%', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker='h', label='DO 70%', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='D0 70%')
ax3.scatter(KaW, XKa, c=LDR, marker='h', label='D0 70%', cmap='jet',
            vmin=-35, vmax=-15)

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataBJ2)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker='.', label='BJ dry', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker='.', label='BJ dry', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='BJ dry')
ax3.scatter(KaW, XKa, c=LDR, marker='.', label='BJ dry', cmap='jet',
            vmin=-35, vmax=-15)

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataBJ2_40)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker='x', label='BJ 40%', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker='x', label='BJ 40%', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='BJ 40%')
ax3.scatter(KaW, XKa, c=LDR, marker='x', label='BJ 40%', cmap='jet',
            vmin=-35, vmax=-15)

Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataBJ2_70)
print(max(LDR))
s = ax[0].scatter(XKa, LDR, c=D0s, marker='v', label='BJ 70%', vmin=vmin,
                  vmax=vmax, cmap='jet')
s = ax[1].scatter(KaW, LDR, c=D0s, marker='v', label='BJ 70%', vmin=vmin,
                  vmax=vmax, cmap='jet')
ax2.plot(LDR, D0s, label='BJ 70%')
ss = ax3.scatter(KaW, XKa, c=LDR, marker='v', label='BJ 70%', cmap='jet',
                 vmin=-35, vmax=-15)

ax[0].grid()
ax[1].grid()
ax[0].legend()
#colorbar = plt.colorbar(mappable=s, ax=ax[0])
#colorbar.set_label(cbarlabel)
colorbar = plt.colorbar(mappable=s, ax=ax[1])
colorbar.set_label(cbarlabel)
colorbar = plt.colorbar(mappable=ss, ax=ax3)

ax[1].set_xlabel('DWR$_{Ka,W}$   [dBZ]')
ax[0].set_xlabel('DWR$_{X,Ka}$   [dBZ]')
ax[0].set_ylabel('LDR    [dB]')
ax[0].set_xlim([0, 25])
ax[1].set_xlim([0, 20])
fig.suptitle('DWR evolution in the melting layer')
#fig.tight_layout()
fig.savefig('LDR_DWR.png',dpi=300)

ax2.set_ylim([0.0, 16])
ax2.invert_yaxis()
ax2.legend()
ax2.grid()
ax2.set_ylabel('D$_0$   [mm]')
ax2.set_xlabel('LDR   [dB]')
ax2.set_title('LDR profile')
fig2.savefig('LDRprofile.pdf')
fig2.savefig('LDRprofile.png',dpi=300)

ax3.set_ylim([0.0, 20])
ax3.legend()
ax3.grid()
ax3.set_ylabel('DWR$_{XKa}$   [dB]')
ax3.set_xlabel('DWR$_{KaW}$   [dB]')
ax3.set_title('3-frequency - LDR')
fig3.savefig('3fLDR.pdf')
fig3.savefig('3fLDR.png',dpi=300)

# %%
fig, ax = plt.subplots(1, 3, sharey=True, figsize=(14, 4))
ax[0].plot(dataDO['Dmax'], dataDO['xpolKa'], label='DO_00', linestyle='-.')
ax[0].plot(meltedDO30['Dmax'], meltedDO30['xpolKa'], label='DO_30',
           linestyle='-.')
ax[0].plot(meltedDO70['Dmax'], meltedDO70['xpolKa'], label='DO_70',
           linestyle='-.')
ax[0].plot(dataBJ2['Dmax'], dataBJ2['xpolKa'], label='BJ_00')
ax[0].plot(dataBJ2_30['Dmax'], dataBJ2_30['xpolKa'], label='BJ_30')
ax[0].plot(dataBJ2_70['Dmax'], dataBJ2_70['xpolKa'], label='BJ_70')
# ax[0].legend()
ax[0].set_title('Xpol Xsec Ka   [mm$^2$]')
ax[0].set_yscale('log')
ax[0].set_xlabel('Dmax     [mm]')

ax[1].plot(dataDO['Dmax'], dataDO['CscaKa'], label='DO_00', linestyle='-.')
ax[1].plot(meltedDO30['Dmax'], meltedDO30['CscaKa'], label='DO_30',
           linestyle='-.')
ax[1].plot(meltedDO70['Dmax'], meltedDO70['CscaKa'], label='DO_70',
           linestyle='-.')
ax[1].plot(dataBJ2['Dmax'], dataBJ2['CscaKa'], label='BJ_00')
ax[1].plot(dataBJ2_30['Dmax'], dataBJ2_30['CscaKa'], label='BJ_30')
ax[1].plot(dataBJ2_70['Dmax'], dataBJ2_70['CscaKa'], label='BJ_70')

ax[2].plot(dataDO['Dmax'], dataDO['Ka'], label='DO_00', linestyle='-.')
ax[2].plot(meltedDO30['Dmax'], meltedDO30['Ka'], label='DO_30', linestyle='-.')
ax[2].plot(meltedDO70['Dmax'], meltedDO70['Ka'], label='DO_70', linestyle='-.')
ax[2].plot(dataBJ2['Dmax'], dataBJ2['Ka'], label='BJ_00')
ax[2].plot(dataBJ2_30['Dmax'], dataBJ2_30['Ka'], label='BJ_30')
ax[2].plot(dataBJ2_70['Dmax'], dataBJ2_70['Ka'], label='BJ_70')
ax[2].set_title('Copol Xsec Ka   [mm$^2$]')
ax[2].set_xlabel('Dmax     [mm]')
ax[2].legend()
ax[2].set_yscale('log')
fig.savefig('DO_BJ_comp_Xpol.png')
# %%
plt.close('all')
reff = lambda mass: np.cbrt(3.0*mass/(4.0*np.pi*917.0))
xeff = lambda r, l: 2.0*np.pi*r/l
qeff = lambda C, r: C/(np.pi*r**2.0)


def mass2x(mass, wl):
    return 2.0*np.pi*np.cbrt(3.0*mass/(4.0*np.pi*917.0))/wl

def reff2rhoBr07(r, ar):
    mass = 4.0*np.pi*r**3.0*916.0/3.0
    size = 0.001*(mass*1000.0/8.9e-5)**(1.0/2.1)
#    print(mass,size)
    vol = np.pi*size**3.0*ar/6.0
    den = mass/vol
    for i, d in enumerate(den):
        if d > 916.0:
            den[i] = 916.0
            size[i] = np.cbrt(6.0*(mass[i]/den[i])/(np.pi*ar))
    return size, den


def reff2rhoBF95(r, ar):
    mass = 4.0*np.pi*r**3.0*916.0/3.0
    size = 0.001*(mass*1000.0/2.4e-5)**(1.0/1.9)
#    print(mass)
#    print(size)
    vol = np.pi*size**3.0*ar/6.0
    den = mass/vol
    for i, d in enumerate(den):
        if d > 916.0:
            den[i] = 916.0
            size[i] = np.cbrt(6.0*(mass[i]/den[i])/(np.pi*ar))
    return size, den


sys.path.append('/work/DBs/scattnlay')
try:
    from scattnlay import scattnlay
    # from scattnlay import fieldnlay
    from pytmatrix import tmatrix, radar, scatter
    import refractiveIndex
except:
    print('not python2')

xmie = 10**np.linspace(-2, 0.8, 400)
xmie = xmie.reshape([len(xmie), 1])
frequency = 94.6e9
wavelength = 299792458.0/frequency
r = scattnlay(xmie, np.ones_like(xmie)*refractiveIndex.ice.n(270.0, frequency))
terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = r

rsoft = xmie*wavelength/(2.0*np.pi)
size, rho = reff2rhoBF95(rsoft, 1.0)
ref_soft = refractiveIndex.snow.n(270.0, frequency, rho)
xsoft = xeff(size*0.5, wavelength)
terms, Qes, Qss, Qas, Qbs, Qpr, gs, Albedo, S1, S2 = scattnlay(xsoft, ref_soft)
corr = ((0.5*size/rsoft)**2.0).reshape(Qas.shape)
Qes = Qes*corr
Qss = Qss*corr
Qas = Qas*corr
Qbs = Qbs*corr

# dmaxes = rsoft*2.0
# masses = 1.0e-3*dataJL.mkg.sort_values().values  # am*(dmaxes*1000.0)**bm
masses = 8.0*10**np.linspace(-13, -5, 1000)
dmaxes = (masses/am)**(1.0/bm)
mm = refractiveIndex.ice.n(270., frequency)
Qbssrg = 0.0*masses
Cbssrg = 0.0*masses
Cassrg = 0.0*masses
Csssrg = 0.0*masses
Qessrg = 0.0*masses
xssrg = xeff(reff(masses), wavelength)
for i, [d, m] in enumerate(zip(dmaxes, masses)):
    Cbssrg[i], Cassrg[i], Csssrg[i], dum, mk = backscattering(frequency, d, mm,
                                                              table='leinonen',
                                                              mass=m)

    Qessrg[i] = (Cassrg[i]+Csssrg[i])/(np.pi*reff(m)**2.0)
    Qessrg[i] = qeff(Cassrg[i]+Csssrg[i], reff(m))
    Qbssrg[i] = Cbssrg[i]/(np.pi*reff(m)**2.0)

limit = 399
axr = 0.6
sizeell, rhoell = reff2rhoBF95(rsoft[:limit], axr)
ref_ell = refractiveIndex.snow.n(270.0, frequency, rhoell)
Qee = 0.0*sizeell
Qse = 0.0*sizeell
Qae = 0.0*sizeell
Qbe = 0.0*sizeell
ge = 0.0*sizeell
for i, [rr, mm] in enumerate(zip(0.5*sizeell, ref_ell)):
    print(i)
    scatt = tmatrix.Scatterer(radius=rr,
                              radius_type=tmatrix.Scatterer.RADIUS_MAXIMUM,
                              wavelength=wavelength, m=mm, axis_ratio=1.0/axr)
    scatt.set_geometry((0.0, 180.0,   0.0, 0.0, 0.0, 0.0))
    Qee[i] = scatter.ext_xsect(scatt)
    # Qse[i] = scatter.sca_xsect(scatt)
    # ge[i] = scatter.asym(scatt)
    scatt.set_geometry((0.0, 180.0, 0.0, 0.0, 0.0, 0.0))
    Qbe[i] = radar.radar_xsect(scatt)

Aeff = np.pi*rsoft[:limit]**2.0
Qee = Qee/Aeff
Qse = Qse/Aeff
Qbe = Qbe/Aeff
Qae = Qee-Qse

plt.figure()
plt.plot(1000.0*dmaxes, Cbssrg, c='r')
plt.scatter(dataJL['Dmax'], dataJL['Wb'])
plt.yscale('log')

plt.figure()
plt.plot(1000.0*dmaxes, Cassrg, c='r')
plt.scatter(dataJL['Dmax'], dataJL['Wa'])
plt.yscale('log')

plt.figure()
plt.plot(1000.0*dmaxes, Csssrg, c='r')
plt.scatter(dataJL['Dmax'], dataJL['Ws'])
plt.yscale('log')

plt.figure()
plt.plot(1000.0*masses, Csssrg, c='r')
plt.scatter(dataJL['mkg'], dataJL['Ws'])
plt.yscale('log')

# %%

fig, ax = plt.subplots(1, 2, figsize=(9, 4), sharex=True)

aeff = reff(dataJL['mkg']*0.001)
ax[0].scatter(xeff(aeff, lamw*0.001), qeff(dataJL['Wb'], aeff), c='0.5', s=1)

aeff = reff(dataRimed['mkg']*0.001)
ax[0].scatter(xeff(aeff, lamw*0.001), qeff(dataRimed['Wb'], aeff), c='k', s=1)

# aeff = reff(1.0e-3*dataDOssrg['mkg'])
# ax[0].plot(xeff(aeff, lamw*0.001), qeff(1.0e-6*dataDOssrg['W'], aeff),
#               c='g')
# ax[0, 0].plot(xeff(aeff, lamx*0.001), qeff(1.0*dataDOssrg['X'], aeff),
#              c='r')
ax[0].plot(xssrg, Qbssrg, c='0.5', lw=3.0)

ax[0].plot(xmie, Qbk, c='k')
ax[0].plot(xmie, Qbs, '--', c='k')
ax[0].plot(xmie[:limit], Qbe, ':', c='k')

mg = 10.0e-6
ax[0].vlines(mass2x(mg, lamx*1.0e-3), ymin=1e-8, ymax=1e2)
ax[0].vlines(mass2x(mg, lama*1.0e-3), ymin=1e-8, ymax=1e2)
ax[0].vlines(mass2x(mg, lamw*1.0e-3), ymin=1e-8, ymax=1e2)
ax[0].text(2.0e-1, 3.0e0, 'X')
ax[0].text(7.0e-1, 3.0e0, 'Ka')
ax[0].text(2.0e0, 3.0e0, 'W')

ax[0].set_xlim([1e-1, 1e1])
ax[0].set_ylim([1e-8, 1e2])
ax[0].set_xscale('log')
ax[0].set_yscale('log')

###############################################################################
aeff = reff(dataJL['mkg']*0.001)
ax[1].scatter(xeff(aeff, lamw*0.001),
              qeff(dataJL['Ws']+dataJL['Wa'], aeff),
              c='0.5', s=1, label='Unrimed aggregates')

aeff = reff(dataRimed['mkg']*0.001)
ax[1].scatter(xeff(aeff, lamw*0.001),
              qeff(dataRimed['Ws']+dataRimed['Wa'], aeff),
              c='k', s=1, label='Rimed aggregates')

ax[1].plot(xmie, Qsca+Qabs, c='k', label='solid-sphere')
ax[1].plot(xmie, Qss+Qas, '--', c='k', label='soft-sphere')
ax[1].plot(xmie[:limit], Qee, ':', c='k', label='soft-spheroid $a_r=0.6$')
ax[1].plot(xssrg, Qessrg, c='0.5', lw=3.0, label='ssrg unrimed')

ax[1].vlines(mass2x(mg, lamx*1.0e-3), ymin=1e-4, ymax=1e1)
ax[1].vlines(mass2x(mg, lama*1.0e-3), ymin=1e-4, ymax=1e1)
ax[1].vlines(mass2x(mg, lamw*1.0e-3), ymin=1e-4, ymax=1e1)
ax[1].text(2.0e-1, 5.0e0, 'X')
ax[1].text(7.0e-1, 5.0e0, 'Ka')
ax[1].text(2.0e0, 5.0e0, 'W')

ax[1].set_xlim([1e-1, 1e1])
ax[1].set_ylim([1e-4, 1e1])
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].legend(loc=4)

ax[0].set_xlabel('$x=2\pi r_{eff}/\lambda$')
ax[1].set_xlabel('$x=2\pi r_{eff}/\lambda$')
ax[0].set_ylabel('$Q_{bk}$')
ax[1].set_ylabel('$Q_{ext}$')
fig.tight_layout()
fig.savefig('aggregates_scattering.png', dpi=300)
fig.savefig('aggregates_scattering.pdf', dpi=300)

# %%
sizes = np.linspace(0.00001, 0.04, 1000)
masses = am*sizes**bm
fx = 9.6e9
mx = refractiveIndex.ice.n(270., fx)
fku = 13.6e9
mku = refractiveIndex.ice.n(270., fx)
fka = 35.6e9
mka = refractiveIndex.ice.n(270., fka)
fw = 94e9
mw = refractiveIndex.ice.n(270., fw)
Cx = 0.0*masses
Cku = 0.0*masses
Cka = 0.0*masses
Cw = 0.0*masses
for i, [d, m] in enumerate(zip(sizes, masses)):
    Cx[i], dum, dum, dum, dum = backscattering(fx, d, mx,
                                               table='leinonen', mass=m)
    Cku[i], dum, dum, dum, dum = backscattering(fku, d, mku,
                                                table='leinonen', mass=m)
    Cka[i], dum, dum, dum, dum = backscattering(fka, d, mka,
                                                table='leinonen', mass=m)
    Cw[i], dum, dum, dum, dum = backscattering(fw, d, mw,
                                               table='leinonen', mass=m)
cols = ['Dmax', 'X', 'Ku', 'Ka', 'W', 'mkg', 'xpolKa']
dataJLssrg = pd.DataFrame(index=sizes, columns=cols)
dataJLssrg['Dmax'] = sizes*1000.0
dataJLssrg['mkg'] = masses*1000.0
dataJLssrg['X'] = Cx
dataJLssrg['Ku'] = Cku
dataJLssrg['Ka'] = Cka
dataJLssrg['W'] = Cw

D0s = np.linspace(0.1, 20., 30)
vmin = 0
vmax = 20
cbarlabel = '$D_0$'

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

mu = 0.0
Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataJLssrg)
s = ax.scatter(KaW, XKa, c=D0s, marker=',',
               label='Leinonen $\mu = 0$', vmin=vmin,
               vmax=vmax, cmap='magma')
mu = 4.0
Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataJLssrg)
s = ax.scatter(KaW, XKa, c=D0s, marker='+',
               label='Leinonen $\mu = 4$', vmin=vmin,
               vmax=vmax, cmap='magma')
mu = 0.0
Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataDOssrg)
s = ax.scatter(KaW, XKa, c=D0s, marker='h',
               label='Ori $\mu = 0$', vmin=vmin,
               vmax=vmax, cmap='magma')
mu = 4.0
Zx, Za, Zw, XKa, KaW, LDR, IWC = integratePSD(dataDOssrg)
s = ax.scatter(KaW, XKa, c=D0s, marker='v',
               label='Ori $\mu = 4$', vmin=vmin,
               vmax=vmax, cmap='magma')

ax.grid()
ax.legend()
colorbar = plt.colorbar(mappable=s, ax=ax)
colorbar.set_label(cbarlabel)

ax.set_xlabel('DWR$_{Ka,W}$   [dB]')
ax.set_ylabel('DWR$_{X,Ka}$   [dB]')
# ax.set_xlim([0, 25])

fig.suptitle('Unrimed aggregates')
fig.savefig('3f_plot_AMS.png',dpi=600)

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.plot(1000.*sizes, 10.0*np.log10(10.0e6*Cx*lamx**4/np.pi**5), label="X-band")
ax.plot(1000.*sizes, 10.0*np.log10(10.0e6*Cka*lama**4/np.pi**5), label="Ka-band")
ax.plot(1000.*sizes, 10.0*np.log10(10.0e6*Cw*lamw**4/np.pi**5), label="W-band")
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim([0.1, 100])
ax.set_ylim([-50, 30])
ax.set_xlabel('Particle size    [mm]')
ax.set_ylabel('Reflectivity    [dBZ]')
ax.set_title('Unrimed aggregates')
ax.legend()
ax.grid()
fig.savefig('XKaW_reflectivity.pdf')
fig.savefig('XKaW_reflectivity.png', dpi=600)

# %%
fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
ax0.scatter(meltedDO10.Dmax, meltedDO10.Ku, vmin=0, vmax=70, label='DO  2014',
            c=10.0*np.ones(meltedDO10.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO20.Dmax, meltedDO20.Ku, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(meltedDO20.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO30.Dmax, meltedDO30.Ku, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(meltedDO30.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO40.Dmax, meltedDO40.Ku, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(meltedDO40.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO50.Dmax, meltedDO50.Ku, vmin=0, vmax=70, label=None,
            c=50.0*np.ones(meltedDO50.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO60.Dmax, meltedDO60.Ku, vmin=0, vmax=70, label=None,
            c=60.0*np.ones(meltedDO60.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(meltedDO70.Dmax, meltedDO70.Ku, vmin=0, vmax=70, label=None,
            c=70.0*np.ones(meltedDO70.Dmax.shape), marker='+', cmap='jet')
ax0.scatter(dataBJ2_10.Dmax, dataBJ2_10.Ku, vmin=0, vmax=70, label=None,
            c=10.0*np.ones(dataBJ2_10.Dmax.shape), marker='h', cmap='jet')
ax0.scatter(dataBJ2_20.Dmax, dataBJ2_20.Ku, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(dataBJ2_20.Dmax.shape), marker='h', cmap='jet')
ax0.scatter(dataBJ2_30.Dmax, dataBJ2_30.Ku, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(dataBJ2_30.Dmax.shape), marker='h', cmap='jet')
ax0.scatter(dataBJ2_40.Dmax, dataBJ2_40.Ku, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(dataBJ2_40.Dmax.shape), marker='h', cmap='jet')
ax0.scatter(dataBJ2_49.Dmax, dataBJ2_49.Ku, vmin=0, vmax=70, label='BJ  2016',
            c=49.0*np.ones(dataBJ2_49.Dmax.shape), marker='h', cmap='jet')
s = ax0.scatter(dataBJ2_70.Dmax, dataBJ2_70.Ku, vmin=0, vmax=70, label=None,
                c=70.0*np.ones(dataBJ2_70.Dmax.shape), marker='h', cmap='jet')

ax1.scatter(meltedDO10.Dmax, meltedDO10.Ka, vmin=0, vmax=70, label='DO  2014',
            c=10.0*np.ones(meltedDO10.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO20.Dmax, meltedDO20.Ka, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(meltedDO20.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO30.Dmax, meltedDO30.Ka, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(meltedDO30.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO40.Dmax, meltedDO40.Ka, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(meltedDO40.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO50.Dmax, meltedDO50.Ka, vmin=0, vmax=70, label=None,
            c=50.0*np.ones(meltedDO50.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO60.Dmax, meltedDO60.Ka, vmin=0, vmax=70, label=None,
            c=60.0*np.ones(meltedDO60.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(meltedDO70.Dmax, meltedDO70.Ka, vmin=0, vmax=70, label=None,
            c=70.0*np.ones(meltedDO70.Dmax.shape), marker='+', cmap='jet')
ax1.scatter(dataBJ2_10.Dmax, dataBJ2_10.Ka, vmin=0, vmax=70, label=None,
            c=10.0*np.ones(dataBJ2_10.Dmax.shape), marker='h', cmap='jet')
ax1.scatter(dataBJ2_20.Dmax, dataBJ2_20.Ka, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(dataBJ2_20.Dmax.shape), marker='h', cmap='jet')
ax1.scatter(dataBJ2_30.Dmax, dataBJ2_30.Ka, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(dataBJ2_30.Dmax.shape), marker='h', cmap='jet')
ax1.scatter(dataBJ2_40.Dmax, dataBJ2_40.Ka, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(dataBJ2_40.Dmax.shape), marker='h', cmap='jet')
ax1.scatter(dataBJ2_49.Dmax, dataBJ2_49.Ka, vmin=0, vmax=70, label='BJ  2016',
            c=49.0*np.ones(dataBJ2_49.Dmax.shape), marker='h', cmap='jet')
s = ax1.scatter(dataBJ2_70.Dmax, dataBJ2_70.Ka, vmin=0, vmax=70, label=None,
                c=70.0*np.ones(dataBJ2_70.Dmax.shape), marker='h', cmap='jet')

ax2.scatter(meltedDO10.Dmax, meltedDO10.W, vmin=0, vmax=70,
            label='Ori et al.  (2014)',
            c=10.0*np.ones(meltedDO10.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO20.Dmax, meltedDO20.W, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(meltedDO20.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO30.Dmax, meltedDO30.W, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(meltedDO30.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO40.Dmax, meltedDO40.W, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(meltedDO40.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO50.Dmax, meltedDO50.W, vmin=0, vmax=70, label=None,
            c=50.0*np.ones(meltedDO50.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO60.Dmax, meltedDO60.W, vmin=0, vmax=70, label=None,
            c=60.0*np.ones(meltedDO60.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(meltedDO70.Dmax, meltedDO70.W, vmin=0, vmax=70, label=None,
            c=70.0*np.ones(meltedDO70.Dmax.shape), marker='+', cmap='jet')
ax2.scatter(dataBJ2_10.Dmax, dataBJ2_10.W, vmin=0, vmax=70, label=None,
            c=10.0*np.ones(dataBJ2_10.Dmax.shape), marker='h', cmap='jet')
ax2.scatter(dataBJ2_20.Dmax, dataBJ2_20.W, vmin=0, vmax=70, label=None,
            c=20.0*np.ones(dataBJ2_20.Dmax.shape), marker='h', cmap='jet')
ax2.scatter(dataBJ2_30.Dmax, dataBJ2_30.W, vmin=0, vmax=70, label=None,
            c=30.0*np.ones(dataBJ2_30.Dmax.shape), marker='h', cmap='jet')
ax2.scatter(dataBJ2_40.Dmax, dataBJ2_40.W, vmin=0, vmax=70, label=None,
            c=40.0*np.ones(dataBJ2_40.Dmax.shape), marker='h', cmap='jet')
ax2.scatter(dataBJ2_49.Dmax, dataBJ2_49.W, vmin=0, vmax=70,
            label='Johnson et al. (2016) aggregate 2',
            c=49.0*np.ones(dataBJ2_49.Dmax.shape), marker='h', cmap='jet')
s = ax2.scatter(dataBJ2_70.Dmax, dataBJ2_70.W, vmin=0, vmax=70, label=None,
                c=70.0*np.ones(dataBJ2_70.Dmax.shape), marker='h', cmap='jet')

ax0.set_ylim([1e-8, 1e1])
ax1.set_ylim([1e-8, 1e1])
ax2.set_ylim([1e-8, 1e1])
ax2.set_xlim([0, 16])
ax0.grid()
ax1.grid()
ax2.grid()
ax2.set_xlabel('maximum dimension    [mm]')
ax0.set_ylabel('radar backscattering cross section    [mm$^2$]')
ax1.set_ylabel('radar backscattering cross section    [mm$^2$]')
ax2.set_ylabel('radar backscattering cross section    [mm$^2$]')
ax0.set_yscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax0.set_title('Ku band    13.6 GHz')
ax1.set_title('Ka band    35.6 GHz')
ax2.set_title('W band    94 GHz')
ax2.legend(loc=4)
colorbar = plt.colorbar(mappable=s, ax=ax0, label='melted fraction   [%]')
colorbar = plt.colorbar(mappable=s, ax=ax1, label='melted fraction   [%]')
colorbar = plt.colorbar(mappable=s, ax=ax2, label='melted fraction   [%]')
fig.tight_layout()
fig.savefig('melted_comp_Dmax_backscatt.pdf')
fig.savefig('melted_comp_Dmax_backscatt.png',dpi=300)

# %%

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.CextX.values, c='k',
           label='melted 10%')
ax[0].plot(1e3*reff(1e-3*meltedDO50.mkg), meltedDO50.CextX.values, c='0.5',
           label='melted 50%')
ax[0].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.CextKa.values, c='k', ls='--')
ax[0].plot(1e3*reff(1e-3*meltedDO50.mkg), meltedDO50.CextKa.values, c='0.5', ls='--')
ax[0].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.CextW.values, c='k', ls='-.')
ax[0].plot(1e3*reff(1e-3*meltedDO50.mkg), meltedDO50.CextW.values, c='0.5', ls='-.')

# ax.set_xlim(0.3,16)
# ax.set_ylim(1.0e-6,1.0e0)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Equivolume radius [mm]')
ax[0].set_ylabel('Extinction cross section  [mm$^2$]')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.X, c='k')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO50.X.values, c='0.5')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.Ka, c='k', ls='--')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO50.Ka.values, c='0.5', ls='--')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO10.W, c='k', ls='-.')
ax[1].plot(1e3*reff(1e-3*meltedDO10.mkg), meltedDO50.W.values, c='0.5', ls='-.')
ax[1].set_xlim(0, 2)
ax[1].set_ylim(1.0e-6, 1.0e2)
ax[1].set_xticks([0.0,0.5,1.0,1.5,2.0])
# ax.set_xscale('log')
ax[1].set_yscale('log')
ax[1].grid()
ax[1].legend()
ax[1].set_xlabel('Equivolume radius [mm]')
ax[1].set_ylabel('Backscattering cross section  [mm$^2$]')
fig.tight_layout()
fig.savefig('book_melted.pdf')
fig.savefig('book_melted.png',dpi=600)


# %%

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.scatter(meltedDO10.Dmax, meltedDO10.mkg, label='Ori et al. (2014)')
ax.scatter(dataBJ2.Dmax, 1.0e6*dataBJ2.mkg,
           label='Johnson et al. (2016) - agg2')
ax.plot(meltedDO10.Dmax, Br07(meltedDO10.Dmax), label='Brandes et al. (2007)')
ax.plot(meltedDO10.Dmax, BF95(meltedDO10.Dmax),
        label='Brown and Francis (1995)')
ax.set_xlim(0.3, 16)
ax.set_ylim(1.0e-6, 1.0e-1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.legend()
ax.set_xlabel('maximum dimension  [mm]')
ax.set_ylabel('mass    [g]')
fig.savefig('melted_mass_size.pdf')
fig.savefig('melted_mass_size.png',dpi=600)

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.scatter(meltedDO10.mkg, meltedDO10.Ku, vmin=0, vmax=70, label='DO  2014',
           c=10.0*np.ones(meltedDO10.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO20.mkg, meltedDO20.Ku, vmin=0, vmax=70, label=None,
           c=20.0*np.ones(meltedDO20.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO30.mkg, meltedDO30.Ku, vmin=0, vmax=70, label=None,
           c=30.0*np.ones(meltedDO30.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO40.mkg, meltedDO40.Ku, vmin=0, vmax=70, label=None,
           c=40.0*np.ones(meltedDO40.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO50.mkg, meltedDO50.Ku, vmin=0, vmax=70, label=None,
           c=50.0*np.ones(meltedDO50.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO60.mkg, meltedDO60.Ku, vmin=0, vmax=70, label=None,
           c=60.0*np.ones(meltedDO60.mkg.shape), marker='+', cmap='jet')
ax.scatter(meltedDO70.mkg, meltedDO70.Ku, vmin=0, vmax=70, label=None,
           c=70.0*np.ones(meltedDO70.mkg.shape), marker='+', cmap='jet')
ax.scatter(1e6*dataBJ2_10.mkg, dataBJ2_10.Ku, vmin=0, vmax=70, label=None,
           c=10.0*np.ones(dataBJ2_10.mkg.shape), marker='h', cmap='jet')
ax.scatter(1e6*dataBJ2_20.mkg, dataBJ2_20.Ku, vmin=0, vmax=70, label=None,
           c=20.0*np.ones(dataBJ2_20.mkg.shape), marker='h', cmap='jet')
ax.scatter(1e6*dataBJ2_30.mkg, dataBJ2_30.Ku, vmin=0, vmax=70, label=None,
           c=30.0*np.ones(dataBJ2_30.mkg.shape), marker='h', cmap='jet')
ax.scatter(1e6*dataBJ2_40.mkg, dataBJ2_40.Ku, vmin=0, vmax=70, label=None,
           c=40.0*np.ones(dataBJ2_40.mkg.shape), marker='h', cmap='jet')
ax.scatter(1e6*dataBJ2_49.mkg, dataBJ2_49.Ku, vmin=0, vmax=70, label='BJ  2016',
           c=49.0*np.ones(dataBJ2_49.mkg.shape), marker='h', cmap='jet')
s = ax.scatter(1e6*dataBJ2_70.mkg, dataBJ2_70.Ku, vmin=0, vmax=70, label=None,
               c=70.0*np.ones(dataBJ2_70.mkg.shape), marker='h', cmap='jet')
ax.set_ylim([1e-8, 1e1])
# ax.set_xlim([0, 16])
ax.grid()
ax.set_xlabel('mass    [g]')
ax.set_ylabel('radar backscattering cross section    [mm$^2$]')
ax.set_yscale('log')
ax.legend()
colorbar = plt.colorbar(mappable=s, ax=ax, label='melted fraction   [%]')
fig.savefig('melted_comp_mkg_backscatt.pdf')
fig.savefig('melted_comp_mkg_backscatt.png',dpi=600)
