# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:25:00 2017

@author: dori
"""

import pandas as pd
import numpy as np

cols = ['Dmax','Nsample','Niter','coH','coV','crosspol']
fileC  = pd.read_csv('ice_fractal_C.dat',sep=' ',names=cols)
fileKu = pd.read_csv('ice_fractal_Ku.dat',sep=' ',names=cols)
fileKa = pd.read_csv('ice_fractal_Ka.dat',sep=' ',names=cols)
fileS  = pd.read_csv('ice_fractal_S.dat',sep=' ',names=cols)
fileW  = pd.read_csv('ice_fractal_W.dat',sep=' ',names=cols)

newcols = ['Dmax','X','Ku','Ka','W','ldr','xpolKa','mkg']
data = pd.DataFrame(index=fileC.index,columns=newcols)
for i in fileC.index:
    data.loc[i] = fileC.loc[i,'Dmax'],1.0e6*fileC.loc[i,'coH'],1.0e6*fileKu.loc[i,'coH'],1.0e6*fileKa.loc[i,'coH'],1.0e6*fileW.loc[i,'coH'],fileKa.loc[i,'crosspol']/fileKa.loc[i,'coH'],1.0e6*fileKa.loc[i,'crosspol'],np.nan

data = data.astype(float)
data.sort(columns='Dmax',inplace=True)
data.to_csv('JT_fractals.csv')