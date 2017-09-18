# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:41:29 2017

@author: dori
"""

import pandas as pd
from glob import glob
from scattDB import scattering

scattfolder = '/work/DBs/melted_aggregate_scaled_reff_Ku_Ka_W_89_165_183/'
subfolder = 'melt3a_aggregate3_060_20091210_040243_20scl_f000001_AEFF_1500_ROT535_13.4/'

avgfiles = glob(scattfolder+subfolder+'*.avg')
scatt = scattering.ScattDDSCAT(avgfiles[0])

scattfolder = '/work/DBs/0/'
subfolder = '400/run000_read_g18m1_775/'
sfld = scattfolder + subfolder
scatt = scattering.ScattADDA(logfile=sfld+'log',muellerfile=sfld+'mueller',
                             csfile=sfld+'CrossSec')