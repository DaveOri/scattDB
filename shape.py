# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:23:18 2017
@author: dori

Copyright (C) 2017 Davide Ori, University of Cologne

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


class shape(object):
    """ General shape of a particle
    still do not know precisely what to put here ...
    """

    def __init__(self):
        print('A shape object has been created')
        
class shapeDDA(shape):
    """ This is a specific DDA shapefile
    it will be flexible enough to read at least ADDA and DDSCAT formats TODO:versions?
    """
    def __init__(self,filename=None):
        """ at the moment we assume the shapefile is from ddscat7.3
        """
        if filename is not None:
            shapefile = open(filename,"r")
            lines = shapefile.readlines()
            #Nlines = len(lines)
            self.Ndipoles = int(lines[-1].split()[0])
            geometry_start_line = len(lines) - self.Ndipoles
            self.shape = pd.read_csv(filename,
                                     sep='\s+',
                                     skiprows=geometry_start_line,
                                     header=None,
                                     names=['N','X','Y','Z','CX','CY','CZ'])
            self.substances = self.shape[['CX','CY','CZ']].drop_duplicates()
            
             
    def draw(self):
        """ Handy function to draw particle shape """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(len(self.substances)):
            currshape=self.shape[(self.shape['CX']==self.substances['CX'][i])&
                                 (self.shape['CY']==self.substances['CY'][i])&
                                 (self.shape['CZ']==self.substances['CZ'][i])]
            xs = currshape.X.values
            ys = currshape.Y.values
            zs = currshape.Z.values
            ax.scatter(xs, ys, zs)#, s=20, c=None, depthshade=True)
            
    def find_dmax(self):
        """ Find maximum dimension """
        
    def find_dsphere(self):
        """ Find diameter of the minimum enclosing sphere """
        
    
        
    