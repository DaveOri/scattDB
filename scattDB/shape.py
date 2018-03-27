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
from scipy.spatial import ConvexHull


class shape(object):
    """ General shape of a particle
    still do not know precisely what to put here ...
    """

    def __init__(self):
        print('A shape object has been created')
        self.Ndipoles = None
        self.shape = None
        self.substances = None
        self.hull = None
        self.dmax = None
        self.dsphere = None

    def find_dmax(self,aeff=None):
        """ Find maximum dimension """
        if self.hull is None:
            self.hull = ConvexHull(self.shape[['X','Y','Z']])
        dmax = 0
        vtx = self.hull.points[self.hull.vertices]
        for i in range(len(vtx)-1):
            for j in range(i+1,len(vtx)):
                d = ((vtx[i]-vtx[j])**2).sum()
                if d > dmax:
                    dmax = d
        if aeff is not None:
            self.d = aeff*np.cbrt(4.*np.pi/(self.Ndipoles*3.))
            self.dmax = self.d*dmax**0.5
        else:
            self.dmax = dmax**0.5
        return self.dmax
        
#    def find_dsphere(self):
#        """ Find diameter of the minimum enclosing sphere """
#        if self.hull is None:
#            self.hull = ConvexHull(self.shape[['X','Y','Z']])
#        print('Method find_dsphere is not implemented yet')
            
    def get_melted_fraction(self,idx=0):
        """ Here I assume there are two substances, one is water, the other ice.
            Not sure about how to distinguish them ...
        """
        if len(self.substances) == 2:
            Ndip1 = len(self.shape[(self.shape['CX']==self.substances['CX'].iloc[idx])&
                                   (self.shape['CY']==self.substances['CY'].iloc[idx])&
                                   (self.shape['CZ']==self.substances['CZ'].iloc[idx])])
            return float(Ndip1)/float(self.Ndipoles)
        elif len(self.substances) == 1:
            return 0.0 # assumed to be pure ice
        else:
            print('This shapefile have ', len(self.substances), ' substances')
            return None
    
    def find_mass(self,aeff=None,d=None):
        if aeff is not None:
            self.d = aeff*np.cbrt(4.*np.pi/(self.Ndipoles*3.))
        else:
            self.d = d
        self.mass = 1e-9*0.917*self.Ndipoles*self.d**3 # if I got microns, this should return milligrams
        return self.mass
        
    def draw(self,ax=None):
        """ Handy function to draw particle shape """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for i in range(len(self.substances)):
            currshape=self.shape[(self.shape['CX']==self.substances['CX'].iloc[i])&
                                 (self.shape['CY']==self.substances['CY'].iloc[i])&
                                 (self.shape['CZ']==self.substances['CZ'].iloc[i])]
            xs = currshape.X.values
            ys = currshape.Y.values
            zs = currshape.Z.values
            ax.scatter(xs, ys, zs, s=1)#, s=20, c=None, depthshade=True)
        
class shapeDDA(shape):
    """ This is a specific DDA shapefile
    it will be flexible enough to read at least ADDA and DDSCAT formats TODO:versions?
    """
    def __init__(self,filename=None):
        """ at the moment we assume the shapefile is from ddscat7.3, or 7.2
        """
        if filename is not None:
            shapefile = open(filename,"r")
            lines = shapefile.readlines()
            #Nlines = len(lines)
            self.Ndipoles = int(float(lines[-1].split()[0]))
            if self.Ndipoles > 1e6: # TODO quick fix for my particles
                self.Ndipoles = int(lines[1].split()[0])
            geometry_start_line = len(lines) - self.Ndipoles
            self.shape = pd.read_csv(filename,
                                     sep='\s+',
                                     skiprows=geometry_start_line,
                                     header=None,
                                     names=['N','X','Y','Z','CX','CY','CZ'])
            self.substances = self.shape[['CX','CY','CZ']].drop_duplicates()
            self.hull = None
            self.dmax = None
            self.dsphere = None

shp400=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/400/shape.adda')
#shp400.draw()

shp1200=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/1200/shape.adda')
#shp1200.draw()

shp1500=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/1500/shape.adda')
#shp1200.draw()

shp1600=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/1600/shape.adda')
#shp1200.draw()

shp2400=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/2400/shape.adda')
#shp2400.draw()

shp5100=shapeDDA('/data/optimice/scattering_databases/DavideOri_2014/10/5100/shape.adda')
#shp5100.draw()
    
