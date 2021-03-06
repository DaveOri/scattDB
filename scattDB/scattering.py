# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:48:16 2017

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
import numpy as np
from scipy import interpolate

def get_line(lines,string):
    return [x for x in lines if string in x][0]

class Scatterer(object):
    def __init__(self):
        print('A scattering object has been created')
        self.wl = None
        self.D = None 
        self.sig_bk = None
        self.sig_abs = None
        self.sig_ext = None
        self.sig_sca = None
        self.S = None
        self.Z = None
        
    def eff2xsect(self,eff):
        return eff*(np.pi*self.aeff**2.)            
    
    @property        
    def radar_xsect(self, h_pol=True):
        """Radar cross section for the current setup.    
    
        Args:
            h_pol: If True (default), use horizontal polarization.
            If False, use vertical polarization.
    
        Returns:
            The radar cross section.
        """        
        Z = self.Z
        if h_pol:
            return 2.*np.pi*(Z[0,0]-Z[0,1]-Z[1,0]+Z[1,1])
        else:
            return 2.*np.pi*(Z[0,0]-Z[0,1]-Z[1,0]+Z[1,1])
        
    @property
    def ldr(self, h_pol=True):
        """
        Linear depolarizarion ratio (LDR) for the current setup.
    
        Args:
            h_pol: If True (default), return LDR_h.
            If False, return LDR_v.
    
        Returns:
           The LDR.
        """
        Z = self.Z
        if h_pol:
            return (Z[0,0]-Z[0,1]+Z[1,0]-Z[1,1])/(Z[0,0]-Z[0,1]-Z[1,0]+Z[1,1])
        else:
            return (Z[0,0]+Z[0,1]-Z[1,0]-Z[1,1])/(Z[0,0]+Z[0,1]+Z[1,0]+Z[1,1])
            
    @property
    def xpol_xsect(self, h_pol=True):
        """
        Cross polarized radar cross section for the current setup
        Args:
            h_pol: If True (default), use horizontal polarization.
            If False, use vertical polarization.
    
        Returns:
            The radar cross section.
        """
        Z = self.Z
        if h_pol:
            return 2.*np.pi*(Z[0,0]-Z[0,1]+Z[1,0]-Z[1,1])
        else:
            return 2.*np.pi*(Z[0,0]+Z[0,1]-Z[1,0]-Z[1,1])

class ScattDDSCAT(Scatterer):
    def __init__(self, filename=None, D=None, melt=0.0, mass=None):
        if filename is not None:
            self.D = D # Size passed from outside, scattering is agnostic of shape and scales unitless
            self.melt = melt
            self.mass = mass
            
            scattfile = open(filename,'r')
            lines = scattfile.readlines()
            self.Ndipoles = int(lines[5].split()[0])
            self.d = float(lines[7].split()[0])
            self.aeff = float(lines[8].split()[1])*1e-3 # convert to millimeters
            self.wl = float(lines[9].split()[1])*1e-3 # convert to millimeters
            self.k = 2*np.pi/self.wl
            qs = lines[31].split()
            self.Qe = float(qs[1])
            self.Qa = float(qs[2])
            self.Qs = float(qs[3])
            self.Qbk = float(qs[6])
            self.g = float(qs[4])
            self.phase = pd.read_csv(filename,sep='\s+',header=37)
            self.sig_bk = self.eff2xsect(self.Qbk)
            self.sig_abs = self.eff2xsect(self.Qa)
            self.sig_ext = self.eff2xsect(self.Qe)
            self.sig_sca = self.eff2xsect(self.Qs)
            idx_bk = self.phase[(self.phase.theta == 180.0) & (self.phase.phi == 0.0)].index[0]
            back = self.phase.iloc[idx_bk]

            self.S = np.array([[back['S_11'],back['S_12']],[back['S_21'],back['S_22']]]) # TODO extend to the full 4x4
            self.Z = self.S/self.k**2.
            
class ScattADDA(Scatterer):
    def __init__(self, logfile=None, muellerfile=None, csfile=None, D=None, melt=0.0, mass=None):
        #if filename is not None:
        self.D = D
        self.melt = melt
        self.mass = mass
        
        logf = open(logfile,'r')
        csf  = open(csfile,'r')
        self.mueller = pd.read_csv(muellerfile,sep=' ',index_col='theta')
        back = self.mueller.loc[180.0]
        #self.mueller.set_index('theta',drop=True,inplace=True)
        lines = logf.readlines()
        self.lines=lines
        CSlines = csf.readlines()
        
        Nstring = get_line(lines,'Total number of occupied dipoles')
        self.Ndipoles = int(Nstring.split()[-1])
        self.wl = float(get_line(lines,'lambda:').split()[-1])*1e-3 # convert to millimeters
        self.dpl= float(get_line(lines,'Dipoles/lambda:').split()[-1])
        self.d = self.wl/self.dpl
        self.mass = 0.917e-3*self.Ndipoles*self.d**3.0
        self.k = 2*np.pi/self.wl

        self.Qe = float(CSlines[1].split()[-1])
        self.Qa = float(CSlines[3].split()[-1])
        self.Qs = self.Qe - self.Qa
        #self.Qbk = float(qs[6])
        #self.g = float(qs[4])
        #self.phase = pd.read_csv(filename,sep='\s+',header=37)
        self.sig_bk  = 2*np.pi*(back.s11+back.s22+2*back.s12)/self.k**2.
        self.sig_abs = float(CSlines[2].split()[-1])*1e-6 # convert to millimeters
        self.sig_ext = float(CSlines[0].split()[-1])*1e-6 # convert to millimeters
        self.sig_sca = self.sig_ext - self.sig_abs 
        self.aeff = np.sqrt(self.sig_sca/(np.pi*self.Qs))
        #idx_bk = self.phase[(self.phase.theta == 180.0) & (self.phase.phi == 0.0)].index[0]
        #back = self.phase.iloc[idx_bk]

        self.S = np.array([[back.s11,back.s12],[back.s21,back.s22]]) # TODO extend to the full 4x4
        self.Z = self.S/self.k**2.
        
#    @property
#    def ldr(self, h_pol=True):
#        """
#        Linear depolarizarion ratio (LDR) for ADDA particles???
#    
#        Args:
#            scatterer: a Scatterer instance.
#            h_pol: If True (default), return LDR_h.
#            If False, return LDR_v.
#    
#        Returns:
#           The LDR.
#        """
#        Z = self.Z
#        if h_pol:
#            return (Z[0,0]-Z[1,1] + Z[1,0]-Z[0,1])/(Z[0,0]+Z[1,1] - Z[1,0]-Z[0,1])
#        else:
#            return (Z[0,0]-Z[1,1] - Z[1,0]+Z[0,1])/(Z[0,0]+Z[1,1] + Z[1,0]+Z[0,1])
            
class ScattDist(object):
    """
    This is a distribution of particles.
    It is initialized as an empty object and filled with particles using the
    add_scatterer() method.
    """
    def __init__(self):
        self.distro = []
        
    def __call__(self, D, quantity):
        """ Calling the distribution should return a value for a particular size"""
        return self.get_value(D, quantity)
        
    def add_scatterer(self,scatterer):
        """ This function simply add a scatterer to the list """
        self.distro.append(scatterer) # TODO check for availability of size
    
#    def sort(self):
#        """ this function should sort the distro according to some size """
        
    def get_value(self, D, quantity):
        """
        This method returns a value for the for the requested quantity at size D
        by linearly interpolating the two closest values (in log space?)
        It also extrapolate assuming boundary values as perpetual
        """
        dist = self.get_distro(quantity)
        dist.sort(axis=0)
        x = dist[:,0]
        y = dist[:,1:]
        interp = interpolate.interp1d(x, y, axis=0, bounds_error=False,
                                      fill_value=(y[0],y[-1]),
                                      assume_sorted=True)
        return interp(D)
               
    def get_distro(self, quantity):
        """
        This function return the numpy array of the quantities asched in the list quantity.
        The first value is always the size of the scatterer
            ---> part_0.D, part_0.quantity[0], part_0.quantity[1], ...
            ---> part_1.D, part_1.quantity[0], part_1.quantity[1], ...
            ---> ...
        """
        if isinstance(quantity, str):
            return np.array([[scat.D,getattr(scat,quantity)] for scat in self.distro])
        else:
            return np.array([[scat.D]+[getattr(scat, quant) for quant in quantity] for scat in self.distro])
    