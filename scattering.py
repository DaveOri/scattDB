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

class scattering(object):
    def __init__(self):
        print('A scattering object has been created')
        self.sigbk = None
        self.sigabs = None
        self.sigext = None
        self.sigsca = None

class scattDDSCAT(scattering):
    def __init__(self,filename=None):
        if filename is not None:
            scattfile = open(filename,'r')
            lines = scattfile.readlines()
            self.lines=lines
            self.Ndipoles = int(lines[5].split()[0])
            self.d = float(lines[7].split()[0])
            self.aeff = float(lines[8].split()[1])
            self.wl = float(lines[9].split()[1])
            qs = lines[31].split()
            self.Qe = float(qs[1])
            self.Qa = float(qs[2])
            self.Qs = float(qs[3])
            self.Qbk = float(qs[6])
            self.g = float(qs[4])
            self.phase = pd.read_csv(filename,sep='\s+',header=37)
            
            
            