# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 12:52:09 2017

@author: dori

Code adapted from psd module of pytmatrix software of Jussi Leinonen

Copyright (C) 2009-2015 Jussi Leinonen, Finnish Meteorological Institute, 
California Institute of Technology

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

import numpy as np


class PSD(object):
    def __call__(self, D):
        if np.shape(D) == ():
            return 0.0
        else:
            return np.zeros_like(D)

    def __eq__(self, other):
        return False


class ExponentialPSD(PSD):
    """Exponential particle size distribution (PSD).
    
    Callable class to provide an exponential PSD with the given 
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = N0 * exp(-Lambda*D)

    Attributes:
        N0: the intercept parameter.
        Lambda: the inverse scale parameter        
        D_max: the maximum diameter to consider (defaults to 11/Lambda,
            i.e. approx. 3*D0, if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """

    def __init__(self, N0=1.0, Lambda=1.0, D_max=None):
        self.N0 = float(N0)
        self.Lambda = float(Lambda)
        self.D_max = 11.0/Lambda if D_max is None else D_max

    def __call__(self, D):
        psd = self.N0 * np.exp(-self.Lambda*D)
        if np.shape(D) == ():
            if D > self.D_max:
                return 0.0
        else:
            psd[D > self.D_max] = 0.0
        return psd

    def __eq__(self, other):
        try:
            return isinstance(other, ExponentialPSD) and \
                (self.N0 == other.N0) and (self.Lambda == other.Lambda) and \
                (self.D_max == other.D_max)
        except AttributeError:
            return False
