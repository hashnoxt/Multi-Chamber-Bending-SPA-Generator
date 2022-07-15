#! /usr/bin/env python
# The MIT License (MIT)
#
# Copyright (c) 2015, EPFL Reconfigurable Robotics Laboratory,
#                     Philip Moseley, philip.moseley@gmail.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np


#--------------------------------------------------------------------------------
# Material model name.
#--------------------------------------------------------------------------------
def name():   return 'ogden3'
def pname():  return 'Ogden-3'
def params(): return 'u1 a1 u2 a2 u3 a3'
def descr():  return 'Ogden Model with order 3 (modified form).'


# NOTE - this is the Abaqus form of the functions. Ogden2004 is similar, but they
#        show these functions as being multiplied by (a[i]/2.0)

#--------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
#--------------------------------------------------------------------------------
def stressU(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-0.5*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-0.5*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-0.5*a3-1.0)) / a3
    return S1+S2+S3


#--------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
#--------------------------------------------------------------------------------
def stressB(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-2.0*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-2.0*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-2.0*a3-1.0)) / a3
    return S1+S2+S3


#--------------------------------------------------------------------------------
# Function defining the planar stress given strain.
#--------------------------------------------------------------------------------
def stressP(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-a3-1.0)) / a3
    return S1+S2+S3


#--------------------------------------------------------------------------------
# Calculate the Ds
#--------------------------------------------------------------------------------
def compressibility(v, u1, a1, u2, a2, u3, a3):
    # This sum is what's in the ABQ manual (and what ABQ calculates with the data).
    # We get an error message which implies that u1 is what ABQ actually expects.
    # I believe the error message to be incorrect; setting u0=u1 typically results
    # in a much less compressible material, even though the error goes away.
    #  u0 = u1
    u0 = u1+u2+u3
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1,0.0,0.0]
