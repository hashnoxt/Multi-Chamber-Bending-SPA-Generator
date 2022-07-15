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
def name():   return 'ogden9'
def pname():  return 'Ogden-9'
def params(): return 'u1 a1 u2 a2 u3 a3 u4 a4 u5 a5 u6 a6 u7 a7 u8 a8 u9 a9'
def descr():  return 'Ogden Model with order 9 (modified form).'

# NOTE - this is the Abaqus form of the functions. Ogden2004 is similar, but they
#        show these functions as being multiplied by (a[i]/2.0)

#--------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
#--------------------------------------------------------------------------------
def stressU(x, u1, a1, u2, a2, u3, a3, u4, a4, u5, a5, u6, a6, u7, a7, u8, a8, u9, a9):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-0.5*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-0.5*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-0.5*a3-1.0)) / a3
    S4 = 2.0 * u4 * (np.power(L,a4-1.0) - np.power(L,-0.5*a4-1.0)) / a4
    S5 = 2.0 * u5 * (np.power(L,a5-1.0) - np.power(L,-0.5*a5-1.0)) / a5
    S6 = 2.0 * u6 * (np.power(L,a6-1.0) - np.power(L,-0.5*a6-1.0)) / a6
    S7 = 2.0 * u7 * (np.power(L,a7-1.0) - np.power(L,-0.5*a7-1.0)) / a7
    S8 = 2.0 * u8 * (np.power(L,a8-1.0) - np.power(L,-0.5*a8-1.0)) / a8
    S9 = 2.0 * u9 * (np.power(L,a9-1.0) - np.power(L,-0.5*a9-1.0)) / a9
    return S1+S2+S3+S4+S5+S6


#--------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
#--------------------------------------------------------------------------------
def stressB(x, u1, a1, u2, a2, u3, a3, u4, a4, u5, a5, u6, a6, u7, a7, u8, a8, u9, a9):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-2.0*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-2.0*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-2.0*a3-1.0)) / a3
    S4 = 2.0 * u4 * (np.power(L,a4-1.0) - np.power(L,-2.0*a4-1.0)) / a4
    S5 = 2.0 * u5 * (np.power(L,a5-1.0) - np.power(L,-2.0*a5-1.0)) / a5
    S6 = 2.0 * u6 * (np.power(L,a6-1.0) - np.power(L,-2.0*a6-1.0)) / a6
    S7 = 2.0 * u7 * (np.power(L,a7-1.0) - np.power(L,-2.0*a7-1.0)) / a7
    S8 = 2.0 * u8 * (np.power(L,a8-1.0) - np.power(L,-2.0*a8-1.0)) / a8
    S9 = 2.0 * u9 * (np.power(L,a9-1.0) - np.power(L,-2.0*a9-1.0)) / a9
    return S1+S2+S3+S4+S5+S6


#--------------------------------------------------------------------------------
# Function defining the planar stress given strain.
#--------------------------------------------------------------------------------
def stressP(x, u1, a1, u2, a2, u3, a3, u4, a4, u5, a5, u6, a6, u7, a7, u8, a8, u9, a9):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L,a1-1.0) - np.power(L,-a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L,a2-1.0) - np.power(L,-a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L,a3-1.0) - np.power(L,-a3-1.0)) / a3
    S4 = 2.0 * u4 * (np.power(L,a4-1.0) - np.power(L,-a4-1.0)) / a4
    S5 = 2.0 * u5 * (np.power(L,a5-1.0) - np.power(L,-a5-1.0)) / a5
    S6 = 2.0 * u6 * (np.power(L,a6-1.0) - np.power(L,-a6-1.0)) / a6
    S7 = 2.0 * u7 * (np.power(L,a7-1.0) - np.power(L,-a7-1.0)) / a7
    S8 = 2.0 * u8 * (np.power(L,a8-1.0) - np.power(L,-a8-1.0)) / a8
    S9 = 2.0 * u9 * (np.power(L,a9-1.0) - np.power(L,-a9-1.0)) / a9
    return S1+S2+S3+S4+S5+S6


#--------------------------------------------------------------------------------
# Calculate the Ds
#--------------------------------------------------------------------------------
def compressibility(v, u1, u2, u3, u4, u5, u6, u7, u8, u9, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    # This sum is what's in the ABQ manual (and what ABQ calculates with the data).
    # We get an error message which implies that u1 is what ABQ actually expects.
    # I believe the error message to be incorrect; setting u0=u1 typically results
    # in a much less compressible material, even though the error goes away.
    #  u0 = u1
    u0 = u1+u2+u3+u4+u5+u6+u7+u8+u9
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
