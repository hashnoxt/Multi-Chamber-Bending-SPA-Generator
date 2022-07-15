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
def name():   return 'ab'
def pname():  return 'Arruda-Boyce'
def params(): return 'mu lambda_m'
def descr():  return 'Arruda-Boyce Model.'


#--------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
#--------------------------------------------------------------------------------
def stressU(x, u, Lm):
    L = 1.0+x
    I1 = np.power(L,2.0) + 2.0*np.power(L,-1.0)
    C = [0.5, 0.05, 11.0/1050.0, 19.0/7000.0, 519.0/673750.0]
    s  = C[0] / Lm
    s += I1 * 2.0*C[1] / np.power(Lm,2.0)
    s += np.power(I1,2) * 3.0*C[2] / np.power(Lm,4.0)
    s += np.power(I1,3) * 4.0*C[3] / np.power(Lm,6.0)
    s += np.power(I1,4) * 5.0*C[4] / np.power(Lm,8.0)
    return 2.0*u*(L-np.power(L,-2.0)) * s


#--------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
#--------------------------------------------------------------------------------
def stressB(x, u, Lm):
    L = 1.0+x
    I1 = 2.0*np.power(L,2.0) + np.power(L,-4.0)
    C = [0.5, 0.05, 11.0/1050.0, 19.0/7000.0, 519.0/673750.0]
    s  = C[0] / Lm
    s += I1 * 2.0*C[1] / np.power(Lm,2.0)
    s += np.power(I1,2) * 3.0*C[2] / np.power(Lm,4.0)
    s += np.power(I1,3) * 4.0*C[3] / np.power(Lm,6.0)
    s += np.power(I1,4) * 5.0*C[4] / np.power(Lm,8.0)
    return 2.0*u*(L-np.power(L,-5.0)) * s


#--------------------------------------------------------------------------------
# Function defining the planar stress given strain.
#--------------------------------------------------------------------------------
def stressP(x, u, Lm):
    L = 1.0+x
    I1 = np.power(L,2.0)+np.power(L,-2.0) + 1.0
    C = [0.5, 0.05, 11.0/1050.0, 19.0/7000.0, 519.0/673750.0]
    s  = C[0] / Lm
    s += I1 * 2.0*C[1] / np.power(Lm,2.0)
    s += np.power(I1,2) * 3.0*C[2] / np.power(Lm,4.0)
    s += np.power(I1,3) * 4.0*C[3] / np.power(Lm,6.0)
    s += np.power(I1,4) * 5.0*C[4] / np.power(Lm,8.0)
    return 2.0*u*(L-np.power(L,-3.0)) * s


#--------------------------------------------------------------------------------
# Calculate the Ds
#--------------------------------------------------------------------------------
def compressibility(v, u, Lm):
    u0 = u * (1.0 + 3.0/(5.0*Lm**2) + 99.0/(175.0*Lm**4) + 513.0/(875.0*Lm**6) + 42039.0/(67375.0*Lm**8))
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1]
