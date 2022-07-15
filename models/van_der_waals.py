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
def name():   return 'vdw'
def pname():  return 'Van der Waals'
def params(): return 'mu lambda_m alpha beta'
def descr():  return 'Van der Waals Model.'


#--------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
#--------------------------------------------------------------------------------
def stressU(x, u, Lm, a, B):
    L = 1.0+x
    I1 = np.power(L,2.0) + 2.0*np.power(L,-1.0)
    I2 = np.power(L,-2.0) + 2.0*L
    I = (1.0-B)*I1 + B*I2
    n = np.sqrt((I-3.0)/(np.power(Lm,2.0)-3.0))
    t1 = (1.0/(1.0-n)) - a * np.sqrt(0.5*(I-3.0))
    t2 = L*(1.0-B) + B
    return u*(1.0-np.power(L,-3.0)) * t1 * t2


#--------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
#--------------------------------------------------------------------------------
def stressB(x, u, Lm, a, B):
    L = 1.0+x
    I1 = 2.0*np.power(L,2.0) + np.power(L,-4.0)
    I2 = 2.0*np.power(L,-2.0) + np.power(L,4.0)
    I = (1.0-B)*I1 + B*I2
    n = np.sqrt((I-3.0)/(np.power(Lm,2.0)-3.0))
    t1 = (1.0/(1.0-n)) - a * np.sqrt(0.5*(I-3.0))
    t2 = 1.0 - B + B*np.power(L,2.0)
    return u*(L-np.power(L,-5.0)) * t1 * t2


#--------------------------------------------------------------------------------
# Function defining the planar stress given strain.
#--------------------------------------------------------------------------------
def stressP(x, u, Lm, a, B):
    L = 1.0+x
    I1 = np.power(L,2.0)+np.power(L,-2.0) + 1.0
    I2 = I1
    I = (1.0-B)*I1 + B*I2
    n = np.sqrt((I-3.0)/(np.power(Lm,2.0)-3.0))
    t1 = (1.0/(1.0-n)) - a * np.sqrt(0.5*(I-3.0))
    return u*(L-np.power(L,-3.0)) * t1


#--------------------------------------------------------------------------------
# Calculate the Ds
#--------------------------------------------------------------------------------
def compressibility(v, u, Lm, a, B):
    u0 = u
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1]
