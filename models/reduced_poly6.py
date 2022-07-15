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
def name():   return 'rpoly6'
def pname():  return 'Reduced-Polynomial-6'
def params(): return 'C10 C20 C30 C40 C50 C60'
def descr():  return 'Reduced Polynomial Model with order 6.'


#--------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
#--------------------------------------------------------------------------------
def stressU(x, C10, C20, C30, C40, C50, C60):
    L = 1.0+x
    I1 = np.power(L,2.0) + 2.0*np.power(L,-1.0)
    s =  C10 + 2.0*C20*(I1-3.0)
    s += 3.0*C30*np.power(I1-3.0,2.0)
    s += 4.0*C40*np.power(I1-3.0,3.0)
    s += 5.0*C40*np.power(I1-3.0,4.0)
    s += 6.0*C40*np.power(I1-3.0,5.0)
    return 2.0*(L-np.power(L,-2.0)) * s


#--------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
#--------------------------------------------------------------------------------
def stressB(x, C10, C20, C30, C40, C50, C60):
    L = 1.0+x
    I1 = 2.0*np.power(L,2.0) + np.power(L,-4.0)
    s =  C10 + 2.0*C20*(I1-3.0)
    s += 3.0*C30*np.power(I1-3.0,2.0)
    s += 4.0*C40*np.power(I1-3.0,3.0)
    s += 5.0*C40*np.power(I1-3.0,4.0)
    s += 6.0*C40*np.power(I1-3.0,5.0)
    return 2.0*(L-np.power(L,-5.0)) * s


#--------------------------------------------------------------------------------
# Function defining the planar stress given strain.
#--------------------------------------------------------------------------------
def stressP(x, C10, C20, C30, C40, C50, C60):
    L = 1.0+x
    I1 = np.power(L,2.0)+np.power(L,-2.0) + 1.0
    s =  C10 + 2.0*C20*(I1-3.0)
    s += 3.0*C30*np.power(I1-3.0,2.0)
    s += 4.0*C40*np.power(I1-3.0,3.0)
    s += 5.0*C40*np.power(I1-3.0,4.0)
    s += 6.0*C40*np.power(I1-3.0,5.0)
    return 2.0*(L-np.power(L,-3.0)) * s


#--------------------------------------------------------------------------------
# Calculate the Ds
#--------------------------------------------------------------------------------
def compressibility(v, C10, C20, C30, C40, C50, C60):
    u0 = 2.0*C10
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1,0.0,0.0,0.0,0.0,0.0]
