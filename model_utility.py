#! /usr/bin/python
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

# Utility functions.
from pprint import pprint
import os.path as path
import numpy as np
import re
import utility
from models import *


#--------------------------------------------------------------------------------
# If new models are added, they must be added to this list.
#--------------------------------------------------------------------------------
def model_list():
    return [arruda_boyce, van_der_waals,
            mooney_rivlin, polynomial2,         # Poly models.
            neo_hookean, yeoh, reduced_poly6,   # Reduced poly models.
            ogden3, ogden6, ogden9]             # Ogden models.

#--------------------------------------------------------------------------------
# Get a model module from a string name.
#--------------------------------------------------------------------------------
def get_model(mname):
    for model in model_list():
        if model.name()==mname: break
    else:
        utility.print_error('Invalid material model: '+mname,False)
        print 'Options are:'
        for model in model_list():
            print '\t',model.name(),'\t',model.pname(),' - ',model.descr()
        print '\n'
        exit(1)
    return model

#--------------------------------------------------------------------------------
# Get name from type
#--------------------------------------------------------------------------------
def get_name(t):
    name = str()
    if t.find('u')>=0: name+='Uniaxial'
    if t.find('b')>=0: name+='Biaxial'
    if t.find('p')>=0: name+='Planar'
    if len(name)==0:
        utility.print_error('Invalid type in get_name: '+t,True)
    return name


#--------------------------------------------------------------------------------
# Get function from type
#--------------------------------------------------------------------------------
def get_function(model,t):
    if t=='u': return model.stressU
    if t=='b': return model.stressB
    if t=='p': return model.stressP
    utility.print_error('Invalid type in get_function: '+t,True)


#--------------------------------------------------------------------------------
# Function to write the parameters to a file.
#--------------------------------------------------------------------------------
def write_matfile(model,descr,t,params,D=[],density=0.0):
    fname = descr+'--'+model.name()+'--'+get_name(t)+'.mat'
    name = model.name()
    order = re.search('(\d+)$',name)
    with open(fname,'w') as f:
        f.write('# Material file for '+model.pname()+' fit to '+get_name(t)+' data.\n')
        f.write('# Note: matname and density fields will be used while creating Abaqus input files.\n')
        f.write('# Parameters are: '+model.params()+'\n')
        f.write('# Units should probably be N/mm^2 (depending on the model). Some parameters may be unitless.\n')
        f.write('matname: ???\n')     # Material name, not model name.
        f.write('density: '+str(density)+'\n')
        if order:
            f.write('model:   '+name[:-len(order.group(0))]+'\n')
            f.write('order:   '+order.group(0)+'\n')
        else:
            f.write('model:   '+name+'\n')
        f.write('params:  '+' '.join(map(str,params))+'\n')
        if D==[]: f.write('nu:      0.48\n')
        else:     f.write('D:       '+' '.join(map(str,D))+'\n')
    return fname


#--------------------------------------------------------------------------------
# Function to write viscoelastic parameters to a file.
#--------------------------------------------------------------------------------
def write_viscomatfile(descr,G0,params,dimensionless=False):
    terms = len(params)/2
    fname = descr+'--'+str(terms)+'terms.mat'
    prony = np.reshape(params,(terms,2)) * 1.0          # 1.0 ensures we create a copy of the data.
    if not dimensionless: prony[:,0] = prony[:,0] / G0  # Non-dimensionalize.
    prony = prony[prony[:,1].argsort()]                 # Sort by time terms.
    with open(fname,'w') as f:
        f.write('# Material file for viscoelastic fit with '+str(terms)+' terms.\n')
        f.write('# Values are given in dimensionless form, except G0 which is in N/mm2.\n')
        f.write('G0:  '+str(G0)+'\n')
        f.write('g:   '+' '.join(map(str,prony[:,0]))+'\n')
        f.write('tau: '+' '.join(map(str,prony[:,1]))+'\n')
    return fname


#--------------------------------------------------------------------------------
# Read the material definition from a file.
#--------------------------------------------------------------------------------
def read_matfile(fname):
    mat = dict()
    mat['file'] = fname
    try:
        with open(fname,'r') as f:
            for line in f:
                vals = (line[:line.rfind('#')]).split()
                if(len(vals)==0): continue
                if  (vals[0]=='#'): continue
                elif(vals[0]=='matname:'): mat['matname'] = vals[1]
                elif(vals[0]=='density:'): mat['density'] = float(vals[1])
                elif(vals[0]=='model:'):   mat['model'] = vals[1]
                elif(vals[0]=='order:'):   mat['order'] = int(vals[1])
                elif(vals[0]=='params:'):  mat['params'] = [float(v) for v in vals[1:]]
                elif(vals[0]=='D:'):       mat['D'] = [float(v) for v in vals[1:]]
                elif(vals[0]=='nu:'):      mat['nu'] = float(vals[1])
                elif(vals[0]=='hyst_S:'):  mat['hyst_S'] = float(vals[1])
                elif(vals[0]=='hyst_A:'):  mat['hyst_A'] = float(vals[1])
                elif(vals[0]=='hyst_m:'):  mat['hyst_m'] = float(vals[1])
                elif(vals[0]=='hyst_C:'):  mat['hyst_C'] = float(vals[1])
                else:
                    print 'WARNING: unrecognized lines in material DB.'
                    print 'model_utility.py may need to be updated for the new format.'
                    print 'LINE: \"'+line+'\"'
    except Exception:
        raise IOError('Material file not found: '+fname)
    if 'nu' in mat and not 'D' in mat:
        if 'order' in mat: model = get_model(mat['model']+str(mat['order']))
        else: model = get_model(mat['model'])
        mat['D'] = model.compressibility(mat['nu'],*mat['params'])
    # pprint(mat)
    return mat


#--------------------------------------------------------------------------------
# Read the viscoelastic definition from a file.
#--------------------------------------------------------------------------------
def read_viscomatfile(fname):
    mat = dict()
    mat['file'] = fname
    try:
        with open(fname,'r') as f:
            for line in f:
                vals = (line[:line.rfind('#')]).split()
                if(len(vals)==0): continue
                if  (vals[0]=='#'): continue
                elif(vals[0]=='G0:'):  mat['G0'] = float(vals[1])
                elif(vals[0]=='g:'):   mat['g'] = [float(v) for v in vals[1:]]
                elif(vals[0]=='k:'):   mat['k'] = [float(v) for v in vals[1:]]
                elif(vals[0]=='tau:'): mat['tau'] = [float(v) for v in vals[1:]]
                else:
                    print 'WARNING: unrecognized lines in viscomaterial DB.'
                    print 'model_utility.py may need to be updated for the new format.'
                    print 'LINE: \"'+line+'\"'
    except Exception:
        raise IOError('Material file not found: '+fname)
    if ('g' in mat and 'tau' in mat) and (not len(mat['g'])==len(mat['tau'])):
        print 'ERROR: viscoelastic parameters must have the same number of terms.'
        exit(1)
    if ('k' in mat and 'tau' in mat) and (not len(mat['g'])==len(mat['tau'])):
        print 'ERROR: viscoelastic parameters must have the same number of terms.'
        exit(1)
    return mat


#--------------------------------------------------------------------------------
# Open and import a dataset.
#   datafile = name of file.dat to import.
#   data     = dictionary to contain imported data.
#   key      = keyname to insert into data.
#--------------------------------------------------------------------------------
def import_dataset(datafile,data,key):
    # fliplr swaps the stress and strain columns. Data is now in "strain,stress" format.
    if not datafile.lower()=='none':
        if not path.isfile(datafile):
            utility.print_error('File does not exist: '+datafile,True)
        data[key] = np.fliplr(np.loadtxt(datafile,comments='#',delimiter=',',dtype=np.longdouble))
        data[key] = np.vstack(([0.0,0.0],data[key]))    # Add a (0,0) datapoint.
        print '  Imported',data[key].shape[0],'datapoints for',get_name(key),'data.'


