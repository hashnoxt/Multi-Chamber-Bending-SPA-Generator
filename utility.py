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
import subprocess as sub
import sys,os,shutil,time


#--------------------------------------------------------------------------------
# Print a failure warning message.
#--------------------------------------------------------------------------------
def print_error(msg,qBOOL):
    print '\n********************************************************************************'
    print '** ERROR: ',msg
    print '********************************************************************************'
    if qBOOL==True: exit(1)


#------------------------------------------------------------------------------
# Create a link, error if source does not exist.
#------------------------------------------------------------------------------
def safe_hardlink(source,target):
    if(os.path.exists(source)):
        try:        # Linux.
            os.link(source,target)
        except:     # Windows.
            shutil.copyfile(source,target)
        return 1
    else:
        print_error('Not found: '+source,False)
        return 0


#------------------------------------------------------------------------------
# Print time elapsed.
#------------------------------------------------------------------------------
def time_elapsed(tstart):
    tfin = time.time()
    gmtime = time.gmtime(tfin-tstart)
    day = int(time.strftime('%d',gmtime))-1
    return time.strftime(str(day)+'-%H:%M:%S d-hh:mm:ss',gmtime)


#------------------------------------------------------------------------------
# Find a file by searching successively higher directories.
#------------------------------------------------------------------------------
def find_file(name):
    new_name = name
    for i in range(5):
        if(os.path.exists(new_name)): break
        new_name = '../'+new_name
    else:
        print_error('Could not find file: '+name,False)
        print '\tCurrent directory: '+os.getcwd()
        print '\tSearched up to:    '+new_name
        exit(1)
    return new_name


#------------------------------------------------------------------------------
# Run a terminal command with error checking, and return results.
#------------------------------------------------------------------------------
def run_cmd(cmd):
    try:        # Linux.
        p = sub.Popen(cmd,stdout=sub.PIPE,stderr=sub.PIPE)
    except:
        try:    # Windows.
            p = sub.Popen(cmd,stdout=sub.PIPE,stderr=sub.PIPE,shell=True)
        except:
            print_error('ERROR: failed to RUN command (not command failed), with message:',False)
            print sys.exc_info()[0]
            print 'Command was:\n\t'+list_to_str(cmd)+'\n'
            exit(1)
    (stdout,stderr) = p.communicate()
    if(p.returncode != 0):
        print('\nERROR: command failed with returncode '+str(p.returncode)+
              '. Command was:\n\t'+list_to_str(cmd)+'\n'+stderr)
        exit(1)
    return stdout


#------------------------------------------------------------------------------
# Run a terminal command with error checking and print output to screen.
#------------------------------------------------------------------------------
def run_cmd_screen(cmd):
    try:        # Linux.
        sub.check_call(cmd)
    except:     # Windows.
        sub.check_call(cmd,shell=True)


#------------------------------------------------------------------------------
# Run a terminal command and print output to file.
#------------------------------------------------------------------------------
def run_cmd_log(cmd,fname):
    f = open(fname,'w')
    try:        # Linux.
        sub.call(cmd,stdout=f,stderr=sub.STDOUT)
    except:     # Windows.
        sub.call(cmd,stdout=f,stderr=sub.STDOUT,shell=True)
    f.close()


#------------------------------------------------------------------------------
# Convert a command list into a prettier string.
#------------------------------------------------------------------------------
def list_to_str(cmd):
    out = '"'
    for i in cmd: out += i+' '
    return out+'"'


#------------------------------------------------------------------------------
# Create and chdir to a new directory. Optionally delete it if it exists.
#------------------------------------------------------------------------------
def open_dir(dir,rm=False):
    if(rm and os.path.exists(dir)): shutil.rmtree(dir)
    if(not os.path.exists(dir)): os.makedirs(dir)
    os.chdir(dir)


#------------------------------------------------------------------------------
# Set common plot options.
# If the x-axis label gets cropped out, adjust the pad_inches parameter.
#------------------------------------------------------------------------------
def plot_opts(fname,xlabel,ylabel,title,location='best'):
    import matplotlib.pyplot as plt
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if location=='outside right':
        plt.legend(loc=2,bbox_to_anchor=(1,1))
    elif not location=='none':
        plt.legend(loc=location,frameon=False,framealpha=0)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(fname)


#------------------------------------------------------------------------------
# Kronecker delta function.
#------------------------------------------------------------------------------
def kd(i,j):
    if(i==j): return 1.0
    return 0.0
