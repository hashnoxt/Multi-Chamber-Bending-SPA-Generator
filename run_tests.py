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

# Run a test or series of tests, with optional variations.
import matplotlib
matplotlib.use('Agg')                   # Matplotlib uses the png renderer for the images. Needed for cluster.
import argparse, sys, os, shutil, time
import os.path as path
import utility
ABAQUS='abaqus'                # Assume Abaqus is on the path.



#--------------------------------------------------------------------------------
# Cluster properties.
#--------------------------------------------------------------------------------
def cluster_properties(name):
    cp = dict()
    if name=='bellatrix':
        cp['name']='bellatrix'
        cp['np']=16
        cp['mem']=32000
        cp['modules']='''
# NOTE: a different queue is required for jobs on Bellatrix longer than 3 days.
module list
# module purge
# module load intelmpi/5.0.1
# module load abaqus/6.14-1/intel-14.0.1
# module load python/2.7.9/intel-15.0.0
'''
    elif name=='aries':
        cp['name']='aries'
        cp['np']=48
        cp['mem']=192000
        cp['modules']='''
eval `modulecmd bash purge`
eval `modulecmd bash load intel/14.0.1`
eval `modulecmd bash load intelmpi/4.1.3`
eval `modulecmd bash load abaqus/6.14-1/intel-14.0.1`
eval `modulecmd bash load python/2.7.8`
'''
    elif name=='deneb':
        cp['name']='deneb'
        cp['np']=16
        cp['mem']=64000
        cp['modules']='''
module list
# module purge
# module load intelmpi/5.0.1
# module load abaqus/6.14-1/intel-14.0.1
# module load python/2.7.9/intel-15.0.0
'''
    elif name=='deneb_himem':
         cp['name']='deneb_himem'
         cp['np']=16
         cp['mem']=256000
         cp['modules']='''
 module list
 # module purge
 # module load intelmpi/5.0.1
 # module load abaqus/6.14-1/intel-14.0.1
 # module load python/2.7.9/intel-15.0.0
 '''
    else:
        utility.print_error('Unrecognized cluster.',True)
    return cp


#--------------------------------------------------------------------------------
# Create the abaqus_v6 environment file.
#   np     = Number of processors to use per node.
#   wdir   = working directory to create file.
#   nn     = Number of nodes to use.
#--------------------------------------------------------------------------------
def create_env(np,wdir,nn=1):
    # Note that this file will be edited for cluster usage with the list of hosts.
    with open(path.join(wdir,'abaqus_v6.env'),'w') as f:
        f.write('import os\n'
                'memory=\"'+str(98/np)+' %\"\n'     # This is memory per process.
                'cpus='+str(np*nn)+'\n')


#--------------------------------------------------------------------------------
# Run an Abaqus problem.
#   inp = Abaqus INP file, including path.
#   outdir = Directory to create/run the problem in.
#   np     = Number of processors to use.
#   descr  = Optional description to add to outdir name.
#--------------------------------------------------------------------------------
def run_abq(inp,outdir,np,descr):
    wd = path.join(outdir,path.splitext(path.basename(inp))[0])
    if descr!="": wd+='-'+descr
    print '\n'
    print '--------------------------------------------------------------------------------'
    print '-- RUNNING ABAQUS:',inp,descr,'('+str(np)+' processors)'
    print '--             IN:',wd
    print '--------------------------------------------------------------------------------'
    if path.exists(wd): shutil.rmtree(wd)
    if not inp.endswith('.inp'):
        utility.print_error('Input file must be a .inp Abaqus run file.',False)
        return 0
    # Prepare the directories.
    os.makedirs(wd)
    inpNEW = path.join(wd,'job.inp')
    if not utility.safe_hardlink(inp,inpNEW):
        shutil.rmtree(wd)
        return 0
    create_env(np,wd)
    # Run the simulation.
    rootdir = os.getcwd()
    os.chdir(wd)
    tstart = time.time()
    cmd=[ABAQUS,'job=job','input=job.inp','interactive']
    try: utility.run_cmd_screen(cmd)
    except Exception:
        utility.print_error(sys.exc_info()[0],False)
        return 0
    finally:
        print 'TIME ELAPSED: ',utility.time_elapsed(tstart)
        os.chdir(rootdir)
    return 1



#--------------------------------------------------------------------------------
# Postprocess Abaqus results.
#   inp = Abaqus INP file, including path.
#   outdir = Directory to create/run the problem in.
#   descr  = Optional description to add to outdir name.
#--------------------------------------------------------------------------------
def run_post(inp,outdir,descr):
    wd = path.join(outdir,path.splitext(path.basename(inp))[0])
    if descr!="": wd+='-'+descr
    print '\n'
    print '--------------------------------------------------------------------------------'
    print '-- RUNNING POSTPROCESSOR:',path.join(wd,'job.odb')
    print '--------------------------------------------------------------------------------'
    # Run the extraction.
    rootdir = os.getcwd()
    try:
        os.chdir(wd)
    except Exception:
        utility.print_error(sys.exc_info()[0],False)
        return 0
    tstart = time.time()
    script = utility.find_file('python/abq_extract_data.py')
    cmd=[ABAQUS,'cae','noGUI='+script,'--','job.odb','data.rpt']
    try:
        utility.run_cmd_screen(cmd)
        with open('data.tmp','w') as of, open('data.rpt','r') as f:
            for line in f:
                line = line.rstrip()    # Remove newline/carriage-return.
                if line.find('NoValue') < 0: of.write(line+'\n')
        shutil.move('data.tmp','data.rpt')
    except Exception:
        utility.print_error(sys.exc_info()[0],False)
        print 'TIME ELAPSED: ',utility.time_elapsed(tstart)
        os.chdir(rootdir)
        return 0
    try:
        import create_result_plots, numpy
        from matplotlib import pyplot
        pyplot.rc('mathtext',default='regular') # Don't use italics for mathmode.
        fig,ax = pyplot.subplots()
        create_result_plots.plot_fem_data(ax,'data.rpt')
        lgd = ax.legend(loc='best',frameon=False,framealpha=0)
        pyplot.savefig('plot.png',bbox_extra_artists=(lgd,), bbox_inches='tight')
        pyplot.close()
        if path.exists('data-energy.rpt'):
            fig,ax = pyplot.subplots()
            create_result_plots.plot_fem_energy(ax,'data-energy.rpt')
            lgd = ax.legend(loc='best',frameon=False,framealpha=0)
            pyplot.savefig('plot-energy.png',bbox_extra_artists=(lgd,), bbox_inches='tight')
            pyplot.close()
        if path.exists('data-strain.csv'):
            fig,ax = pyplot.subplots()
            data = numpy.loadtxt('data-strain.csv',delimiter=',')
            ax.plot(data[:,0],data[:,1],'o-',label='Nominal Strain')
            ax.set_xlabel('Time (meaningless units)')
            ax.set_ylabel('Nominal Strain')
            ax.grid()
            pyplot.savefig('plot-strain.png', bbox_inches='tight')
            pyplot.close()
    except Exception:
        utility.print_error('Failed to create result plots.',False)
    print 'TIME ELAPSED: ',utility.time_elapsed(tstart)
    os.chdir(rootdir)
    return 1


#--------------------------------------------------------------------------------
# Create and submit slurm jobs on the cluster.
#   inp    = Abaqus INP file, including path.
#   cp     = dictionary of cluster properties.
#   time   = Runtime to request.
#   outdir = Directory to create/run the problem in.
#   nn     = Number of nodes to use.
#   np     = Number of processors to use per node.
#   descr  = Optional description to add to outdir name.
#--------------------------------------------------------------------------------
def run_slurm(inp,cp,time,outdir,nn,np,descr):
    outdir = path.abspath(outdir)
    dirname = path.splitext(path.basename(inp))[0]
    if descr!="": dirname+='-'+descr
    wd = path.join(outdir,dirname)
    if wd.startswith('/home'):
        print 'ERROR: /home is readonly during execution, use /scratch.'
        return 0,''
    if path.exists(wd): shutil.rmtree(wd)
    if not inp.endswith('.inp'):
        utility.print_error('Input file must be a .inp Abaqus run file.',False)
        return 0,''
    # Prepare the directories.
    os.makedirs(wd)
    shutil.copyfile(inp, path.join(wd,'job.inp'))
    create_env(np,wd,nn)
    # Create the submission script.
    extfile = path.abspath(utility.find_file('python/abq_extract_data.py'))
    jobfile = path.join(wd,dirname+'.job')
    runcmd = ABAQUS+' job=job input=job.inp interactive'
    postcmd = ABAQUS+' cae noGUI='+extfile+' -- job.odb data.rpt'
    # TODO - need ability to set email address here.
    with open(jobfile,'w') as f:
        f.write('''#!/bin/bash
# mem directive is per node.
#SBATCH --nodes '''+str(nn)+'''
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node '''+str(np)+'''
#SBATCH --mem '''+str(cp['mem']*np/cp['np'])+'''
#SBATCH --time '''+time+'''
#SBATCH --workdir '''+wd+'''
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=philip.moseley@epfl.ch
#SBATCH --share
''')
        f.write(cp['modules'])
        f.write('''echo '**********************************************************************'
echo 'Starting execution at' `date`
echo 'Job script:' $0
''')
        if nn!=1: f.write('''
scontrol show hostname $SLURM_JOB_NODELIST | paste -d -s > hostlist
HLIST="["
for i in $(cat hostlist)
do
HLIST="$HLIST['$i','''+str(np)+'''],"
done
HLIST=`echo $HLIST | sed -e "s/,$/]/"`
echo "mp_host_list=$HLIST" >> abaqus_v6.env
echo "Hostlist: $HLIST"
''')
        f.write('''echo '**********************************************************************'
'''+runcmd+'''
sleep 5
echo
echo '**********************************************************************'
echo 'Running postprocessor.'
echo '**********************************************************************'
'''+postcmd+'''
echo
echo 'Finished execution at' `date`
echo " ****** END OF JOB ******"''')
    return 1,jobfile


#--------------------------------------------------------------------------------
# Create and submit an array job on the cluster.
#   outdir = Directory to create/run the problem in.
#   cp     = dictionary of cluster properties.
#   time   = Runtime to request.
#   nn     = Number of nodes to use.
#   np     = Number of processors to use per node.
#   descr  = Optional description to add to outdir name.
#--------------------------------------------------------------------------------
def run_array(outdir,jobs,NS,cp,time,nn,np,descr):
    job = path.join(outdir,'array'+descr+'.job')
    with open(job,'w') as f:
        f.write('''#!/bin/bash
# mem directive is per node.
#SBATCH --nodes '''+str(nn)+'''
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node '''+str(np)+'''
#SBATCH --mem '''+str(cp['mem']*np/cp['np'])+'''
#SBATCH --time '''+time+'''
#SBATCH --workdir '''+outdir+'''
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=philip.moseley@epfl.ch
#SBATCH --share
''')
        if NS>0: f.write('#SBATCH --array=0-'+str(len(jobs)-1)+'%'+str(NS)+'\n')
        else:    f.write('#SBATCH --array=0-'+str(len(jobs)-1)+'\n')
        f.write('case "$SLURM_ARRAY_TASK_ID" in\n')
        for i,j in enumerate(jobs):
            wd = path.dirname(j)
            f.write('\t'+str(i)+')\n')
            f.write('\t\tcd '+wd+'\n')
            f.write('\t\tchmod +x '+j+'\n')
            f.write('\t\t'+j+'\n')
            f.write('\t\t;;\n')
        f.write('\t*) ;;\n')
        f.write('esac\n')
    print '\n\nSubmitting array job with',len(jobs),'jobs...'
    utility.run_cmd_screen(['sbatch',job])


#--------------------------------------------------------------------------------
# Main.
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    tinit = time.time()
    # Handle user input.
    parser = argparse.ArgumentParser(description='Run a test or series of tests, with optional variations.',
                                     epilog='Example: run_tests.py -a --sort . ../inputs/cube/*')
    parser.add_argument('outdir',help='Base directory to create output.')
    parser.add_argument('inp',nargs='+',help='Abaqus INP file or list of INP files.')
    parser.add_argument('-c','--cluster',choices=['bellatrix','aries','deneb','deneb_himem'],help='Create and submit Slurm jobs on the cluster.')
    parser.add_argument('-a','--all',action='store_true',help='Run all the actions (eg, --run and --post).')
    parser.add_argument('-r','--run',action='store_true',help='Run an Abaqus simulation.')
    parser.add_argument('-p','--post',action='store_true',help='Postprocess results.')
    parser.add_argument('-nn','--nodes',metavar='NN',default=1,type=int,help='Number of nodes to use (default=1).')
    parser.add_argument('-np','--procs',metavar='NP',type=int,help='Number of processors to use per node (default=8 for non-cluster).')
    parser.add_argument('--sort',action='store_true',help='Numerically sort inputs from smallest to largest.')
    parser.add_argument('--descr',help='Optional description to add to directory name.',default='')
    parser.add_argument('--qtime',default='00:30:00',metavar='hh:mm:ss',help='Requested queue time on cluster (default=30min).')
    parser.add_argument('--array',metavar='NS',type=int,default=-1,help='Submit the jobs as an array, with <=NS running simultaneous.'
                            'Good for avoiding using all the Abaqus licenses. NS=0 runs all array jobs simultaneously.')
    args = parser.parse_args()

    tC = 0
    tA = 0
    # Handle inputs.
    if args.sort:
        try:
            from natsort import natsort
            args.inp = natsort(args.inp)
        except:
            print 'WARNING: no natsort module found, sorting not available.'
    if args.cluster: cp=cluster_properties(args.cluster)
    if not args.procs:
        if args.cluster: np=cp['np']
        else: np=8
    else: np=args.procs

    # Print a summary.
    print 'Tasks to run:',len(args.inp)
    for inp in args.inp:
        print '\t',inp

    # Run tasks.
    jobs = list()
    for inp in args.inp:
        if args.cluster:
            r,jobfile = run_slurm(inp,cp,args.qtime,args.outdir,args.nodes,np,args.descr)
            tC += r
            tA += 1
            # Submit the job.
            if args.array<0:
                # log = utility.run_cmd(['sbatch','--qos=debug',jobfile])
                log = utility.run_cmd(['sbatch',jobfile])
                print log.rstrip(),'--',jobfile
            else: jobs.append(jobfile)
        else:
            if args.all or args.run:  tC+=run_abq(inp,args.outdir,np,args.descr); tA+=1
            if args.all or args.post: tC+=run_post(inp,args.outdir,args.descr); tA+=1

    # Submit array job if requested.
    if args.array>=0: run_array(args.outdir,jobs,args.array,cp,args.qtime,args.nodes,np,args.descr)

    print '\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print 'TASKS ATTEMPTED:    ',tA
    print 'TASKS COMPLETED:    ',tC
    print 'TOTAL TIME ELAPSED: ',utility.time_elapsed(tinit)
    if tA!=tC: print '\t!!!!!!!!!! ',tA-tC,'tasks failed  !!!!!!!!!!'
    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

