#!/usr/bin/python2

import sys
import os
import commands
from commands import getstatusoutput
from commands import getoutput
import datetime

from parser_utils import *

def getProxy():
    stat,out = commands.getstatusoutput("voms-proxy-info -e --valid 5:00")
    if stat:
        raise Exception,"voms proxy not found or validity  less than 5 hours:\n%s" % out
    stat,out = commands.getstatusoutput("voms-proxy-info -p")
    out = out.strip().split("\n")[-1] # remove spurious java info at ic.ac.uk
    if stat:
        raise Exception,"Unable to voms proxy:\n%s" % out
    proxy = out.strip("\n")
    return proxy


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def lxbatchSubmitJob(run, path, cfg, outdir, queue, job_dir, dryrun):
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'.sh'
    jobtar = job_dir+'/job.tar'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('cp '+jobtar+' ./ \n')
    f.write ('tar -xf job.tar \n')
    f.write ('source scripts/setup.sh \n')
    f.write ('make -j 2 \n')
    f.write ('cp '+path+'/'+cfg+' job.cfg \n\n')
    f.write ('bin/H4Reco job.cfg '+run+'\n\n')
    f.write ('cp ntuples/*'+run+'.root '+outdir+'\n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; bsub -q ' + queue + ' ' + '-u ' + os.environ['USER'] + '@cern.ch ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def htcondorSubmitJob(runs, path, cfg, outdir, queue, job_dir, dryrun):
    jobname = job_dir+'/H4Reco_condor'
    jobtar = job_dir+'/job.tar'
    #---H4Reco script
    fsh = open (jobname+'.sh', 'w')
    fsh.write ('#!/bin/sh' + '\n\n')
    fsh.write ('declare -a runs=('+' '.join([str(run) for run in runs])+')\n\n')
    fsh.write ('cp '+jobtar+' ./ \n')
    fsh.write ('tar -xf job.tar \n')
    fsh.write ('source scripts/setup.sh \n')
    fsh.write ('make -j 2 \n')
    fsh.write ('cp '+path+'/'+cfg+' job.cfg \n\n')
    fsh.write ('bin/H4Reco job.cfg ${runs[${1}]}\n\n')
    fsh.write ('cp ntuples/*${runs[${1}]}.root '+outdir+'\n')
    fsh.close ()
    #---HTCondor submit file
    fsub = open (jobname+'.sub', 'w')    
    fsub.write('+JobFlavour = "'+queue+'"\n\n')
    fsub.write('executable  = '+jobname+'.sh\n')
    fsub.write('arguments   = $(ProcId)\n')
    fsub.write('output      = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).out\n')
    fsub.write('error       = '+job_dir+'/output/h4reco.$(ClusterId).$(ProcId).err\n')
    fsub.write('log         = '+job_dir+'/log/h4reco.$(ClusterId).log\n\n')
    fsub.write('x509userproxy = '+getProxy()+' \n\n')
    fsub.write('max_retries = 3\n')
    fsub.write('queue '+str(len(runs))+'\n')
    fsub.close()
    #---submit
    getstatusoutput('chmod 755 '+jobname+ '*')
    getstatusoutput('mkdir -p '+job_dir+'/output')
    getstatusoutput('mkdir -p '+job_dir+'/log')    
    if not dryrun:
        ret = getstatusoutput('condor_submit '+jobname+'.sub')
        print(ret)
        
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def herculesSubmitJob(run, path, cfg, outdir, queue, job_dir, dryrun):
    jobname = job_dir+'/H4Reco_'+queue+'_'+run+'.sh'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('mkdir -p /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('cd /gwpool/users/$USER/pool/'+run+'/ \n')
    f.write ('wget https://github.com/simonepigazzini/H4Analysis/archive/master.zip \n')
    f.write ('unzip master.zip \n')
    f.write ('cd H4Analysis-master/ \n\n')
    f.write ('cp '+path+'/ntuples/Template*.root ./ntuples/ \n')
    f.write ('cp '+path+cfg+' job.cfg \n')
    f.write ('source scripts/setup.sh \n')
    f.write ('make -j \n\n')
    f.write ('bin/H4Reco job.cfg '+run+'\n\n')
    f.write ('cp ntuples/*'+run+'.root '+outdir+'\n')
    f.write ('cd /gwpool/users/$USER/pool/ \n')
    f.write ('rm -r '+run+' \n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+job_dir+'; qsub -q ' + queue + ' ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def resubmitFailed(runs, outdir):
    
    successful_runs = getoutput("ls -lhrt "+outdir+"| sed 's:^.*_::g' | sed 's:.root::g'")
    
    successful_runs = set(successful_runs.split('\n'))
    
    return [run for run in runs if run not in successful_runs]    
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'submit H4Reco step to lsf')
    parser.add_argument('-r', '--runs' , action=customAction, help='run to be processed, either list or file')
    parser.add_argument('-q', '--queue' , default = 'longlunch', help='batch queue/condor flavour (def: longlunch)')
    parser.add_argument('-s', '--storage' , default = '/store/group/dpg_ecal/alca_ecalcalib/ECALTB_H4_Fall2015/', help='storage path')
    parser.add_argument('-v', '--version' , default = 'v1', help='production version')
    parser.add_argument('-c', '--cfg' , default = '../cfg/H4DAQ_base.cfg', help='production version')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    parser.add_argument('--batch' , default='condor', help='batch system to use')
    parser.add_argument('--resub' , action="store_true", default=False, help='resubmit failed jobs')
    
    args = parser.parse_args ()

    ## check ntuple version
    stageOutDir = args.storage+'/ntuples_'+args.version+'/'
    stageOutDir = stageOutDir.replace('//', '/')
    
    if args.batch == 'lxbatch':
        if getoutput('ls '+stageOutDir) == "":
            print "ntuples version "+args.version+" directory on eos already exist! no jobs created."
            exit(0)
    getstatusoutput('mkdir -p '+stageOutDir)    
    
    ## job setup
    local_path = getoutput('pwd')
    date = datetime.datetime.now().strftime("%d-%m-%Y")
    job_dir = local_path+"/"+date+"_ntuples_"+args.version
    getstatusoutput('mkdir -p '+job_dir)
    if local_path.endswith('scripts'):
        local_path = local_path[:-len('scripts')]

    if args.cfg.startswith('../'):
        args.cfg = args.cfg[len('../'):]
    
    if len(args.runs) == 1 and os.path.isfile(args.runs[0]):
        runs_file = open(args.runs[0], 'r')
        args.runs  = []
        if runs_file:
            for run in runs_file:
                args.runs.append(run.rstrip())

    ## resubmit failed
    if args.resub:
        args.runs = resubmitFailed(args.runs, stageOutDir)

    ## create jobs
    getstatusoutput('tar --exclude-vcs --exclude="*ntuples*" -cjf '+job_dir+'/job.tar -C '+local_path+' .')
    print 'submitting', len(args.runs), 'jobs to queue', args.queue
    if args.batch == 'condor':
        htcondorSubmitJob(args.runs, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun) 
    else:
        for run in args.runs:
            print 'submitting run: ', run
            if args.batch == 'lxbatch':
                lxbatchSubmitJob(run, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun) 
            if args.batch == 'hercules':
                herculesSubmitJob(run, local_path, args.cfg, stageOutDir, args.queue, job_dir, args.dryrun) 
