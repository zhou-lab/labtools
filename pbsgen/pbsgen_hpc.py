#!/RIS/HPC_apps/AMD/python/python-2.7.2/bin/python

import subprocess
import sys
from optparse import OptionParser

template="""
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N {self.jobname}
#PBS -d {self.working_dir}
#PBS -e {self.stderr}
#PBS -o {self.stdout}
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.ppn},walltime={self.time}
#PBS -V

{self.commands}
        
"""

class Job():
    
    pbsdir    = '/scratch/bcb/wzhou1/pbs/'
    stdoutdir = '/scratch/bcb/wzhou1/stdout'
    stderrdir = '/scratch/bcb/wzhou1/stderr'
    def __init__(self):

        self.jobroot     = "WandingJob"
        self.working_dir = "/scratch/bcb/wzhou1/"
        self.batch_cnt   = 1
        self.ppn         = "24"
        self.queue       = "medium"
        self.hours       = 24

    def gen(self):
        
        self.jobname = self.jobroot+"_batch_%d" % self.batch_cnt
        self.stdout  = "%s/%s.stdout" % (Job.stdoutdir, self.jobname)
        self.stderr  = "%s/%s.stderr" % (Job.stderrdir, self.jobname)
        self.time    = "%d:00:00" % self.hours
        if self.hours <= 1:
            self.queue = "short"
        with open('%s/%s.pbs' % (Job.pbsdir, self.jobname), 'w') as fout:
            fout.write(template.format(self=self))
        self.batch_cnt += 1

    def clean(self):
        import os
        map(os.unlink, [os.path.join(Job.pbsdir, _) for _ in os.listdir(Job.pbsdir) if _.endswith(".pbs")])
        map(os.unlink, [os.path.join(Job.stdoutdir, _) for _ in os.listdir(Job.stdoutdir) if _.endswith(".stdout")])
        map(os.unlink, [os.path.join(Job.stderrdir, _) for _ in os.listdir(Job.stderrdir) if _.endswith(".stderr")])
        self.batch_cnt = 1

    def set_memory(self, memory_in_G):
        self.ppn = int(memory_in_G / 1.3) + 1
        if self.ppn >= 32:
            raise Exception("memory exceeded")
        
    def set_hour(self, hours):
        self.hours = hours


def main(opts, args):
    job = Job()
    job.clean()
    
    if opts.arg_fns:
        command_args = zip(*[[line.strip() for line in open(_)] for _ in opts.arg_fns])
        if not command_args:
            return
        n_batches = (len(command_args)-1) / opts.bsize + 1
        for i in xrange(n_batches):
            command = ""
            for command_arg in command_args[i*opts.bsize:(i+1)*opts.bsize]:
                command += (args[0] % command_arg)+"\n"
            job.commands = command
            job.set_memory(opts.memory)
            job.set_hour(opts.wtime)
            job.gen()
            
        print "generated %d pbs files." % n_batches
    else:
        job.commands = args[0]
        job.set_memory(opts.memory)
        job.set_hour(opts.wtime)
        job.gen()
                
if __name__ == "__main__":

    parser = OptionParser("Usage: %prog [options] command")
    parser.add_option("-m", "--memory",
                      dest="memory",
                      type="int",
                      default=2,
                      help="maximum memory in Gigabytes [2]",
                      metavar="memory")
    parser.add_option("-t", "--wtime",
                      dest="wtime",
                      type="int",
                      default=24,
                      help="wall time limit in hours [24]",
                      metavar="wtime")
    parser.add_option("-b", "--batchsize",
                      dest="bsize",
                      type="int",
                      default=1,
                      help="batch size",
                      metavar="batchsize")
    parser.add_option("-a", "--args",
                      dest="arg_fns",
                      type="string",
                      action="append",
                      help="argument list file names, will be used to substitute %s in command following the order of appearance. (Note that multiple argument must be specified by \"-a arg1 -a arg2\")",
                      metavar="args")
    
    
    options, args = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit()
        
    main(options, args)
    
    
