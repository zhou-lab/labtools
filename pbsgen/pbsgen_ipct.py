#!/usr/bin/env python

import subprocess
import sys
import os
from optparse import OptionParser

template="""
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N {self.jobname}
#PBS -d {self.working_dir}
#PBS -e {self.stderr}
#PBS -o {self.stdout}
#PBS -l nodes=1:ppn={self.ppn},walltime={self.time}
#PBS -V

{self.commands}
        
"""

class Job():
    
    pbsdir    = '/home/wzhou1/pbs'
    stdoutdir = '/home/wzhou1/stdout'
    stderrdir = '/home/wzhou1/stderr'
    def __init__(self):

        self.jobroot     = "WandingJob"
        self.working_dir = "/home/wzhou1/"
        self.batch_cnt   = 1
        self.ppn         = "1"
        self.hours       = 200
        if not os.path.exists(Job.stdoutdir):
            raise Exception("stdout dir %s not existent" % Job.stdoutdir)
        if not os.path.exists(Job.stderrdir):
            raise Exception("stderr dir %s not existent" % Job.stderrdir)

    def set_batch_cnt(self):
        self.batch_cnt = 0
        while (True):
            if not os.path.exists('%s/%s_batch_%d.pbs' % (Job.pbsdir, self.jobroot, self.batch_cnt)):
                break
            self.batch_cnt += 1

    def gen(self):

        self.set_batch_cnt()
        self.jobname = self.jobroot+"_batch_%d" % self.batch_cnt
        self.stdout  = "%s/%s.stdout" % (Job.stdoutdir, self.jobname)
        self.stderr  = "%s/%s.stderr" % (Job.stderrdir, self.jobname)
        self.time    = "%d:00:00" % self.hours
        if self.hours <= 4:
            self.queue = "short"
        with open('%s/%s.pbs' % (Job.pbsdir, self.jobname), 'w') as fout:
            fout.write(template.format(self=self))

    def clean(self):
        import os
        map(os.unlink, [os.path.join(Job.pbsdir, _) for _ in os.listdir(Job.pbsdir) if _.endswith(".pbs")])
        map(os.unlink, [os.path.join(Job.stdoutdir, _) for _ in os.listdir(Job.stdoutdir) if _.endswith(".stdout")])
        map(os.unlink, [os.path.join(Job.stderrdir, _) for _ in os.listdir(Job.stderrdir) if _.endswith(".stderr")])
        self.batch_cnt = 1

    def set_hour(self, hours):
        self.hours = hours


def main(opts, args):
    job = Job()
    if not opts.nc:
        job.clean()

    if opts.arg_fns:
        command_args = zip(*[[line.strip() for line in open(_)] for _ in opts.arg_fns])
        if not command_args:
            return
        n_batches = (len(command_args)-1) / opts.bsize + 1
        for i in xrange(n_batches):
            command = "echo `date`\n\n"
            for command_arg in command_args[i*opts.bsize:(i+1)*opts.bsize]:
                command += (args[0] % command_arg)+"\n"
            command += "\n\necho `date` Done.\n"
            job.commands = command
            job.set_hour(opts.wtime)
            job.gen()

        print "generated %d pbs files." % n_batches

    else:
        job.commands = args[0]
        job.set_hour(opts.wtime)
        job.gen()
    
if __name__ == "__main__":
    parser = OptionParser("Usage: %prog [options] command")
    parser.add_option("-n", "--noclean",
                      dest="nc", action="store_true", help="no cleaning of pbs directory")
    parser.add_option("-t", "--wtime",
                      dest="wtime",
                      type="int",
                      default=96,
                      help="wall time limit in hours [96]",
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
