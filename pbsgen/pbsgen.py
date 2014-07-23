#!/usr/bin/env python

import subprocess
import sys
import os
import argparse
import types

template="""
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N {self.jobname}
#PBS -e {self.stderr}
#PBS -o {self.stdout}
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.ppn},walltime={self.time}
#PBS -V

{self.commands}
        
"""

class Job():
    
    
    def __init__(self, args):

        self.pbs        = args.pbs
        self.jobroot    = "WandingJob"
        self.ppn        = "1"
        self.hour       = args.hour
        self.stderrdir  = args.stderr
        self.stdoutdir  = args.stdout
        self.set_queue = types.MethodType(args.set_queue, self)
        if not os.path.exists(self.stdoutdir):
            raise Exception("stdout dir %s not existent" % self.stdoutdir)
        if not os.path.exists(self.stderrdir):
            raise Exception("stderr dir %s not existent" % self.stderrdir)

    def next_index_auto(self):
        batch_index = 0
        while (True):
            if not os.path.exists('%s/%s_batch_%d.pbs' %
                                  (self.pbs, self.jobroot, batch_index)):
                return batch_index
            batch_index += 1

    def gen(self, batch_index):

        if batch_index < 0:
            batch_index = self.next_index_auto()

        self.jobname = self.jobroot+"_batch_%d" % batch_index
        self.stdout = "%s/%s.stdout" % (self.stdoutdir, self.jobname)
        self.stderr = "%s/%s.stderr" % (self.stderrdir, self.jobname)
        self.time = "%d:00:00" % self.hour
        self.set_queue()
        with open('%s/%s.pbs' % (self.pbs, self.jobname), 'w') as fout:
            fout.write(template.format(self=self))

    def clean(self):
        import os
        map(os.unlink, [os.path.join(self.pbs, _) for _ in os.listdir(self.pbs) if _.endswith(".pbs")])
        map(os.unlink, [os.path.join(self.stdoutdir, _) for _ in os.listdir(self.stdoutdir) if _.endswith(".stdout")])
        map(os.unlink, [os.path.join(self.stderrdir, _) for _ in os.listdir(self.stderrdir) if _.endswith(".stderr")])


class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    rgs = indstr.split(',')
    indices = Indices()
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices

def main_clean(args):

    job = Job(args)
    job.clean()

def main_one(args):

    job = Job(args)
    job.commands = "echo `date`\n\n"
    job.commands += args.command
    job.commands += "\n\necho `date` Done.\n"
    job.gen(args.index)

def main_batch(args):

    job = Job(args)
    job.clean()

    command_args_t = []
    for argpair in args.args.split(','):
        argkey, argfn, argcol = argpair.split(":")
        argcol = int(argcol) - 1
        file_args = []
        for line in open(argfn):
            fields = line.strip().split(args.delim)
            file_args.append((argkey, fields[argcol]))
        command_args_t.append(file_args)
    command_args = zip(*command_args_t)

    n_batches = (len(command_args)-1) / args.bsize + 1
    for i in xrange(n_batches):
        command = "echo `date`\n\n"
        for command_arg in command_args[i*args.bsize:(i+1)*args.bsize]:
            command += args.command.format(**dict(command_arg))+"\n"
        command += "\n\necho `date` Done.\n"
        job.commands = command
        job.gen(i)

    print "generated %d pbs files." % n_batches

def main_submit(args):

    job = Job(args)
    beg,end = map(int, args.srange.split(','))
    pbs_fns = [fn for fn in os.listdir(job.pbs) if fn.endswith('.pbs')]
    pbs_fns.sort()
    for fn in pbs_fns[beg-1:end]:
        print job.pbs+'/'+fn
        subprocess.check_call(['qsub', job.pbs+'/'+fn])

def add_default_settings(parser, d):

    """ default setting """

    parser.add_argument("-pbs", default=d.scriptdir, 
                        help="pbs directory [%s]" % d.scriptdir)
    parser.add_argument("-hour", default=d.hour, 
                        help="wall time in hour [%d]" % d.hour)
    parser.add_argument("-stdout", default=d.stdoutdir,
                        help="stdout [%s]" % d.stdoutdir)
    parser.add_argument('-stderr', default=d.stderrdir,
                        help='stderr [%s]' % d.stderrdir)
    parser.add_argument("-ppn", default=d.ppn, 
                        help="ppn [%d]" % d.ppn)
    parser.add_argument("-jobroot", default="WandingJob", 
                        help="Job root [WandingJob]")
    parser.add_argument("-memG", default=d.memG, 
                        help="memory in G [%d]" % d.memG)

def pbsgen_main(setting, set_queue):
    
    parser = argparse.ArgumentParser(description="generate pbs files")
    subparsers = parser.add_subparsers()

    ##### clean ####
    parser_clean = subparsers.add_parser("clean", help="clean pbs files")
    add_default_settings(parser_clean, setting)
    parser_clean.set_defaults(func=main_clean)

    ###### one #####
    parser_one = subparsers.add_parser("one", help="generate one pbs script")
    parser_one.add_argument('command', help='command to run')
    parser_one.add_argument("-index", type=int, default=-1, help="index of next pbs script [inferred from existing file names]")
    add_default_settings(parser_one, setting)
    parser_one.set_defaults(func=main_one)

    ###### batch #####
    parser_batch = subparsers.add_parser("batch", help="batch generate pbs scripts")
    parser_batch.add_argument('command', help='command to run')
    parser_batch.add_argument('-args', help="mappable argument, will be used to substitute {symbol} structure in command. (e.g., -args symbol1:argfile1:col1,symbol2:argfile2:col2)")
    parser_batch.add_argument('-delim', default='\t', help='delimiter in the argument table')
    parser_batch.add_argument('-bsize', type=int, default=1, help="batch size")
    add_default_settings(parser_batch, setting)
    parser_batch.set_defaults(func=main_batch)

    ###### submit #####
    parser_submit = subparsers.add_parser('submit', help='submit jobs')
    parser_submit.add_argument('srange', help='range in the sorted list of scripts, 1-based, e.g., 1,100')
    add_default_settings(parser_submit, setting)
    parser_submit.set_defaults(func=main_submit)

    args = parser.parse_args()
    args.set_queue = set_queue
    args.func(args)

class Default:

    def __init__(self):

        self.pbsdir = "/home/wzhou1/pbs03"
        self.scriptdir = '%s/pbs/' % self.pbsdir
        self.stdoutdir = '%s/stdout/' % self.pbsdir
        self.stderrdir = '%s/stderr/' % self.pbsdir
        self.hour = 96
        self.ppn = 1
        self.memG = 2


if __name__ == "__main__":

    def set_queue_htc(job):
        if job.hour <= 4:
            job.queue = "short"
        else:
            job.queue = "medium"

    def set_queue_ipct(job):
        job.queue = 'batch'

    def set_queue_hpc(job):
        if job.hour <= 1:
            job.queue = "short"
        else:
            job.queue = "medium"

    pbsgen_main(Default(), set_queue_htc)
