#!/usr/bin/env python

from pbsgen import pbsgen_main, Default

def set_queue(job):
    job.queue = 'batch'


default = Default()
default.pbsdir = "/home/wzhou1/pbs03"
default.scriptdir = '%s/pbs/' % default.pbsdir
default.stdoutdir = '%s/stdout/' % default.pbsdir
default.stderrdir = '%s/stderr/' % default.pbsdir
default.hour = 96
default.ppn = 1
default.memG = 2


pbsgen_main(default, set_queue)
