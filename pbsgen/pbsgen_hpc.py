#!/usr/bin/env python

from pbsgen import pbsgen_main, Default

def set_queue(job):
    if job.hour <= 1:
        job.queue = "short"
    elif job.hour > 24:
        job.queue = "long"
    else:
        job.queue = "medium"

default = Default()
default.pbsdir = "/scratch/bcb/wzhou1/pbs"
default.scriptdir = '%s/pbs/' % default.pbsdir
default.stdoutdir = '%s/stdout/' % default.pbsdir
default.stderrdir = '%s/stderr/' % default.pbsdir
default.hour = 24
default.ppn = 1
default.memG = 2

pbsgen_main(default, set_queue)
