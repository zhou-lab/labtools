#!/usr/bin/env python

from pbsgen import pbsgen_main, Default

def set_queue(job):
    job.queue = 'default'

default = Default()
default.pbsdir = '/primary/home/wanding.zhou/pbs'
default.scriptdir = '%s/pbs/' % default.pbsdir
default.stdoutdir = '%s/stdout/' % default.pbsdir
default.stderrdir = '%s/stderr/' % default.pbsdir
default.hour = 12
default.ppn = 1
default.memG = 2

pbsgen_main(default, set_queue)
