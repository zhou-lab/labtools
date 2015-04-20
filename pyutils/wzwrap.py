#!/usr/bin/env python

import os, sys, re

lw = 80
for line in sys.stdin:

    line = line.strip()
    fields = line.split('\t')[:4]
    line = '\t'.join(fields[4:])
    
    while len(line) > 0:
        k = lw-3
        print('   '+line[:k])
        line = line[k:]
    
