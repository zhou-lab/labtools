#!/usr/bin/env python

import subprocess
import os, sys, re

for line in sys.stdin:
    ar = re.split(r'[\'"]', line.strip())

    B = []
    for i, e in enumerate(ar):
        if i % 2 == 0:
            if len(e.strip())>0:
                B.extend(e.strip().split())
        else:
            B.append(e.strip())
    # print B
    subprocess.check_call(B)
    # subprocess.call(B)
    # subprocess.call(['transvar', 'panno', '--ccds', '-i', 'PIK3CA:E545K'])

