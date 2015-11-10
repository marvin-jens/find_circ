#!/usr/bin/env python
import sys
max_l = float(sys.argv[1])

for line in sys.stdin:
    passed = False
    if line.startswith("#"):
        passed = True
    else:
        cols = line.rstrip().split("\t")
        start,end = float(cols[1]),float(cols[2])
        if end - start <= max_l:
            passed = True
        
    if passed:
        print line.rstrip()


