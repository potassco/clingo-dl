#!/usr/bin/env python

from subprocess import call
from sys import exit
calls =  [("clingo-dl examples/taskassignment/encoding-dl.lp examples/taskassignment/tai4_4_1.lp -c n=132 0 --outf=3",30)]
rc = 0
for (i, rc_expected) in calls:
    rc_call = call(i, shell=True)
    if (rc_call != rc_expected):
        rc = 1
        print("Failed: '{}' with error code ".format(i, rc_call))
if rc == 0:
    print ("Success: {} tests passed".format(len(calls)))
exit(rc)
