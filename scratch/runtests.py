from subprocess import call
calls =  [("clingoDL examples/taskassignment/encoding-dl.lp examples/taskassignment/tai4_4_1.lp -c n=132 0 --outf=3",30)]
for (i,r) in calls:
    return_code = call(i, shell=True)
    if (return_code != r):
        print("Failed: " + i + " with error code ", return_code)
