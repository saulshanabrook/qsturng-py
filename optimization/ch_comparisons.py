# Copyright (c) 2011, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""used to estimate RMSD compared to the C-H algorithm for parameter
optimization"""

import math
import time
import numpy as np
from qsturng import qsturng, v_keys
from make_tbls import R

t0 = time.time()
inf = float('inf')

def read_ch(fname):
    with open(fname) as f:
        lines = f.readlines()

    ps,qs = zip(*[L.split(',') for L in lines])

    return map(float, ps), map(float, qs)

ps = np.linspace(.5, .999, 100)
rms = 0.
##for r in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
##          22,23,24,25,26,27,28,29,30,35,40,50,60,70,80,90,100,200]:
i=0
nan_count = 0
hmm = 0
error = []
for r in [r for r in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
          25,30,35,40,50,60,70,80] if r not in R.keys()]:
    for v in [v for v in [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
              20,22,24,26,30,35,40,50,60,90,120,240,480,inf] if v not in v_keys]:
        ps_ch, qs_ch = read_ch('CH_r=%i,v=%.0f.dat'%(r,v))

        # There are cases where R's qtukey algorithm fails to
        # converge and results in NaN
        for ch, gl, p in zip(qs_ch[:-10],
                             [qsturng(p,r,v) for p in ps[:-10]],
                             ps[:-10]):

##            print ch, gl
            if not math.isnan(ch):
                
                rms += (ch-gl)**2
                ne = (gl-ch)/ch

                if abs(ne) > .05:
                    print ne,r,v, ch,p
        
                error.append(ne)
                
            else:
                nan_count += 1
        i += 1
##        print r,v,rms

rms = math.sqrt(rms / (i*100.))

print rms

print time.time() - t0
import pylab
import numpy as np
from random import shuffle
error = np.array(error)
sq_err = error**2

indices = range(len(error))
sem=[]
for i in range(100):
    shuffle(indices)
    sem.append(np.std(error[indices[:100]])/10.)

pylab.figure()
pylab.title(r'$\mathrm{Bootstrap} \; \mathrm{SEM}$')
pylab.hist(sem)
xticks=pylab.xticks()[0]
pylab.xticks(xticks,[r'$%.1e$'%t for t in xticks])
yticks=pylab.yticks()[0]
pylab.yticks(yticks,[r'$%i$'%t for t in yticks])
pylab.savefig('interp_bootstrap_sd.png')
pylab.close()

mu=[]
for i in range(100):
    shuffle(indices)
    mu.append(math.sqrt(np.sum(sq_err[indices[:100]])/100))

pylab.figure()
pylab.title(r'$\mathrm{Bootstrap} \; \mathrm{RMSD}$')
pylab.hist(sem)
xticks=pylab.xticks()[0]
pylab.xticks(xticks,[r'$%.1e$'%t for t in xticks])
yticks=pylab.yticks()[0]
pylab.yticks(yticks,[r'$%i$'%t for t in yticks])
pylab.savefig('interp_bootstrap_rms.png')
pylab.close()
