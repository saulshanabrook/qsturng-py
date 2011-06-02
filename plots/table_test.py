# Copyright (c) 2011 BSD, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""compares values obtained by qsturng to the values in the T table"""

from qsturng import qsturng,p_keys,v_keys
from qsturng.make_tbls import T,R
from collections import Counter

##error = []
##n=0
##for p in T:
##    for v in T[p]:
##        for r in R.keys():
##            
##            actual = T[p][v][R[r]]
##            estimate = qsturng(p, r, v)
##            error.append((estimate - actual)/actual)
##            n+=1
##
##print n

ps,vs,rs = [],[],[]
n=0
for p in T:
    for v in T[p]:
        for r in R.keys():
            ps.append(p)
            vs.append(v)
            rs.append(r)
            n += 1


indices = range(n)

from random import shuffle
shuffle(indices)
for i in indices[:100]:
    print ps[i],rs[i],vs[i],qsturng(ps[i],rs[i],vs[i])

            
import pylab

pylab.figure(figsize=(6,8))
pylab.subplots_adjust(top=.9, hspace=.07)

pylab.subplot(2,1,1)
pylab.hist(error, bins=100)
pylab.ylim([70.,7000.])
xticks=pylab.xticks()[0]
pylab.xticks(xticks, ['' for t in xticks])
yticks=pylab.yticks()[0]
pylab.yticks(yticks, [r'$%i$'%t for t in yticks])

pylab.text(0,7500,r'$\frac{\hat{q_p}(r,v)-q_p(r,v)}{q_p(r,v)}$',
           fontsize=16, horizontalalignment='center')

pylab.subplot(2,1,2)
pylab.hist(error, bins=100)
pylab.ylim([0.,70.])
xticks=pylab.xticks()[0]
pylab.xticks(xticks, [r'$%.2f$'%t for t in xticks])

yticks=pylab.yticks()[0]
pylab.yticks(yticks, [r'$%i$'%t for t in yticks])

pylab.savefig('error.png')
pylab.close()

import numpy as np
from random import shuffle

error = np.array(error)
indices = range(n)

sem=[]
for i in range(100):
    shuffle(indices)
    sem.append(np.std(error[indices[:100]])/10.)

pylab.figure()
pylab.title(r'$\mathrm{Bootstrap} \; \mathrm{SEM}$')
##pylab.text(2.5e-4, 16, r'$\mathrm{100 samples} \\ \mathrm{100 samplings}$')
pylab.hist(sem)
xticks=pylab.xticks()[0]
pylab.xticks(xticks,[r'$%.1e$'%t for t in xticks])
yticks=pylab.yticks()[0]
pylab.yticks(yticks,[r'$%i$'%t for t in yticks])
pylab.savefig('bootstrap_sd.png')
pylab.close()
