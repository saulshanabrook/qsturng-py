# Copyright (c) 2011 BSD, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""compares values obtained by qsturng to the random values obtained by R's
qtukey function"""

import math
import time
import numpy as np
import pylab

from random import shuffle
from qsturng import qsturng, v_keys
from qsturng.make_tbls import R

def read_ch(fname):
    with open(fname) as f:
        lines = f.readlines()
    ps,rs,vs,qs = zip(*[L.split(',') for L in lines])
    return map(float, ps), map(float, rs),map(float, vs), map(float, qs)

ps, rs, vs, actuals = read_ch('../tests/bootleg.dat')
estimates = qsturng(ps, rs, vs)
actuals = np.array(actuals)
estimates = np.array(estimates)
errors = (estimates - actuals) / actuals
indices = range(1000)

shuffle(indices)

from pprint import pprint
pprint([(ps[i],rs[i],vs[i],qsturng(ps[i],rs[i],vs[i])) for i in indices[:100]])


means, medians, sems = [], [], []
for i in xrange(100):
    print i
    shuffle(indices)
    means.append(np.mean(errors[indices[:100]]))
    medians.append(np.median(errors[indices[:100]]))
    sems.append(np.std(errors[indices[:100]])/10.)
    
pylab.figure()
pylab.hist(errors, bins=100)
yticks = pylab.yticks()[0]
pylab.yticks(yticks, [r'$%i$'%t for t in yticks])
##pylab.xlim([-.004, .004])
xticks = pylab.xticks()[0]
pylab.xticks(xticks, [r'$%.3f$'%t for t in xticks])
pylab.text(0, 475,
           r'$\frac{\mathrm{qsturng}(p,r,v)-'
           '\mathrm{qtukey}(p,r,v)}'
           '{\mathrm{qtukey}(p,r,v)}$',
           fontsize=16, horizontalalignment='center')
pylab.savefig('non-tbl.png')
pylab.close()

##pylab.figure()
##pylab.hist(errors, bins=100)
##pylab.ylim([0.,20.])
##yticks = pylab.yticks()[0]
##pylab.yticks(yticks, [r'$%i$'%t for t in yticks])
##xticks = pylab.xticks()[0]
##pylab.xticks(xticks, [r'$%.3f$'%t for t in xticks])
##pylab.text(0, 1900, r'$\frac{\hat{q_p}(r,v)-q_p(r,v)}{q_p(r,v)}$',
##           fontsize=16, horizontalalignment='center')
##pylab.savefig('non-tbl_detail.png')
##pylab.close()

pylab.figure()
pylab.hist(means, bins=10)
pylab.savefig('bootstrap_means.png')
pylab.close()

pylab.figure()
pylab.hist(medians, bins=10)
pylab.savefig('bootstrap_medians.png')
pylab.close()

pylab.figure()
pylab.hist(sems, bins=10)
pylab.savefig('bootstrap_sems.png')
pylab.close()
