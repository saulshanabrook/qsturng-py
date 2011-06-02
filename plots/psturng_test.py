# Copyright (c) 2011 BSD, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""generates q-values using qsturng, and tests to make psturng finds
the right p-value"""

import math
import time
import numpy as np
import pylab

from qsturng import qsturng, psturng, v_keys
from qsturng.make_tbls import R

n = 100
ps = np.random.random(n)*(.999 - .1) + .1
rs = np.random.random_integers(2, 100, n)
vs = np.random.random(n)*998. + 2.
qs = qsturng(ps, rs, vs)
t0=time.time()
estimates = psturng(qs, rs, vs)
import pprint
pprint.pprint([(p,r,v,q) for p,r,v,q in zip(ps,rs,vs,qs)])

print time.time()-t0

actuals = 1. - ps
errors = estimates - actuals

pylab.figure()
pylab.hist(errors, bins=100)
yticks = pylab.yticks()[0]
pylab.yticks(yticks, [r'$%i$'%t for t in yticks])
xticks = pylab.xticks()[0]
pylab.xticks(xticks, [r'$%.0e$'%t for t in xticks])
pylab.text(0, 475,
           r'$\frac{\mathrm{psturng}(q,r,v)-p}{p}$',
           fontsize=16, horizontalalignment='center')
pylab.savefig('non-tbl_psturng_.png')
pylab.close()
