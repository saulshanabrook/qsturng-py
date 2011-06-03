# Copyright (c) 2011, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""this script makes the Figure 1. on the google code website. It provides a
gross visual verification of the algorithm"""

import numpy as np
##from qsturng import qsturng,p_keys,v_keys
from qsturng import p_keys,v_keys
from qsturng import qsturng
from qsturng.make_tbls import T,R
import pylab

inf = float('inf')

vs= [2,18,inf]
ps = np.linspace(.1, .999, 500)
pylab.figure(figsize=(9, 12))
pylab.subplots_adjust(bottom=.08, left=.07, right=.96, top=.97, hspace=.09)
for i, v in enumerate(vs):
    pylab.subplot(3, 1, i+1)
    for r in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,60,80,100]:
        
        print v,r
        qs = [qsturng(p, r, v) for p in ps]
        if r in [20, 30]:
            if r == 20:
                pylab.plot(ps, qs, 'k', alpha=.7,
                          label=r'$p \, \mathrm{interpolated}$')
            else:
                pylab.plot(ps, qs, 'k', alpha=.7)
        else:
            pylab.plot(ps, qs, 'k', alpha=.3)

        if r == 2:
            pylab.scatter([p for p in p_keys],
                          [T[p][(v,1e38)[v==inf]][R[r]] for p in p_keys],
                          marker='x',
                          label=r'$\mathrm{from} \, T \, \mathrm{table}$')
        else:
            pylab.scatter([p for p in p_keys],
                          [T[p][(v,1e38)[v==inf]][R[r]] for p in p_keys],
                          marker='x')
    for r in range(21,31):
        qs = [qsturng(p, r, v) for p in ps]
        if r == 21:
            pylab.plot(ps, qs, 'b', alpha=.4,
                      label=r'$p \, \mathrm{and} \, r \, \mathrm{interpolated}$')
        else:
            pylab.plot(ps, qs, 'b', alpha=.4)
                
    pylab.xlim([.09,1.01])
    print pylab.ylim()
    pylab.ylim([0.,pylab.ylim()[1]])
    yticks = pylab.yticks()[0]
    pylab.yticks(yticks, ['$%i$'%t for t in yticks])
    
    if i == 0:
        pylab.legend(loc=2)
        pylab.xticks([x for x in p_keys],
                     ['' for x in p_keys])
        pylab.title(r'$v = %i$'%v)
    elif i == 1:
        pylab.ylabel(r'$\mathrm{Studentized \, Range} \, (q) \, \mathrm{Approximation}$',
                     fontsize=14)
        pylab.xticks([x for x in p_keys],
                     ['' for x in p_keys])
        pylab.text(.55, .75, '$r \, = \, 2$')
        pylab.text(.55, 3.5, '$r \, = \, 20$')
        pylab.text(.55, 4.5, '$r \, = \, 30$')
        pylab.text(.55, 5.5, '$r \, = \, 100$')
        pylab.title(r'$v = %i$'%v)
    else:
        pylab.xlabel(r'$p$', fontsize=14)
        pylab.xticks([x for x in p_keys],
                     [('$%.3f$'%x, '')[x in [.990,.995]] for x in p_keys],
                     rotation='vertical')
        pylab.text(.55, .75, '$r \, = \, 2$')
        pylab.text(.55, 3.5, '$r \, = \, 20$')
        pylab.text(.55, 4.5, '$r \, = \, 30$')
        pylab.text(.55, 5.5, '$r \, = \, 100$')
        pylab.title(r'$v = \infty$')
            
pylab.savefig('qsturng,dpi=100.png', dpi=100)
pylab.savefig('qsturng,dpi=300.png', dpi=300)
pylab.close()
