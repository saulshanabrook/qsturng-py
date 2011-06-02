# Copyright (c) 2011, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

"""this script uses basic evolutionary algorithm to find parameter estimates
for qsturng0 when p < .85. qsturng0 is used to increase the performance of
the scalar minimization used by psturng.

In the end, eyeball fitting p<.85 seemed to work out better."""

import csv
import getopt
import math
import pylab
import sys

from copy import deepcopy
from math import sqrt
import numpy as np
from numpy import array
from numpy import e
from random import choice
from numpy.random import random
from numpy.random import randint
from pprint import pprint
from string import Template
from time import time

from qsturng import _qsturng0ea, _qsturng, _qsturng0
from make_tbls import T, R

class ind:
    def __init__(self,func,p=30):
        self.f = func['f']
        self.xmin = func['xmin']
        self.xmax = func['xmax']
        self.range = func['xmax']-func['xmin']
        self.p=p
        
        # initialize with random initial states
        self.X = array([random()*self.range+self.xmin for i in xrange(5)]+[random()*100+100])
        self.set_fitness()

    def mutate(self, bump):
        for j in xrange(self.p):
            val=self.xmin-1.
##            while val<self.xmin or val>self.xmax:
            val = self.X[j] + (random()-.5)*2.*bump*self.range
            self.X[j]=val

        self.set_fitness()
            
    def crossover_2pt(self, other):
        new1,new2=deepcopy(self),deepcopy(other)
        
        [i1,i2]=sorted([randint(self.p),randint(self.p)]) # i1 <= i2
        new1.X[i1:i2],new2.X[i1:i2]=other.X[i1:i2],self.X[i1:i2]
        new1.set_fitness(); new2.set_fitness()

        return [new1,new2]

    def crossover_arithmetic(self, other, alpha=.35):
        new1,new2=deepcopy(self),deepcopy(other)

        new1.X = alpha*self.X + (1.-alpha)*other.X
        new2.X = alpha*other.X + (1.-alpha)*self.X        
        new1.set_fitness(); new2.set_fitness()

        return [new1,new2]
    
    def crossover_uniform(self, other, prob=.35):
        new1,new2=deepcopy(self),deepcopy(other)

        indices=[i for i in range(self.p) if random()<prob]
        new1.X[indices],new2.X[indices]=other.X[indices],self.X[indices]
        new1.set_fitness(); new2.set_fitness()

        return [new1,new2]
    
    def set_fitness(self):
        self.fitness=self.f(self.X)

    def __str__(self):
        return str(self.X)

    def __cmp__(self, other):
        return cmp(self.fitness,other.fitness)
    
def generational_GA(popsize,tournament_size,imax,func_dict,crossover,bump,k,p):
    # start clock
    start_time=time()

    f=func_dict['f']
    c=crossovers[crossover]
    
    # build pop
    inds = [ind(func_dict,p=p) for i in xrange(popsize)]
    print 'pop built'

    for i in xrange(imax):
        parent1=sorted([inds[randint(popsize)] for m in xrange(tournament_size)])[0]
        parent2=sorted([inds[randint(popsize)] for m in xrange(tournament_size)])[0]

        inds.extend(c(parent1,parent2))
        if random() < .8 : inds[-1].mutate(bump)
        if random() < .8 : inds[-2].mutate(bump)
        inds=sorted(inds)[:popsize]
##        print inds[0].fitness
        
    return {'func'    : f,
            'X'       : inds[0].X,
            'f(X)'    : inds[0].fitness,
            'runtime' : time() - start_time}

crossovers={
            '2point'        : ind.crossover_2pt,
            'uniform.2'     : lambda self, other : ind.crossover_uniform(self, other, .2),
            'uniform.35'    : lambda self, other : ind.crossover_uniform(self, other, .35),
            'uniform.5'     : lambda self, other : ind.crossover_uniform(self, other, .5),
            'arithmetic.2'  : lambda self, other : ind.crossover_arithmetic(self, other, .2),
            'arithmetic.35' : lambda self, other : ind.crossover_arithmetic(self, other, .35),
            'arithmetic.5'  : lambda self, other : ind.crossover_arithmetic(self, other, .5)
           }

inf = float('inf')
# p values that are defined in the A table
p_keys = [.1,.5,.675,.75,.8,.85,.9,.95,.975,.99,.995,.999]

# v values that are defined in the A table
v_keys = range(2, 21) + [24, 30, 40, 60, 120, inf]
print len(v_keys)

def err_func(c):
    # pick 100 random combinations of p, r, v

    e = 0.
    for i in xrange(30):
        p = .1
        r = choice(R.keys())
        v = choice(v_keys)
        
        actual = T[p][(v,1e38)[v==inf]][R[r]]
##        estimate = _qsturng0(p, r, v)
        estimate = _qsturng0ea(c, p, r, v)
        e += (abs(actual-estimate)/actual)**2.
        
    for i in xrange(70):
        p = choice([.5,.675,.75,.8,.85])
        r = choice(R.keys())
        v = choice(v_keys)
        
        actual = T[p][(v,1e38)[v==inf]][R[r]]
##        estimate = _qsturng0(p, r, v)
        estimate = _qsturng0ea(c, p, r, v)
        e += (abs(actual-estimate)/actual)**2.

    return e

##func_dict = { 'f': err_func,
##              'xmin':-5.12,
##              'xmax': 5.12}
##
##for i in xrange(8):
##    print generational_GA(200,5,10000,func_dict,'2point',.001,0,6)

from scipy.optimize import fminbound

print fminbound(full_output=1)
