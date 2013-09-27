#! /usr/bin/env python

import sys
from time import time
from random import random

from UnionFind import *

N=int(sys.argv[1])

a=UF(N)

t1=time()
for i in range(N):
    for j in range(i+1,N):
        if random()>0.6:
            a.merge(i,j)
        if j>i+N/1000: break
    if i%10000==0:
        print >> sys.stderr, "Mergind for (%d/%d)\r"%(i,N),

print >> sys.stderr,"Cost: %.2f.              "%(time()-t1)
print >> sys.stderr,"Clusters:",a.count()
