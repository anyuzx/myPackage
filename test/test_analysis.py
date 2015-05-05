import numpy as np
import sys
import scipy.spatial.distance as distance

from .. import core
from .. import analysis

# ==================================================
# test pdist and cdist module
a = np.random.rand(10,3)
dist = core.pdist(a)
dist_scipy = np.tril(distance.squareform(distance.pdist(a)))

if np.all(dist - dist_scipy < 0.0001):
	sys.stdout.write('**** Test core.pdist PASS ****\n')
else:
	sys.stdout.write('**** Failure: Test core.pdist NOT PASS ****\n')

b = np.random.rand(8,3)
dist = core.cdist(a,b)
dist_scipy = distance.cdist(a,b)

if np.all(dist - dist_scipy < 0.0001):
	sys.stdout.write('**** Test core.pdist PASS ****\n')
else:
	sys.stdout.write('**** Failure: Test core.cdist NOT PASS ****\n')

# ==================================================
# test subchaindist tool
a = np.random.rand(10,3)
dist = core.pdist(a)
test1 = analysis.polymertool.subchaindist(dist)

test2 = np.zeros(len(dist)-1)
for i in range(len(dist)-1):
	for j in range(i+1,len(dist)):
		test2[j-i-1] += dist[j,i]/(len(dist)-(j-i))

if np.all(test1 - test2 < 0.0001):
	sys.stdout.write('**** Test analysis.polymertool.subchaindist PASS ****\n')
else:
	sys.stdout.write('**** Failure: Test analysis.polymertool.subchaindist NOT PASS ****\n')