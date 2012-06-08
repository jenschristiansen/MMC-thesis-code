
## Chang & Yang numerical example 2 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

import CY2solver
import numpy as np
import matplotlib.pyplot as plt

from time import time

# xx1 = first order velocity
# xx2 = second order velocity

dims = [4,8,12,16,24]#,32]
dims = [4,8,16,20]
#dims = [2,4,8,16]
methods = ['MF','CY1','CY2','QNE1','QNE2','QNE1h','QNE2h']
methods = ['MF','CY1','QNE2','QNE1h']
#methods = ['MF','QNE2']

# output is 'results', a 3-dim array: dims x vars x methods, vars = output variables
t1 = time()
results = CY2solver.driver(dims,methods)
t2 = time()
print 'Time elapsed: %0.3f seconds.' % (t2-t1)

DictOut = {}
DictOut['dims']=dims
DictOut['methods']=methods
DictOut['results']=results

import pickle

output = open('CY2output004.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(DictOut, output)

output.close()

### Notes ###

# Saved files:

# CY2output001.npy ~17 min
# dims = [4,24]
# methods = ['CY2','QNE2','QNE1h']

# CY2output003.npy ~26 min
#dims = [4,8,12,16]
#methods = ['MF','CY1','CY2','QNE1','QNE2','QNE1h','QNE2h']


