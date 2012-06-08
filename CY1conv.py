
## Chang & Yang numerical example 1 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

import CY1solver
import numpy as np
import matplotlib.pyplot as plt

from time import time

# xx1 = first order velocity
# xx2 = second order velocity

dims = [10,20,40,80,160]
methods = ['MF','CY1','CY2','QNE1','QNE2','QNE1h','QNE2h']
#dims = [40]
methods = ['MF','CY1','QNE2','QNE1h']
methods = ['CY2','QNE1','QNE2h']

# output is 'results', a 3-dim array: dims x vars x methods, vars = output variables
t1 = time()
results = CY1solver.driver(dims,methods)
t2 = time()
print 'Time elapsed: %0.3f seconds.' % (t2-t1)

DictOut = {}
DictOut['dims']=dims
DictOut['methods']=methods
DictOut['results']=results
DictOut['model']=CY1solver.model

import pickle

output = open('CY1_BG01_3_bad_methods.pkl', 'wb') # set output file name here

# Pickle dictionary using protocol 0.
pickle.dump(DictOut, output)

# Pickle the list using the highest protocol available.
#pickle.dump(selfref_list, output, -1)

output.close()


### Notes ###

# Saved files:

# CY1output001.npy
# dims = [10,20,40,80,160]
# methods = ['MF','CY1','CY2','QNE1','QNE2','QNE1h','QNE2h']


