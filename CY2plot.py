
## Chang & Yang numerical example 2 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

import numpy as np
import matplotlib.pyplot as plt

import pickle,pprint

pkl_file = open('CY2output004.pkl', 'rb')

DictIn = pickle.load(pkl_file)
#pprint.pprint(DictIn)

pkl_file.close()

dims = DictIn['dims']
methods = DictIn['methods']
results = DictIn['results']

lm = len(methods)
II = len(dims)-1

# print velocity L2 error on largest grid
print '\nVelocity L2 error:'
for i,method in enumerate(methods):
	print method,results[II,3,i]

# print vorticity L2 error on largest grid
print '\nVorticity L2 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,4,i]

# print pressure L2 error on largest grid
print '\nPressure L2 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,5,i]

# print velocity H1 error on largest grid
print '\nVelocity H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,9,i]

# print vorticity H1 error on largest grid
print '\nVorticity H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,10,i]

# print pressure H1 error on largest grid
print '\nPressure H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,11,i]

# print total energy relative error on largest grid
print '\nTotal energy rel. error:'
for i,method in enumerate(methods):
	print method,results[II,15,i]

# print residuals on largest grid
print '\nResiduals:'
for i,method in enumerate(methods):
	print method,results[II,20:23,i]

# print solver time on largest grid
print '\nSolver time:'
for i,method in enumerate(methods):
	print method,results[II,23,i]

var = 3 # output variable of interest, 3=Velocity L2 error

x = np.array(dims)
arrConv = np.zeros([lm,2,4]) # arrConv = methods x L2/H1 x velocity/vorticity/pressure/TERE
# Calc convergence order
A = np.vstack([np.log(x), np.ones(len(x))]).T
print '\nConvergence order for velocity (L^2 norm):'
for i,method in enumerate(methods):
	m,c = np.linalg.lstsq(A,np.log(results[:,var,i]))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,0,0]=abs(m)

print '\nConvergence order for velocity (H^1 norm):'
for i,method in enumerate(methods):
	# calc H1err from H1 semi-norm
	H1err = np.sqrt(pow(results[:,9,i],2)+pow(results[:,var,i],2))
	m,c = np.linalg.lstsq(A,np.log(H1err))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,1,0]=abs(m)
var = 4
print '\nConvergence order for vorticity (L^2 norm):'
for i,method in enumerate(methods):
	if method<>"MF":
		m,c = np.linalg.lstsq(A,np.log(results[:,var,i]))[0]
		print method,"%.3g" % abs(m)
		arrConv[i,0,1]=abs(m)

print '\nConvergence order for vorticity (H^1 norm):'
for i,method in enumerate(methods):
	# calc H1err from H1 semi-norm
	if method<>"MF":
		H1err = np.sqrt(pow(results[:,10,i],2)+pow(results[:,var,i],2))
		m,c = np.linalg.lstsq(A,np.log(H1err))[0]
		print method,"%.3g" % abs(m)
		arrConv[i,1,1]=abs(m)
var = 5
print '\nConvergence order for pressure (L^2 norm):'
for i,method in enumerate(methods):
	m,c = np.linalg.lstsq(A,np.log(abs(results[:,var,i])))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,0,2]=abs(m)

print '\nConvergence order for pressure (H^1 norm):'
for i,method in enumerate(methods):
	# calc H1err from H1 semi-norm
	H1err = np.sqrt(pow(results[:,11,i],2)+pow(results[:,var,i],2))
	m,c = np.linalg.lstsq(A,np.log(H1err))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,1,2]=abs(m)

print '\nConvergence order for TERE:'
for i,method in enumerate(methods):
	m,c = np.linalg.lstsq(A,np.log(abs(results[:,15,i])))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,0,3]=abs(m)

print '\nUnknowns, nonzero entries, sparsity:'
for i,method in enumerate(methods):
	print method,"%.3e, %.3e, %.3e" % (results[II,29,i],results[II,30,i],results[II,31,i])

print '\nSolver iterations:'
for i,method in enumerate(methods):
	print method,results[II,32,i] # returns three values with index 26,27,28...

# Write output for Latex table
f = open('CY2_latex1.txt', 'w')
for i,method in enumerate(methods):
	s1 = [method,results[II,3,i],results[II,4,i],results[II,5,i],results[II,15,i],\
arrConv[i,0,0],arrConv[i,1,0],results[II,23,i]]
	s = "%s & %.3e & %.3e & %.3e & %.3g & %.3g & %.3g & %.1f\\\\\n" % \
(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6],s1[7])
	f.write(s)
f.close()

# Write output for Latex table (convergence rates)
f = open('CY2_latex2.txt', 'w')
for i,method in enumerate(methods):
	s1 = [method,arrConv[i,0,0],arrConv[i,1,0],arrConv[i,0,1],arrConv[i,1,1],arrConv[i,0,2],arrConv[i,1,2],arrConv[i,0,3]]
	s = "%s & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g \\\\\n" % \
(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6],s1[7])
	f.write(s)
f.close()

if lm<=2:
	plotDim = lm*100+10
elif lm<=4:
	plotDim = 220
elif lm<=6:
	plotDim = 320
elif lm<=8:
	plotDim = 420
elif lm==9:
	plotDim = 330
elif lm<=12:
	plotDim = 430

# Make loglog plot of var
var = 3 # Velocity L2 error
plt.figure(1)
plt.subplots_adjust(hspace=0.4)
a = plotDim
for i,method in enumerate(methods):
	a += 1
	plt.subplot(a)
	plt.loglog(dims, results[:,var,i])
	plt.xticks(dims,dims)
	plt.grid(True)
	plt.title(method)
	plt.xlim(dims[0],dims[II])

plt.savefig('CY2_loglog1.pdf')
#plt.show()

# Make loglog plot of var
var = 15 # Total energy relative error
plt.figure(2)
plt.subplots_adjust(hspace=0.4)
a = plotDim
for i,method in enumerate(methods):
	a += 1
	plt.subplot(a)
	plt.loglog(x, abs(results[:,var,i]))
	plt.xticks(dims,dims)
	plt.grid(True)
	plt.title(method)
	plt.xlim(dims[0],dims[II])

plt.savefig('CY2_loglog2.pdf')
#plt.show()

if 0:'''
var = 23 # Time to solve
# make regular plot of var
plt.figure(3)
plt.subplots_adjust(hspace=0.4)
a = plotDim
for i,method in enumerate(methods):
	a += 1
	plt.subplot(a)
	plt.plot(x, results[:,var,i])
	plt.xticks(dims,dims)
	plt.grid(True)
	plt.title(method)
	plt.xlim(dims[0],dims[II])

plt.savefig('CY2_plot1.pdf')
#plt.show()
'''

plt.figure(4)
plt.subplots_adjust(hspace=0.4)
plt.subplot(221)
plt.title('Time to solve')
lines = plt.plot(x,results[:,23,:]) # time to solve
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'QNE2', 'QNE1h'),
        'upper left')
plt.subplot(222)
plt.title('Solver iterations')
plt.plot(x,results[:,32,:]) # krylov iterations
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.subplot(223)
plt.title('Unknowns')
plt.plot(x,results[:,29,:]) # unknowns
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.subplot(224)
plt.title('Nonzero elements')
plt.plot(x,results[:,30,:]) # nonzero elements
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.savefig('CY2_plot1.pdf')
plt.show()



