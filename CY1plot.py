
## Chang & Yang numerical example 1 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

import numpy as np
import matplotlib.pyplot as plt

import pickle,pprint

pkl_file = open('CY1_CY01_4_methods.pkl', 'rb')

DictIn = pickle.load(pkl_file)
#pprint.pprint(DictIn)

pkl_file.close()

dims = DictIn['dims']
methods = DictIn['methods']
results = DictIn['results']
model = DictIn['model']

lm = len(methods)
II = len(dims)-1

print 'Model problem: %s' %(model)
print 'Dimensions: ',dims

# vvp L2 error, tot enrgy rel error, convergence order,

# print velocity L2 error on largest grid
print '\nVelocity L2 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,2,i]

# print vorticity L2 error on largest grid
print '\nVorticity L2 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,3,i]

# print pressure L2 error on largest grid
print '\nPressure L2 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,4,i]

# print velocity H1 error on largest grid
print '\nVelocity H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,7,i]

# print vorticity H1 error on largest grid
print '\nVorticity H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,8,i]

# print pressure H1 error on largest grid
print '\nPressure H1 error:'
for i,method in enumerate(methods):
	print method,"%.3e" % results[II,9,i]

# print total energy relative error on largest grid
print '\nTotal energy rel. error:'
for i,method in enumerate(methods):
	print method,"%.3g %%" % results[II,13,i]

# print residuals on largest grid
print '\nResiduals:'
for i,method in enumerate(methods):
	print method,"%.3e, %.3e, %.3e" % (results[II,18,i],results[II,19,i],results[II,20,i])

# print solver time on largest grid
print '\nSolver time:'
for i,method in enumerate(methods):
	print method,"%.1f" % results[II,21,i]

var = 2 # output variable of interest, 2=Velocity L2 error

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
	H1err = np.sqrt(pow(results[:,7,i],2)+pow(results[:,var,i],2))
	m,c = np.linalg.lstsq(A,np.log(H1err))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,1,0]=abs(m)
var = 3
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
		H1err = np.sqrt(pow(results[:,8,i],2)+pow(results[:,var,i],2))
		m,c = np.linalg.lstsq(A,np.log(H1err))[0]
		print method,"%.3g" % abs(m)
		arrConv[i,1,1]=abs(m)
var = 4
print '\nConvergence order for pressure (L^2 norm):'
for i,method in enumerate(methods):
	m,c = np.linalg.lstsq(A,np.log(results[:,var,i]))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,0,2]=abs(m)

print '\nConvergence order for pressure (H^1 norm):'
for i,method in enumerate(methods):
	# calc H1err from H1 semi-norm
	H1err = np.sqrt(pow(results[:,9,i],2)+pow(results[:,var,i],2))
	m,c = np.linalg.lstsq(A,np.log(H1err))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,1,2]=abs(m)

print '\nConvergence order for TERE:'
for i,method in enumerate(methods):
	m,c = np.linalg.lstsq(A,np.log(abs(results[:,13,i])))[0]
	print method,"%.3g" % abs(m)
	arrConv[i,0,3]=abs(m)

print '\nUnknowns, nonzero entries, sparsity:'
for i,method in enumerate(methods):
	print method,"%.3e, %.3e, %.3e" % (results[II,27,i],results[II,28,i],results[II,29,i])

# print total energy relative error on largest grid
print '\nTotal energy rel. error:'
for i,method in enumerate(methods):
	print method, results[:,13,i]

# Write output for Latex table (mixed stuff)
f = open('CY1_latex1.txt', 'w')
for i,method in enumerate(methods):
	s1 = [method,results[II,2,i],results[II,3,i],results[II,4,i],results[II,13,i],\
arrConv[i,0,0],arrConv[i,1,0],results[II,21,i]]
	s = "%s & %.3e & %.3e & %.3e & %.3g & %.3g & %.3g & %.1f\\\\\n" % \
(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6],s1[7])
	f.write(s)
f.close()

# Write output for Latex table (convergence rates)
f = open('CY1_latex2.txt', 'w')
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

hspace = 0.6

# Make loglog plot of var
var = 2 # Velocity L2 error
plt.figure(1)
plt.subplots_adjust(hspace=hspace)
a = plotDim
for i,method in enumerate(methods):
	a += 1
	plt.subplot(a)
	plt.loglog(dims, results[:,var,i])
	plt.xticks(dims,dims)
	plt.grid(True)
	plt.title(method)
	plt.xlim(dims[0],dims[II])
plt.savefig('CY1_loglog1.pdf')
#plt.show()

# Make loglog plot of var
var = 13 # Total energy relative error
var = 7 # H1 velocity error
plt.figure(2)
plt.subplots_adjust(hspace=hspace)
a = plotDim
for i,method in enumerate(methods):
	a += 1
	plt.subplot(a)
	plt.loglog(x, abs(results[:,var,i]))
	plt.xticks(dims,dims)
	plt.grid(True)
	plt.title(method)
	plt.xlim(dims[0],dims[II])
plt.savefig('CY1_loglog2.pdf')
#plt.show()

var = 21 # Time to solve
# make regular plot of var
plt.figure(3, figsize=(7, 11), dpi=80)
plt.subplots_adjust(hspace=0.3)
plt.subplot(311)
plt.title('Time to solve (sec)')
lines = plt.plot(x,results[:,21,:]) # time to solve
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'CY2', 'QNE1', 'QNE2', 'QNE1h', 'QNE2h'),
        'upper left')
plt.subplot(312)
plt.title('Unknowns')
lines = plt.plot(x,results[:,27,:]) # unknowns
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'CY2', 'QNE1', 'QNE2', 'QNE1h', 'QNE2h'),
        'upper left')
plt.subplot(313)
plt.ticklabel_format(axis='y',style='sci')
plt.title('Nonzero elements')
plt.plot(x,results[:,28,:]) # nonzero elements
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'CY2', 'QNE1', 'QNE2', 'QNE1h', 'QNE2h'),
        'upper left')
plt.savefig('CY1_plot1.pdf')
plt.show()

if 0:'''
var = 21 # Time to solve
# make regular plot of var
plt.figure(3, figsize=(7, 3), dpi=80)
plt.title('Time to solve (sec)')
lines = plt.plot(x,results[:,21,:]) # time to solve
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'CY2', 'QNE1', 'QNE2', 'QNE1h', 'QNE2h'),
        'upper left')
plt.savefig('CY1_plot1.pdf')
#plt.show()

plt.figure(4, figsize=(7, 7), dpi=80)
plt.subplots_adjust(hspace=0.3)
plt.subplot(211)
plt.title('Unknowns')
lines = plt.plot(x,results[:,27,:]) # unknowns
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.legend( lines,
        ('MF', 'CY1', 'CY2', 'QNE1', 'QNE2', 'QNE1h', 'QNE2h'),
        'upper left')
plt.subplot(212)
plt.ticklabel_format(axis='y',style='sci')
plt.title('Nonzero elements')
plt.plot(x,results[:,28,:]) # nonzero elements
plt.xticks(dims,dims)
plt.grid(True)
plt.xlim(dims[0],dims[II])
plt.savefig('CY1_plot2.pdf')
plt.show()
'''

