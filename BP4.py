
## Replication of Borrvall & Petersson (2003) numerical examples

# Created: 	15/04/12
# Last edited: 	27/05/12
# Author: 	Jens V. Christiansen

# Solution procedure

if 0: '''
0) Make an initial guess for the design function \alpha(\rho)
1) Solve Stokes equations in some domain using \alpha(\rho) 
2) Compute energy functional J_\alpha (u)
3) Find variation of energy functional \delta J
4) Update \alpha(\rho) using some algorithm (MMA)
5) Compute some measure of the quality of the solution
6) If the solution is still not good enough, go to step 1. Otherwise stop.
'''

#----------------------------------------#
# Numerical example 4: Double slit	 #
#----------------------------------------#

import numpy as np
import ksmma
from BP4import import *

from time import time
t1 = time()

# Set parameters (set mesh dim and solver method in BP1InitDolfin) #

mma_tol = 0.05 # stop MMA when objective function changes less than this (in %)
max_ite = 125
gamma = 1/3.0 # allowed fluid volume

### Call PDE solver ###

# set initial guess for rho vector
rho = gamma*np.ones(N)
xval = rho
#xval = np.load('BP4iniGuess.npy') 


print 'Number of elements: ',N

from time import time
t2 = time()
f0val,fval,df0dx,dfdx,U,qval = SolveStokes(xval,useLS,SSargs,0)
t3 = time()
f0ini = f0val

# Save initial solution for plot
U0 = Function(W)
U0.vector()[:] = U.vector()[:]

### Call MMA ###

fmax = gamma*np.ones(M)

print fval
print 'Initial guess:'
print 'x-values:'
print xval
print 'Objective function f0=%.4f' % (f0val)
print 'Constraint value: %.4f' % (fval)
print 'Time taken to solve PDE: %.3f seconds.\n' % (t3-t2) 

# optimization iteration
it = 0
cont = 1
f0old = f0val

info = np.zeros([200,6])

#for k in range(iterations) :
while cont:
	it += 1
	# make MMA step
	t2 = time()
	ksmma.mmasub(it,M,N,GEPS,iyfree,xval,xmma,\
		 xmin,xmax,xold1,xold2,xlow,xupp, \
		 alfa,beta,a,b,c,y,z,ulam, \
		 f0val,fval,fmax,df0dx,dfdx, \
		 p,q,p0,q0,uu,gradf,dsrch,hessf)
	ksmma.xupdat(N,it,xmma,xval,xold1,xold2)
	t3 = time()

	# evaluate everything again
	t4 = time()
	f0val,fval,df0dx,dfdx,U,qval = SolveStokes(xval,useLS,SSargs,it)
	t5 = time()

	f0chg = f0val-f0old
	f0chg_pc = 100*(f0val-f0old)/f0old
	f0old = f0val
	vio   = 100*(fval-gamma)/gamma
	print 'Iteration: %.2g' % it
	print 'q-value: %.3g' % qval
	print 'x-values: ', xval
	print 'Objective function f0=%.4f' % f0val
	print 'Change in f0 = %.4f or %.4f%%' % (f0chg,f0chg_pc)
	print 'Constraint value: %.4f' % fval
	print 'Constraint violation: %.3f %%' % vio
	print 'Time taken to run MMA: %.3f' % (t3-t2)
	print 'Time taken to solve PDE: %.3f\n' % (t5-t4)
	
	info[it,:] = [qval,f0val,f0chg,f0chg_pc,fval,vio]
	if (it>50 and abs(f0chg_pc)<mma_tol and abs(vio)<mma_tol) or it == max_ite:
		cont = 0

### Write to log ###
from datetime import datetime
now = datetime.now()
f = open('../BP_log.txt', 'a')
s1 = "Run finished at: %s\nProblem: BP4\nMethod: %s\ndim: %.3g\nTime taken: %.3g min %.3f sec\nIterations: %.3g\nUnknowns: %.6g\n" % \
(now.strftime("%Y-%m-%d %H:%M:%S"),method,dim,int((t5-t1)/60),(t5-t1) % 60,it,unkn)
s2 = "Gamma: %.3g\ndelta: %.3g\nmma_tol: %.3g\nmax_ite: %.3g\n" % (gamma,delta,mma_tol,max_ite)
s3 = "Obj.func. f0\n    initial: %.4g\n    final: %.4g\n    change: %.3g %%\nConstraint\n    fval: %.4g\n    violation: %.4g %%\n\n" % (f0ini,f0val,100*(f0val-f0ini)/f0ini,fval,vio) 
f.write(s1)
f.write(s2)
f.write(s3)
f.close()

### write iteration info to file ###
filename = 'BP4info'+now.strftime("%Y-%m-%d_%H-%M-%S")+'.txt'
f = open(filename,'w')
f.write("iteration,qval,f0val,f0chg,f0chg_pc,fval,vio\n")
j = 1
while info[j,0]<>0:
	s1 = "%.3g," % j
	for i in range(0,6):
		s1 += "%.4g," % info[j,i]
	s1 += "\n"
	f.write(s1)
	j += 1
f.close()

### Display results ###

print 'Solution done.\nTotal time taken: %.3f seconds on a %.3g x %.3g grid.\n' % (t5-t1,int(delta*dim),dim)
print 'Solver: %s, mma_tol: %.3g, max_ite: %.3g, gamma: %.3g' % (method,mma_tol,max_ite,gamma)

# Get sub-functions
if useLS:
	u, w, p = U.split()
	u0, w0, p0 = U0.split()
else:
	u, p = U.split()
	u0, p0 = U0.split()

rho = Function(G)
rho.vector()[:] = xval

#np.save('BP4iniGuess', xval)

# Save solution in VTK format
ufile_pvd = File("results/BP4_velocity.pvd")
ufile_pvd << u
pfile_pvd = File("results/BP4_pressure.pvd")
pfile_pvd << p
rhofile_pvd = File("results/BP4_rho.pvd")
rhofile_pvd << rho

if useLS:
	wfile_pvd = File("results/BP4_vorticity.pvd")
	wfile_pvd << w

# Plot solution
plot(u,mesh = mesh,wireframe = True,axes = True)
if useLS:
	plot(w,mesh = mesh,wireframe = True,axes = True)
plot(p,mesh = mesh,wireframe = True,axes = True)
plot(rho,mesh = mesh,wireframe = True,axes = True)
#plot(diff,mesh = mesh,wireframe = True,axes = True)
#plot(u0,mesh = mesh,wireframe = True,axes = True)
#plot(p0,mesh = mesh,wireframe = True,axes = True)
#plot(mesh,axes = True)
interactive()










