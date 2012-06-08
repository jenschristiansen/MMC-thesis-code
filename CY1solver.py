
## Chang & Yang numerical example 1 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

from dolfin import *
import numpy as np

cols = 34 # number of output variables

def driver(dims,methods):
	results = np.zeros((len(dims),cols,len(methods)))
	for j,method in enumerate(methods):
		output = np.zeros((len(dims),cols))
		for i,dim in enumerate(dims):
			output[i,:] = calc(method,dim)
			#if (i==(len(dims)-1) and j==(len(methods)-1)):
			#	it = 0
			#	for key,value in ED.iteritems():
			#		print it,key,value
			#		it += 1
		results[:,:,j] = output
	return results

## Shared stuff ##

model = "BG"

# True solution
if model=="CY":
	u_x = "pow(x[0],2)*pow(1.0-x[0],2)*(2.0*x[1]-6.0*pow(x[1],2)+4.0*pow(x[1],3))"
	u_y = "pow(x[1],2)*pow(1.0-x[1],2)*(-2.0*x[0]+6.0*pow(x[0],2)-4.0*pow(x[0],3))"
	u_t = Expression((u_x,u_y))
	w_t = Expression("pow(x[0],2)*pow(1.0-x[0],2)*(-2.0+12.0*x[1]-12.0*pow(x[1],2)) \
	+pow(x[1],2)*pow(1.0-x[1],2)*(-2.0+12.0*x[0]-12.0*pow(x[0],2))")
	p_t = Expression("pow(x[0],2)+pow(x[1],2)-20.0/3.0*x[0]*x[1]+x[0]+x[1]")
	f1 = "1.0 + 2.0*x[0] - 20.0/3.0*x[1] \
	-4.0*(1.0-6.0*x[0]+6.0*x[0]*x[0])*pow(-1.0+x[1],2)*x[1] \
	-4.0*(1.0-6.0*x[0]+6.0*x[0]*x[0])*(-1.0+x[1])*x[1]*x[1] \
	-12.0*pow(-1.0+x[0],2)*x[0]*x[0]*(-1.0+2.0*x[1])"
	f2 = "1.0 - 20.0/3.0*x[0] + 2.0*x[1] \
	+12.0*(-1.0+2.0*x[0])*pow(-1.0+x[1],2)*x[1]*x[1] \
	+4.0*(1.0-6.0*x[1]+6.0*x[1]*x[1])*pow(-1.0+x[0],2)*x[0] \
	+4.0*(1.0-6.0*x[1]+6.0*x[1]*x[1])*(-1.0+x[0])*x[0]*x[0]"
	f = Expression((f1,f2))
	Du00 = "0.2e1*x[0]*pow(0.1e1-x[0],0.2e1)*(0.2e1*x[1]-0.6e1*x[1]*x[1]+0.4e1*pow(x[1],0.3e1))-0.2e1*x[0]*x[0]*(0.1e1-x[0])*(0.2e1*x[1]-0.6e1*x[1]*x[1]+0.4e1*pow(x[1],0.3e1))"
	Du01 = "x[0]*x[0]*pow(0.1e1-x[0],0.2e1)*(0.2e1-0.12e2*x[1]+0.12e2*x[1]*x[1])"
	Du10 = "x[1]*x[1]*pow(0.1e1-x[1],0.2e1)*(-0.2e1+0.12e2*x[0]-0.12e2*x[0]*x[0])"
	Du11 = "0.2e1*x[1]*pow(0.1e1-x[1],0.2e1)*(-0.2e1*x[0]+0.6e1*x[0]*x[0]-0.4e1*pow(x[0],0.3e1))-0.2e1*x[1]*x[1]*(0.1e1-x[1])*(-0.2e1*x[0]+0.6e1*x[0]*x[0]-0.4e1*pow(x[0],0.3e1))"
	Du = Expression(((Du00,Du01),(Du10,Du11))) # true grad(u)
	Du1 = Expression((Du00,Du01))
	Du2 = Expression((Du10,Du11))
	Dw = Expression(("0.2e1*x[0]*pow(0.1e1-x[0],0.2e1)*(-0.2e1+0.12e2*x[1]-0.12e2*x[1]*x[1])-0.2e1*x[0]*x[0]*(0.1e1-x[0])*(-0.2e1+0.12e2*x[1]-0.12e2*x[1]*x[1])+x[1]*x[1]*pow(0.1e1-x[1],0.2e1)*(0.12e2-0.24e2*x[0])","x[0]*x[0]*pow(0.1e1-x[0],0.2e1)*(0.12e2-0.24e2*x[1])+0.2e1*x[1]*pow(0.1e1-x[1],0.2e1)*(-0.2e1+0.12e2*x[0]-0.12e2*x[0]*x[0])-0.2e1*x[1]*x[1]*(0.1e1-x[1])*(-0.2e1+0.12e2*x[0]-0.12e2*x[0]*x[0])"))
	Dp = Expression(("0.2e1*x[0]-0.20e2/0.3e1*x[1]+0.1e1","0.2e1*x[1]-0.20e2/0.3e1*x[0]+0.1e1"))
	
elif model=="BG":
	u_x = "exp(x[0]) * cos(x[1]) + sin(x[1])"
	u_y = "-exp(x[0]) * sin(x[1]) + 0.1e1 - pow(x[0], 0.3e1)"
	u_t = Expression((u_x,u_y))
	w_t = Expression("-0.3e1 * x[0] * x[0] - cos(x[1])")
	p_t = Expression("sin(x[1])*cos(x[0])+x[0]*x[1]*x[1]-0.1e1/0.6e1-sin(0.1e1)*(0.1e1-cos(0.1e1))")
	f1 = "sin(x[1])-sin(x[1])*sin(x[0])+x[1]*x[1]"
	f2 = "0.6e1*x[0]+cos(x[1])*cos(x[0])+0.2e1*x[0]*x[1]"
	f = Expression((f1,f2))
	Du00 = "exp(x[0])*cos(x[1])"
	Du01 = "-exp(x[0])*sin(x[1])+cos(x[1])"
	Du10 = "-exp(x[0])*sin(x[1])-0.3e1*x[0]*x[0]"
	Du11 = "-exp(x[0])*cos(x[1])"
	Du = Expression(((Du00,Du01),(Du10,Du11))) # true grad(u)
	Du1 = Expression((Du00,Du01))
	Du2 = Expression((Du10,Du11))
	Dw = Expression(("-0.6e1*x[0]","sin(x[1])"))
	Dp = Expression(("-sin(x[1])*sin(x[0])+x[1]*x[1]","cos(x[1])*cos(x[0])+0.2e1*x[0]*x[1]"))

# Boundaries
def boundary(x, on_boundary):
	return x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or x[0] > 1-DOLFIN_EPS \
or x[1] > 1-DOLFIN_EPS
def corner(x, on_boundary):
    return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS
noslip = Constant((0.0, 0.0))

def calc(method,dim):
	print 'Solving for method: ',method
	print 'Mesh size: ',dim
	from time import time
	t1 = time()

	ERR = np.zeros(cols)
	# create mesh
	mesh = UnitSquare(dim, dim) # domain is unit square [0,1]x[0,1]
	h=1.0/dim
	h2=h*h

	if method=='MF':
		V1 = VectorFunctionSpace(mesh, "CG", 2)
		Q = FunctionSpace(mesh, "CG", 1) 
		W = V1*Q
	elif (method=='CY1' or method=='QNE1' or method=='QNE1h'):
		V1 = VectorFunctionSpace(mesh, "CG", 1)
		V2 = FunctionSpace(mesh, "CG", 1)
		Q = FunctionSpace(mesh, "CG", 1)
		W = MixedFunctionSpace([V1,V2,Q])
	elif (method=='CY2' or method=='QNE2' or method=='QNE2h'):
		V1 = VectorFunctionSpace(mesh, "CG", 2)
		V2 = FunctionSpace(mesh, "CG", 1)
		Q = FunctionSpace(mesh, "CG", 1)
		W = MixedFunctionSpace([V1,V2,Q])
	
	if model=="CY":
		bc0 = DirichletBC(W.sub(0), noslip, boundary)
	elif model=="BG":
		bc0 = DirichletBC(W.sub(0), u_t, boundary)
		
	if method=='MF':
		bc1 = DirichletBC(W.sub(1), p_t, corner, "pointwise")
	else:
		bc1 = DirichletBC(W.sub(2), p_t, corner, "pointwise")
	bcs = [bc0, bc1]
	
	if method=='MF':
		(u, p) = TrialFunctions(W)
		(v, q) = TestFunctions(W)
	else:
		(u, w, p) = TrialFunctions(W)
		(v, phi, q) = TestFunctions(W)

	if method=='MF':
		a = (inner(grad(u),grad(v))-div(v)*p-q*div(u))*dx
		L = inner(f,v)*dx
	elif (method=='CY1' or method=='CY2'):
		a = (inner(curl(w)+grad(p),curl(phi)+grad(q))+inner(curl(u)-w,curl(v)-phi)+inner(div(u),div(v)))*dx
		L = inner(f,curl(phi)+grad(q))*dx
	elif (method=='QNE1' or method=='QNE2'):
		a = (h2*inner(curl(w)+grad(p),curl(phi)+grad(q))+inner(curl(u)-w,curl(v)-phi) \
+ inner(div(u),div(v)))*dx
		L = h2*inner(f,curl(phi)+grad(q))*dx
	elif (method=='QNE1h' or method=='QNE2h'):
		h = Circumradius(u.cell())
		a = (inner((curl(w)+grad(p)),h*(curl(phi)+grad(q)))+inner(curl(u)-w,curl(v)-phi) \
+ inner(div(u),div(v)))*dx  + inner(h*curl(w),h*curl(phi))*dx
		L = inner(f,h*(curl(phi)+grad(q)))*dx
	
	U = Function(W)
	
	t2=time()
	solve(a==L,U,bcs)
	t3=time()
	print 'Solver took %0.3f seconds.' % (t3-t2)

	if method=='MF':
		u,p = U.split()
	else:
		u,w,p = U.split()

	# Calculate errors
	if dim <= 120:
		Vp = VectorFunctionSpace(mesh,"CG",5) # memory intensive!
		Tp = TensorFunctionSpace(mesh,"CG",3)
		Qp = FunctionSpace(mesh,"CG",5) # memory intensive!
	else:
		Vp = VectorFunctionSpace(mesh,"CG",2)
		Tp = TensorFunctionSpace(mesh,"CG",2)
		Qp = FunctionSpace(mesh,"CG",2)
		
	u1=Expression(u_x)
	u2=Expression(u_y)
	u_tp = project(u_t,Vp)
	u1_tp = project(u1,Qp)
	u2_tp = project(u2,Qp)
	w_tp = project(w_t,Qp)
	p_tp = project(p_t,Qp)
	f_tp = project(f,Vp)
	Du_tp = project(Du,Tp)
	Du1_tp = project(Du1,Vp)
	Du2_tp = project(Du2,Vp)
	Dw_tp = project(Dw,Vp)
	Dp_tp = project(Dp,Vp)
	ED = {} # dictionary

	# L2 errors
	ERR[0]=ED['u1err']= sqrt(abs(assemble(inner(u[0]-u1,u[0]-u1)*dx)))
	ERR[1]=ED['u2err']= sqrt(abs(assemble(inner(u[1]-u2,u[1]-u2)*dx)))
	ERR[2]=ED['u_err']= sqrt(abs(assemble(inner(u-u_t,u-u_t)*dx)))
	print 'Vel. L2 err: ',ED['u_err']
	if method<>'MF':
		ERR[3]= ED['w_err'] = sqrt(assemble((w-w_t)**2*dx))
		print 'Vor. L2 err: ',ED['w_err']
	ERR[4]=ED['p_err']= sqrt(assemble((p-p_t)**2*dx))
	
	# H1 (semi-norm) errors
	ERR[5]=ED['u1errH1']= sqrt(assemble(inner(grad(u[0])-Du1_tp,grad(u[0])-Du1_tp)*dx))
	ERR[6]=ED['u2errH1']= sqrt(assemble(inner(grad(u[1])-Du2_tp,grad(u[1])-Du2_tp)*dx))
	ERR[7]=ED['u_errH1']= sqrt(assemble(inner(grad(u)-Du_tp,grad(u)-Du_tp)*dx))
	if method<>'MF':
		ERR[8]=ED['w_errH1']= sqrt(assemble(inner(grad(w)-Dw_tp,grad(w)-Dw_tp)*dx))	
	ERR[9]=ED['p_errH1']= sqrt(assemble(inner(grad(p)-Dp_tp,grad(p)-Dp_tp)*dx))

	# Infinity norm errors
	#print dir(u)
	#print u[0],dir(u[0])
	#print u.vector()
	#print u.vector().array()
	#print u.vector().array()[0]
	#print w.vector()
	#ED['u1errInf'] = abs(u.vector().array()[0]-project(u1,Q).vector().array()).max()
	#ED['u2errInf'] = abs(u.vector().array()[1]-project(u2,Q).vector().array()).max()
	#ED['w_errInf'] = abs(w.vector().array()-project(w_true,Q).vector().array()).max()
	#ED['p_errInf'] = abs(p.vector().array()-project(p_true,Q).vector().array()).max()
	#raw_input("Press Enter to continue...")
	#quit()

	# Energy
	ERR[10]=ED['E'] = assemble(0.5*inner(grad(u),grad(u))*dx-inner(f,u)*dx)
	ERR[11]=ED['E_tr'] = assemble(0.5*inner(Du_tp,Du_tp)*dx-inner(f,u_t)*dx)
	ERR[12]=ED['E_err'] = ED['E']-ED['E_tr']
	ERR[13]=ED['E_re'] = 100*ED['E_err']/ED['E_tr']

	ERR[30]=ED['E2'] = ED['E'] + assemble(0.5*inner(u,u)*dx)
	ERR[31]=ED['E_tr2'] = ED['E_tr'] + assemble(0.5*inner(u_tp,u_tp)*dx)
	ERR[32]=ED['E_err2'] = ED['E2']-ED['E_tr2']
	ERR[33]=ED['E_re2'] = 100*ED['E_err2']/ED['E_tr2']

	ERR[14]=ED['F'] = assemble(inner(f_tp,f_tp)*dx) # total force

	ERR[15]=ED['gradU_tot'] = sqrt(assemble(inner(Du_tp,Du_tp)*dx))
	ERR[16]=ED['gradU_err'] = sqrt(assemble(inner(grad(u)-Du_tp,grad(u)-Du_tp)*dx))
	ERR[17]=ED['gradU_re'] = 100*ED['gradU_err']/ED['gradU_tot']

	# Residuals
	if method<>'MF':
		ERR[18]=ED['r1']=assemble(inner(curl(w)+grad(p)-f,curl(w)+grad(p)-f)*dx)
		ERR[19]=ED['r2']=assemble(inner(curl(u)-w,curl(u)-w)*dx)
		ERR[20]=ED['r3']=assemble(inner(div(u),div(u))*dx)

	t4 = time()
	print 'Error calculation took %0.3f seconds.' % (t4-t3)
	
	cnt = 0
	t5 = time()
	#A, bb = assemble_system(a, L, bcs)
	#N = A.size(0)
	#for i in range(N):
	#	row = A.getrow(i)
	#	row = row[1] # get second array, the one with actual values
	#	cnt += sum(abs(row)>DOLFIN_EPS)
	t6 = time()

	ERR[21]=ED['tSolve']= t3-t2 # time to solve
	ERR[22]=ED['tSetupSolve']= t3-t1 # time to setup and solve
	ERR[23]=ED['tErrors']= t4-t3 # time to calc errors
	ERR[24]=ED['tSparsity']= t6-t5 # time to calc sparsity
	ERR[25]=ED['tTotal']= t6-t1 # total time elapsed
	
	ERR[26]=ED['dim']=dim
	ERR[27]=ED['unknowns']=len(U.vector().array())
	#ERR[28]=ED['nonzero']=cnt
	#ERR[29]=ED['sparsity']=float(cnt)/float(N*N)
	
	return ERR













