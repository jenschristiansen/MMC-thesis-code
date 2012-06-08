
### Dolfin initialization ###

from dolfin import *
import numpy as np

method = 'QNE2' # methods=['MF','CY1','QNE2','QNE1h']
if method<>'MF':
	useLS = 1
else:
	useLS = 0
usePrBC = 0 # prescribe p=0 at outlet

if method<>'MF':
	usePrBC = 0 # LS requires velocity BCs

# create mesh
dim = 100
mesh = UnitSquare(dim, dim) # domain is unit square [0,1]x[0,1]
h = 1.0/dim
h2 = h*h

# Define function spaces
if method=='CY1':
	ww = 1
	V1 = VectorFunctionSpace(mesh, "CG", 1) # velocity
elif method=='QNE1h':
	ww = h
	V1 = VectorFunctionSpace(mesh, "CG", 1) # velocity
elif method=='QNE2':
	ww = h2
	V1 = VectorFunctionSpace(mesh, "CG", 2) # velocity
if useLS:
	V2 = FunctionSpace(mesh, "CG", 1) # vorticity
	Q = FunctionSpace(mesh, "CG", 1) # pressure
	W = MixedFunctionSpace([V1,V2,Q]) # see p. 109 in FEniCS manual
else:
	V = VectorFunctionSpace(mesh, "CG", 2)
	Q = FunctionSpace(mesh, "CG", 1)
	W = V * Q

G = FunctionSpace(mesh,"DG",0)

unkn = len(Function(W).vector().array())
print 'Number of unknowns: ',unkn

# Boundaries
def topAndRight(x, on_boundary):
	return x[1] > 1-DOLFIN_EPS or x[0] > 1-DOLFIN_EPS
def left_hole(x, on_boundary):
	return x[0] < DOLFIN_EPS and (x[1] >= 0.7 and x[1] <= 0.9)
def left_wall(x, on_boundary):
	return x[0] < DOLFIN_EPS and (x[1] < 0.7 or x[1] > 0.9)
def bottom_hole(x, on_boundary):
	return x[1] < DOLFIN_EPS and (x[0] >= 0.7 and x[0] <= 0.9)
def bottom_wall(x, on_boundary):
	return x[1] < DOLFIN_EPS and (x[0] < 0.7 or x[0] > 0.9)
def corner(x, on_boundary): 
	return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0))
bc01 = DirichletBC(W.sub(0), noslip, topAndRight)
bc02 = DirichletBC(W.sub(0), noslip, left_wall)
bc03 = DirichletBC(W.sub(0), noslip, bottom_wall)

# Inflow / outflow BCs for velocity
inflow = Expression(("1.0-pow(10.0*x[1]-8.0,2)", "0.0"))
outflow = Expression(("0.0","-1.0+pow(10.0*x[0]-8.0,2)"))

bc11 = DirichletBC(W.sub(0), inflow, left_hole)
bc12 = DirichletBC(W.sub(0), outflow, bottom_hole)

# Boundary condition for pressure at outflow
zero = Constant(0.0)
bc21 = DirichletBC(W.sub(1), zero, bottom_hole)
bc22 = DirichletBC(W.sub(1), zero, corner, 'pointwise')

bc31 = DirichletBC(W.sub(1), zero, corner, "pointwise")

if useLS:
	bc31 = DirichletBC(W.sub(2), zero, corner, "pointwise")

# Collect boundary conditions
if usePrBC:
	bcs = [bc01, bc02, bc03, bc11, bc21] # pressure + velocity BCs
else:
	bcs = [bc01, bc02, bc03, bc11, bc12, bc22] # pin pressure, velocity BCs
	#bcs = [bc01, bc02, bc03, bc11, bc12] # plain velocity BCs
if useLS:
	bcs = [bc01, bc02, bc03, bc11, bc12, bc31] # pin pressure, velocity BCs

# Initialize piecewise constant coefficient vectors
alpha = Function(G)
grad_alpha = Function(G)
dummy = Function(G)
test = Function(G)
qqq = TestFunction(G)
N = len(alpha.vector().array())

# Define constants
mu = 1
alpha_low = 2.5*mu/100**2
alpha_high = 2.5*mu/0.01**2

# alpha >= 0 is inverse permeability: alpha close to zero => fluid flow
# rho \in [0,1]. rho = 0 => no flow, rho = 1 => free flow

# q = qval > 0 parameter determines the level of 'grey' in the optimal design
qval = 0.1

def alpha_q(rho):
	N = len(rho)
	alpha = np.zeros(N)
	for i in range(N):
		alpha[i] = alpha_high+(alpha_low-alpha_high)*rho[i]*(1+qval)/(rho[i]+qval)
	return alpha

def grad_alpha_q(rho):
	N = len(rho)
	grad_alpha = np.zeros(N)
	for i in range(N):
		grad_alpha[i] = -(-alpha_low+alpha_high)*(1+qval)*qval/(rho[i]+qval)**2
	return grad_alpha

# Define variational problem

f = Constant((0.0, 0.0))

if useLS:
	(u, w, p) = TrialFunctions(W) # w = omega
	(v, phi, q) = TestFunctions(W)
	L = ww*inner(f,curl(phi)+grad(q))*dx
	SSargs = [u,p,v,q,w,phi]
else:
	(u, p) = TrialFunctions(W)
	(v, q) = TestFunctions(W)
	L = inner(f,v)*dx
	SSargs = [u,p,v,q]

U = Function(W)

### End of Dolfin initialization ### 

### MMA initialization ###

#    The problem is assumed to be on the following form:
#
#    minimize  f_0(x) + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
#
#    subject to  f_i(x) - a_i*z - y_i <= fmax_i ,   i=1,..,M
#                       xmin_j <= x_j <= xmax_j ,   j=1,..,N
#                                 y_i >= 0 ,        i=1,..,M
#                                   z >= 0 .
	
# DEFINE PROBLEM DIMENSIONS & FUNCTIONS
# number of constraints
M = 1
# number of variables
GEPS = 1.0E-07
xmin = 0  *np.ones(N)
xmax = 1 *np.ones(N)

a    = 0.0  *np.ones(M)
c    = 1000.*np.ones(M)

# ALLOCATE OTHER STORAGE
xold1 =  np.zeros(N)
xold2 =  np.zeros(N)
xmma  =  np.zeros(N)
xlow  =  np.zeros(N) # lower bound for x
xupp  =  np.zeros(N) # upper bound for x
alfa  =  np.zeros(N)
beta  =  np.zeros(N)
df0dx =  np.zeros(N) # gradient of objective function
p0    =  np.zeros(N)
q0    =  np.zeros(N)

uu    =  np.zeros(M)
ulam  =  np.zeros(M)
b     =  np.zeros(M)
y     =  np.zeros(M)
gradf =  np.zeros(M)
dsrch =  np.zeros(M)
hessf =  np.zeros(M*(M+1)/2)
iyfree=  np.zeros(M,dtype=np.int32)

dfdx  =  np.zeros(N*M) # Jacobian of constraint function
p     =  np.zeros(N*M)
q     =  np.zeros(N*M)

z     =  0.

### End of MMA initialization ###


### Begin Stokes solver ###

def SolveStokes(rho,useLS,ssArgs):
	u = SSargs[0]
	p = SSargs[1]
	v = SSargs[2]
	q = SSargs[3]
	if useLS:
		w = SSargs[4]
		phi = SSargs[5]

	# build alpha(rho) coefficient vector
	alpha_vec = alpha_q(rho)
	alpha.vector()[:] = alpha_vec
	grad_alpha_vec = grad_alpha_q(rho)
	grad_alpha.vector()[:] = grad_alpha_vec

	# Update a() bilinear form
	if useLS:
		a=(ww*inner(curl(w)+grad(p)+alpha*u,curl(phi)+grad(q)+alpha*v) \
+inner(curl(u)-w,curl(v)-phi)+inner(div(u),div(v)))*dx
	else:
		a=alpha*inner(u,v)*dx+mu*inner(grad(u),grad(v))*dx-p*div(v)*dx-q*div(u)*dx

	# Solve
	from time import time
	t2=time()
	#solver.solve(U.vector(), bb)
	solve(a==L,U,bcs)
	t3=time()

	# Get sub-functions
	if useLS:
		u, w, p = U.split()
	else:
		u, p = U.split()
	
	# Calculate energy functional and its gradient (vector)
	energy = 0.5*(alpha*inner(u,u)+mu*inner(grad(u),grad(u)))*dx-inner(f,u)*dx
	E = assemble(energy)
	Lqq = 0.5 * grad_alpha * inner(u,u) * qqq *dx
	grad_E = assemble(Lqq)

	# Calculate fluid volume
	dummy.vector()[:] = rho
	fluid_vol = assemble(dummy*dx)
	dummy.vector()[:] = np.ones(N)
	dom_vol = assemble(dummy*dx)

	# Calculate constraint
	fval = fluid_vol/dom_vol # = constraint value
	dfdx = np.zeros(N)
	for c in cells(mesh):
		dfdx[c.index()] = c.volume()
	
	# Calculate residuals
	if useLS:
		r1 = assemble(ww*inner(curl(w)+grad(p)+alpha*u-f,curl(w)+grad(p)+alpha*u-f)*dx)
		r2 = assemble(inner(curl(u)-w,curl(u)-w)*dx)
		r3 = assemble(inner(div(u),div(u))*dx)

	#print 'Energy = ',E
	print 'Energy gradient: '
	print grad_E.array()
	#print 'Domain volume = ',dom_vol
	#print 'Fluid volume = ',fluid_vol
	#print 'Ratio = ',fval
	if useLS:
		print 'Residuals:',r1,r2,r3

	f0val = E
	
	df0dx = grad_E
	
	return f0val,fval,df0dx,dfdx,U

### End Stokes solver ###



