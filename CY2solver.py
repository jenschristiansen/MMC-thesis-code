
## Chang & Yang numerical example 2 ##

# Convergence analysis

# Created on April 18, 2012
# Last revised on May 30, 2012
# Author: Jens V. Christiansen

#------------------------------------#

# Anton:
# You can do even further speedups, by not reallocating/reassembling the sparsity pattern after you have # assembled the system once.
# Basically look at p. 98 in Fenics book excerpt: 
#
# assemble(a, tensor=A,reset_sparsity=false)
#

from dolfin import *
import numpy as np

cols = 37 # number of output variables

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

# t solution
u_x = "2*x[0]*x[0]*pow((1-x[0]),2)*x[1]*pow((1-x[1]),2)*x[2]*(1-x[2])- \
2*x[0]*x[0]*pow((1-x[0]),2)*x[1]*x[1]*(1-x[1])*x[2]*(1-x[2])- \
2*x[0]*x[0]*pow((1-x[0]),2)*x[1]*(1-x[1])*x[2]*pow((1-x[2]),2)+ \
2*x[0]*x[0]*pow((1-x[0]),2)*x[1]*(1-x[1])*x[2]*x[2]*(1-x[2])"
u_y = "2*x[0]*(1-x[0])*x[1]*x[1]*pow((1-x[1]),2)*x[2]*pow((1-x[2]),2) \
-2*x[0]*(1-x[0])*x[1]*x[1]*pow((1-x[1]),2)*x[2]*x[2]*(1-x[2])-2*x[0]* \
pow((1-x[0]),2)*x[1]*x[1]*pow((1-x[1]),2)*x[2]*(1-x[2])+2*x[0]*x[0]* \
(1-x[0])*x[1]*x[1]*pow((1-x[1]),2)*x[2]*(1-x[2])"
u_z = "2*x[0]*pow((1-x[0]),2)*x[1]*(1-x[1])*x[2]*x[2]*pow((1-x[2]),2) \
-2*x[0]*x[0]*(1-x[0])*x[1]*(1-x[1])*x[2]*x[2]*pow((1-x[2]),2)-2*x[0]* \
(1-x[0])*x[1]*pow((1-x[1]),2)*x[2]*x[2]*pow((1-x[2]),2)+2*x[0]*(1-x[0]) \
*x[1]*x[1]*(1-x[1])*x[2]*x[2]*pow((1-x[2]),2)"
u_t = Expression((u_x,u_y,u_z))

w1="((0.4e1-0.8e1*x[2])*pow(x[1],0.4e1)+(0.16e2*x[2]-0.8e1)*pow(x[1],0.3e1)+(0.4e1-0.8e1*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1]+0.4e1*pow(x[2],0.4e1)-0.8e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*pow(x[0],0.3e1)+((-0.4e1+0.12e2*x[2]*x[2])*pow(x[1],0.4e1)+(0.8e1-0.24e2*x[2]*x[2])*pow(x[1],0.3e1)+(-0.4e1+0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[1]*x[1]-0.4e1*pow(x[2],0.4e1)+0.8e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[0]*x[0]+((-0.12e2*x[2]*x[2]+0.8e1*x[2])*pow(x[1],0.4e1)+(0.24e2*x[2]*x[2]-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1])*x[0]"
w2="((0.4e1-0.8e1*x[2])*pow(x[1],0.3e1)+(-0.4e1+0.12e2*x[2]*x[2])*x[1]*x[1]+(-0.12e2*x[2]*x[2]+0.8e1*x[2])*x[1])*pow(x[0],0.4e1)+((0.16e2*x[2]-0.8e1)*pow(x[1],0.3e1)+(0.8e1-0.24e2*x[2]*x[2])*x[1]*x[1]+(0.24e2*x[2]*x[2]-0.16e2*x[2])*x[1])*pow(x[0],0.3e1)+((0.4e1-0.8e1*x[2])*pow(x[1],0.3e1)+(-0.4e1+0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[1]*x[1]+(-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[1])*x[0]*x[0]+((-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1])*x[0]+(0.4e1*pow(x[2],0.4e1)-0.8e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*pow(x[1],0.3e1)+(-0.4e1*pow(x[2],0.4e1)+0.8e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[1]*x[1]"
w3="((0.12e2*x[2]*x[2]-0.12e2*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]+0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*pow(x[0],0.4e1)+((-0.24e2*x[2]*x[2]+0.24e2*x[2])*x[1]*x[1]+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*x[1]-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*pow(x[0],0.3e1)+((0.12e2*x[2]*x[2]-0.12e2*x[2])*pow(x[1],0.4e1)+(-0.24e2*x[2]*x[2]+0.24e2*x[2])*pow(x[1],0.3e1)+(0.24e2*x[2]*x[2]-0.24e2*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]+0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[0]*x[0]+((-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*pow(x[1],0.4e1)+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]*x[1])*x[0]+(0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*pow(x[1],0.4e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[1]*x[1]"
w_t = Expression((w1,w2,w3))

p_t = Expression("-2.0*x[0]*x[1]*x[2]+x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[0]*x[1] \
+x[0]*x[2]+x[1]*x[2]-x[0]-x[1]-x[2]")

f1="(0.8e1*pow(x[1],0.3e1)-0.24e2*x[1]*x[1]*x[2]+(-0.8e1+0.24e2*x[2]*x[2])*x[1]-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*pow(x[0],0.4e1)+(-0.16e2*pow(x[1],0.3e1)+0.48e2*x[1]*x[1]*x[2]+(0.16e2-0.48e2*x[2]*x[2])*x[1]+0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*pow(x[0],0.3e1)+((0.8e1+0.48e2*x[2]*x[2]-0.48e2*x[2])*pow(x[1],0.3e1)+(-0.48e2*pow(x[2],0.3e1)+0.24e2*x[2])*x[1]*x[1]+(-0.8e1+0.48e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2])*x[1]-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[0]*x[0]+(0.2e1+(-0.48e2*x[2]*x[2]+0.48e2*x[2])*pow(x[1],0.3e1)+(0.48e2*pow(x[2],0.3e1)-0.48e2*x[2])*x[1]*x[1]+(-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2])*x[1])*x[0]+(0.8e1*x[2]*x[2]-0.8e1*x[2])*pow(x[1],0.3e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]*x[1]+(0.1e1+0.8e1*pow(x[2],0.3e1)-0.8e1*x[2]*x[2]-0.2e1*x[2])*x[1]-0.1e1+x[2]"
f2="(-0.8e1*pow(x[1],0.4e1)+0.16e2*pow(x[1],0.3e1)+(-0.8e1-0.48e2*x[2]*x[2]+0.48e2*x[2])*x[1]*x[1]+(0.48e2*x[2]*x[2]-0.48e2*x[2])*x[1]-0.8e1*x[2]*x[2]+0.8e1*x[2])*pow(x[0],0.3e1)+(0.24e2*pow(x[1],0.4e1)*x[2]-0.48e2*pow(x[1],0.3e1)*x[2]+(0.48e2*pow(x[2],0.3e1)-0.24e2*x[2])*x[1]*x[1]+(-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2])*x[1]+0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*x[0]*x[0]+((0.8e1-0.24e2*x[2]*x[2])*pow(x[1],0.4e1)+(-0.16e2+0.48e2*x[2]*x[2])*pow(x[1],0.3e1)+(0.8e1-0.48e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[1]*x[1]+(0.48e2*pow(x[2],0.3e1)-0.48e2*x[2]*x[2])*x[1]+0.1e1-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2]*x[2]-0.2e1*x[2])*x[0]+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*pow(x[1],0.4e1)+(-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2])*pow(x[1],0.3e1)+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*x[1]*x[1]+0.2e1*x[1]-0.1e1+x[2]"
f3="((0.8e1+0.48e2*x[2]*x[2]-0.48e2*x[2])*x[1]*x[1]+(-0.8e1-0.48e2*x[2]*x[2]+0.48e2*x[2])*x[1]+0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*pow(x[0],0.3e1)+((-0.8e1-0.48e2*x[2]*x[2]+0.48e2*x[2])*pow(x[1],0.3e1)+(0.8e1-0.24e2*pow(x[2],0.4e1)+0.48e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.48e2*x[2])*x[1])*x[0]*x[0]+((0.8e1+0.48e2*x[2]*x[2]-0.48e2*x[2])*pow(x[1],0.3e1)+(-0.8e1+0.24e2*pow(x[2],0.4e1)-0.48e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.48e2*x[2])*x[1]*x[1]-0.2e1*x[1]+0.1e1-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[0]+(-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.1e1+0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1]-0.1e1+0.2e1*x[2]"
f = Expression((f1,f2,f3))

Du00="((-0.16e2*x[2]*x[2]+0.16e2*x[2])*pow(x[1],0.3e1)+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*x[1]*x[1]+(-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2]*x[2])*x[1])*pow(x[0],0.3e1)+((0.24e2*x[2]*x[2]-0.24e2*x[2])*pow(x[1],0.3e1)+(-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2])*x[1]*x[1]+(0.24e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2])*x[1])*x[0]*x[0]+((-0.8e1*x[2]*x[2]+0.8e1*x[2])*pow(x[1],0.3e1)+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1])*x[0]"
Du01="((-0.12e2*x[2]*x[2]+0.12e2*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*x[1]-0.4e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*pow(x[0],0.4e1)+((0.24e2*x[2]*x[2]-0.24e2*x[2])*x[1]*x[1]+(-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2])*x[1]+0.8e1*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*pow(x[0],0.3e1)+((-0.12e2*x[2]*x[2]+0.12e2*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2])*x[1]-0.4e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*x[0]*x[0]"
Du02="((-0.8e1*x[2]+0.4e1)*pow(x[1],0.3e1)+(-0.4e1+0.12e2*x[2]*x[2])*x[1]*x[1]+(-0.12e2*x[2]*x[2]+0.8e1*x[2])*x[1])*pow(x[0],0.4e1)+((0.16e2*x[2]-0.8e1)*pow(x[1],0.3e1)+(0.8e1-0.24e2*x[2]*x[2])*x[1]*x[1]+(0.24e2*x[2]*x[2]-0.16e2*x[2])*x[1])*pow(x[0],0.3e1)+((-0.8e1*x[2]+0.4e1)*pow(x[1],0.3e1)+(-0.4e1+0.12e2*x[2]*x[2])*x[1]*x[1]+(-0.12e2*x[2]*x[2]+0.8e1*x[2])*x[1])*x[0]*x[0]"
Du10="((0.12e2*x[2]*x[2]-0.12e2*x[2])*pow(x[1],0.4e1)+(-0.24e2*x[2]*x[2]+0.24e2*x[2])*pow(x[1],0.3e1)+(0.12e2*x[2]*x[2]-0.12e2*x[2])*x[1]*x[1])*x[0]*x[0]+((-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*pow(x[1],0.4e1)+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]*x[1])*x[0]+(0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*pow(x[1],0.4e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.4e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[1]*x[1]"
Du11="((0.16e2*x[2]*x[2]-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.24e2*x[2]*x[2]+0.24e2*x[2])*x[1]*x[1]+(0.8e1*x[2]*x[2]-0.8e1*x[2])*x[1])*pow(x[0],0.3e1)+((-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2])*pow(x[1],0.3e1)+(0.24e2*pow(x[2],0.3e1)-0.24e2*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1])*x[0]*x[0]+((0.16e2*pow(x[2],0.3e1)-0.16e2*x[2]*x[2])*pow(x[1],0.3e1)+(-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1])*x[0]"
Du12="((0.8e1*x[2]-0.4e1)*pow(x[1],0.4e1)+(-0.16e2*x[2]+0.8e1)*pow(x[1],0.3e1)+(0.8e1*x[2]-0.4e1)*x[1]*x[1])*pow(x[0],0.3e1)+((0.4e1-0.12e2*x[2]*x[2])*pow(x[1],0.4e1)+(-0.8e1+0.24e2*x[2]*x[2])*pow(x[1],0.3e1)+(0.4e1-0.12e2*x[2]*x[2])*x[1]*x[1])*x[0]*x[0]+((0.12e2*x[2]*x[2]-0.8e1*x[2])*pow(x[1],0.4e1)+(-0.24e2*x[2]*x[2]+0.16e2*x[2])*pow(x[1],0.3e1)+(0.12e2*x[2]*x[2]-0.8e1*x[2])*x[1]*x[1])*x[0]"
Du20="((-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.12e2*x[2]*x[2])*x[1]*x[1]+(0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.12e2*x[2]*x[2])*x[1])*x[0]*x[0]+((0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1])*x[0]+(-0.4e1*pow(x[2],0.4e1)+0.8e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.4e1*pow(x[2],0.4e1)-0.8e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*x[1]*x[1]"
Du21="((-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1]+0.4e1*pow(x[2],0.4e1)-0.8e1*pow(x[2],0.3e1)+0.4e1*x[2]*x[2])*pow(x[0],0.3e1)+((0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.12e2*x[2]*x[2])*x[1]*x[1]-0.4e1*pow(x[2],0.4e1)+0.8e1*pow(x[2],0.3e1)-0.4e1*x[2]*x[2])*x[0]*x[0]+((-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.12e2*x[2]*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1])*x[0]"
Du22="((-0.16e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.8e1*x[2])*x[1]*x[1]+(0.16e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[1])*pow(x[0],0.3e1)+((0.16e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*pow(x[1],0.3e1)+(-0.16e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.8e1*x[2])*x[1])*x[0]*x[0]+((-0.16e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.8e1*x[2])*pow(x[1],0.3e1)+(0.16e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[1]*x[1])*x[0]"
Du=Expression(((Du00,Du01,Du02),(Du10,Du11,Du12),(Du20,Du21,Du22)))
Du1=Expression((Du00,Du01,Du02))
Du2=Expression((Du10,Du11,Du12))
Du3=Expression((Du20,Du21,Du22))

Dw00="((-0.24e2*x[2]+0.12e2)*pow(x[1],0.4e1)+(0.48e2*x[2]-0.24e2)*pow(x[1],0.3e1)+(-0.24e2*x[2]+0.12e2)*x[1]*x[1]+(-0.24e2*pow(x[2],0.4e1)+0.48e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2])*x[1]+0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.12e2*x[2]*x[2])*x[0]*x[0]+((-0.8e1+0.24e2*x[2]*x[2])*pow(x[1],0.4e1)+(0.16e2-0.48e2*x[2]*x[2])*pow(x[1],0.3e1)+(-0.8e1+0.24e2*pow(x[2],0.4e1)-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2])*x[1]*x[1]-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[0]+(-0.12e2*x[2]*x[2]+0.8e1*x[2])*pow(x[1],0.4e1)+(0.24e2*x[2]*x[2]-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1]"
Dw01="((-0.32e2*x[2]+0.16e2)*pow(x[1],0.3e1)+(0.48e2*x[2]-0.24e2)*x[1]*x[1]+(-0.16e2*x[2]+0.8e1)*x[1]-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*pow(x[0],0.3e1)+((0.48e2*x[2]*x[2]-0.16e2)*pow(x[1],0.3e1)+(0.24e2-0.72e2*x[2]*x[2])*x[1]*x[1]+(-0.8e1+0.24e2*pow(x[2],0.4e1)-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2])*x[1])*x[0]*x[0]+((-0.48e2*x[2]*x[2]+0.32e2*x[2])*pow(x[1],0.3e1)+(0.72e2*x[2]*x[2]-0.48e2*x[2])*x[1]*x[1]+(-0.24e2*pow(x[2],0.4e1)+0.48e2*pow(x[2],0.3e1)-0.48e2*x[2]*x[2]+0.16e2*x[2])*x[1]+0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[0]"
Dw02="(-0.8e1*pow(x[1],0.4e1)+0.16e2*pow(x[1],0.3e1)-0.8e1*x[1]*x[1]+(-0.32e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2]-0.16e2*x[2])*x[1]+0.16e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*pow(x[0],0.3e1)+(0.24e2*pow(x[1],0.4e1)*x[2]-0.48e2*pow(x[1],0.3e1)*x[2]+(0.48e2*pow(x[2],0.3e1)-0.72e2*x[2]*x[2]+0.48e2*x[2])*x[1]*x[1]-0.16e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.8e1*x[2])*x[0]*x[0]+((-0.24e2*x[2]+0.8e1)*pow(x[1],0.4e1)+(0.48e2*x[2]-0.16e2)*pow(x[1],0.3e1)+(-0.48e2*pow(x[2],0.3e1)+0.72e2*x[2]*x[2]-0.48e2*x[2]+0.8e1)*x[1]*x[1]+(0.32e2*pow(x[2],0.3e1)-0.48e2*x[2]*x[2]+0.16e2*x[2])*x[1])*x[0]"
Dw10="((-0.32e2*x[2]+0.16e2)*pow(x[1],0.3e1)+(0.48e2*x[2]*x[2]-0.16e2)*x[1]*x[1]+(-0.48e2*x[2]*x[2]+0.32e2*x[2])*x[1])*pow(x[0],0.3e1)+((0.48e2*x[2]-0.24e2)*pow(x[1],0.3e1)+(0.24e2-0.72e2*x[2]*x[2])*x[1]*x[1]+(0.72e2*x[2]*x[2]-0.48e2*x[2])*x[1])*x[0]*x[0]+((-0.16e2*x[2]+0.8e1)*pow(x[1],0.3e1)+(-0.8e1+0.24e2*pow(x[2],0.4e1)-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2])*x[1]*x[1]+(-0.24e2*pow(x[2],0.4e1)+0.48e2*pow(x[2],0.3e1)-0.48e2*x[2]*x[2]+0.16e2*x[2])*x[1])*x[0]+(-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*pow(x[1],0.3e1)+(0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[1]"
Dw11="((-0.24e2*x[2]+0.12e2)*x[1]*x[1]+(-0.8e1+0.24e2*x[2]*x[2])*x[1]-0.12e2*x[2]*x[2]+0.8e1*x[2])*pow(x[0],0.4e1)+((0.48e2*x[2]-0.24e2)*x[1]*x[1]+(0.16e2-0.48e2*x[2]*x[2])*x[1]+0.24e2*x[2]*x[2]-0.16e2*x[2])*pow(x[0],0.3e1)+((-0.24e2*x[2]+0.12e2)*x[1]*x[1]+(-0.8e1+0.24e2*pow(x[2],0.4e1)-0.48e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2])*x[1]-0.12e2*pow(x[2],0.4e1)+0.24e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*x[0]*x[0]+((-0.24e2*pow(x[2],0.4e1)+0.48e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2])*x[1]*x[1]+0.8e1*pow(x[2],0.4e1)-0.16e2*pow(x[2],0.3e1)+0.8e1*x[2]*x[2])*x[0]+(0.12e2*pow(x[2],0.4e1)-0.24e2*pow(x[2],0.3e1)+0.12e2*x[2]*x[2])*x[1]*x[1]+(-0.8e1*pow(x[2],0.4e1)+0.16e2*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1]"
Dw12="(-0.8e1*pow(x[1],0.3e1)+0.24e2*x[1]*x[1]*x[2]+(-0.24e2*x[2]+0.8e1)*x[1])*pow(x[0],0.4e1)+(0.16e2*pow(x[1],0.3e1)-0.48e2*x[1]*x[1]*x[2]+(0.48e2*x[2]-0.16e2)*x[1])*pow(x[0],0.3e1)+(-0.8e1*pow(x[1],0.3e1)+(0.48e2*pow(x[2],0.3e1)-0.72e2*x[2]*x[2]+0.48e2*x[2])*x[1]*x[1]+(-0.48e2*pow(x[2],0.3e1)+0.72e2*x[2]*x[2]-0.48e2*x[2]+0.8e1)*x[1])*x[0]*x[0]+((-0.32e2*pow(x[2],0.3e1)+0.48e2*x[2]*x[2]-0.16e2*x[2])*pow(x[1],0.3e1)+(0.32e2*pow(x[2],0.3e1)-0.48e2*x[2]*x[2]+0.16e2*x[2])*x[1])*x[0]+(0.16e2*pow(x[2],0.3e1)-0.24e2*x[2]*x[2]+0.8e1*x[2])*pow(x[1],0.3e1)+(-0.16e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2]-0.8e1*x[2])*x[1]*x[1]"
Dw20="((0.48e2*x[2]*x[2]-0.48e2*x[2])*x[1]*x[1]+(-0.32e2*pow(x[2],0.3e1)+0.32e2*x[2])*x[1]+0.16e2*pow(x[2],0.3e1)-0.16e2*x[2]*x[2])*pow(x[0],0.3e1)+((-0.72e2*x[2]*x[2]+0.72e2*x[2])*x[1]*x[1]+(0.48e2*pow(x[2],0.3e1)-0.48e2*x[2])*x[1]-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[0]*x[0]+((0.24e2*x[2]*x[2]-0.24e2*x[2])*pow(x[1],0.4e1)+(-0.48e2*x[2]*x[2]+0.48e2*x[2])*pow(x[1],0.3e1)+(0.48e2*x[2]*x[2]-0.48e2*x[2])*x[1]*x[1]+(-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2])*x[1]+0.8e1*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[0]+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*pow(x[1],0.4e1)+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*pow(x[1],0.3e1)+(-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[1]*x[1]"
Dw21="((0.24e2*x[2]*x[2]-0.24e2*x[2])*x[1]-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*pow(x[0],0.4e1)+((-0.48e2*x[2]*x[2]+0.48e2*x[2])*x[1]+0.16e2*pow(x[2],0.3e1)-0.16e2*x[2])*pow(x[0],0.3e1)+((0.48e2*x[2]*x[2]-0.48e2*x[2])*pow(x[1],0.3e1)+(-0.72e2*x[2]*x[2]+0.72e2*x[2])*x[1]*x[1]+(0.48e2*x[2]*x[2]-0.48e2*x[2])*x[1]-0.8e1*pow(x[2],0.3e1)+0.8e1*x[2])*x[0]*x[0]+((-0.32e2*pow(x[2],0.3e1)+0.32e2*x[2])*pow(x[1],0.3e1)+(0.48e2*pow(x[2],0.3e1)-0.48e2*x[2])*x[1]*x[1]+(-0.16e2*pow(x[2],0.3e1)+0.16e2*x[2])*x[1])*x[0]+(0.16e2*pow(x[2],0.3e1)-0.16e2*x[2]*x[2])*pow(x[1],0.3e1)+(-0.24e2*pow(x[2],0.3e1)+0.24e2*x[2]*x[2])*x[1]*x[1]+(0.8e1*pow(x[2],0.3e1)-0.8e1*x[2]*x[2])*x[1]"
Dw22="((0.24e2*x[2]-0.12e2)*x[1]*x[1]+(0.8e1-0.24e2*x[2]*x[2])*x[1]+0.12e2*x[2]*x[2]-0.8e1*x[2])*pow(x[0],0.4e1)+((-0.48e2*x[2]+0.24e2)*x[1]*x[1]+(0.48e2*x[2]*x[2]-0.16e2)*x[1]-0.24e2*x[2]*x[2]+0.16e2*x[2])*pow(x[0],0.3e1)+((0.24e2*x[2]-0.12e2)*pow(x[1],0.4e1)+(-0.48e2*x[2]+0.24e2)*pow(x[1],0.3e1)+(0.48e2*x[2]-0.24e2)*x[1]*x[1]+(0.8e1-0.24e2*x[2]*x[2])*x[1]+0.12e2*x[2]*x[2]-0.8e1*x[2])*x[0]*x[0]+((0.8e1-0.24e2*x[2]*x[2])*pow(x[1],0.4e1)+(0.48e2*x[2]*x[2]-0.16e2)*pow(x[1],0.3e1)+(0.8e1-0.24e2*x[2]*x[2])*x[1]*x[1])*x[0]+(0.12e2*x[2]*x[2]-0.8e1*x[2])*pow(x[1],0.4e1)+(-0.24e2*x[2]*x[2]+0.16e2*x[2])*pow(x[1],0.3e1)+(0.12e2*x[2]*x[2]-0.8e1*x[2])*x[1]*x[1]"

Dw_t = Expression(((Dw00,Dw01,Dw02),(Dw10,Dw11,Dw12),(Dw20,Dw21,Dw22)))
Dp_t = Expression(("-0.2e1*x[1]*x[2]+0.2e1*x[0]+x[1]+x[2]-0.1e1","-0.2e1*x[0]*x[2]+0.2e1*x[1]+x[0]+x[2]-0.1e1","-0.2e1*x[0]*x[1]+0.2e1*x[2]+x[0]+x[1]-0.1e1"))

# Boundaries
def boundary(x, on_boundary):
	return x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or x[2] < DOLFIN_EPS \
or x[0] > 1-DOLFIN_EPS or x[1] > 1-DOLFIN_EPS or x[2] > 1-DOLFIN_EPS
def corner(x, on_boundary):
    return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS and x[2] < DOLFIN_EPS

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0, 0.0))

def calc(method,dim):
	useKrylovSolver = 0
	SpecifyPrecond = 1
	print 'Solving for method: ',method
	print 'Mesh size: ',dim
	from time import time
	t1 = time()

	ERR = np.zeros(cols)

	if (method == 'MF' and dim > 20):
		return ERR # because of MPI out-of-memory error

	# create mesh
	mesh = UnitCube(dim, dim, dim)
	h=1.0/dim
	h2=h*h

	if method=='MF':
		V1 = VectorFunctionSpace(mesh, "CG", 2)
		Q = FunctionSpace(mesh, "CG", 1) 
		W = V1*Q
	elif (method=='CY1' or method=='QNE1' or method=='QNE1h'):
		V1 = VectorFunctionSpace(mesh, "CG", 1)
		V2 = VectorFunctionSpace(mesh, "CG", 1)
		Q = FunctionSpace(mesh, "CG", 1)
		W = MixedFunctionSpace([V1,V2,Q])
	elif (method=='CY2' or method=='QNE2' or method=='QNE2h'):
		V1 = VectorFunctionSpace(mesh, "CG", 2)
		V2 = VectorFunctionSpace(mesh, "CG", 1)
		Q = FunctionSpace(mesh, "CG", 1)
		W = MixedFunctionSpace([V1,V2,Q])
	
	bc0 = DirichletBC(W.sub(0), noslip, boundary)
	if method=='MF':
		bc1 = DirichletBC(W.sub(1), Constant(0.0), corner, "pointwise")
	else:
		bc1 = DirichletBC(W.sub(2), Constant(0.0), corner, "pointwise")
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

	ED = {} # dictionary
	ED['unknowns']=len(U.vector().array())

	print 'Number of unknowns: ',ED['unknowns']
	
	if dim >= 12:
		useKrylovSolver = 1

	t2=time()
	if useKrylovSolver:
		A, bb = assemble_system(a, L, bcs)
		if method<>'MF':
			solver = KrylovSolver("cg","jacobi")
			solver.set_operator(A)
		else:
			solver = KrylovSolver("tfqmr","amg")
			if SpecifyPrecond:
				b = inner(grad(u),grad(v))*dx + p*q*dx
				P, btmp = assemble_system(b, L, bcs)
				solver.set_operators(A,P)
			else:
				solver.set_operator(A)
				
		solver.parameters["absolute_tolerance"] = 1.0E-10
		solver.parameters["relative_tolerance"] = 1.0E-06
		solver.parameters["maximum_iterations"] = 10000
		solver.parameters["report"] = True
		solver.parameters["monitor_convergence"] = True
		ED['krylov_ite'] = solver.solve(U.vector(), bb)
	else:
		solve(a==L,U,bcs)
		ED['krylov_ite'] = 0
	t3=time()
	print 'Solver took %0.3f seconds.' % (t3-t2)

	if method=='MF':
		u,p = U.split()
	else:
		u,w,p = U.split()

	# Calculate errors
	Ve = VectorFunctionSpace(mesh, "CG", 2)
	Qe = FunctionSpace(mesh, "CG", 2)
	Te = TensorFunctionSpace(mesh, "CG", 2)

	u1=Expression(u_x)
	u2=Expression(u_y)
	u3=Expression(u_z)
	if dim<=12:
		Du_tp = project(Du,Te) # very memory intensive
	u_tp = project(u_t,Ve)
	u1_tp = project(u1,Qe)
	u2_tp = project(u2,Qe)
	u3_tp = project(u3,Qe)
	w_tp = project(w_t,Ve)
	p_tp = project(p_t,Qe)
	f_p = project(f,Ve)
	Du1_tp = project(Du1,Ve)
	Du2_tp = project(Du2,Ve)
	Du3_tp = project(Du3,Ve)
	if dim<=12:
		Dw_tp = project(Dw_t,Te)
	Dp_tp = project(Dp_t,Ve)
	
	# L2 errors
	ERR[0]=ED['u1err']= sqrt(assemble(inner(u[0]-u1,u[0]-u1)*dx))
	ERR[1]=ED['u2err']= sqrt(assemble(inner(u[1]-u2,u[1]-u2)*dx))
	ERR[2]=ED['u3err']= sqrt(assemble(inner(u[2]-u3,u[2]-u3)*dx))
	ERR[3]=ED['u_err']= sqrt(assemble(inner(u-u_t,u-u_t)*dx))
	print 'Vel. L2 err: ',ED['u_err']
	if method<>'MF':
		ERR[4]= ED['w_err'] = sqrt(assemble(inner(w-w_t,w-w_t)*dx))
		print 'Vor. L2 err: ',ED['w_err']
	ERR[5]=ED['p_err']= sqrt(assemble((p-p_t)**2*dx))
	
	# H1 (semi-norm) errors
	ERR[6]=ED['u1errH1']= sqrt(assemble(inner(grad(u[0])-Du1_tp,grad(u[0])-Du1_tp)*dx))
	ERR[7]=ED['u2errH1']= sqrt(assemble(inner(grad(u[1])-Du2_tp,grad(u[1])-Du2_tp)*dx))
	ERR[8]=ED['u3errH1']= sqrt(assemble(inner(grad(u[2])-Du3_tp,grad(u[2])-Du3_tp)*dx))

	if dim <= 12:	
		ERR[9]=ED['u_errH1']= sqrt(assemble(inner(grad(u)-Du_tp,grad(u)-Du_tp)*dx))
	else:
		ERR[9]=ED['u_errH1']= sqrt(assemble(inner(grad(u-u_tp),grad(u-u_tp))*dx))
	if method<>'MF':
		if dim <= 12:
			ERR[10]=ED['w_errH1']=sqrt(assemble(inner(grad(w)-Dw_tp,grad(w)-Dw_tp)*dx))
		else:		
			ERR[10]=ED['w_errH1']=sqrt(assemble(inner(grad(w-w_tp),grad(w-w_tp))*dx))
	ERR[11]=ED['p_errH1']= sqrt(assemble(inner(grad(p)-Dp_tp,grad(p)-Dp_tp)*dx))

	# Infinity norm errors
	#print dir(u)
	#print u[0],dir(u[0])
	#print u.vector()
	#print u.vector().array()
	#print u.vector().array()[0]
	#print w.vector()
	#ED['u1errInf'] = abs(u.vector().array()[0]-project(u1,Q).vector().array()).max()
	#ED['u2errInf'] = abs(u.vector().array()[1]-project(u2,Q).vector().array()).max()
	#ED['w_errInf'] = abs(w.vector().array()-project(w_t,Q).vector().array()).max()
	#ED['p_errInf'] = abs(p.vector().array()-project(p_t,Q).vector().array()).max()
	#raw_input("Press Enter to continue...")
	#quit()

	# Energy
	ERR[12]=ED['E'] = assemble(0.5*inner(grad(u),grad(u))*dx-inner(f,u)*dx)

	if dim <= 12:
		ERR[13]=ED['E_tr'] = assemble(0.5*inner(Du_tp,Du_tp)*dx-inner(f,u_t)*dx)
	else:
		ERR[13]=ED['E_tr'] = assemble(0.5*inner(grad(u_tp),grad(u_tp))*dx-inner(f,u_t)*dx)	
	ERR[14]=ED['E_err'] = ED['E']-ED['E_tr']
	ERR[15]=ED['E_re'] = 100*ED['E_err']/ED['E_tr']

	ERR[33]=ED['E2'] = ED['E'] + assemble(0.5*inner(u,u)*dx)
	ERR[34]=ED['E_tr2'] = ED['E_tr'] + assemble(0.5*inner(u_tp,u_tp)*dx)
	ERR[35]=ED['E_err2'] = ED['E2']-ED['E_tr2']
	ERR[36]=ED['E_re2'] = 100*ED['E_err2']/ED['E_tr2']

	ERR[16]=ED['F'] = assemble(inner(f_p,f_p)*dx) # total force
	if dim <= 12:
		ERR[17]=ED['gradU_tot'] = sqrt(assemble(inner(Du_tp,Du_tp)*dx))
		ERR[18]=ED['gradU_err'] = ED['u_errH1']
	else:
		ERR[17]=ED['gradU_tot'] = sqrt(assemble(inner(grad(u_tp),grad(u_tp))*dx))
		ERR[18]=ED['gradU_err'] = sqrt(assemble(inner(grad(u)-grad(u_tp),grad(u)-grad(u_tp))*dx))
	ERR[19]=ED['gradU_re'] = 100*ED['gradU_err']/ED['gradU_tot']

	# Residuals
	if method<>'MF':
		ERR[20]=ED['r1']=assemble(inner(curl(w)+grad(p)-f,curl(w)+grad(p)-f)*dx)
		ERR[21]=ED['r2']=assemble(inner(curl(u)-w,curl(u)-w)*dx)
		ERR[22]=ED['r3']=assemble(inner(div(u),div(u))*dx)

	t4 = time()
	print 'Error calculation took %0.3f seconds.' % (t4-t3)
	
	if 0:'''
	# count non-zero entries in A
	if not useKrylovSolver:
		A, bb = assemble_system(a, L, bcs)
	cnt = 0
	t5 = time()
	N = A.size(0)
	
	for i in range(N):
		row = A.getrow(i)
		row = row[1] # get second array, the one with actual values
		cnt += sum(abs(row)>DOLFIN_EPS)
	'''
	t6 = time()
	
	#print 'Counting non-zero elements took %0.3f seconds.' % (t6-t5)

	ERR[23]=ED['tSolve']= t3-t2 # time to solve
	ERR[24]=ED['tSetupSolve']= t3-t1 # time to setup and solve
	ERR[25]=ED['tErrors']= t4-t3 # time to calc errors
	#ERR[26]=ED['tSparsity']= t6-t5 # time to calc sparsity
	ERR[27]=ED['tTotal']= t6-t1 # total time elapsed
	
	ERR[28]=ED['dim']=dim
	#ERR[29]=ED['unknowns']=len(U.vector().array())
	#ERR[30]=ED['nonzero']=cnt
	#ERR[31]=ED['sparsity']=float(cnt)/float(N*N)
	ERR[32]=ED['krylov_ite']
		
	return ERR













