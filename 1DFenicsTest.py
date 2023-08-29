import os
import sys

import dolfin
from fenics import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from matplotlib.ticker import FixedLocator, FixedFormatter
from mpl_toolkits import mplot3d


## Početne konstante, vrijeme, korak vremenski, parametri modela

T = 30.0
num_steps = 1000
dt = T / num_steps

## Predefinirani parametri bistabilnog slučaja

# alpha_1  = 2.0
# beta_1 	 = 8.0
# tau_1 	 = 1.0
# lambda_1 = 1.0
# d_1 	 = 0.001
# alpha_2  = 7.0
# tau_2 	 = 13.0/2.0
# lambda_2 = 26.0
# d_2 	 = 0.001

## Predefinirani parametri monostabilnog slučaja

alpha_1  = 2.0
beta_1 	 = 8.0
tau_1 	 = 1.0
lambda_1 = 1.0
d_1 	 = 0.001
alpha_2  = 1.0
tau_2 	 = 42.0/43.0
lambda_2 = 1.0
d_2 	 = 0.001

print("alpha_1/lambda_1")
print(alpha_1/lambda_1)
print("")

print("lambda_2/alpha_2")
print(lambda_2/alpha_2)
print("")

print("beta_1*tau_1/lambda_1")
print(beta_1*tau_1/lambda_1)
print("")

print("lambda_2/(alpha_2*tau_2)")
print(lambda_2/(alpha_2*tau_2))
print("")

sqrt1 = sqrt(beta_1/lambda_1 - alpha_1/(tau_1*lambda_2))
sqrt2 = sqrt(lambda_2/(alpha_2*tau_1) - alpha_1/(tau_1*lambda_1))
print("(sqrt1 - sqrt2)**2")
print((sqrt1 - sqrt2)**2)
print("")


input("Press to continue...")

## Mesh i funkcijski prostor

## koliko puta koliko elemenata po bridu
nx = 1000

## Jedinični kvadrat
mesh = UnitIntervalMesh(nx)

## Elementi
P1 = FiniteElement("P", mesh.ufl_cell(), 1)

## Prostor, Mixed u ovom slučaju označava činjenicu da je sustav više nepoznatih funkcija
V = FunctionSpace(mesh, MixedElement([P1, P1]))

## Nepoznanice
## sol   	- rješenje
## sol_n 	- "prethodno rješenje", bitno za evoluciju, u vremenu t=0 je početni uvjet
## M, A     - Pojedinačne funkcije koje želimo dobiti
## M_t, A_t - Njihovi "test funkcija" ekvivalenti
sol = Function(V)
sol_n = Function(V)
M, A = split(sol)
M_t, A_t = TestFunctions(V)

## Početni uvjet, prvo ih definiramo, dodijelimo ih "prethodnom rješenju" i nakon toga 
## 	zapravo uzmemo nazad varijable koje koristimo u formulaciji

M_0 = project(Expression("2.0+3.0*(0.45<=x[0])*(x[0]<=0.55)", degree=1), V.sub(0).collapse())
A_0 = project(Expression("0.0+0.5*(0.45<=x[0])*(x[0]<=0.55)", degree=1), V.sub(1).collapse())
#A_0 = project(Expression("0.0+20.50*(0.35<=x[0])*(x[0]<=0.65)", degree=1), V.sub(1).collapse())

assign(sol_n, [M_0, A_0])
M_n, A_n = split(sol_n)


## Forme

a = ((1/dt)*(M - M_n)*M_t + lambda_1*M*M_t + d_1*dot(grad(M),grad(M_t)) \
	+ (1/dt)*(A - A_n)*A_t +  lambda_2*A*A_t + d_2*dot(grad(A),grad(A_t)))*dx 


f =  ((alpha_1/(1.0+A/tau_1)+beta_1*A/(1.0+A/tau_1))*M_t + (alpha_2*A/(1.0+A/tau_2)*M)*A_t)*dx



## Početak evolucije
t = 0

vtkfile = File('./1D/solution.pvd')

M_n, A_n = sol_n.split()

np.set_printoptions(threshold = np.inf)


M_ARRY = []
M_ARRY.append(np.array(M_n.vector())[::2])
A_ARRY = []
A_ARRY.append(np.array(A_n.vector())[1::2])


#plot(M_n)
#plt.show()

vtkfile << (M_n, t)

for n in range(num_steps):

	t+= dt

	solve(a - f== 0, sol)

	sol_n.assign(sol)


	M_n, A_n = sol_n.split()
	#if(t <= 4*dt):
		#plot(M_n)
		#plt.show()


	M_ARRY.append(np.array(M_n.vector())[::2])
	A_ARRY.append(np.array(A_n.vector())[1::2])

	vtkfile << (M_n, t)

	

	# # Greška
	# u_e = interpolate(u_D, V)
	# err = np.abs(u_e.vector().array() - u.vector().array()).max()
	# print('t = %.2f: error = %.3g' % (t, error))

	# input("Press enter to continue...")
	

NNN = (np.array(M_n.vector()))[::2]




#print(np.shape(M_ARRY))
#print(M_ARRY)

## Konstrukcija mesha za plotanje
X = (np.linspace(0, 1, (nx + 1)))
Y = (np.linspace(0, T, (num_steps + 1)))
Xx, Yy = np.meshgrid(X,Y)


## Sve moguće za plot M_n rješenja kroz vrijeme

M_ARRY = np.array(M_ARRY)

cmap = cm.get_cmap('viridis')


hfm = plt.figure()
hm = hfm.add_subplot(111, projection='3d')

vminM = np.min(M_ARRY)
vmaxM = np.max(M_ARRY)

surface_M = hm.plot_surface(Xx, Yy, M_ARRY, cmap=cmap, vmin=vminM, vmax=vmaxM)


cbar = hfm.colorbar(surface_M, ax=hm, ticks=np.linspace(vminM, vmaxM, 10))
cbar.set_label("M(x,t)")

hm.set_xlabel("x")
hm.set_ylabel("vrijeme")
hm.set_zlabel("M(x,t)")
hm.set_title("1D - M(x,t) - Monostabilan.")

plt.show()

## Sve moguće za plot A_n rješenja kroz vrijeme, potpuno isti postupak kao i za M_n

cmap = cm.get_cmap('viridis')

A_ARRY = np.array(A_ARRY)

hfa = plt.figure()
ha = hfa.add_subplot(111, projection='3d')

vminA = np.min(A_ARRY)
vmaxA = np.max(A_ARRY)

surface_A = ha.plot_surface(Xx, Yy, A_ARRY, cmap=cmap, vmin=vminA, vmax=vmaxA)

cbar = hfa.colorbar(surface_A, ax=ha, ticks=np.linspace(vminA, vmaxA, 10))
cbar.set_label("A(x,t)")

ha.set_xlabel("x")
ha.set_ylabel("vrijeme")
ha.set_zlabel("A(x,t)")
ha.set_title("1D - A(x,t) - Monostabilan.")

plt.show()


#print(np.array(sol_n.vector()))