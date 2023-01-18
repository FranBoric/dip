import os
import sys

import dolfin
from fenics import *
import numpy as np
import matplotlib.pyplot as plt


## Početne konstante, vrijeme, korak vremenski, parametri modela

T = 5.0
num_steps = 100
dt = T / num_steps

alpha_1  = 2.0
beta_1 	 = 8.0
tau_1 	 = 1.0
lambda_1 = 1.0
d_1 	 = 1.0
alpha_2  = 7.0
tau_2 	 = 13.0/2.0
lambda_2 = 26.0
d_2 	 = 1.0




## Mesh i funkcijski prostor

## koliko puta koliko elemenata po bridu
nx = 50

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

M_0 = project(Expression("2.0+7.0*(0.45<=x[0])*(x[0]<=0.55)", degree=1), V.sub(0).collapse())
A_0 = project(Constant(0), V.sub(1).collapse())
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
vtkfile << (M_n, t)

for n in range(num_steps):

	t+= dt

	solve(a - f== 0, sol)


	sol_n.assign(sol)
	M_n, A_n = sol_n.split()
	vtkfile << (M_n, t)

	

	# Greška
	# u_e = interpolate(u_D, V)
	# err = np.abs(u_e.vector().array() - u.vector().array()).max()
	# print('t = %.2f: error = %.3g' % (t, error))

	#input("Press enter to continue...")
	
