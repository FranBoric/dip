import os
import sys

import dolfin
from fenics import *
import numpy as np
import matplotlib.pyplot as plt

os.system("rm ./2D/*")
#input("Press enter to continue..")


## Početne konstante, vrijeme, korak vremenski, parametri modela

T = 10.0
num_steps = 500
dt = T / num_steps

## Predefinirani parametri bistabilnog slučaja

# alpha_1  = 2.0
# beta_1 	 = 8.0
# tau_1 	 = 1.0
# lambda_1 = 1.0
# d_1 	 = 1.0
# alpha_2  = 7.0
# tau_2 	 = 13.0/2.0
# lambda_2 = 26.0
# d_2 	 = 1.0

## Predefinirani parametri monostabilnog slučaja

alpha_1  = 2.0
beta_1 	 = 8.0
tau_1 	 = 1.0
lambda_1 = 1.0
d_1 	 = 1.0
alpha_2  = 1.0
tau_2 	 = 42.0/43.0
lambda_2 = 1.0
d_2 	 = 1.0



## Mesh i funkcijski prostor

## koliko puta koliko elemenata po bridu
nx = 100
ny = 10

dolje_lijevo_x = 0.0
dolje_lijevo_y = 0.0
gore_desno_x = 1.0
gore_desno_y = 0.001

## Jedinični kvadrat
#mesh = UnitSquareMesh(nx,ny)
mesh = RectangleMesh(Point(dolje_lijevo_x , dolje_lijevo_y), Point(gore_desno_x, gore_desno_y), nx, ny, "right")

## Elementi
P1 = FiniteElement("P", triangle, 1)

## Prostor, Mixed u ovom slučaju označava činjenicu da je sustav više nepoznatih funkcija
V = FunctionSpace(mesh, MixedElement([P1, P1]))


## Definiranje ruba, default je homogeni Neumann
## Ovaj poziv će nam "dati rub"
bdr = MeshFunction('size_t', mesh, mesh.topology().dim()-1)

## slijedećih nekoliko linija u suštini kaže gdje je boundary (oko y = 1, tolerancija 1e-8)
## i nakon toga ga označi kao rub 2
class BoundaryY1(SubDomain):
	tol = 1E-8
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1], gore_desno_y, 1E-8)
bx3 = BoundaryY1()
bx3.mark(bdr, 2)

## Redefiniranje integracije po rubu ds kako bi mogli koristiti ds(2) kao integriranje po rubu označenom s 2
ds = Measure('ds', domain=mesh, subdomain_data=bdr)






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

#M_0 = project(Expression("2.0+7.0*(0.4<=x[0])*(x[0]<=0.6)", degree=1), V.sub(0).collapse())
#A_0 = project(Constant(0), V.sub(1).collapse())

M_0 = project(Expression("2.0+7.0*(0.4<=x[0])*(x[0]<=0.6)", degree=1), V.sub(0).collapse())
A_0 = project(Expression("0.0+0.1*(0.4<=x[0])*(x[0]<=0.6)", degree=1), V.sub(1).collapse())

assign(sol_n, [M_0, A_0])
M_n, A_n = split(sol_n)


## Forme

a = ((1/dt)*(M - M_n)*M_t + lambda_1*M*M_t + d_1*dot(grad(M),grad(M_t)) \
	+ (1/dt)*(A - A_n)*A_t +  lambda_2*A*A_t + d_2*dot(grad(A),grad(A_t)))*dx 


f =  (alpha_2*A/(1.0+A/tau_2)*M)*A_t*dx


## Rubni uvjet, u ovom slučaju jedan jedini Neumannov na granici 2

g = ((gore_desno_y*(alpha_1/(1.0+A/tau_1)+beta_1*A/(1.0+A/tau_1))*M_t)/d_1)*ds(2)

## Početak evolucije
t = 0

vtkfile = File('./2D/solution.pvd')

M_n, A_n = sol_n.split()

np.set_printoptions(threshold = np.inf)

# input("Press enter to continue..")

plt.figure()

#plot(M_n)
#plt.show()

vtkfile << (M_n, t)

for n in range(num_steps):

	t+= dt

	solve(a - f - g== 0, sol, solver_parameters={"newton_solver":{"relative_tolerance":1E-8},"newton_solver":{"maximum_iterations":50}})



	sol_n.assign(sol)
	M_n, A_n = sol_n.split()


	vtkfile << (M_n, t)

	

	# Greška
	# u_e = interpolate(u_D, V)
	# err = np.abs(u_e.vector().array() - u.vector().array()).max()
	# print('t = %.2f: error = %.3g' % (t, error))

# 	#input("Press enter to continue...")
	
