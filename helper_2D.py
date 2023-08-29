import os
import sys

import dolfin
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import *
from matplotlib import cm


## koliko puta koliko elemenata po bridu
nx = 100
ny = 10

dolje_lijevo_x = 0.0
dolje_lijevo_y = 0.0
gore_desno_x = 1.0
gore_desno_y = 0.0001

T = 10.0
num_steps = 500
dt = T / num_steps


np.set_printoptions(threshold = np.inf)

for i in range(0, num_steps+int(num_steps/10), int(num_steps/10)):

	## Konstrukcija mesha za plotanje
	X = (np.linspace(0, gore_desno_x, (nx + 1)))
	Y = (np.linspace(0, gore_desno_y, (ny + 1)))
	Xx, Yy = np.meshgrid(X,Y)


	data = read_csv("./2D/test_" + str(i) + ".csv")

	dejta = data["f_10-0"].tolist()

	M_ARRY = np.array(dejta)
	M_ARRY = M_ARRY.reshape((ny + 1),(nx + 1))

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
	hm.set_title("Vrijeme " + str(i*dt))

	plt.show()