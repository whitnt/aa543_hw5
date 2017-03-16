""" 
AA543 - CFD, HW#5, Winter 2017

2D Contour plots of simulation variables

!!!! NOTE: will throw error if no variation in field (nothing to controur)!!!! 
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

# Read grid point coordinate files
x = np.loadtxt('output/cell_x.dat')
y = np.loadtxt('output/cell_y.dat')

# Read variable files
rho     = np.loadtxt('output/rho.dat')
rho_u   = np.loadtxt('output/rho_u.dat')
rho_v   = np.loadtxt('output/rho_v.dat')
rho_E   = np.loadtxt('output/rho_E.dat')

print(len(rho), len(x), len(y))

# Plot on single figure
f, ((sub1, sub2),(sub3, sub4)) = plt.subplots(2,2, sharex=True, sharey=True, num=2)
plot1 = sub1.tricontourf(x, y, rho, levels=np.linspace(rho.min(),rho.max(),20))
#~ sub1.title("rho")
plt.colorbar(plot1)
plot2 = sub2.tricontourf(x, y, rho_u, levels=np.linspace(rho_u.min(),rho_u.max(),20))
plt.colorbar(plot2)
plot3 = sub3.tricontourf(x, y, rho_v, levels=np.linspace(rho_v.min(),rho_v.max(),20))
plt.colorbar(plot3)
plot4 = sub4.tricontourf(x, y, rho_E, levels=np.linspace(rho_E.min(),rho_E.max(),20))
plt.colorbar(plot4)
plt.show()
