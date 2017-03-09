""" 
AA543 - CFD, HW#3, Winter 2017

Plots grids
"""

from matplotlib import pyplot as plt
import numpy as np


# Import airfoil profile
I_max  = 128
J_max = 64
num = []
x = np.zeros((I_max, J_max)) 
y = np.zeros((I_max, J_max)) 

i = 0
j=0
for line in open('cell_x.dat', 'r'):
    if j < 63 :
        x[i,j] = float(line)
        j = j +1
    else:
        x[i,j] = float(line)
        #~ print(x[0,j], ", ", x[-1,j]
        j = 0
        i = i + 1
        
i = 0
j = 0
for line in open('cell_y.dat', 'r'):
    if j < 63 :
        y[i,j] = float(line)
        j = j +1
    else:
        y[i,j] = float(line)
        #~ print(y[0,j], ", ", y[-1,j])
        j = 0
        i = i + 1


# Plot explicit vs analytic for time step specified
for r in range(len(x)):
    plt.plot(x[r,:], y[r,:], 'g')
for c in range(len(x[0])):
    plt.plot(x[:,c], y[:,c], 'r')
plt.xlabel('x')
plt.ylabel('y')
#~ plt.xlim([-1,2])
#~ plt.ylim([-1.5,1.5])
plt.show()

# Plot all three against each other
