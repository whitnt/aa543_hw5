#!/usr/bin/python
# Unpacks the txt files for x and y data

import numpy as np
import os # For input/output

# File names
fnamexc = 'converged_grid_x.txt'
fnameyc = 'converged_grid_y.txt'

# Load them
dataxc = np.loadtxt(fnamexc, skiprows=0)
datayc = np.loadtxt(fnameyc, skiprows=0)

# Restack them
dataxc = hstack(dataxc)
datayc = hstack(datayc)

# Write to files
np.savetxt('grid_x.dat',dataxc)
np.savetxt('grid_y.dat',datayc)

    
# Plot like this (if desired)
# figure()
# for i in range(mj+1):
#     plot(dataxc[i, :], datayc[i, :], color='k')

# for i in range(mi):
#     plot(dataxc[:, i], datayc[:, i], color='k')
        
# xlabel('Position x')
# ylabel('Position y')
# axis('equal')
