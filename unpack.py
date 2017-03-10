#!/usr/bin/python
# Unpacks the txt files for x and y data

from pylab import * # For plotting and math
import os # For input/output

# File names
fnamexc = 'converged_grid_x.txt'
fnameyc = 'converged_grid_y.txt'

# Load them
dataxc = loadtxt(fnamexc, skiprows=0)
datayc = loadtxt(fnameyc, skiprows=0)

# Restack them
dataxc = hstack((dataxc, dataxc[:, 0].reshape(-1, 1)))
datayc = hstack((datayc, datayc[:, 0].reshape(-1, 1)))

# Plot like this (if desired)
# figure()
# for i in range(mj+1):
#     plot(dataxc[i, :], datayc[i, :], color='k')

# for i in range(mi):
#     plot(dataxc[:, i], datayc[:, i], color='k')
        
# xlabel('Position x')
# ylabel('Position y')
# axis('equal')
