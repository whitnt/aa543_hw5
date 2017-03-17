#!/usr/bin/python
# Script to plot the c++ output data
# Daniel Crews 3/16/2017

from pylab import * # For plotting and math
import os # For input/output

mi = 129
mj = 64

##### Import Data
## File names
xname = 'converged_grid_x.txt'
yname = 'converged_grid_y.txt'
xcenter_name = 'cell_x.dat'
ycenter_name = 'cell_y.dat'

#xnormname = 'dsj_x.dat'
#ynormname = 'dsj_y.dat'
xnormname = 'dsi_x.dat'
ynormname = 'dsi_y.dat'

##### Import data
## GRID
gridx = loadtxt(xname, skiprows=0)
gridy = loadtxt(yname, skiprows=0)
dataxc = hstack((gridx, gridx[:, 0].reshape(-1, 1)))
datayc = hstack((gridy, gridy[:, 0].reshape(-1, 1)))

## NORMALS
xnorm = loadtxt(xnormname, skiprows=0)
ynorm = loadtxt(ynormname, skiprows=0)

## CENTERS
cellx = loadtxt(xcenter_name, skiprows=0)
celly = loadtxt(ycenter_name, skiprows=0)

#print(xnorm.shape)
#print(dataxc.shape)

## Reshape
xnorm = reshape(xnorm, (65, 128), order='F')
ynorm = reshape(ynorm, (65, 128), order='F')

cellx = reshape(cellx, (64, 128), order='F')
celly = reshape(celly, (64, 128), order='F')

#xnorm = reshape(xnorm, (64, 129), order='F')
#ynorm = reshape(ynorm, (64, 129), order='F')

halfgridx = zeros(dataxc.shape)
halfgridy = zeros(dataxc.shape)

#print(halfgridx.shape)
#print(xnorm.shape)

##### Calculate x and y positions of normal
for i in range(mj):
    halfgridx[i, :] = dataxc[i, :] + (dataxc[i+1, :] - dataxc[i, :])/2.
    halfgridy[i, :] = datayc[i, :] + (datayc[i+1, :] - datayc[i, :])/2.

for j in range(mi):
    halfgridx[:, j] = dataxc[:, j] + (dataxc[:, j+1] - dataxc[:, j])/2.
    halfgridy[:, j] = datayc[:, j] + (datayc[:, j+1] - datayc[:, j])/2.

#print(halfgridx[0:65][0])
#print(dataxc[0:65][0])

## Plot the grid
figure()
for i in range(mj+1):
    plot(dataxc[i, :], datayc[i, :], color='k')

for i in range(mi):
    plot(dataxc[:, i], datayc[:, i], color='k')

for i in range(127):
    scatter(cellx[:, i], celly[:, i], color='g')

for i in range(65):
    for j in range(128):
        arrow(halfgridx[i, j], halfgridy[i, j],
              xnorm[i, j], ynorm[i, j], fc="k", ec="k",
              head_width=0.0005, head_length=0.0001)

xlabel('Position x')
ylabel('Position y')
axis('equal')

#figure()
#for i in range(mj+1):
#    plot(halfgridx[i, :], halfgridy[i, :], color='k')

show()
