#! /usr/bin/env python

''' @brief Generates surface points and normals for a capsule-shaped immersed boundary.
Only works when there's up-down symmetry.
'''

from math import sqrt, pi
import sys
import numpy as np
from matplotlib import pyplot as plt

# the four vertices
x1 = 6.0; y1 = -0.25
x2 = 6.0; y2 = 0.25
x3 = 4.0; y3 = 1.0
x4 = 4.0; y4 = -1.0

# number of points in the 4 sides
np1 = 100; np2 = 150; np3 = 50; np4 = 150
ncirc1 = 50; ncirc2 = 400

# output file
ofile = "ib-capsule.dat"

# centres of end circles
m2 = (y2-y3)/(x2-x3); m4 = (y1-y4)/(x1-x4)

c1x = m2*m4/(m4-m2)* (y2- y1 + x2/m2 - x1/m4)
c1y = (y1+y2)/2.0
c2x = m2*m4/(m4-m2)* (y3-y4 + x3/m2 - x4/m4)
c2y = (y3 + y4)/2.0

# radii of circles at the ends
rad1 = sqrt((c1x-x2)**2 + (c1y-y2)**2)
rad2 = sqrt((c2x-x3)**2 + (c2y-y3)**2)

## points and normals of quad

t1 = np.linspace(0,1.0,np1, endpoint = False, dtype = np.float64)
t2 = np.linspace(0,1.0,np2, endpoint = False, dtype = np.float64)
t3 = np.linspace(0,1.0,np3, endpoint = False, dtype = np.float64)
t4 = np.linspace(0,1.0,np4, endpoint = False, dtype = np.float64)

pq1 = np.zeros((np1,2), dtype = np.float64)
pq2 = np.zeros((np2,2), dtype = np.float64)
pq3 = np.zeros((np3,2), dtype = np.float64)
pq4 = np.zeros((np4,2), dtype = np.float64)
nq1 = np.zeros((np1,2), dtype = np.float64)
nq2 = np.zeros((np2,2), dtype = np.float64)
nq3 = np.zeros((np3,2), dtype = np.float64)
nq4 = np.zeros((np4,2), dtype = np.float64)

pq1[:,0] = x1 + (x2-x1)*t1; pq1[:,1] = y1 + (y2-y1)*t1
pq2[:,0] = x2 + (x3-x2)*t2; pq2[:,1] = y2 + (y3-y2)*t2
pq3[:,0] = x3 + (x4-x3)*t3; pq3[:,1] = y3 + (y4-y3)*t3
pq4[:,0] = x4 + (x1-x4)*t4; pq4[:,1] = y4 + (y1-y4)*t4

mag2 = sqrt((y3-y2)**2 + (x3-x2)**2)
mag4 = sqrt((y1-y4)**2 + (x1-x4)**2)

nq1[:,0] = 1.0; nq1[:,1] = 0.0
nq2[:,0] = (y3-y2)/mag2; nq2[:,1] = -(x3-x2)/mag2
nq3[:,0] = -1.0; nq3[:,1] = 0.0
nq4[:,0] = (y1-y4)/mag4; nq4[:,1] = -(x1-x4)/mag4

nq1[0,0] = y1-y4; nq1[0,1] = -(x1-x4)
nq3[0,0] = y3-y2; nq3[0,1] = -(x3-x2)

## points and normals of circles

tc1 = np.linspace(0,2.0*pi,ncirc1, endpoint = False, dtype = np.float64)
tc2 = np.linspace(0,2.0*pi,ncirc2, endpoint = False, dtype = np.float64)

pcirc1 = np.zeros((ncirc1,2), dtype = np.float64)
pcirc1[:,0] = c1x*np.ones(ncirc1) + rad1*np.cos(tc1)
pcirc1[:,1] = c1y*np.ones(ncirc1) + rad1*np.sin(tc1)

pcirc2 = np.zeros((ncirc2,2), dtype = np.float64)
pcirc2[:,0] = c2x*np.ones(ncirc2) + rad2*np.cos(tc2)
pcirc2[:,1] = c2y*np.ones(ncirc2) + rad2*np.sin(tc2)

norcirc1 = np.zeros((ncirc1,2), dtype = np.float64)
norcirc1[:,0] = np.cos(tc1)
norcirc1[:,1] = np.sin(tc1)

norcirc2 = np.zeros((ncirc2,2), dtype = np.float64)
norcirc2[:,0] = np.cos(tc2)
norcirc2[:,1] = np.sin(tc2)

# plot
'''
plt.scatter(pq1[:,0], pq1[:,1])
plt.scatter(pq2[:,0], pq2[:,1])
plt.scatter(pq3[:,0], pq3[:,1])
plt.scatter(pq4[:,0], pq4[:,1])
plt.scatter(pcirc1[:,0], pcirc1[:,1])
plt.scatter(pcirc2[:,0], pcirc2[:,1])
plt.show() '''

## write 
fout = open(ofile, 'w')
fout.write(str(3)+'\n\n')
fout.write(str(np1+np2+np3+np4) + '\n' + str(ncirc1) + '\n' + str(ncirc2) + '\n\n')
fout.close()
fout = open(ofile, 'ab')
np.savetxt(fout, pq1)
np.savetxt(fout, pq2)
np.savetxt(fout, pq3)
np.savetxt(fout, pq4)
fout.write(b'\n')
np.savetxt(fout, nq1)
np.savetxt(fout, nq2)
np.savetxt(fout, nq3)
np.savetxt(fout, nq4)
fout.write(b'\n\n')
np.savetxt(fout, pcirc1)
fout.write(b'\n')
np.savetxt(fout, norcirc1)
fout.write(b'\n\n')
np.savetxt(fout, pcirc2)
fout.write(b'\n')
np.savetxt(fout, norcirc2)
fout.close()