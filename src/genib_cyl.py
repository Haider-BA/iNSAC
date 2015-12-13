#! /usr/bin/env python

''' @brief Generates surface points and normals for a disk-shaped immersed boundary.
'''

from math import sqrt, pi
import sys
import numpy as np

if(len(sys.argv) < 6):
	print("Not enough arguments! Need five arguments:")
	print("\tCenter x, center y, radius, number of points, output file name.")
	sys.exit()

cx = float(sys.argv[1])
cy = float(sys.argv[2])
rad = float(sys.argv[3])
npoin = int(sys.argv[4])
ofile = sys.argv[5]

t = np.linspace(0,2.0*pi,npoin, dtype = np.float64)

pts = np.zeros((npoin,2), dtype = np.float64)
pts[:,0] = cx*np.ones(npoin) + rad*np.cos(t)
pts[:,1] = cy*np.ones(npoin) + rad*np.sin(t)

nors = np.zeros((npoin,2), dtype = np.float64)
nors[:,0] = np.cos(t)
nors[:,1] = np.sin(t)

fout = open(ofile, 'w')
fout.write(str(1)+'\n\n')
fout.write(str(npoin) + '\n\n')
fout.close()
fout = open(ofile, 'ab')
np.savetxt(fout, pts)
fout.write(b'\n')
np.savetxt(fout, nors)
fout.close()