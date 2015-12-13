#! /usr/bin/env python

''' @brief Generates surface points and normals for an ellipse-shaped immersed boundary.
The parametric equation is
r(t) = r_0 + a cos(t) i + b sin(t) j.
So the normal is
n(t) = b cos(t) i + a sin(t) j.
'''

from math import sqrt, pi
import sys
import numpy as np

if(len(sys.argv) < 7):
	print("Not enough arguments! Need five arguments:")
	print("\tCenter x, center y, semi-major axis, semi-minor axis, number of points, output file name.")
	sys.exit()

# read data from command line
cx = float(sys.argv[1])
cy = float(sys.argv[2])
a = float(sys.argv[3])
b = float(sys.argv[4])
npoin = int(sys.argv[5])
ofile = sys.argv[6]

# compute points and normals
t = np.linspace(0,2.0*pi,npoin, endpoint = False, dtype = np.float64)

pts = np.zeros((npoin,2), dtype = np.float64)
pts[:,0] = cx*np.ones(npoin) + a*np.cos(t)
pts[:,1] = cy*np.ones(npoin) + b*np.sin(t)

normmag = sqrt(a**2 + b**2)

nors = np.zeros((npoin,2), dtype = np.float64)
nors[:,0] = b/normmag * np.cos(t)
nors[:,1] = a/normmag * np.sin(t)

# write data to file
fout = open(ofile, 'w')
fout.write(str(1)+'\n\n')
fout.write(str(npoin) + '\n\n')
fout.close()
fout = open(ofile, 'ab')
np.savetxt(fout, pts)
fout.write(b'\n')
np.savetxt(fout, nors)
fout.close()