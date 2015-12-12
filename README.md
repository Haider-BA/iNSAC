iNSAC
=====

Aditya Kashi, 2015

An 2nd-order finite volume solver for the 2D incompressible Navier Stokes equations on structured grids. The artificial compressibility approximation of the iNS equations is used.

Questionable parts
------------------

- Is beta calculated every time step or not?
- What is the reference velocity used for beta calculation?
- What is being done about corner ghost cells? We need it for parallel CV gradients.
	Ans: One could set them according to left and right boundaries.
