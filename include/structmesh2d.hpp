/** 
\file structmesh2d.hpp
\brief This file provides a handler class for structured meshes.
\author Aditya Kashi
*/

// Include for array storage and operations
#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

#define __STRUCTMESH2D_H 1

using namespace amat;
using namespace std;

namespace acfd {

/**	\brief This class encapsulates the 2D structured mesh.

	It reads from a mesh file and computes coordinates of cell centres, cell volumes and area vectors of each face.
	NOTE that the index of a cell is same as that of its lower-left corner node, as opposed to the convention in class.
*/
class Structmesh2d
{
	Array2d<double> x;
	Array2d<double> y;
	int imx;			///< Number of points along x
	int jmx;			///< Number of points along y

	Array2d<double> xc;	///< x-coordinates of cell centres
	Array2d<double> yc;	///< y-coordinates of cell centres
	
	/**	\brief Area vectors of faces 1 and 2 for each face.
	
		del[0](i,j) and del[1](i,j) are x- and y-components of area of face 1; 
		del[2](i,j) and del[3](i,j) are x- and y-components of area of face 2.
		Note that the vectors point along +i direction for the for the i-faces and the +j direction for the j-faces (rather than -i and -j).
	*/
	Array2d<double>* del;
	int ndelcomp;			///< Number of components of del
	
	bool allocdel;			///< Stores whether or not del has been allocated.

	Array2d<double> vol;		///< Contains the measure (area in 2D) of each cell.

public:
	Structmesh2d();
	Structmesh2d(int nxpoin, int nypoin);

	void setup(int nxpoin, int nypoin);

	~Structmesh2d();

	void readmesh(string fname);

	void writemesh_vtk(string fname) const;

	/** \brief Calculates cell centers, volumes and area vectors and store in xc, yc, vol  and del respectively. */
	void preprocess();

	int gimx() const;						///< Returns the number of (real) grid points in i-direction.
	int gjmx() const;						///< Returns the number of (real) grid points in j-direction.
	double gx(int i, int j) const;			///< Returns x-coordinates of grid points.
	double gy(int i, int j) const;			///< Returns y-coordinates of grid points.
	double gxc(int i, int j) const;			///< Returns x-coordinates of cell centers.
	double gyc(int i, int j) const;			///< Returns y-coordinates of cell centers.
	double gdel(int i, int j, int idat) const;		///< Returns del[idat](i,j); refer to [del](@ref del)
	double gvol(int i, int j) const;		///< Returns 'volune' (area) of a cell.
};

Structmesh2d::Structmesh2d()
{
	allocdel = false;
	ndelcomp = 4;
}

Structmesh2d::Structmesh2d(int nxpoin, int nypoin)
{
	imx = nxpoin;
	jmx = nypoin;
	/// Allocate space for point coordinates. Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);
	/// Allocate space for cell centers. Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);
	vol.setup(imx+1,jmx+1);

	ndelcomp = 4;
	del = new Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
}

Structmesh2d::~Structmesh2d()
{
	if(allocdel)
	{
		delete [] del;
		allocdel = false;
	}
}

void Structmesh2d::setup(int nxpoin, int nypoin)
{
	imx = nxpoin;
	jmx = nypoin;
	/// Allocate space for point coordinates. Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);
	/// Allocate space for cell centers. Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);
	vol.setup(imx+1,jmx+1);

	ndelcomp = 4;
	del = new Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
}

void Structmesh2d::readmesh(string meshname)
{
	ifstream fin(meshname);
	fin >> imx >> jmx;

	/// Allocate space for point coordinates. Two extra in each direction for ghost cells.
	x.setup(imx+2,jmx+2);
	y.setup(imx+2,jmx+2);
	/// Allocate space for cell centers. Space required = (imx-1) real cells + 2 ghost cells; similarly for j direction.
	xc.setup(imx+1,jmx+1);
	yc.setup(imx+1,jmx+1);
	
	if(allocdel)
		delete [] del;

	del = new Array2d<double>[ndelcomp];
	for(int i = 0; i < ndelcomp; i++)
		del[i].setup(imx+1,jmx+1);
	allocdel = true;
	vol.setup(imx+1,jmx+1);

	for(int j = 1; j <= jmx; j++)
		for(int i = 1; i <= imx; i++)
			fin >> x(i,j) >> y(i,j);
	
	std::cout << "Structmesh2d: readmesh(): Mesh read. Points in i-dir: " << imx << ", points in j-dir: " << jmx << "." << endl;
}

void Structmesh2d::preprocess()
{
	/// Add ghost points. The ghost point is such that the boundary point is the arithmetic mean of the ghost point and the first interior point.
	for(int j = 1; j <= jmx; j++)
	{
		x(0,j) = 2*x(1,j) - x(2,j);
		y(0,j) = 2*y(1,j) - y(2,j);
		x(imx+1,j) = 2*x(imx,j) - x(imx-1,j);
		y(imx+1,j) = 2*y(imx,j) - y(imx-1,j);
	}
	for(int i = 1; i <= imx; i++)
	{
		x(i,0) = 2*x(i,1) - x(i,2);
		y(i,0) = 2*y(i,1) - y(i,2);
		x(i,jmx+1) = 2*x(i,jmx) - x(i,jmx-1);
		y(i,jmx+1) = 2*y(i,jmx) - y(i,jmx-1);
	}
	// corner points:
	x(0,0) = 2*x(1,0) - x(2,0);
	y(0,0) = 2*y(1,0) - y(2,0);
	x(0,jmx+1) = 2*x(1,jmx+1) - x(2,jmx+1);
	y(0,jmx+1) = 2*y(1,jmx+1) - y(2,jmx+1);
	x(imx+1,0) = 2*x(imx,0) - x(imx-1,0);
	y(imx+1,0) = 2*y(imx,0) - y(imx-1,0);
	x(imx+1,jmx+1) = 2*x(imx,jmx+1) - x(imx-1,jmx+1);
	y(imx+1,jmx+1) = 2*y(imx,jmx+1) - y(imx-1,jmx+1);

	/// Calculate cell centers, volumes and face areas. We calculate volumes by Heron's formula.
	cout << "Structmesh2d: preprocess(): Calculating cell centers, volume and del" << endl;
	double a,b,c,s;
	for(int i = 0; i <= imx; i++)
		for(int j = 0; j <= jmx; j++)
		{
			if((i>0 || j>0) && (i < imx || j < imx))	// ignore corner ghost cells 0,0 and imx,jmx
			{
				// set cell centers
				xc(i,j) = 0.25*(x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1));
				yc(i,j) = 0.25*(y(i,j) + y(i+1,j) + y(i,j+1) + y(i+1,j+1));
				// next, calculate volumes
				a = sqrt( (x(i,j)-x(i,j+1))*(x(i,j)-x(i,j+1)) + (y(i,j)-y(i,j+1))*(y(i,j)-y(i,j+1)) );
				b = sqrt( (x(i,j+1)-x(i+1,j+1))*(x(i,j+1)-x(i+1,j+1)) + (y(i,j+1)-y(i+1,j+1))*(y(i,j+1)-y(i+1,j+1)) );
				c = sqrt( (x(i,j)-x(i+1,j+1))*(x(i,j)-x(i+1,j+1)) + (y(i,j)-y(i+1,j+1))*(y(i,j)-y(i+1,j+1)) );
				s = (a+b+c)*0.5;
				vol(i,j) = sqrt(s*(s-a)*(s-b)*(s-c));
				
				a = sqrt( (x(i+1,j)-x(i+1,j+1))*(x(i+1,j)-x(i+1,j+1)) + (y(i+1,j)-y(i+1,j+1))*(y(i+1,j)-y(i+1,j+1)) );
				b = sqrt( (x(i,j)-x(i+1,j))*(x(i,j)-x(i+1,j)) + (y(i,j)-y(i+1,j))*(y(i,j)-y(i+1,j)) );
				s = (a+b+c)*0.5;
				vol(i,j) += sqrt(s*(s-a)*(s-b)*(s-c));
			}

			// Compute face area vectors
			// Note that face normals point in the positive i or positive j directions.
			del[0](i,j) = y(i+1,j+1) - y(i+1,j);
			del[1](i,j) = -(x(i+1,j+1) - x(i+1,j));
			del[2](i,j) = -(y(i+1,j+1) - y(i,j+1));
			del[3](i,j) = x(i+1,j+1) - x(i,j+1);
		}
}

int Structmesh2d::gimx() const { return imx; }
int Structmesh2d::gjmx() const { return jmx; }

double Structmesh2d::gx(int i, int j) const
{ return x.get(i,j); }

double Structmesh2d::gy(int i, int j) const
{ return y.get(i,j); }

double Structmesh2d::gxc(int i, int j) const
{ return xc.get(i,j); }

double Structmesh2d::gyc(int i, int j) const
{ return yc.get(i,j); }

double Structmesh2d::gdel(int i, int j, int idat) const
{ return del[idat].get(i,j); }

double Structmesh2d::gvol(int i, int j) const
{ return vol.get(i,j); }

} // end namespace
