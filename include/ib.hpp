/** @file ib.hpp
 * @brief Basic infrastructure for immersed boundary (IB) capability in 2D structured grid finite volume code.
 * @author Aditya Kashi
 * @date Dec 11, 2015
 */

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _ACONSTANTS_H
#include <aconstants.h>
#endif

#ifndef __STRUCTMESH2D_H
#include <structmesh2d.hpp>
#endif

#ifndef __IB_H
#define __IB_H 1

using amat::Array2d;
using namespace acfd;
using namespace std;

namespace acfd {

/// Class to handle immersed boundary (IB) in finite volume code.
/** Implements
 * - surface detection, ie., classification of cells into field cells, band cells and interior cells
 * - imposition of boundary condition on IB
 */
class ImmersedBoundary
{
	Structmesh2d* m;				///< mesh associated with the simulation
	
	static int ncellvnbd;			///< number of vertex-neighbors of each cell in the mesh; 8 for a structured quad mesh
	static vector<int> idif;		///< stores i-offsets for a local numbering of neighboring cells
	static vector<int> jdif;		///< stores j-offsets for a local numbering of neighboring cells
	
	int nobj;						///< number of immersed objects
	vector<int> nspobj;				///< number of given surface points for each object
	
	/// array to hold x-coordinates of surface points for each immersed object
	/** We choose a vector of vector<double> over an Array2d<double> as each immersed object can have a different number of surface points. */
	vector<vector<double>> spx;

	/// array to hold y-coordinates of surface points for each immersed object
	/** We choose a vector of vector<double> over an Array2d<double> as each immersed object can have a different number of surface points. */
	vector<vector<double>> spy;

	vector<vector<double>> snx;		///< x-coordinates of surface normals at points corresponding to spx and spy
	vector<vector<double>> sny;		///< y-coordinates of surface normals at points corresponding to spx and spy

	vector<Array2d<int>> itag;		///< itag[iobj](i,j) points to that surface point number of immersed object 'iobj' that is closest to mesh cell (i,j).
	Array2d<int> lpri;				///< contains the number of the closest immersed object for each mesh cell
	vector<Array2d<double>> dist;	///< stores distance to closest surface point of each immersed object for each mesh cell
	Array2d<double> distg;			///< stores the distance of each mesh cell to the closest surface point of the *closest* immersed object - 'global dist'
	Array2d<double> hh;				///< 'Heavyside function' for each cell; it is 1.0 if the cell is a band cell or an interior cell, and zero otherwise.

public:

	/// Reads immersed object data from a file and allocates class variables
	void setup(Structmesh2d* m, string objectfile);

	/// set up neighborhood data
	static void setnbd();

	/// Classify cells into field, band and interior cells
	/** Compute distances, signed distances and global distances, and the [heavyside function](@ref hh) */
	void classify();
};

void ImmersedBoundary::setup(Structmesh2d* mesh, string objectfile)
{
	m = mesh;
	
	cout << "ImmersedBoundary: setup(): Opening immersed objects file " << objectfile << endl;
	ifstream fin(objectfile);

	fin >> nobj;
	cout << "ImmersedBoundary: setup(): Number of objects = " << nobj << endl;

	nspobj.resize(nobj);
	spx.resize(nobj);
	spy.resize(nobj);
	snx.resize(nobj);
	sny.resize(nobj);
	itag.resize(nobj);
	dist.resize(nobj);

	for(int i = 0; i < nobj; i++)
	{
		fin >> nspobj[i];
		spx[i].resize(nspobj[i]);
		spy[i].resize(nspobj[i]);
		snx[i].resize(nspobj[i]);
		sny[i].resize(nspobj[i]);
		itag[i].setup(m->gimx()+1,m->gjmx()+1);
		distg[i].setup(m->gimx()+1;m->gjmx()+1);
	}
	lpri.setup(m->gimx()+1,m->gjmx()+1);
	distg.setup(m->gimx()+1,m->gjmx()+1);
	hh.setup(m->gimx()+1,m->gjmx()+1);

	for(int iobj = 0; iobj < nobj; iobj++)
	{
		// read surface point coords
		for(int ip = 0; ip < nspobj[iobj]; ip++)
			fin >> spx[iobj][ip] >> spy[iobj][ip];
		// read surface normal components
		for(int ip = 0; ip < nspobj[iobj]; ip++)
			fin >> snx[iobj][ip] >> sny[iobj][ip];
	}

	fin.close();
}

static void ImmersedBoundary::setnbd()
{
	ncellvnbd = 8;
	idif.reserve(ncellvnbd);
	jdif.reserve(ncellvnbd);

	idif[0] = -1;	jdif[0] = 1;
	idif[1] = 0;	jdif[1] = 1;
	idif[2] = 1;	jdif[2] = 1;
	idif[3] = 1;	jdif[3] = 0;
	idif[4] = 1;	jdif[4] = -1;
	idif[5] = 0;	jdif[5] = -1;
	idif[6] = -1;	jdif[6] = -1;
	idif[7] = -1;	jdif[7] = 0;
}

void ImmersedBoundary::classify()
{
	/** The union of the IB surfaces is a zero-level set of a signed distance funtion. Thus we need a signed distance function from each grid point to the nearest IB surface point.*/
	int i,j,k, iob;
	double dd, dotp;

	/// first, we compute the signed distance function dist
	for(int iobj = 0; iobj < nobj; iobj++)
	{
		for(j = 0; j <= m->gjmx(); j++)
			for(i = 0; i <= m->gimx(); i++)
			{
				dist[iobj](i,j) = BIG_NUMBER;
				for(k = 0; k < nspobj[iobj]; k++)
				{
					dd = pow(m->gxc(i,j)-spx[iobj][k],2) + pow(m->gyc(i,j)-spy[iobj][k], 2)
					if(dd < dist[iobj](i,j))
					{
						dist[iobj](i,j) = dd;
						itag[iobj](i,j) = k;
					}
				}
				dist[iobj](i,j) = sqrt(dist[iobj](i,j));

				/** To get the sign, we compute the dot product of the vector from the closest surface point to the grid point in question, and the surface norma vector at that point. If this dot product is positive, the grid point is outside the immersed object, else it's inside. */
				
				dotp = ( m->gxc(i,j)-spx[iobj][itag[iobj](i,j)] ) * snx[iobj][itag[iobj](i,j)] + ( m->gyc(i,j)-spy[iobj][itag[iobj](i,j)] ) * sny[iobj][itag[iobj](i,j)];
				if(fabs(dotp) < ZERO_TOL)
					dist[iobj](i,j) = 0;
				else
					dist[iobj](i,j) *= dotp/fabs(dotp);
			}
	}
	
	/** Next, we compute the global signed distance function distg . If distg is positive or zero for a cell, it is treated as a field/band cell; if distg is negative, it's an interior cell.
	*/
	for(j = 0; j <= m->gjmx(); j++)
		for(i = 0; i <= m->gimx(); i++)
		{
			distg(i,j) = BIG_NUMBER;
			lpri(i,j) = 1;				// some initial assumption for closest immersed object

			for(iob = 0; iob < nobj; iob++)
				if(dist[iob](i,j) < distg(i,j))
				{
					distg(i,j) = dist[iob](i,j);
					lpri(i,j) = iob;
				}
		}
	
	/** Finally, we compute the heavyside function hh to complete the classification. */
}

} // end namespace
#endif
