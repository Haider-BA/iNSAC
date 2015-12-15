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
	/// mesh associated with the simulation
	Structmesh2d* m;
	
	/// number of vertex-neighbors of each cell in the mesh; 8 for a structured quad mesh
	int ncellvnbd;
	///< stores i-offsets for a local numbering of neighboring cells
	vector<int> idif;
	///< stores j-offsets for a local numbering of neighboring cells
	vector<int> jdif;
	
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
	Array2d<double> hh;				///< 'Heavyside function' for each cell; it is 1.0 if the cell is a band cell or an interior cell, and zero for field cells.

	/// Weight of each surrounding cell in the interpolation point, for every band cell
	vector<Array2d<double>> weights;

	/// Exponent for the tangential velocity interpolation function
	/** 1.0 for laminar flow, 1/7 for turbulent flow */
	double p;

	/// data for interpolation of tangential velocity to the interpolation point (which is not a grid point, in general)
	/** Store, for each band cell, the value of \$ (\frac{d_b}{d_I}})^p \$, where \$ d_b \$ is the distance of the band cell center from the closest immersed surface, and \$d_I \$ is the ditance of the interpolation point from the immersed surface.
	 */
	Array2d<double> acoef;

	/// data for computation of normal velocity at interpolation point for each band cell
	Array2d<double> bcoef;

public:

	/// Reads immersed object data from a file and allocates class variables
	void setup(Structmesh2d* m, string objectfile);

	/// set up neighborhood data
	void setnbd();

	/// Classify cells into field, band and interior cells
	/** Compute distances, signed distances and global distances, and the [heavyside function](@ref hh) */
	void classify();

	/// Compute weights of surrounding cells in the interpolation point for each band cell, and acoef and bcoef
	void compute_interpolation_data();
	
	int gnobj() const;
	int gnspobj(int iobj) const;
	int gncellvnbd() const;
	int gidif(int k) const;
	int gjdif(int k) const;
	double gspx(int n_object, int ipoin) const;
	double gspy(int n_object, int ipoin) const;
	double gsnx(int n_object, int ipoin) const;
	double gsny(int n_object, int ipoin) const;
	double gdistg(int i, int j) const;
	double gitag(int i, int j, int n_object) const;
	double glpri(int i, int j) const;
	double ghh(int i, int j) const;
	double gacoef(int i, int j) const;
	double gbcoef(int i, int j) const;
	double gweights(int i, int j, int n_object) const;

	/// Can be used to update the IB geometry.
	/** The arguments are added to spx, spy, snx and sny respectively.
	 */
	void update_ib(const vector<vector<double>>& dspx, const vector<vector<double>>& dspy, const vector<vector<double>>& dsnx, const vector<vector<double>>& dsny);
	
	/// Used to update the position of the IB without changing the normals
	void displace_ib(const vector<vector<double>>& dspx, const vector<vector<double>>& dspy);
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
		dist[i].setup(m->gimx()+1,m->gjmx()+1);
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

	setnbd();
}

void ImmersedBoundary::setnbd()
{
	ncellvnbd = 8;
	idif.resize(ncellvnbd);
	jdif.resize(ncellvnbd);

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
		for(i = 1; i <= m->gimx()-1; i++)
			for(j = 1; j <= m->gjmx()-1; j++)
			{
				dist[iobj](i,j) = BIG_NUMBER;
				for(k = 0; k < nspobj[iobj]; k++)
				{
					dd = pow(m->gxc(i,j)-spx[iobj][k],2) + pow(m->gyc(i,j)-spy[iobj][k], 2);
					if(dd < dist[iobj](i,j))
					{
						dist[iobj](i,j) = dd;
						itag[iobj](i,j) = k;
					}
				}
				dist[iobj](i,j) = sqrt(dist[iobj](i,j));

				/** To get the sign, we compute the dot product of the vector from the closest surface point to the grid point in question, and the surface normal vector at that point. If this dot product is positive, the grid point is outside the immersed object, else it's inside. */
				
				dotp = ( m->gxc(i,j)-spx[iobj][itag[iobj](i,j)] ) * snx[iobj][itag[iobj](i,j)] + ( m->gyc(i,j)-spy[iobj][itag[iobj](i,j)] ) * sny[iobj][itag[iobj](i,j)];
				if(fabs(dotp) < ZERO_TOL)
					dist[iobj](i,j) = 0;
				else
					dist[iobj](i,j) *= dotp/fabs(dotp);
			}
	}
	
	/** Next, we compute the global signed distance function distg . If distg is positive or zero for a cell, it is treated as a field/band cell; if distg is negative, it's an interior cell.
	*/
	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
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
	int icflag, iq, jq;
	
	// iterate over real cells first
	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			hh(i,j) = 0;
			icflag = 0;
			for(k = 0; k < ncellvnbd; k++)
			{
				iq = i + idif[k];
				jq = j + jdif[k];
				if(distg(i,j) >= 0 && distg(iq,jq) < 0) icflag = 1;
			}
			if(icflag == 1)
				hh(i,j) = 1.0;
			if(distg(i,j) < 0) hh(i,j) = 1.0;
		}
	// so now, hh is 1.0 for interior cells (cells inside the immersed objects) and band cells
	
	// copy hh values to ghost cells
	for(i = 1; i <= m->gimx()-1; i++)
	{
		hh(i,0) = hh(i,1);
		hh(i,m->gjmx()) = hh(i,m->gjmx()-1);
	}
	for(j = 1; j <= m->gjmx()-1; j++)
	{
		hh(0,j) = hh(1,j);
		hh(m->gimx(),j) = hh(m->gimx()-1,j);
	}
	// corner cells
	hh(0,0) = hh(1,1);
	hh(0,m->gjmx()) = hh(1,m->gjmx()-1);
	hh(m->gimx(),0) = hh(m->gimx()-1,1);
	hh(m->gimx(),m->gjmx()) = hh(m->gimx()-1,m->gjmx()-1);

	//lpri.mprint();
	//itag[0].mprint();
}

void ImmersedBoundary::compute_interpolation_data()
{
	int i,j,k, iq,jq;

	weights.resize(ncellvnbd);
	for(i = 0; i < ncellvnbd; i++)
		weights[i].setup(m->gimx()+1,m->gjmx()+1);
	acoef.setup(m->gimx()+1,m->gjmx()+1);
	bcoef.setup(m->gimx()+1,m->gjmx()+1);

	// interpolation function exponent for laminar flow
	p = 1.0;

	// vector from band cell center to neighboring cell center, called d_k in the notes
	double dx, dy;
	// dot product of d_k and n (normal to immersed surface passing through current band cell center)
	double deld;
	// some values needed
	double deld1 = 0, deld2 = 0;
	// perpendicular distance from a neighboring cell center to the surface normal drawn from the surface point nearest to the current band cell
	double dcross;
	// normal, at a band cell, of the nearest immersed surface
	double nxd, nyd;
	// band cell centers
	double xcd, ycd;
	// centers of cells surrounding the a band cell
	double xcp, ycp;

	double dist_plus, dist_minus, dratio, term;
	
	double eps = 1e-8;

	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			// check if cell i,j is a band cell; CAREFUL HERE - checking equality of doubles
			if(hh(i,j)==1.0 && distg(i,j) >= 0.0)
			{
				for(k = 0; k < ncellvnbd; k++)
					weights[k](i,j) = 0.0;
				// get surface normal passing through band cell center and the coords of band cell center
				nxd = snx [ lpri(i,j) ] [ itag[lpri(i,j)](i,j) ];
				nyd = sny [ lpri(i,j) ] [ itag[lpri(i,j)](i,j) ];
				xcd = m->gxc(i,j);
				ycd = m->gyc(i,j);

				for(k = 0; k < ncellvnbd; k++)
				{
					iq = i + idif[k];
					jq = j + jdif[k];
					
					xcp = m->gxc(iq,jq); ycp = m->gyc(iq,jq);
					dx = xcp - xcd; dy = ycp - ycd;
					
					deld = dx*nxd + dy*nyd;
					dcross = sqrt(dx*dx + dy*dy - deld*deld);
					
					// due to some reason, check if neighboring cell center is further away from the surface than the band cell center
					// check if the neighboring cell is a field cell
					// again CAREFUL, testing equality of doubles
					if(deld > 0 && hh(iq,jq) == 0)
					{
						weights[k](i,j) = 1.0/(dcross + eps);
						deld1 += weights[k](i,j);
						deld2 += weights[k](i,j)*deld;
					}
				}

				if(fabs(deld1) < ZERO_TOL)
					for(k = 0; k < ncellvnbd; k++)
					{
						iq = i + idif[k];
						jq = j + jdif[k];
						
						xcp = m->gxc(iq,jq); ycp = m->gyc(iq,jq);
						dx = xcp - xcd; dy = ycp - ycd;
						
						deld = dx*nxd + dy*nyd;
						dcross = sqrt(dx*dx + dy*dy - deld*deld);
						
						// due to some reason, check if neighboring cell center is further away from the surface than the band cell center
						// slight difference from Dr Edwards' code here, but should amount to the same thing since hh(i,j) >= 0 is always true.
						if(deld > 0)
						{
							weights[k](i,j) = 1.0/(dcross + eps);
							deld1 += weights[k](i,j);
							deld2 += weights[k](i,j)*deld;
						}

					}
				
				// normalize the weights, and finally get deld as the distance from the band cell center to the interpolation point (d_IB)
				if(fabs(deld1) > ZERO_TOL)
				{
					for(k = 0; k < ncellvnbd; k++)
						weights[k](i,j) /= deld1;
					deld = deld2/deld1;
				}
				else
				{
					cout << "ImmersedBoundary: compute_interpolation_data(): No interpolation point found for band cell (" << i << "," << j << ")" << endl;
					deld = 0;
				}
				
				// finally get coefficients for tangential (acoef) and normal (bcoef) velocity interpolation
				//
				// dist_plus is distance (from surface) of a point midway between B (band cell center) and I (interpolation point)
				// dist_minus is half the distance from surface to B
				// dratio is the rato of d_IB to d_B

				dratio = deld / distg(i,j);
				acoef(i,j) = pow(1.0 / (1.0+dratio), p);
				dist_plus = distg(i,j) + 0.5*deld;
				dist_minus = 0.5*distg(i,j);
				term = pow(dist_minus/dist_plus, p) / dratio;
				bcoef(i,j) = term/(1.0 + term);					// for more complicated normal velocity that attempts to satisfy continuity
				// bcoef(i,j) = distg(i,j)/(deld+distg(i,j))	// for linear interpolation, same as tangential component
			}
		}
}

// accessor functions
int ImmersedBoundary::gnobj() const
{ return nobj; }

int ImmersedBoundary::gnspobj(int iobj) const
{ return nspobj.at(iobj); }

int ImmersedBoundary::gncellvnbd() const
{ return ncellvnbd; }

int ImmersedBoundary::gidif(int k) const
{ return idif.at(k); }

int ImmersedBoundary::gjdif(int k) const
{ return jdif.at(k); }

double ImmersedBoundary::gspx(int n_object, int ipoin) const
{ return spx.at(n_object).at(ipoin); }

double ImmersedBoundary::gspy(int n_object, int ipoin) const
{ return spy[n_object][ipoin]; }

double ImmersedBoundary::gsnx(int n_object, int ipoin) const
{ return snx.at(n_object).at(ipoin); }

double ImmersedBoundary::gsny(int n_object, int ipoin) const
{ return sny.at(n_object).at(ipoin); }

double ImmersedBoundary::gdistg(int i, int j) const
{ return distg.get(i,j); }

double ImmersedBoundary::gitag(int i, int j, int n_object) const
{ return itag[n_object].get(i,j); }

double ImmersedBoundary::glpri(int i, int j) const
{ return lpri.get(i,j); }

double ImmersedBoundary::ghh(int i, int j) const
{ return hh.get(i,j); }

double ImmersedBoundary::gacoef(int i, int j) const
{ return acoef.get(i,j); }

double ImmersedBoundary::gbcoef(int i, int j) const
{ return bcoef.get(i,j); }

double ImmersedBoundary::gweights(int i, int j, int n_nbr) const
{ return weights.at(n_nbr).get(i,j); }

void ImmersedBoundary::update_ib(const vector<vector<double>>& dspx, const vector<vector<double>>& dspy, const vector<vector<double>>& dsnx, const vector<vector<double>>& dsny)
{
	for(int iob = 0; iob < nobj; iob++)
	{
		for(int i = 0; i < nspobj[iob]; i++)
		{
			spx[iob][i] += dspx.at(iob).at(i);
			spy[iob][i] += dspy.at(iob).at(i);
			snx[iob][i] += dsnx.at(iob).at(i);
			sny[iob][i] += dsny.at(iob).at(i);
		}
	}
}

void ImmersedBoundary::displace_ib(const vector<vector<double>>& dspx, const vector<vector<double>>& dspy)
{
	for(int iob = 0; iob < nobj; iob++)
	{
		for(int i = 0; i < nspobj[iob]; i++)
		{
			spx[iob][i] += dspx.at(iob).at(i);
			spy[iob][i] += dspy.at(iob).at(i);
		}
	}
}

} // end namespace
#endif
