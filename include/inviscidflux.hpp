/** This file contains inviscid flux schemes for the artificial compressiblity (AC) formulation of the incompressible Navier-Stokes equations.
*/

#ifndef __AARRAY2D_H
#include <aarray2d.hpp>
#endif

#ifndef __STRUCTMESH2D_H
#include <structmesh2d.hpp>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#define __INVISCIDFLUX_H 1

namespace acfd {

double minmod_avg(double a, double b)
{
	double signa, signb;
	if(fabs(a) <= ZERO_TOL) signa = 0;
		else signa = a/fabs(a);
	if(fabs(b) <= ZERO_TOL) signb = 0;
		else signb = b/fabs(b);

	if(fabs(a) <= fabs(b))
		return (signa+signb)*0.5*fabs(a);
	else
		return (signa+signb)*0.5*fabs(b);
}

/// Abstract class to compute pressure difference dp across each face as (p_L - p_R)
class PressureReconstruction
{
protected:
	Structmesh2d* m;
	Array2d<double>* u;
	Array2d<double>* dp;
public:
	void setup(Structmesh2d* mesh, Array2d<double>* unknowns, Array2d<double>* delp);
	virtual void compute_pressure_difference() = 0;
};

void PressureReconstruction::setup(Structmesh2d* mesh, Array2d<double>* unknowns, Array2d<double>* delp)
{
	cout << "PressureReconstruction: Contructor called." << endl;
	m = mesh;
	u = unknowns;
	dp = delp;
}

/** \brief Implements the basic first-order pressure reconstruction needed for the mass flux in the Rhie-Chow method.
*
* Computes pressure difference across each face in the mesh (including faces shared between ghost cells).
*/
class BasicPR : public PressureReconstruction
{
public:
	void compute_pressure_difference();
};

void BasicPR::compute_pressure_difference()
{
	// we include the i=0 and j=0 ghost cells
	for(int j = 0; j <= m->gjmx()-1; j++)
		for(int i = 0; i <= m->gimx()-1; i++)
		{
			dp[0](i,j) = u[0].get(i,j) - u[0].get(i+1,j);		// p_L - p_R
			dp[1](i,j) = u[0].get(i,j) - u[0].get(i,j+1);		// p_L - p_R in j-direction, so it's more like p_down - p_up
		}
}

/** \brief Implements a second-order accurate TVD pressure reconstruction for the mass flux needed for the Rhie-Chow scheme for colocated grid.
*/
class TVDPR : public PressureReconstruction
{
public:
	void compute_pressure_difference();
};

void TVDPR::compute_pressure_difference()
{
	int i,j;
	for(j = 1; j <= m->gjmx()-2; j++)
		for(i = 1; i <= m->gimx()-2; i++)
		{
			dp[0](i,j) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i+1,j)-u[0](i,j), u[0].get(i,j)-u[0].get(i-1,j));				// p_L
			dp[0](i,j) -= u[0].get(i+1,j) - 0.5*minmod_avg(u[0].get(i+1,j)-u[0].get(i,j), u[0].get(i+2,j)-u[0].get(i+1,j));		// p_R
			dp[1](i,j) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i,j+1)-u[0](i,j), u[0].get(i,j)-u[0].get(i,j-1));				// p_down
			dp[1](i,j) -= u[0].get(i,j+1) - 0.5*minmod_avg(u[0].get(i,j+1)-u[0].get(i,j), u[0].get(i,j+2)-u[0].get(i,j+1));		// p_up
		}

	// Now values for remaining cells (real or ghost) - just copy values from adjacent interior cells.
	// We do not need values for cells i = imx and j = jmx (ghost cells). We also probably don't need dp for corner cells.
	for(j = 1; j <= m->gjmx()-2; j++)
	{
		dp[0](0,j) = dp[0](1,j);
		dp[0](m->gimx()-1,j) = dp[0](m->gimx()-2,j);
		dp[1](0,j) = u[0].get(0,j) +0.5*minmod_avg(u[0].get(0,j+1)-u[0](0,j), u[0].get(0,j)-u[0].get(0,j-1));				// p_down
		dp[1](0,j) -= u[0].get(0,j+1) - 0.5*minmod_avg(u[0].get(0,j+1)-u[0].get(0,j), u[0].get(0,j+2)-u[0].get(0,j+1));		// p_up
		dp[1](m->gimx()-1,j) = u[0].get(m->gimx()-1,j) +0.5*minmod_avg(u[0].get(m->gimx()-1,j+1)-u[0](m->gimx()-1,j), u[0].get(m->gimx()-1,j)-u[0].get(m->gimx()-1,j-1));
		dp[1](m->gimx()-1,j) -= u[0].get(m->gimx()-1,j+1) -0.5*minmod_avg(u[0].get(m->gimx()-1,j+1)-u[0].get(m->gimx()-1,j), u[0].get(m->gimx()-1,j+2)-u[0].get(m->gimx()-1,j+1));
	}
	for(i = 1; i <= m->gimx()-2; i++)
	{
		dp[1](i,0) = dp[1](i,1);
		dp[1](i,m->gjmx()-1) = dp[1](i,m->gjmx()-2);
		dp[0](i,0) = u[0].get(i,0) +0.5*minmod_avg(u[0].get(i+1,0)-u[0](i,0), u[0].get(i,0)-u[0].get(i-1,0));				// p_L
		dp[0](i,0) -= u[0].get(i+1,0) - 0.5*minmod_avg(u[0].get(i+1,0)-u[0].get(i,0), u[0].get(i+2,0)-u[0].get(i+1,0));		// p_R
		dp[0](i,m->gjmx()-1) = u[0].get(i,j) +0.5*minmod_avg(u[0].get(i+1,j)-u[0](i,j), u[0].get(i,j)-u[0].get(i-1,j));
		dp[0](i,m->gjmx()-1) -= u[0].get(i+1,m->gjmx()-1) - 0.5*minmod_avg(u[0].get(i+1,m->gjmx()-1)-u[0].get(i,m->gjmx()-1), u[0].get(i+2,m->gjmx()-1)-u[0].get(i+1,m->gjmx()-1));
	}
}

/// Class to add contributions of inviscid fluxes to the residual.
/** Implements a Rhie-Chow stabilization in the mass flux, for which pressure reconstruction is needed.
* @see PressureReconstruction
*/
class InviscidFlux
{
	Structmesh2d* m;
	Array2d<double>* u;			///< Contains one (imx+1) x (jmx+1) array each for pressure, x-velocity and y-velocity.
	Array2d<double>* res;		///< Contains residuals corresponding to [u](@ref u).
	Array2d<double>* beta;		///< Artificial compressiblity factor for each cell.
	double rho;					///< Density of fluid.
	PressureReconstruction* pr;

	/// Pressure difference across each face.
	/** dp[0] contains pressure difference across i-faces. dp[1] contains pressure differences across j-faces. */
	Array2d<double>* dp;	
	
	double c;					///< Rhie-Chow constant.
	int nvar;					///< Number of unknowns.
	bool isallocdp;

	/** \brief Integers indicating the type of boundary for each of the 4 boundaries.
	*
	*	0. velocity inflow
	*	1. pressure outflow
	*	2. no-slip wall
	*	3. slip wall
	*/
	vector<int> bcflags;
	
	/// Values of prescribed quantities corresponding to the flags in bcflags.
	/** The value for each boundary consists of two numbers - this is needed for inflow velocity.
	* For other types of boundary, we just read the first value corresponding to that boundary.
	*/
	vector<vector<double>> bvalues;

	/// one flag for each of the 4 boundaries, indicating whether or not that boundary is a wall.
	vector<double> wall;

public:
	void setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residuals, Array2d<double>* _beta, double _rho, string pressurereconstruction, vector<int> _bcflag, vector<vector<double>> _bvalues);
	~InviscidFlux();
	/// Add the inviscid flux contribution to the [residual](@ref res)
	void compute_fluxes();
};

void InviscidFlux::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residuals, Array2d<double>* _beta, double _rho, string pressurereconstruction, vector<int> _bcflag, vector<vector<double>> _bvalues)
{
	cout << "InviscidFlux: setup(): Setting up inviscid flux calculator" << endl;
	m = mesh;
	u = unknown;
	res = residuals;
	beta = _beta;
	rho = _rho;
	bvalues = _bvalues;
	bcflags = _bcflag;
	dp = new Array2d<double>[2];
	for(int i = 0; i<2; i++)
		dp[i].setup(m->gimx()+1, m->gjmx()+1);
	isallocdp = true;
	
	/// Sets the pressure reconstruction scheme to be used based on the last argument - "basic" or "TVD".
	if(pressurereconstruction=="basic")
	{
		pr = new BasicPR;
		cout << "InviscidFlux: setup(): BasicPR selected.\n";
	}
	else if(pressurereconstruction=="tvd")
		pr = new TVDPR;
	else {
		cout << "InviscidFlux: setup(): Pressure reconstruction scheme requested does not exist. Choosing basic first-order scheme." << endl;
		pr = new BasicPR;
	}
	pr->setup(m, unknown, dp);
	
	/// Rhie-Chow constant [c](@ref c) is maximum 0.5.
	c = 0.5;
	
	nvar = 3;
	
	/*cout << "BC flags ";
	for(int i = 0; i < 4; i++)
		cout << bcflags[i] << " ";
	cout << "\nB values:\n";
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 2; j++)
			cout << bvalues[i][j] << " ";
		cout << endl;
	}*/
	
	// account for solid walls; ie bcflag values of 2 or 3
	wall.resize(4);		// one flag for each of the 4 boundaries.
	for(int i = 0; i < 4; i++)
		if(bcflags[i] == 2 || bcflags[i] == 3)
			wall[i] = 0.0;
		else
			wall[i] = 1.0;
	
	/*cout << "Wall:";
	for(int i = 0; i < 4; i++)
		cout << " " << wall[i];
	cout << endl;*/
}

InviscidFlux::~InviscidFlux()
{
	if(isallocdp) {
		delete [] dp;
		delete pr;
	}
}

void InviscidFlux::compute_fluxes()
{
	pr->compute_pressure_difference();

	// add inviscid flux contribution to residuals
	//cout << "InviscidFlux: compute_flux(): Computing inviscid fluxes now..." << endl;
	int i, j, k;
	double area, nx, ny, eigen, vdotn;
	double bhalf2; vector<double> uhalf(nvar);			// interface values for each unknown
	Array2d<double> g(m->gimx(),nvar);
	
	for(j = 1; j <= m->gjmx()-1; j++)
	{
		for(i = 0; i <= m->gimx()-1; i++)
		{
			area = sqrt(m->gdel(i,j,0)*m->gdel(i,j,0) + m->gdel(i,j,1)*m->gdel(i,j,1));
			nx = m->gdel(i,j,0)/area;
			ny = m->gdel(i,j,1)/area;
			for(k = 0; k < nvar; k++)
				uhalf[k] = 0.5*(u[k](i,j) + u[k](i+1,j));

			// get average of beta^2
			bhalf2 = 0.5*(beta->get(i,j)*beta->get(i,j) + beta->get(i+1,j)*beta->get(i+1,j));

			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			
			// now get eigenvalue for Rhie-Chow
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4.0*bhalf2) );

			// Boundary faces need to be treated differently to account for wall BCs.
			if(i == 0)
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2)*wall[2];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[2] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[2] + uhalf[0]*ny);
			}
			else if(i == m->gimx()-1)
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2)*wall[0];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[0] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[0] + uhalf[0]*ny);
			}
			else
			{
				g(i,0) = area*(rho*vdotn + 0.5*c*eigen*dp[0](i,j)/bhalf2);
				g(i,1) = area*(rho*vdotn*uhalf[1] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2] + uhalf[0]*ny);
			}
		}
		for(i = 1; i <= m->gimx()-1; i++)
			for(k = 0; k < nvar; k++)
				res[k](i,j) += g(i,k) - g(i-1,k);
	}

	// now we add contribution of j-fluxes
	Array2d<double> h(m->gjmx(),nvar);

	for(i = 1; i <= m->gimx()-1; i++)
	{
		for(j = 0; j <= m->gjmx()-1; j++)
		{
			area = sqrt(m->gdel(i,j,2)*m->gdel(i,j,2) + m->gdel(i,j,3)*m->gdel(i,j,3));
			nx = m->gdel(i,j,2)/area;
			ny = m->gdel(i,j,3)/area;
			for(k = 0; k < nvar; k++)
				uhalf[k] = 0.5*(u[k](i,j) + u[k](i,j+1));
			bhalf2 = 0.5*(beta->get(i,j)*beta->get(i,j) + beta->get(i,j+1)*beta->get(i,j+1));
			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4.0*bhalf2) );
			
			if(j == 0)
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2)*wall[3];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[3] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[3] + uhalf[0]*ny); 
			}
			else if(j == m->gjmx()-1)
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2)*wall[1];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[1] + uhalf[0]*ny); 
			}
			else
			{
				h(j,0) = area*(rho*vdotn + 0.5*c*eigen*dp[1](i,j)/bhalf2);
				h(j,1) = area*(rho*vdotn*uhalf[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2] + uhalf[0]*ny); 
			}
		}
		for(j = 1; j <= m->gjmx()-1; j++)
			for(k = 0; k < nvar; k++)
				res[k](i,j) += h(j,k) - h(j-1,k);
	}
}

} // end namespace
