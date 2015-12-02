/** This file contains inviscid flux schemes for the artificial compressiblity (AC) formulation of the incompressible Navier-Stokes equations.
*/

#ifndef __AARRAY2D_H
#include <aarray2d.hpp>
#endif

#ifndef __STRUCTMESH2D_H
#include <structmesh.h>
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

/// Compute pressure difference dp across each face as (p_L - p_R)
class PressureReconstruction
{
	Structmesh2d* m;
	Array2d<double>* u;
	Array2d<double>* dp;
public:
	void setup(Structmesh2d* mesh, Array2d<double>* unknowns, Array2d<double>* delp);
	virtual void compute_pressure_difference() = 0;
};

void PressureReconstrucion::setup(Structmesh2d* mesh, Array2d<double>* unknowns, Array2d<double>* delp)
{
	m = mesh;
	u = unknowns;
	dp = delp;
}

/** \brief Implements the basic first-order pressure reconstruction for the mass flux according to Rhie-Chow method.
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

/** \brief Implements a second-order accurate TVD pressure reconstruction for the mass flux according to Rhie-Chow scheme for colocated grid.
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
		dp[1](m->gimx()-1,j) = u[0].get(,j) +0.5*minmod_avg(u[0].get(,j+1)-u[0](,j), u[0].get(,j)-u[0].get(,j-1));
		dp[1](m->gimx()-1,j) -= u[0].get(m->gimx()-1,j+1) - 0.5*minmod_avg(u[0].get(m->gimx()-1,j+1)-u[0].get(m->gimx()-1,j), u[0].get(m->gimx()-1,j+2)-u[0].get(m->gimx()-1,j+1));
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
	Array2d<double>* dp;		///< Pressure difference across each face.
	double c;					///< Rhie-Chow constant.
	int nvar;					///< Number of unknowns.
	bool isallocdp;

	/** \brief Integers indicating the tyoe of boundary for each of the 4 boundaries.
	*
	*	1. velocity inflow
	*	2. pressure outflow
	*	3. no-slip wall
	*	4. slip wall
	*/
	vector<int> bcflags;
	
	/// Values of prescribed quantities corresponding to the flags in bcflags.
	/** The value for each boundary consists of two numbers - this is needed for inflow velocity.
	* For other types of boundary, we just read the first value corresponding to that boundary.
	*/
	vector<vector<double>> bvalues;
public:
	void setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residuals, Array2d<double>* _beta, double _rho, string pressurereconstruction, vector<int> _bcflag, vector<vector<double> _bvalues);
	~InviscidFlux();
	/// Add the inviscid flux contribution to the [residual](@ref res)
	void compute_flux();
};

void InviscidFlux::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residuals, Array2d<double>* _beta, double _rho, string pressurereconstruction, vector<int> _bcflag, vector<vector<double> _bvalues)
{
	m = mesh;
	u = unknown;
	res = residuals;
	beta = _beta;
	rho = _rho;
	bvalues = _bvalues;
	bcflags = _bcflag;
	dp = new Array2d<double>(2);
	for(int i = 0; i<2; i++)
		dp[i].setup(m->gimx()+1, m->gjmx()+1);
	isallocdp = true;
	/// Sets the pressure reconstruction scheme to be used based on the last argument - "basic" or "TVD".
	switch(pressurereconstruction)
	{
		case "basic":
			pr = new BasicPR();
			break;
		case "TVD":
			pr = new TVDPR();
			break;
		default:
			cout << "InviscidFlux: setup(): Pressure reconstruction scheme requested does not exist. Choosing basic first-order scheme." << endl;
			pr = new BasicPR();
	}
	pr->setup(mesh, unkown, dp);
	/// Rhie-Chow constant [c](@ref c) is set as 0.5.
	c = 0.5;
	nvar = 3;
}

InviscidFlux::~InviscidFlux()
{
	if(isallocdp) {
		delete [] dp;
		delete pr;
	}
}

void InviscidFlux::compute_flux()
{
	pr->compute_pressure_difference();

	// account for solid walls.
	vector<int> wall(4);		// one flag for each of the 4 boundaries.
	for(int i = 0; i < 4; i++)
		if(bcflags[i] == 2 || bcflags[i] == 3)
			wall[i] = 0;
		else
			wall[i] = 1;

	// add inviscid flux contribution to residuals
	cout << "InviscidFlux: compute_flux(): Computing inviscid fluxes now..." << endl;
	// We consider the first ghost cell layer for i, but not the last.
	int i, j, k;
	double area, nx, ny, eigen, vdotn;
	double bhalf; vector<double> uhalf(3);			// interface values for each unknown
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
			bhalf = 0.5*(beta->get(i,j) + beta->get(i+1,j);
			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			// now get eigenvalue for Rhie-Chow
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4*bhalf*bhalf) );

			// Boundary faces need to be treated differently to account for wall BCs.
			if(i == 0)
			{
				g(i,0) = area*(rho*vdotn + c*eigen*dp[0](i,j)/(bhalf*bhalf))*wall[2];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[2] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[2] + uhalf[0]*ny); 
			}
			else if(i == m->gimx()-1)
			{
				g(i,0) = area*(rho*vdotn + c*eigen*dp[0](i,j)/(bhalf*bhalf))*wall[0];
				g(i,1) = area*(rho*vdotn*uhalf[1]*wall[0] + uhalf[0]*nx);
				g(i,2) = area*(rho*vdotn*uhalf[2]*wall[0] + uhalf[0]*ny); 
			}
			else
			{
				g(i,0) = area*(rho*vdotn + c*eigen*dp[0](i,j)/(bhalf*bhalf));
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
			bhalf = 0.5*(beta->get(i,j) + beta(i,j+1));
			vdotn = uhalf[1]*nx + uhalf[2]*ny;
			eigen = 0.5*( fabs(vdotn) + sqrt(vdotn*vdotn + 4*bhalf*bhalf) );
			
			if(j == 0)
			{
				h(j,0) = area*(rho*vdotn + c*eigen*dp[1](i,j)/(bhalf*bhalf))*wall[3];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[3] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[3] + uhalf[0]*ny); 
			}
			else if(j == m->gjmx()-1)
			{
				h(j,0) = area*(rho*vdotn + c*eigen*dp[1](i,j)/(bhalf*bhalf))*wall[1];
				h(j,1) = area*(rho*vdotn*uhalf[1]*wall[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2]*wall[1] + uhalf[0]*ny); 
			}
			else
			{
				h(j,0) = area*(rho*vdotn + c*eigen*dp[1](i,j)/(bhalf*bhalf));
				h(j,1) = area*(rho*vdotn*uhalf[1] + uhalf[0]*nx);
				h(j,2) = area*(rho*vdotn*uhalf[2] + uhalf[0]*ny); 
			}
		}
		for(j = 1; j < m->gjmx()-1; j++)
			for(k = 0; k < nvar; k++)
				res[k](i,j) += h(j,k) - h(j-1,k);
	}
}

}
