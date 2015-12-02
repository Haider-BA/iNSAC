#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

#ifndef __STRUCTMESH2D_H
#include "structmesh2d.hpp"
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#define __GRADIENTSCHEMESINS_H 1

using namespace amat;
using namespace acfd;
using namespace std;

namespace acfd {

/**	\brief Base class for gradient reconstruction schemes for AC formulation of iNS equations. */
/**	NOTE: compute_fluxes() *increments* the residual res. If res contains anything non-zero, contribution by gradients is *added* to it; it is not replaced.
*/
class GradientSchemeIns
{
protected:
	int nvar;				///< Number of variables in u.
	double mu;				///< Viscosity (dynamic)
	Structmesh2d* m;		///< Associated mesh
	Array2d<double>* u;		///< The unknown with which to compute the graient flux
	Array2d<double>* res;	///< The residual containing fluxes at each grid cell
	
	/** \brief LHS of AF scheme corresponding to a particular (compact) gradient model.
	
		- a(i,j,0) corresponds to coeff of u(i-1,j)
		- a(i,j,1) corresponds to coeff of u(i,j-1)
		- a(i,j,2) corresponds to coeff of u(i,j)
		- a(i,j,3) corresponds to coeff of u(i,j+1)
		- a(i,j,4) corresponds to coeff of u(i+1,j),
	 in the row of Ax=b corresponding to the (i,j) cell.
	*/
	Array2d<double>* a;

	/// \brief Volumes of the thin-layer CVs for the gradient.
	///
	/// Two components: first component for +i face and the other for +j face. Arranged like Structmesh2d::del.
	vector<Array2d<double>> dualvol;

public:
	virtual void setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc);
	void compute_CV_volumes();		///< Computes [volumes](@ref dualvol) of thin-layer CVs. To be precomputed just once.

	virtual void compute_s();		///< Required for Normal tangent gradient decomposition scheme.
	
	virtual void compute_fluxes() = 0;
	
	/** \brief Computes LHS arrays corresponding to a particular gradient scheme.
	
		Make sure to execcute [calculate_CV_volumes](@ref calculate_CV_volumes) before calling this function.
	*/
	virtual void compute_lhs() = 0;
};

void GradientSchemeIns::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc)
{
	nvar = 3;
	mu = visc;
	m = mesh;
	u = unknown;
	res = residual;
	a = lhs;
	dualvol.resize(2);
	for(int i = 0; i < 2; i++)
		dualvol[i].setup(m->gimx()+1, m->gjmx()+1);
}

void GradientSchemeIns::compute_CV_volumes()
{
	// iterate over cells
	for(int i = 0; i <= m->gimx()-1; i++)
		for(int j = 0; j <= m->gjmx()-1; j++)
			if((i>0 || j>0) && (i < m->gimx() || j < m->gjmx()))
			{
				dualvol[0](i,j) = (m->gvol(i+1,j) + m->gvol(i,j))*0.5;
				dualvol[1](i,j) = (m->gvol(i,j+1) + m->gvol(i,j))*0.5;
			}
	// We can do the above since volumes of all cells, real and ghost, have been computed in Structmesh2d::preprocess.
}

void GradientSchemeIns::compute_s()
{ }

//----------------- end of base class GradientScheme -----------------------------//


/**	\brief Thin layer gradient reconstruction scheme. 
*
* We consider grad u at a face to be influenced by only the change in normal component of u at the face.
*/
class ThinLayerGradientIns : public GradientSchemeIns
{
public:
	/// Initializes data and computes LHS
	void setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc);
	/// Adds viscous flux contribution to the residual
	void compute_fluxes();
	void compute_lhs();
};

void ThinLayerGradientIns::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs)
{
	GradientSchemeIns::setup(mesh, unknown, residual, lhs, visc);
	GradientSchemeIns::compute_CV_volumes();
	compute_lhs();
}

void ThinLayerGradientIns::compute_lhs()
{
	for(int i = 1; i <= m->gimx()-1; i++)
		for(int j = 1; j <= m->gjmx()-1; j++)
		{
			a[0](i,j) = (m->gdel(i-1,j,0)*m->gdel(i-1,j,0) + m->gdel(i-1,j,1)*m->gdel(i-1,j,1))/dualvol[0](i-1,j);
			a[1](i,j) = (m->gdel(i,j-1,2)*m->gdel(i,j-1,2) + m->gdel(i,j-1,3)*m->gdel(i,j-1,3))/dualvol[1](i,j-1);
			a[3](i,j) = (m->gdel(i,j,2)*m->gdel(i,j,2) + m->gdel(i,j,3)*m->gdel(i,j,3))/dualvol[1](i,j);
			a[4](i,j) = (m->gdel(i,j,0)*m->gdel(i,j,0) + m->gdel(i,j,1)*m->gdel(i,j,1))/dualvol[0](i,j);
			a[2](i,j) = -1.0 * (a[0](i,j) + a[1](i,j) + a[3](i,j) + a[4](i,j));
		}
}

/** \note We *subtract* the viscous contribution from the residual as the viscous term has a negative sign in iNS equations. */
void ThinLayerGradientIns::compute_fluxes()
{
	vector<double> g(nvar), h(nvar);
	for(int i = 1; i <= m->gimx()-1; i++)
		for(int j = 1; j <= m->gjmx()-1; j++)
		{
			for(k = 1; k < nvar; k++)
			{
				g[k] = (u[k].get(i+1,j) - u[k].get(i,j))*a[4](i,j) - (u[k].get(i,j)-u[k].get(i-1,j))*a[0](i,j);
				h[k] = (u[k].get(i,j+1) - u[k].get(i,j))*a[3](i,j) - (u[k].get(i,j)-u[k].get(i,j-1))*a[1](i,j);
				res[k](i,j) -= mu*(g[k]+h[k]);
			}
		}
}

/**	\brief Computes LHS and residual for FVM solution of Poisson equation using Normal-tangent gradient decomposition model.
*
*	Applied properly, the scheme is non-compact and has 13 terms in the LHS for each cell. However, we implement only the 5 compact terms for the LHS. 
*	But we treat the residual fully. This requires specification of more than one layer of ghost cells. This is taken care of without introducing another layer;
*	rather, we compute a `ghost state' on-the-fly using the BCs.
* \todo TODO : This class is not yet ready for iNS:
* - Implement for a vector unknown
* - Multiply by viscosity at appropriate places
* - Remove boundary-condition treatment, and do something about the extended stencil for the real boundary cells.
*/
class NormTanGradientIns : public GradientSchemeIns
{
	/** \brief Stores a unit vector in the direction of the line joining the two cells across each face.
	*
	*  Organized like Structmesh2d::del.
	*/
	vector<Array2d<double>> svect;

	/** \brief Magnitude of the corresponding vectors in [svect](@ref svect).
	*
	* dels[0] refers to magnitude of the vector for +i face and dels[1] refers to magnitude of the vector for +j face.
	*/
	vector<Array2d<double>> dels;

public:
	void compute_s();
	void compute_lhs();
	void compute_fluxes();
};

void NormTanGradientIns::compute_s()
{
	svect.resize(4);
	for(int i = 0; i < 4; i++)
		svect[i].setup(m->gimx()+1, m->gjmx()+1);
	
	dels.resize(2);
	for(int i = 0; i < 2; i++)
		dels[i].setup(m->gimx()+1, m->gjmx()+1);

	/* NOTE: We are not calculating svect and dels for the ghost cells i=imx+1 and j=jmx+1. */
	
	for(int i = 0; i <= m->gimx()-1; i++)
		for(int j = 0; j <= m->gjmx()-1; j++)
		{
			svect[0](i,j) = m->gxc(i+1,j) - m->gxc(i,j);
			svect[1](i,j) = m->gyc(i+1,j) - m->gyc(i,j);
			dels[0](i,j) = sqrt(svect[0](i,j)*svect[0](i,j) + svect[1](i,j)*svect[1](i,j));
			svect[0](i,j) /= dels[0](i,j);
			svect[1](i,j) /= dels[0](i,j);
			
			svect[2](i,j) = m->gxc(i,j+1) - m->gxc(i,j);
			svect[3](i,j) = m->gyc(i,j+1) - m->gyc(i,j);
			dels[1](i,j) = sqrt(svect[2](i,j)*svect[2](i,j) + svect[3](i,j)*svect[3](i,j));
			svect[2](i,j) /= dels[1](i,j);
			svect[3](i,j) /= dels[1](i,j);
		}
}

/** As an approximation, we use the same coefficients as for the Thin Layer gradient model. */
void NormTanGradientIns::compute_lhs()
{
	vector<double> temp(2);
	int k;
	double coeff;
	
	for(int i = 1; i <= m->gimx()-1; i++)
		for(int j = 1; j <= m->gjmx()-1; j++)
		{
			a[0](i,j) = (m->gdel(i-1,j,0)*m->gdel(i-1,j,0) + m->gdel(i-1,j,1)*m->gdel(i-1,j,1))/dualvol[0](i-1,j);
			a[1](i,j) = (m->gdel(i,j-1,2)*m->gdel(i,j-1,2) + m->gdel(i,j-1,3)*m->gdel(i,j-1,3))/dualvol[1](i,j-1);
			a[3](i,j) = (m->gdel(i,j,2)*m->gdel(i,j,2) + m->gdel(i,j,3)*m->gdel(i,j,3))/dualvol[1](i,j);
			a[4](i,j) = (m->gdel(i,j,0)*m->gdel(i,j,0) + m->gdel(i,j,1)*m->gdel(i,j,1))/dualvol[0](i,j);
			a[2](i,j) = -1.0 * (a[0](i,j) + a[1](i,j) + a[3](i,j) + a[4](i,j));
		}
}

void NormTanGradientIns::compute_fluxes()
{
	Array2d<double> cdelu(5,2);		// to store cell-wise gradient values for each of the 5 cells in a loop iteration.
	int k,d;

	// first iterate over interior cells
	for(int i = 2; i <= m->gimx()-2; i++)
		for(int j = 2; j <= m->gjmx()-2; j++)
		{
			for(d = 0; d < 2; d++)
			{
				// i,j-1
				cdelu(0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
					- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+u->get(i,j-2))*m->gdel(i,j-2,2+d));
				//i-1,j
				cdelu(1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
					- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
				//i,j
				cdelu(2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
					- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
				//i+1,j
				cdelu(3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
					- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
				//i,j+1
				cdelu(4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+u->get(i,j+2))*m->gdel(i,j+1,2+d)
					- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
			}
			
			// now calculate contributions from the 4 faces
			(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(3,0))*m->gdel(i,j,0)+(cdelu(2,1)+cdelu(3,1))*m->gdel(i,j,1) 
				- ((cdelu(2,0)+cdelu(3,0))*svect[0](i,j)+(cdelu(2,1)+cdelu(3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
				+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
			(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(1,0))*m->gdel(i-1,j,0)+(cdelu(2,1)+cdelu(1,1))*m->gdel(i-1,j,1) 
				- ((cdelu(2,0)+cdelu(1,0))*svect[0](i-1,j)+(cdelu(2,1)+cdelu(1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
				+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
			(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(0,0))*m->gdel(i,j-1,2)+(cdelu(2,1)+cdelu(0,1))*m->gdel(i,j-1,3) 
				- ((cdelu(2,0)+cdelu(0,0))*svect[2](i,j-1)+(cdelu(2,1)+cdelu(0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
				+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
			(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(4,0))*m->gdel(i,j,2)+(cdelu(2,1)+cdelu(4,1))*m->gdel(i,j,3) 
				- ((cdelu(2,0)+cdelu(4,0))*svect[2](i,j)+(cdelu(2,1)+cdelu(4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
				+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
		}
	
	// boundary cells
	// If it's a homogeneous Neumann boundary, the "second" ghost cell value is the same as the ghost cell value;
	// If it's a Dirichlet boundary, the value is 2*bvalue - u(second interior cell).
	// We take care of corner boundary cells during the treatment of boundaries 1 and 3 (the i=const boundaries).
	// TODO: Implement this second ghost cell value for non-homogeneous Neumann condition.
	
	// boundary 1
	double secondghost;
	double sg1, sg2, sg3, sg4;		// for corner boundary cells
	vector<double> arvec(2);		// del for second ghost cells
	vector<double> crvec(2);		// extra del required in the 4 corner boundary cells

	if(bcflag[3] >= 1)	// Neumann
	{
		sg1 = u->get(m->gimx()-1,0);
		sg3 = u->get(1,0);
	}
	else {
		sg1 = 2*bvalue[3] - u->get(m->gimx()-1,2);
		sg3 = 2*bvalue[3] - u->get(1,2);
	}
	
	if(bcflag[1] >= 1)	// Neumann
	{
		sg2 = u->get(m->gimx()-1,m->gjmx());
		sg4 = u->get(1,m->gjmx());
	}
	else {
		sg2 = 2*bvalue[1] - u->get(m->gimx()-1,m->gjmx()-2);
		sg4 = 2*bvalue[1] - u->get(1,m->gjmx()-2);
	}

	int i = m->gimx()-1;
	crvec[0] = -(m->gy(i+1,0) - m->gy(i,0));
	crvec[1] = m->gx(i+1,0) - m->gx(i,0);
	for(int j = 1; j <= m->gjmx()-1; j++)
	{
		if(bcflag[0] >= 1)
			secondghost = u->get(i+1,j);
		else
			secondghost = 2*bvalue[0] - u->get(i-1,j);
		
		for(d = 0; d < 2; d++)
		{
			// i,j-1
			cdelu(0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
				- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+( j>1 ? u->get(i,j-2):sg1))*( j>1 ? m->gdel(i,j-2,2+d):crvec[d]));
			//i-1,j
			cdelu(1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
				- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
			//i,j
			cdelu(2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
				- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
			//i+1,j
			cdelu(3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+secondghost)*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
				- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
			//i,j+1
			cdelu(4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+( j < m->gjmx()-1 ? u->get(i,j+2):sg2))*m->gdel(i,j+1,2+d)
				- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
		}
		
		// now calculate contributions from the 4 faces
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(3,0))*m->gdel(i,j,0)+(cdelu(2,1)+cdelu(3,1))*m->gdel(i,j,1) 
			- ((cdelu(2,0)+cdelu(3,0))*svect[0](i,j)+(cdelu(2,1)+cdelu(3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
			+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(1,0))*m->gdel(i-1,j,0)+(cdelu(2,1)+cdelu(1,1))*m->gdel(i-1,j,1) 
			- ((cdelu(2,0)+cdelu(1,0))*svect[0](i-1,j)+(cdelu(2,1)+cdelu(1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
			+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(0,0))*m->gdel(i,j-1,2)+(cdelu(2,1)+cdelu(0,1))*m->gdel(i,j-1,3) 
			- ((cdelu(2,0)+cdelu(0,0))*svect[2](i,j-1)+(cdelu(2,1)+cdelu(0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
			+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(4,0))*m->gdel(i,j,2)+(cdelu(2,1)+cdelu(4,1))*m->gdel(i,j,3) 
			- ((cdelu(2,0)+cdelu(4,0))*svect[2](i,j)+(cdelu(2,1)+cdelu(4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
			+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
	}

	// boundary 2
	int j = m->gjmx()-1;
	for(int i = 2; i <= m->gimx()-2; i++)
	{
		if(bcflag[1] >= 1)
			secondghost = u->get(i,j+1);
		else
			secondghost = 2*bvalue[1] - u->get(i,j-1);
		for(d = 0; d < 2; d++)
		{
			// i,j-1
			cdelu(0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
				- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+u->get(i,j-2))*m->gdel(i,j-2,2+d));
			//i-1,j
			cdelu(1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
				- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
			//i,j
			cdelu(2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
				- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
			//i+1,j
			cdelu(3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
				- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
			//i,j+1
			cdelu(4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+secondghost)*m->gdel(i,j+1,2+d)
				- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
		}
		
		// now calculate contributions from the 4 faces
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(3,0))*m->gdel(i,j,0)+(cdelu(2,1)+cdelu(3,1))*m->gdel(i,j,1) 
			- ((cdelu(2,0)+cdelu(3,0))*svect[0](i,j)+(cdelu(2,1)+cdelu(3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
			+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(1,0))*m->gdel(i-1,j,0)+(cdelu(2,1)+cdelu(1,1))*m->gdel(i-1,j,1) 
			- ((cdelu(2,0)+cdelu(1,0))*svect[0](i-1,j)+(cdelu(2,1)+cdelu(1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
			+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(0,0))*m->gdel(i,j-1,2)+(cdelu(2,1)+cdelu(0,1))*m->gdel(i,j-1,3) 
			- ((cdelu(2,0)+cdelu(0,0))*svect[2](i,j-1)+(cdelu(2,1)+cdelu(0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
			+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(4,0))*m->gdel(i,j,2)+(cdelu(2,1)+cdelu(4,1))*m->gdel(i,j,3) 
			- ((cdelu(2,0)+cdelu(4,0))*svect[2](i,j)+(cdelu(2,1)+cdelu(4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
			+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
	}

	// boundary 3
	i = 1;
	crvec[0] = -(m->gy(i+1,0) - m->gy(i,0));
	crvec[1] = m->gx(i+1,0) - m->gx(i,0);
	for(int j = 1; j <= m->gjmx()-1; j++)
	{
		if(bcflag[2] >= 1)
			secondghost = u->get(i-1,j);
		else
			secondghost = 2*bvalue[2] - u->get(i+1,j);
		arvec[0] = m->gy(i-1,j+1) - m->gy(i-1,j);
		arvec[1] = -(m->gx(i-1,j+1) - m->gx(i-1,j));
		for(d = 0; d < 2; d++)
		{
			// i,j-1
			cdelu(0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
				- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+(j>1 ? u->get(i,j-2):sg3))*(j>1 ? m->gdel(i,j-2,2+d):crvec[d]));
			//i-1,j
			cdelu(1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
				- (u->get(i-1,j)+secondghost)*arvec[d] - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
			//i,j
			cdelu(2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
				- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
			//i+1,j
			cdelu(3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
				- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
			//i,j+1
			cdelu(4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+(j<m->gjmx()-1 ? u->get(i,j+2):sg4))*m->gdel(i,j+1,2+d)
				- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
		}
		
		// now calculate contributions from the 4 faces
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(3,0))*m->gdel(i,j,0)+(cdelu(2,1)+cdelu(3,1))*m->gdel(i,j,1) 
			- ((cdelu(2,0)+cdelu(3,0))*svect[0](i,j)+(cdelu(2,1)+cdelu(3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
			+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(1,0))*m->gdel(i-1,j,0)+(cdelu(2,1)+cdelu(1,1))*m->gdel(i-1,j,1) 
			- ((cdelu(2,0)+cdelu(1,0))*svect[0](i-1,j)+(cdelu(2,1)+cdelu(1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
			+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(0,0))*m->gdel(i,j-1,2)+(cdelu(2,1)+cdelu(0,1))*m->gdel(i,j-1,3) 
			- ((cdelu(2,0)+cdelu(0,0))*svect[2](i,j-1)+(cdelu(2,1)+cdelu(0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
			+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(4,0))*m->gdel(i,j,2)+(cdelu(2,1)+cdelu(4,1))*m->gdel(i,j,3) 
			- ((cdelu(2,0)+cdelu(4,0))*svect[2](i,j)+(cdelu(2,1)+cdelu(4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
			+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
	}

	// boundary 4
	j = 1;
	for(int i = 2; i <= m->gjmx()-2; i++)
	{
		if(bcflag[3] >= 1)
			secondghost = u->get(i,j-1);
		else
			secondghost = 2*bvalue[3] - u->get(i,j+1);
		arvec[0] = -(m->gy(i+1,j-1) - m->gy(i,j-1));
		arvec[1] = m->gx(i+1,j-1) - m->gx(i,j-1);
		for(d = 0; d < 2; d++)
		{
			// i,j-1
			cdelu(0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
				- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+secondghost)*arvec[d]);
			//i-1,j
			cdelu(1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
				- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
			//i,j
			cdelu(2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
				- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
			//i+1,j
			cdelu(3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
				- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
			//i,j+1
			cdelu(4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+u->get(i,j+2))*m->gdel(i,j+1,2+d)
				- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
		}
		
		// now calculate contributions from the 4 faces
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(3,0))*m->gdel(i,j,0)+(cdelu(2,1)+cdelu(3,1))*m->gdel(i,j,1) 
			- ((cdelu(2,0)+cdelu(3,0))*svect[0](i,j)+(cdelu(2,1)+cdelu(3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
			+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(1,0))*m->gdel(i-1,j,0)+(cdelu(2,1)+cdelu(1,1))*m->gdel(i-1,j,1) 
			- ((cdelu(2,0)+cdelu(1,0))*svect[0](i-1,j)+(cdelu(2,1)+cdelu(1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
			+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
		(*res)(i,j) -= 0.5*((cdelu(2,0)+cdelu(0,0))*m->gdel(i,j-1,2)+(cdelu(2,1)+cdelu(0,1))*m->gdel(i,j-1,3) 
			- ((cdelu(2,0)+cdelu(0,0))*svect[2](i,j-1)+(cdelu(2,1)+cdelu(0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
			+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
		(*res)(i,j) += 0.5*((cdelu(2,0)+cdelu(4,0))*m->gdel(i,j,2)+(cdelu(2,1)+cdelu(4,1))*m->gdel(i,j,3) 
			- ((cdelu(2,0)+cdelu(4,0))*svect[2](i,j)+(cdelu(2,1)+cdelu(4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
			+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
	}
}

}




