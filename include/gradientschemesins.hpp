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
/**	\note compute_fluxes() *increments* the residual res. If res contains anything non-zero, contribution by gradients is *added* to it; it is not replaced.
*/
class GradientSchemeIns
{
protected:
	int nvar;				///< Number of variables in u.
	int ndim;
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
	ndim = 2;
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

/// Parallel CV model for gradient reconstruction
class ParallelCVGradientIns : public GradientSchemeIns
{
	/// flux across any i-face
	vector<double> g;
	/// flux across any j-face
	vector<double> h;

	/// CV normal
	vector<double> deln;
	/// gradient
	vector<double> grad;
	
public:
	void setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc);
	
	/// This function is a dummy for this class
	void compute_lhs();

	void compute_fluxes();
};

void ParallelCVGradientIns::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc)
{
	cout << "ParallelCVGradientIns: setup() called." << endl;
	GradientSchemeIns::setup(mesh, unknown, residual, lhs, visc);
	GradientSchemeIns::compute_CV_volumes();

	g.resize(m->gimx());
	h.resize(m->gjmx());
	deln.resize(ndim);
	grad.resize(ndim);
}

void ParallelCVGradientIns::compute_lhs()
{ }

void ParallelCVGradientIns::compute_fluxes()
{
	int i,j,k,l;
	double f1, f2, f3, f4;
	
	for(k = 1; k < nvar; k++)
	{
		// i- fluxes
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			for(i = 0; i <= m->gimx()-1; i++)
			{
				// values at the 4 faces of the parallel CV
				f1 = u[k].get(i+1,j);
				f3 = u[k].get(i,j);
				f2 = 0.25*(u[k].get(i,j) + u[k].get(i+1,j) + u[k].get(i,j+1) + u[k].get(i+1,j+1));
				f4 = 0.25*(u[k].get(i,j) + u[k].get(i+1,j) + u[k].get(i,j-1) + u[k].get(i+1,j-1));
				
				// normal to the i-face of CV
				deln[0] = -(m->gyc(i+1,j) - m->gyc(i,j));
				deln[1] = (m->gxc(i+1,j) - m->gxc(i,j));

				for(l = 0; l < ndim; l++)
				{
					grad[l] = (f1-f3)*m->gdel(i,j,l) + (f2-f4)*deln[l];
					grad[l] /= dualvol[0].get(i,j);
				}

				g[i] = grad[0]*m->gdel(i,j,0) + grad[1]*m->gdel(i,j,1);
			}
			for(i = 1; i <= m->gimx()-1; i++)
				res[k](i,j) -= mu*(g[i] - g[i-1]);
		}

		// j- fluxes
		for(i = 1; i <= m->gimx()-1; i++)
		{
			for(j = 0; j <= m->gjmx()-1; j++)
			{
				// values at 4 faces of CV
				f1 = u[k].get(i,j+1);
				f3 = u[k].get(i,j);
				f4 = 0.25* (u[k].get(i,j) + u[k].get(i,j+1) + u[k].get(i-1,j) + u[k].get(i-1,j+1));
				f2 = 0.25* (u[k].get(i,j) + u[k].get(i,j+1) + u[k].get(i+1,j) + u[k].get(i+1,j+1));

				deln[0] = m->gyc(i,j+1) - m->gyc(i,j);
				deln[1] = -(m->gxc(i,j+1) - m->gxc(i,j));

				for(l = 0; l < ndim; l++)
				{
					grad[l] = (f1-f3)*m->gdel(i,j,2+l) + (f2-f4)*deln[l];
					grad[l] /= dualvol[1].get(i,j);
				}

				h[j] = grad[0]*m->gdel(i,j,2) + grad[1]*m->gdel(i,j,3);
			}
			for(j = 1; j <= m->gjmx()-1; j++)
				res[k](i,j) -= mu*(h[j] - h[j-1]);
		}
	}
}

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

void ThinLayerGradientIns::setup(Structmesh2d* mesh, Array2d<double>* unknown, Array2d<double>* residual, Array2d<double>* lhs, double visc)
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
	int i,j,k;
	vector<double> g(nvar), h(nvar);
	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
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
	Array2d<double>* cdelu;		// to store cell-wise gradient values for each of the 5 cells in a loop iteration.
	cdelu = new Array2d<double>[nvar];
	for(int i = 0; i < nvar; i++)
		cdelu[i].setup(5,2);
	int i,j,k,d;

	// first iterate over interior cells
	for(i = 2; i <= m->gimx()-2; i++)
		for(j = 2; j <= m->gjmx()-2; j++)
		{
			for(d = 0; d < 2; d++)
			{
				for(k = 1; k < nvar; k++)
				{
					//Array2d<double>* u = &(this->u[k]);
					// i,j-1
					cdelu[k](0,d) = 0.5/m->gvol(i,j-1)*( (u[k].get(i,j-1)+u[k].get(i+1,j-1))*m->gdel(i,j-1,d) + (u[k].get(i,j-1)+u[k].get(i,j))*m->gdel(i,j-1,2+d)
						- (u[k].get(i,j-1)+u[k].get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u[k].get(i,j-1)+u[k].get(i,j-2))*m->gdel(i,j-2,2+d));
					//i-1,j
					cdelu[k](1,d) = 0.5/m->gvol(i-1,j)*( (u[k].get(i-1,j)+u[k].get(i,j))*m->gdel(i-1,j,d) + (u[k].get(i-1,j)+u[k].get(i-1,j+1))*m->gdel(i-1,j,2+d)
						- (u[k].get(i-1,j)+u[k].get(i-2,j))*m->gdel(i-2,j,d) - (u[k].get(i-1,j)+u[k].get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
					//i,j
					cdelu[k](2,d) = 0.5/m->gvol(i,j)*( (u[k].get(i,j)+u[k].get(i+1,j))*m->gdel(i,j,d) + (u[k].get(i,j)+u[k].get(i,j+1))*m->gdel(i,j,2+d)
						- (u[k].get(i,j)+u[k].get(i-1,j))*m->gdel(i-1,j,d) - (u[k].get(i,j)+u[k].get(i,j-1))*m->gdel(i,j-1,2+d));
					//i+1,j
					cdelu[k](3,d) = 0.5/m->gvol(i+1,j)*( (u[k].get(i+1,j)+u[k].get(i+2,j))*m->gdel(i+1,j,d) + (u[k].get(i+1,j)+u[k].get(i+1,j+1))*m->gdel(i+1,j,2+d)
						- (u[k].get(i+1,j)+u[k].get(i,j))*m->gdel(i,j,d) - (u[k].get(i+1,j)+u[k].get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
					//i,j+1
					cdelu[k](4,d) = 0.5/m->gvol(i,j+1)*( (u[k].get(i,j+1)+u[k].get(i+1,j+1))*m->gdel(i,j+1,d) + (u[k].get(i,j+1)+u[k].get(i,j+2))*m->gdel(i,j+1,2+d)
						- (u[k].get(i,j+1)+u[k].get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u[k].get(i,j+1)+u[k].get(i,j))*m->gdel(i,j,2+d));
				}
			}
			
			// now calculate contributions from the 4 faces
			for(k = 1; k < nvar; k++)
			{
				Array2d<double>* u = &(this->u[k]);
				res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](3,0))*m->gdel(i,j,0)+(cdelu[k](2,1)+cdelu[k](3,1))*m->gdel(i,j,1) 
					- ((cdelu[k](2,0)+cdelu[k](3,0))*svect[0](i,j)+(cdelu[k](2,1)+cdelu[k](3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
					+ (u[k].get(i+1,j)-u[k].get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
				res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](1,0))*m->gdel(i-1,j,0)+(cdelu[k](2,1)+cdelu[k](1,1))*m->gdel(i-1,j,1) 
					- ((cdelu[k](2,0)+cdelu[k](1,0))*svect[0](i-1,j)+(cdelu[k](2,1)+cdelu[k](1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
					+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
				res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](0,0))*m->gdel(i,j-1,2)+(cdelu[k](2,1)+cdelu[k](0,1))*m->gdel(i,j-1,3) 
					- ((cdelu[k](2,0)+cdelu[k](0,0))*svect[2](i,j-1)+(cdelu[k](2,1)+cdelu[k](0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
					+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
				res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](4,0))*m->gdel(i,j,2)+(cdelu[k](2,1)+cdelu[k](4,1))*m->gdel(i,j,3) 
					- ((cdelu[k](2,0)+cdelu[k](4,0))*svect[2](i,j)+(cdelu[k](2,1)+cdelu[k](4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
					+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
			}
		}
	
	/** For boundary cells:
	 * For the "second layer" of ghost cells, we simply copy flow velocities from the first layer of ghost cells.
	 * We take care of corner boundary cells during the treatment of boundaries 1 and 3 (the i=const boundaries).
	 */
	
	// boundary 1
	vector<double> arvec(2);		// del for second ghost cells
	vector<double> crvec(2);

	i = m->gimx()-1;
	crvec[0] = -(m->gy(i+1,0) - m->gy(i,0));
	crvec[1] = m->gx(i+1,0) - m->gx(i,0);
	for(j = 1; j <= m->gjmx()-1; j++)
	{
		
		for(d = 0; d < 2; d++)
		{
			for(k = 1; k < nvar; k++)
			{
				Array2d<double>* u =&(this->u[k]);
				// i,j-1
				cdelu[k](0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
					- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+( j>1 ? u->get(i,j-2):u->get(i,j-1)))*( j>1 ? m->gdel(i,j-2,2+d):m->gdel(i,j-1,2+d)));
				//i-1,j
				cdelu[k](1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
					- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
				//i,j
				cdelu[k](2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
					- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
				//i+1,j
				cdelu[k](3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+1,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
					- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
				//i,j+1
				cdelu[k](4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+(j < m->gjmx()-1 ? u->get(i,j+2):u->get(i,j+1)))*m->gdel(i,j+1,2+d)
					- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
			}
		}
		
		// now calculate contributions from the 4 faces
		for(k = 1; k < nvar; k++)
		{
			Array2d<double>* u =&(this->u[k]);
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](3,0))*m->gdel(i,j,0)+(cdelu[k](2,1)+cdelu[k](3,1))*m->gdel(i,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](3,0))*svect[0](i,j)+(cdelu[k](2,1)+cdelu[k](3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
				+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](1,0))*m->gdel(i-1,j,0)+(cdelu[k](2,1)+cdelu[k](1,1))*m->gdel(i-1,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](1,0))*svect[0](i-1,j)+(cdelu[k](2,1)+cdelu[k](1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
				+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](0,0))*m->gdel(i,j-1,2)+(cdelu[k](2,1)+cdelu[k](0,1))*m->gdel(i,j-1,3) 
				- ((cdelu[k](2,0)+cdelu[k](0,0))*svect[2](i,j-1)+(cdelu[k](2,1)+cdelu[k](0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
				+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](4,0))*m->gdel(i,j,2)+(cdelu[k](2,1)+cdelu[k](4,1))*m->gdel(i,j,3) 
				- ((cdelu[k](2,0)+cdelu[k](4,0))*svect[2](i,j)+(cdelu[k](2,1)+cdelu[k](4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
				+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
		}
	}

	// boundary 2
	j = m->gjmx()-1;
	for(i = 2; i <= m->gimx()-2; i++)
	{
		for(k = 1; k < nvar; k++)
		{
			Array2d<double>* u =&(this->u[k]);
			for(d = 0; d < 2; d++)
			{
				// i,j-1
				cdelu[k](0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
					- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+u->get(i,j-2))*m->gdel(i,j-2,2+d));
				//i-1,j
				cdelu[k](1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
					- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
				//i,j
				cdelu[k](2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
					- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
				//i+1,j
				cdelu[k](3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
					- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
				//i,j+1
				cdelu[k](4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+u->get(i,j+1)*m->gdel(i,j+1,2+d)
					- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d)));
			}
			
			// now calculate contributions from the 4 faces
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](3,0))*m->gdel(i,j,0)+(cdelu[k](2,1)+cdelu[k](3,1))*m->gdel(i,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](3,0))*svect[0](i,j)+(cdelu[k](2,1)+cdelu[k](3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
				+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](1,0))*m->gdel(i-1,j,0)+(cdelu[k](2,1)+cdelu[k](1,1))*m->gdel(i-1,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](1,0))*svect[0](i-1,j)+(cdelu[k](2,1)+cdelu[k](1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
				+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](0,0))*m->gdel(i,j-1,2)+(cdelu[k](2,1)+cdelu[k](0,1))*m->gdel(i,j-1,3) 
				- ((cdelu[k](2,0)+cdelu[k](0,0))*svect[2](i,j-1)+(cdelu[k](2,1)+cdelu[k](0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
				+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](4,0))*m->gdel(i,j,2)+(cdelu[k](2,1)+cdelu[k](4,1))*m->gdel(i,j,3) 
				- ((cdelu[k](2,0)+cdelu[k](4,0))*svect[2](i,j)+(cdelu[k](2,1)+cdelu[k](4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
				+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
		}
	}

	// boundary 3
	i = 1;
	crvec[0] = -(m->gy(i+1,0) - m->gy(i,0));
	crvec[1] = m->gx(i+1,0) - m->gx(i,0);
	for(j = 1; j <= m->gjmx()-1; j++)
	{
		arvec[0] = m->gy(i-1,j+1) - m->gy(i-1,j);
		arvec[1] = -(m->gx(i-1,j+1) - m->gx(i-1,j));
		for(k = 1; k < nvar; k++)
		{
			Array2d<double>* u =&(this->u[k]);
			for(d = 0; d < 2; d++)
			{
				// i,j-1
				cdelu[k](0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
					- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+(j>1 ? u->get(i,j-2):u->get(i,j-1)))*(j>1 ? m->gdel(i,j-2,2+d):m->gdel(i,j-1,2+d)));
				//i-1,j
				cdelu[k](1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
					- (u->get(i-1,j)+u->get(i-1,j))*arvec[d] - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
				//i,j
				cdelu[k](2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
					- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
				//i+1,j
				cdelu[k](3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
					- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
				//i,j+1
				cdelu[k](4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+(j<m->gjmx()-1 ? u->get(i,j+2):u->get(i,j+1)))*m->gdel(i,j+1,2+d)
					- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
			}
			
			// now calculate contributions from the 4 faces
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](3,0))*m->gdel(i,j,0)+(cdelu[k](2,1)+cdelu[k](3,1))*m->gdel(i,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](3,0))*svect[0](i,j)+(cdelu[k](2,1)+cdelu[k](3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
				+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](1,0))*m->gdel(i-1,j,0)+(cdelu[k](2,1)+cdelu[k](1,1))*m->gdel(i-1,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](1,0))*svect[0](i-1,j)+(cdelu[k](2,1)+cdelu[k](1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
				+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](0,0))*m->gdel(i,j-1,2)+(cdelu[k](2,1)+cdelu[k](0,1))*m->gdel(i,j-1,3) 
				- ((cdelu[k](2,0)+cdelu[k](0,0))*svect[2](i,j-1)+(cdelu[k](2,1)+cdelu[k](0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
				+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](4,0))*m->gdel(i,j,2)+(cdelu[k](2,1)+cdelu[k](4,1))*m->gdel(i,j,3) 
				- ((cdelu[k](2,0)+cdelu[k](4,0))*svect[2](i,j)+(cdelu[k](2,1)+cdelu[k](4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
				+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
		}
	}

	// boundary 4
	j = 1;
	for(i = 2; i <= m->gjmx()-2; i++)
	{
		arvec[0] = -(m->gy(i+1,j-1) - m->gy(i,j-1));
		arvec[1] = m->gx(i+1,j-1) - m->gx(i,j-1);
		for(k = 1; k < nvar; k++)
		{
			Array2d<double>* u =&(this->u[k]);
			for(d = 0; d < 2; d++)
			{
				// i,j-1
				cdelu[k](0,d) = 0.5/m->gvol(i,j-1)*( (u->get(i,j-1)+u->get(i+1,j-1))*m->gdel(i,j-1,d) + (u->get(i,j-1)+u->get(i,j))*m->gdel(i,j-1,2+d)
					- (u->get(i,j-1)+u->get(i-1,j-1))*m->gdel(i-1,j-1,d) - (u->get(i,j-1)+u->get(i,j-1))*arvec[d]);
				//i-1,j
				cdelu[k](1,d) = 0.5/m->gvol(i-1,j)*( (u->get(i-1,j)+u->get(i,j))*m->gdel(i-1,j,d) + (u->get(i-1,j)+u->get(i-1,j+1))*m->gdel(i-1,j,2+d)
					- (u->get(i-1,j)+u->get(i-2,j))*m->gdel(i-2,j,d) - (u->get(i-1,j)+u->get(i-1,j-1))*m->gdel(i-1,j-1,2+d));
				//i,j
				cdelu[k](2,d) = 0.5/m->gvol(i,j)*( (u->get(i,j)+u->get(i+1,j))*m->gdel(i,j,d) + (u->get(i,j)+u->get(i,j+1))*m->gdel(i,j,2+d)
					- (u->get(i,j)+u->get(i-1,j))*m->gdel(i-1,j,d) - (u->get(i,j)+u->get(i,j-1))*m->gdel(i,j-1,2+d));
				//i+1,j
				cdelu[k](3,d) = 0.5/m->gvol(i+1,j)*( (u->get(i+1,j)+u->get(i+2,j))*m->gdel(i+1,j,d) + (u->get(i+1,j)+u->get(i+1,j+1))*m->gdel(i+1,j,2+d)
					- (u->get(i+1,j)+u->get(i,j))*m->gdel(i,j,d) - (u->get(i+1,j)+u->get(i+1,j-1))*m->gdel(i+1,j-1,2+d));
				//i,j+1
				cdelu[k](4,d) = 0.5/m->gvol(i,j+1)*( (u->get(i,j+1)+u->get(i+1,j+1))*m->gdel(i,j+1,d) + (u->get(i,j+1)+u->get(i,j+2))*m->gdel(i,j+1,2+d)
					- (u->get(i,j+1)+u->get(i-1,j+1))*m->gdel(i-1,j+1,d) - (u->get(i,j+1)+u->get(i,j))*m->gdel(i,j,2+d));
			}
			
			// now calculate contributions from the 4 faces
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](3,0))*m->gdel(i,j,0)+(cdelu[k](2,1)+cdelu[k](3,1))*m->gdel(i,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](3,0))*svect[0](i,j)+(cdelu[k](2,1)+cdelu[k](3,1))*svect[1](i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1)))
				+ (u->get(i+1,j)-u->get(i,j))*(svect[0](i,j)*m->gdel(i,j,0)+svect[1](i,j)*m->gdel(i,j,1))/dels[0](i,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](1,0))*m->gdel(i-1,j,0)+(cdelu[k](2,1)+cdelu[k](1,1))*m->gdel(i-1,j,1) 
				- ((cdelu[k](2,0)+cdelu[k](1,0))*svect[0](i-1,j)+(cdelu[k](2,1)+cdelu[k](1,1))*svect[1](i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1)))
				+ (u->get(i,j)-u->get(i-1,j))*(svect[0](i-1,j)*m->gdel(i-1,j,0)+svect[1](i-1,j)*m->gdel(i-1,j,1))/dels[0](i-1,j);
			res[k](i,j) -= 0.5*((cdelu[k](2,0)+cdelu[k](0,0))*m->gdel(i,j-1,2)+(cdelu[k](2,1)+cdelu[k](0,1))*m->gdel(i,j-1,3) 
				- ((cdelu[k](2,0)+cdelu[k](0,0))*svect[2](i,j-1)+(cdelu[k](2,1)+cdelu[k](0,1))*svect[3](i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3)))
				+ (u->get(i,j)-u->get(i,j-1))*(svect[2](i,j-1)*m->gdel(i,j-1,2)+svect[3](i,j-1)*m->gdel(i,j-1,3))/dels[1](i,j-1);
			res[k](i,j) += 0.5*((cdelu[k](2,0)+cdelu[k](4,0))*m->gdel(i,j,2)+(cdelu[k](2,1)+cdelu[k](4,1))*m->gdel(i,j,3) 
				- ((cdelu[k](2,0)+cdelu[k](4,0))*svect[2](i,j)+(cdelu[k](2,1)+cdelu[k](4,1))*svect[3](i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3)))
				+ (u->get(i,j+1)-u->get(i,j))*(svect[2](i,j)*m->gdel(i,j,2)+svect[3](i,j)*m->gdel(i,j,3))/dels[1](i,j);
		}
	}
	delete [] cdelu;

	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
			for(k = 1; k < nvar; k++)
				res[k](i,j) *= mu;
}

}




