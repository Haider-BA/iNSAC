#ifndef __INVISCIDFLUX_H
#include <inviscidflux.hpp>
#endif

#ifndef __GRADIENTSCHEMESINS_H
#include <gradientschemesins.hpp>
#endif

using namespace amat;
using namespace acfd;
using namespace std;

/** \brief Implements the main solution loop of stead-state iNS equations in artificial compressibility form.
*/
class Steady_insac
{
	Structmesh2d* m;				///< Pointer to mesh.
	int nvar;						///< Number of unkowns - 3 in this case.
	Array2d<double>* u;				///< Values of the 3 unkowns - index 0 is pressure, 1 is x-velocity and 2 is y-velocity.
	Array2d<double>* res;			///< Residuals for each of the 3 quanitities.
	Array2d<double> beta;			///< Artificial compressibility in each cell.
	Array2d<double> dt;				///< Local time step for each cell.
	double cfl;						///< CFL number
	vector<int> bcflags;			///< Flags indicating the type of boundary for each of the 4 boundaries.
	vector<vector<double>> bvalues;	///< Boundary values for each of the 4 boundaries (2 each, at most)
	double rho;						///< Density
	double mu;						///< Viscosity

	InviscidFlux* invf;				///< Inviscid flux object
	GradientSchemeIns* grad;		///< Gradient reconstruction object for viscous flux
	Array2d<double>* visc_lhs;		///< Required for gradient reconstruction.

	string pressurescheme;			///< Denotes which pressure reconstruction scheme to use in mass flux.
	double uref;					///< Some reference fluid velocity
	int ndim;						///< No. of spatial dimensions (2)
	double tol;						///< Relative residual tolerance to decide convergence to steady-state
	int maxiter;					///< Max number of time steps

	bool isalloc;

	Array2d<double> upoint[3];
	
	Array2d<double>* vel[2];			///< Required for output of velocity

public:
	Steady_insac();

	/// Sets up the iNS problem.
	/** \param _bcflags contains 4 integers denoting the type of boundary for each boundary.
	 - 0 is a velocity inlet
	 - 1 is a pressure outlet
	 - 2 is a no-slip wall
	 - 3 is a slip wall
	\param bvalues contains the corresponding boundary values - at most two per boundary.
	\param gradscheme is string which is either "thinlayer" or "normtan", describing the gradient reconstruction scheme to use for viscous fluxes.
	\param pressure_scheme is a string (either "basic" or "tvd") describing the pressure reconstruction to use for the Rhie-Chow mass flux.
	\param refvel is some reference fluid velocity.
	\param CFL is the C.F.L. number to use.
	\param tolerance is the relative tolerance
	\param maxiters is the maximum number of time steps
	*/
	void setup(Structmesh2d* mesh, double dens, double visc, vector<int> _bcflags, vector<vector<double>> _bvalues, string gradscheme, string pressure_scheme, double refvel, double CFL, double tolerance, int maxiters);

	~Steady_insac();

	/// Computes artificial conmpressibility ([beta](@ref beta) ) for each cell.
	void compute_beta();
	
	/// Sets quantities to ghost cells.
	void setBCs();
	
	void setInitialConditions();
	
	void compute_timesteps();
	
	/// Contains the main solver time loop.
	void solve();
	
	/// Computes solutions at grid nodes by averaging cell values around each node
	void getPointSolution();

	Array2d<double>* getpressure();
	Array2d<double>** getvelocity();
};

Steady_insac::Steady_insac() {
	isalloc = false;
}

void Steady_insac::setup(Structmesh2d* mesh, double dens, double visc, vector<int> _bcflags, vector<vector<double>> _bvalues, string gradscheme, string pressure_scheme, double refvel, double CFL, double tolerance, int maxiters)
{
	ndim = 2;
	m = mesh;
	nvar = 3;
	rho = dens;
	mu = visc;
	bcflags = _bcflags;
	bvalues = _bvalues;
	pressurescheme = pressure_scheme;
	uref = refvel;
	if(gradscheme == "normtan")
		grad = new NormTanGradientIns;
	else if(gradscheme == "parallelcv")
		grad = new ParallelCVGradientIns;
	else
		grad = new ThinLayerGradientIns;
	
	invf = new InviscidFlux;

	u = new Array2d<double>[nvar];
	res = new Array2d<double>[nvar];
	visc_lhs = new Array2d<double>[5];

	isalloc = true;

	for(int i = 0; i < nvar; i++)
	{
		u[i].setup(m->gimx()+1, m->gjmx()+1);
		res[i].setup(m->gimx()+1, m->gjmx()+1);
	}
	for(int i = 0; i < 5; i++)
		visc_lhs[i].setup(m->gimx()+1, m->gjmx()+1);
	beta.setup(m->gimx()+1,m->gjmx()+1);
	dt.setup(m->gimx()+1,m->gjmx()+1);
	cfl = CFL;
	tol = tolerance;
	maxiter = maxiters;

	invf->setup(m, u, res, &beta, rho, pressure_scheme, _bcflags, _bvalues);
	grad->setup(m,u,res,visc_lhs, mu);
	
	for(int k = 0; k < 3; k++)
		upoint[k].setup(m->gimx()+1,m->gjmx()+1);
}

Steady_insac::~Steady_insac()
{
	if(isalloc) {
		delete [] u;
		delete [] res;
		delete [] visc_lhs;
		delete grad;
		delete invf;
	}
}

/** We currently do not consider viscosity in calculating the artificial compressibility [beta](@ref beta).
*/
void Steady_insac::compute_beta()
{
	int i,j; double vmag, uvisc, h;
	for(i = 0; i <= m->gimx(); i++)
		for(j = 0; j <= m->gjmx(); j++)
		{
			vmag = sqrt(u[1].get(i,j)*u[1].get(i,j) + u[2].get(i,j)*u[2].get(i,j));
			if(vmag >= uref)
				beta(i,j) = vmag;
			else
				beta(i,j) = uref;
			
			// account for viscous effect
			/*h = sqrt((m->gx(i+1,j+1)-m->gx(i,j+1))*(m->gx(i+1,j+1)-m->gx(i,j+1)) + (m->gy(i+1,j+1)-m->gy(i,j+1))*(m->gy(i+1,j+1)-m->gy(i,j+1)));
			h += sqrt((m->gx(i+1,j+1)-m->gx(i+1,j))*(m->gx(i+1,j+1)-m->gx(i+1,j)) + (m->gy(i+1,j+1)-m->gy(i+1,j))*(m->gy(i+1,j+1)-m->gy(i+1,j)));
			h /= 2.0;
			uvisc = (mu/rho)/h;
			if(beta(i,j) < uvisc) beta(i,j) = uvisc;*/
		}
}

/** Note that for an inlet boundary, a parabolic profile is imposed with maximum velocity as that given by the [bvalues](@ref bvalues) entries. We assume that the respective boundaries are parallel to x- or y-axis, so the inlets work only for straight boundaries parallel to the axes.
* \note For the time being, inlet is only allowed for boundary 2, ie, for the i=1 boundary.
* \note The corner ghost cells are not directly given values; they are just set as the average of their neighboring ghost cells.
*/
void Steady_insac::setBCs()
{
	int i,j,k;
	double vdotn;

	// boundary 0
	i = m->gimx();
	if(bcflags[0] == 1)		// pressure outlet
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			u[0](i,j) = bvalues[0][0];
			u[1](i,j) = u[1].get(i-1,j);
			u[2](i,j) = u[2].get(i-1,j);
		}
	else if(bcflags[0] == 2)	// no-slip wall
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			u[0](i,j) = u[0](i-1,j);
			u[1](i,j) = 2.0*bvalues[0][0] - u[1](i-1,j);
			u[2](i,j) = 2.0*bvalues[0][1] - u[2](i-1,j);
		}
	else 			// slip-wall
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			u[0](i,j) = u[0](i-1,j);
			vdotn = u[1](i-1,j)*m->gdel(i-1,j,0) + u[2](i-1,j)*m->gdel(i-1,j,1);
			u[1](i,j) = u[1](i-1,j) - 2*vdotn*m->gdel(i-1,j,0);
			u[2](i,j) = u[2](i-1,j) - 2*vdotn*m->gdel(i-1,j,1);
		}

	// boundary 2
	i = 0;
	if(bcflags[2] == 0)			// parabolic velocity inlet
	{
		double a,b,c, rm;
		rm = (m->gy(1,1) + m->gy(1,m->gjmx()))/2.0;		// mid point
		a = bvalues[2][0] / ( rm*rm - 2.0*rm*rm + m->gy(1,1)*m->gy(1,m->gjmx()) );
		b = -a*2.0*rm;
		c = -a*m->gy(1,1)*m->gy(1,1) - b*m->gy(1,1);
		for(j = 1; j <= m->gjmx()-1; j++) {
			u[1](i,j) = a*m->gyc(i,j)*m->gyc(i,j) + b*m->gyc(i,j) + c;
			u[2](i,j) = 0;
			u[0](i,j) = u[0](i+1,j);
		}
	}
	else if(bcflags[2] == 1)		// pressure outlet
	{
		for(j = 1; j <= m->gjmx()-1; j++) 
		{
			u[0](i,j) = bvalues[2][0];
			u[1](i,j) = u[1](i+1,j);
			u[2](i,j) = u[2](i+1,j);
		}
	}
	else if(bcflags[2] == 2)	// no-slip wall
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			u[0](i,j) = u[0](i+1,j);
			u[1](i,j) = 2.0*bvalues[2][0] - u[1](i+1,j);
			u[2](i,j) = 2.0*bvalues[2][1] - u[2](i+1,j);
		}
	else 			// slip-wall
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			u[0](i,j) = u[0](i+1,j);
			// we use negative of the del value as the normal always points in the positive i direction.
			vdotn = u[1](i+1,j)*(-1*m->gdel(i,j,0)) + u[2](i+1,j)*(-1*m->gdel(i,j,1));
			u[1](i,j) = u[1](i+1,j) - 2*vdotn*m->gdel(i+1,j,0);
			u[2](i,j) = u[2](i+1,j) - 2*vdotn*m->gdel(i+1,j,1);
		}

	// boundary 1
	j = m->gjmx();
	if(bcflags[1] == 1)		// pressure outlet
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = bvalues[1][0];
			u[1](i,j) = u[1](i,j-1);
			u[2](i,j) = u[2](i,j-1);
		}
	else if(bcflags[1] == 2)	// no-slip wall
	{
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = u[0](i,j-1);
			u[1](i,j) = 2*bvalues[1][0] - u[1](i,j-1);
			u[2](i,j) = 2*bvalues[1][1] - u[2](i,j-1);
		}
	}
	else if(bcflags[1] == 3)	// slip wall
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = u[0](i,j-1);
			vdotn = u[1](i,j-1)*m->gdel(i,j-1,2) + u[2](i,j-1)*m->gdel(i,j-1,3);
			u[1](i,j) = u[1](i,j-1) - 2.0*vdotn*m->gdel(i,j-1,2);
			u[2](i,j) = u[2](i,j-1) - 2.0*vdotn*m->gdel(i,j-1,3);
		}

	// boundary 3
	j = 0;
	if(bcflags[3] == 1)		// pressure outlet
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = bvalues[3][0];
			u[1](i,j) = u[1](i,j+1);
			u[2](i,j) = u[2](i,j+1);
		}
	else if(bcflags[3] == 2)	// no-slip wall
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = u[0](i,j+1);
			u[1](i,j) = 2*bvalues[3][0] - u[1](i,j+1);
			u[2](i,j) = 2*bvalues[3][1] - u[2](i,j+1);
		}
	else if(bcflags[3] == 3)	// slip wall
		for(i = 1; i <= m->gimx()-1; i++)
		{
			u[0](i,j) = u[0](i,j+1);
			// we use negative of the del values as the normal always points in the positive j-direction.
			vdotn = u[1](i,j+1)*(-1)*m->gdel(i,j,2) + u[2](i,j+1)*(-1)*m->gdel(i,j,3);
			u[1](i,j) = u[1](i,j+1) - 2.0*vdotn*m->gdel(i,j,2);
			u[2](i,j) = u[2](i,j+1) - 2.0*vdotn*m->gdel(i,j,3);
		}

	for(k = 0; k < nvar; k++)
	{
		u[k](0,0) = 0.5*(u[k](1,0)+u[k](0,1));
		u[k](0,m->gjmx()) = 0.5*(u[k](1,m->gjmx())+u[k](0,m->gjmx()-1));
		u[k](m->gimx(),0) = 0.5*(u[k](m->gimx(),1)+u[k](m->gimx()-1,0));
		u[k](m->gimx(),m->gjmx()) = 0.5*(u[k](m->gimx()-1,m->gjmx())+u[k](m->gimx(),m->gjmx()-1));
	}
}

/** Computes local time-step for each cell using characteristics of the system in the normal directions.
* Face normals 'in a cell' are calculated by averaging those of the two faces bounding the cell in each the i- and j-directions.
*/
void Steady_insac::compute_timesteps()
{
	int i,j, dim;
	vector<double> areavi(ndim), areavj(ndim);		// area vectors
	vector<double> inormal(ndim), jnormal(ndim);	// unit normal vectors
	double areai, areaj;							// area magnitudes
	double vdotni, vdotnj, eigeni, eigenj, voldt;
	//cout << "Steady_insac: compute_timesteps(): Now computing time steps for next iteration..." << endl;

	for(i = 1; i <= m->gimx()-1; i++)
		for(j = 1; j <= m->gjmx()-1; j++)
		{
			for(dim = 0; dim < ndim; dim++) {
				areavi[dim] = 0.5*(m->gdel(i,j,dim) + m->gdel(i-1,j,dim));
				areavj[dim] = 0.5*(m->gdel(i,j,2+dim) + m->gdel(i,j-1,2+dim));
			}
			areai = sqrt(areavi[0]*areavi[0] + areavi[1]*areavi[1]);
			areaj = sqrt(areavj[0]*areavj[0] + areavj[1]*areavj[1]);
			for(dim = 0; dim < ndim; dim++)
			{
				inormal[dim] = areavi[dim]/areai;
				jnormal[dim] = areavj[dim]/areaj;
			}
			// we now have face geometrical info "at the cell centers"

			vdotni = u[1](i,j)*inormal[0] * u[2](i,j)*inormal[1];
			vdotnj = u[1](i,j)*jnormal[0] * u[2](i,j)*jnormal[1];

			eigeni = 0.5*( fabs(vdotni) + sqrt(vdotni*vdotni + 4.0*beta(i,j)*beta(i,j)) );
			eigenj = 0.5*( fabs(vdotnj) + sqrt(vdotnj*vdotnj + 4.0*beta(i,j)*beta(i,j)) );
			voldt = eigeni*areai + eigenj*areaj;
			dt(i,j) = m->gvol(i,j)/voldt * cfl;
		}
}

/** Sets all quantities in all cells to zero. */
void Steady_insac::setInitialConditions()
{
	int i,j,k;
	for(i = 0; i <= m->gimx(); i++)
		for(j = 0; j <= m->gjmx(); j++)
		{
			for(k = 0; k < nvar; k++)
			{
				u[k](i,j) = 0;
				res[k](i,j) = 0;
			}
		}
}

/** Make sure the [Steady_insac](@ref Steady_insac) object has been [setup](@ref setup) and initialized with some [initial condition](@ref setInitialConditions). Both the momentum-magnitude residual and the mass flux are taken as convergence criteria; the tolerance for the mass flux is the square-root of the tolerance for the relative momentum-magnitude residual. [tol](@ref tol) is the latter.
*/
void Steady_insac::solve()
{
	setInitialConditions();
	compute_beta();

	int i,j,k;
	double resnorm, resnorm0, massflux, dtv;
	cout << "Steady_insac: solve(): Beginning the time-march." << endl;
	for(int n = 0; n < maxiter; n++)
	{
		// calculate stuff needed for this iteration
		setBCs();
		
		compute_timesteps();

		// compute fluxes
		for(k = 0; k < nvar; k++)
			res[k].zeros();
		invf->compute_fluxes();
		grad->compute_fluxes();

		// check convergence
		resnorm = 0; massflux = 0;
		for(i = 1; i <= m->gimx()-1; i++)
			for(j = 1; j <= m->gjmx()-1; j++)
			{
				resnorm += (res[1].get(i,j)*res[1].get(i,j) + res[2].get(i,j)*res[2].get(i,j))*m->gvol(i,j);
				massflux += res[0].get(i,j);
			}
		resnorm = sqrt(resnorm);
		if(n == 0) resnorm0 = resnorm;
		if(n == 1 || n%10 == 0) {
			cout << "Steady_insac: solve(): Iteration " << n << ": relative L2 norm of residual = " << resnorm/resnorm0 << ", net mass flux = " << massflux << endl;
			cout << "  L2 norm of residual = " << resnorm << endl;
		}
		if(resnorm/resnorm0 < tol && fabs(massflux) < sqrt(tol))
		{
			cout << "Steady_insac: solve(): Converged in " << n << " iterations. Norm of final residual = " << resnorm << ", final net mass flux = " << massflux << endl;
			break;
		}

		// update u
		for(i = 1; i <= m->gimx()-1; i++)
			for(j = 1; j <= m->gjmx()-1; j++)
			{
				dtv = dt.get(i,j)/m->gvol(i,j);
				u[0](i,j) = u[0].get(i,j) - res[0].get(i,j)*beta.get(i,j)*beta.get(i,j)*dtv;
				for(k = 1; k < nvar; k++)
					u[k](i,j) = u[k].get(i,j) - res[k].get(i,j)*dtv/rho;
			}
	}
	
	getPointSolution();
}

Array2d<double>* Steady_insac::getpressure()
{
	return &(upoint[0]);
}

Array2d<double>** Steady_insac::getvelocity()
{
	/// [vel](@ref vel) is an array of 2 [Array2d<double>](@ref Array2d)'s, one each to hold a pointer to a component Array2d of velocity.
	vel[0] = &upoint[1];
	vel[1] = &upoint[2];
	return vel;
}

void Steady_insac::getPointSolution()
{ 
	// interior points
	for(int k = 0; k < 3; k++)
		for(int i = 1; i <= m->gimx(); i++)
			for(int j = 1; j <= m->gjmx(); j++)
			{
				upoint[k](i,j) = (u[k].get(i,j) + u[k].get(i-1,j) + u[k].get(i,j-1) + u[k].get(i-1,j-1))*0.25;
			}
}
