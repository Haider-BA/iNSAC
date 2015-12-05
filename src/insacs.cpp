/** @file insacs.cpp
 * @brief Driver function for stead-state iNS AC solution
 */

#include "ins.hpp"

using namespace amat;
using namespace acfd;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please provide the name of a control file. Quitting.\n";
		return -1;
	}

	string meshfile, outfile, pressurescheme, gradscheme, dum;
	double tol, rho, mu, cfl, refvel;
	int maxiter;
	vector<int> bcflags(4);
	vector<vector<double>> bvalues(4);
	for(int i = 0; i < 4; i++)
		bvalues[i].resize(2);

	ifstream conf(argv[1]);

	// read control file
	conf >> dum; conf >> meshfile;
	conf >> dum; conf >> outfile;
	conf >> dum; conf >> rho;
	conf >> dum; conf >> mu;
	conf >> dum; conf >> refvel;
	conf >> dum; conf >> gradscheme;
	conf >> dum; conf >> pressurescheme;
	conf >> dum; conf >> cfl;
	conf >> dum; conf >> tol;
	conf >> dum; conf >> maxiter;
	conf >> dum;
	for(int i = 0; i < 4; i++)
		conf >> bcflags[i];
	conf >> dum;
	for(int i = 0; i < 4; i++)
	{
		// for a np-slip wall, we can take x- and y-velocities
		// for a velocity inlet, we just take the max velocity of a parabolic profile

		if(bcflags[i] == 2)
			conf >> bvalues[i][0] >> bvalues[i][1];
		else
			conf >> bvalues[i][0];
	}

	conf.close();

	Structmesh2d m;
	m.readmesh(meshfile);
	m.preprocess();

	Steady_insac ins;
	ins.setup(&m, rho, mu, bcflags, bvalues, gradscheme, pressurescheme, refvel, cfl, tol, maxiter);
	ins.solve();

	cout << endl;
	return 0;
}
