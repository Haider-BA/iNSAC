#include <structmesh2d.hpp>
#include <aoutput_struct.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

int main()
{
	Structmesh2d m;
	m.readmesh("../../grid4.grid");

	Array2d<double>* dum;
	Array2d<double>** dum2;
	string* sdum;
	string title = "project2_mesh";
	Structdata2d d(&m, 0, dum, sdum, 0, dum2, sdum, title);
	d.writevtk("grid4.vtk");

	cout << endl;
	return 0;
}
