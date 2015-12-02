/** Creates a rectilinear grid from data passed as command-line arguments.
  Requires 5 arguments:
  - Number of points in i-direction
  - Number of points in j-direction
  - Length of domain in i-direction
  - Length of domain in j-direction
  - Output mesh file name.
*/ 

#include <structmesh2d.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cout << "Need 5 arguments: No. of points in i-dir, number of points in j-dir, x-length, y-length  and output file name." << endl;
		return -1;
	}
	int imx = stoi(argv[1]);
	int jmx = stoi(argv[2]);
	double xlength = stod(argv[3]);
	double ylength = stod(argv[4]);

	ofstream fout(argv[5]);
	fout << setprecision(20);

	fout << imx << " " << jmx << '\n';
	for(int j = 0; j < jmx; j++)
		for(int i = 0; i < imx; i++)
			fout << xlength/(imx-1)*i << " " << ylength/(jmx-1)*j << '\n';
	
	fout.close();
	cout << endl;
	return 0;
}
