/** @file aoutput_struct.hpp
 * @brief Provides a class to manage output of mesh data to VTK-type files.
 * @author Aditya Kashi
*/

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <string>
#endif

#ifndef __AARRAY2D_H
#include "aarray2d.hpp"
#endif

#ifndef __STRUCTMESH2D_H
#include "structmesh2d.hpp"
#endif

#define __AOUTPUT_STRUCT_H 1

using namespace amat;
using namespace acfd;
using namespace std;

namespace acfd {

///	Class for managing output of analysis data for simulations on structured 2D grids.
/**	We assume 1-based arrays for all array-quantities.
* Vectors are taken as an array (for each different vector quantity) of arrays (for each component of vector) of matrices over i,j.
*/
class Structdata2d
{
	int ndim;			///< Dimension of the problem
	int nscalars;		///< Number of scalar quantities
	int nvectors;		///< Number of vector quantities
	Structmesh2d* m;
	Array2d<double>* scalars;
	Array2d<double>** vectors;
	string* scalarnames;
	string* vectornames;
	string title;
public:
	Structdata2d(Structmesh2d* mesh, int n_scalars, Array2d<double>* _scalars, string* scalar_names, int n_vectors, Array2d<double>** _vectors, string* vector_names, string title);

	/// Writes data to a file in legacy VTK format
	void writevtk(string fname);
};

Structdata2d::Structdata2d(Structmesh2d* mesh, int n_scalars, Array2d<double>* _scalars, string* scalar_names, int n_vectors, Array2d<double>** _vectors, string* vector_names, string title)
{
	ndim = 2;
	m = mesh;
	nscalars = n_scalars;
	nvectors = n_vectors;
	scalars = _scalars;
	vectors = _vectors;
	scalarnames = scalar_names;
	vectornames = vector_names;
}

void Structdata2d::writevtk(string fname)
{
	cout << "Structdata2d: writevtk(): Writing data to file " << fname << endl;
	ofstream fout(fname);
	fout << "# vtk DataFile Version 2.0\n";
	fout << title << '\n';
	cout << title << '\n';
	fout << "ASCII\n";
	fout << "DATASET STRUCTURED_GRID\n";
	fout << "DIMENSIONS " << m->gimx() << " " << m->gjmx() << " 1\n";
	fout << "POINTS " << m->gimx()*m->gjmx() << " float\n";
	for(int j = 1; j <= m->gjmx(); j++)
		for(int i = 1; i <= m->gimx(); i++)
			fout << m->gx(i,j) << " " << m->gy(i,j) << " " << 0.0 << '\n';
	// Now output data
	if(nscalars > 0 || nvectors > 0)
	{
		fout << "POINT_DATA " << m->gimx()*m->gjmx() << '\n';
		for(int isca = 0; isca < nscalars; isca++)
		{
			fout << "SCALARS " << scalarnames[isca] << " float 1\n";
			fout << "LOOKUP_TABLE default\n";
			for(int j = 1; j <= m->gjmx(); j++)
				for(int i = 1; i <= m->gimx(); i++)
					fout << scalars[isca].get(i,j) << '\n';
		}
		for(int iv = 0; iv < nvectors; iv++)
		{
			fout << "VECTORS " << vectornames[iv] << " float\n";
			for(int j = 1; j <= m->gjmx(); j++)
				for(int i = 1; i <= m->gimx(); i++)
				{
					for(int idim = 0; idim < ndim; idim++)
						fout << vectors[iv][idim].get(i,j) << " ";
					for(int idim = 3-ndim; idim > 0; idim--)
						fout << "0 ";
					fout << '\n';
				}
		}
	}
}

}
