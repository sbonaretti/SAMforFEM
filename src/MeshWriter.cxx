/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshWriter.h>

namespace mesh{

	// constructor
	MeshWriter::MeshWriter(){

		_meshFlag = false;
		_mesh = vtkPolyData::New();
		_youngModulusElemFlag = false;
		_youngModulusNodeFlag = false;
		_youngModulus = vtkDoubleArray::New();
		_loadFlag = false;
		_loadMagnitude = vtkDoubleArray::New();
		_loadPoints = vtkPoints::New();
		_bcFlag = false;
		_bcType = vtkDoubleArray::New();
		_bcPoints = vtkPoints::New();

	}

	// destructor
	MeshWriter::~MeshWriter(){
	}

}