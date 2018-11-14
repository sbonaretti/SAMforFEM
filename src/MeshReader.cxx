/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshReader.h>

namespace mesh{

	// constructor
	MeshReader::MeshReader(){

		_mesh = vtkPolyData::New();
	}

	// destructor
	MeshReader::~MeshReader(){
		
		_mesh->Delete();
	
	}

}