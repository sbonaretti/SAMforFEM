/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointReader.h>

namespace points{

	// constructor
	PointReader::PointReader(){

		_points = vtkPoints::New();
	}

	// destructor
	PointReader::~PointReader(){
	}

}