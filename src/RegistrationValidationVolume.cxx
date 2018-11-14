/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationValidationVolume.h>

#include <ImageHandler.h>
#include <MeshReaderMorpherVolume.h>
#include <PointWriterXyz.h>

#include <itkImageRegionIterator.h>

#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkSTLWriter.h>


using namespace image;
using namespace mesh;
using namespace points;

namespace registration{

	// constructor
	RegistrationValidationVolume::RegistrationValidationVolume(){
	}

	// destructor
	RegistrationValidationVolume::~RegistrationValidationVolume(){
	}

}