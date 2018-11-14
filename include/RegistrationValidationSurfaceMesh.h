/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTRATIONVALIDATIONSURFACEMESH_H
#define REGISTRATIONVALIDATIONSURFACEMESH_H

#include <RegistrationValidationSurface.h>

namespace registration{


	/**
	* It computes the validation of the Ansys Morpher registration for the outer surface.
	*
	* It compares the original .stl file with the surface extracted from the volume .cdb file
	* (to do that it recalls the class MeshExtractOuterSurface). The comparison is done in two directions:
	* .cdb nodes -> .stl surface and .stl nodes -> .cdb surface. The max distance is considered in both cases
	* and the max bextween the two max is then considered as measure of the quality of the registration. (Hausdorff distance)
	* The average distance is also taken into account in the two cases.
	* 
	* Use SetFlag() to define the kind of mesh. ()
	*/

	class RegistrationValidationSurfaceMesh: public RegistrationValidationSurface{
		
		public:
		
		// constructor/destructor
		RegistrationValidationSurfaceMesh();
		~RegistrationValidationSurfaceMesh();

		// accessors
		/**
		* Sets the volume mesh
		*/
		void SetFlag (int flag) {_flag = flag;}

		// member functions
		/**
		* Computes the surface morphing validation 
		*/
		void Update();

	protected:
		int _flag;

	};
}
#endif //REGISTRATIONVALIDATIONSURFACEMESH_H