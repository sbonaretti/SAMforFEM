/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTRATIONVALIDATIONSURFACEIMAGE_H
#define REGISTRATIONVALIDATIONSURFACEIMAGE_H

#include <RegistrationValidationSurface.h>

namespace registration{

	
	/**
	* It computes the validation of the image registration for the outer surface.
	*
	* It extracts the surface mesh from both the original (rigidly registered) and 
	* the registered images with the same threshold and evaluates the distances between them.
	* Meshes are simplyfy and smoothed in order to have about the same number of nodes 
	* as in the registrationValidatonSurfaceMesh case, where the input meshes have less points than the marching cube ones.
	* The comparison is done in two directions: original nodes -> registered surface and 
	* registered nodes -> original surface. The max distance is considered in both cases
	* and the max bextween the two max is then considered as measure of the quality of the registration. (Hausdorff distance)
	* The average distance is also taken into account in the two cases.
	*/

	class RegistrationValidationSurfaceImage: public RegistrationValidationSurface{

		public:
		
		// constructor/destructor
		RegistrationValidationSurfaceImage();
		~RegistrationValidationSurfaceImage();

		// member functions
		/**
		* Computes the surface morphing validation 
		*/
		void Update();

	};
}
#endif //REGISTRATIONVALIDATIONSURFACEIMAGE_H

