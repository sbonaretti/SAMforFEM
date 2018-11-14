/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTATIONVALIDATIONVOLUMEIMAGE_H
#define REGISTATIONVALIDATIONVOLUMEIMAGE_H

#include <RegistrationValidationVolume.h>

#include <ImageHandler.h>

using namespace image;

namespace registration{

	/**
	* It calculates the DICE coefficient for the bone three layers masks.
	* Use SetReferenceMaskFileName() and  SetMovingMaskFileNames() before Update()
	*/

	class RegistrationValidationVolumeImage: public RegistrationValidationVolume{

	public:

		// constructor/destructor
		RegistrationValidationVolumeImage();
		~RegistrationValidationVolumeImage();

		// member functions
		/**
		* Computes the volume morphing validation 
		*/
		void Update();

	};
}
#endif //REGISTATIONVALIDATIONVOLUMEIMAGE_H
