/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTATIONVALIDATIONVOLUMEMESH_H
#define REGISTATIONVALIDATIONVOLUMEMESH_H

#include <RegistrationValidationVolume.h>

using namespace image;


namespace registration{

	/**
	* It computes the validation of the Ansys Morpher registration for the complete volume.
	*
	* For the reference bone, each node is labelled as part of the cortical bone, trabecular bone or marrow using the segmentation mask (classified image).
	* For each moving mesh, the nodes are compared with the reference one: if they belong to the same part of the bone,
	* they are labelled as "right", otherwise as "wrong". Average and standard deviation in % are given for the right positioned nodes
	* Average, standard deviation and standard error are given for the wrong positioned nodes.
	*/

	class RegistrationValidationVolumeMesh: public RegistrationValidationVolume{

	public:

		// constructor/destructor
		RegistrationValidationVolumeMesh();
		~RegistrationValidationVolumeMesh();

		// member functions
		/**
		* Computes the volume morphing validation 
		*/
		void Update();

	private:
		/**
		* Defines to which part of the bone the voxels belong (reference bone)
		*/
		void ReferenceLabels(vtkPolyData* mesh, ImageHandler::ImageType::Pointer image, 
			vtkIntArray* &cortical, vtkIntArray* &trabecular, vtkIntArray* &marrow);
		/**
		* Evaluates if correspondence among nodes is the same for the moving as for the reference 
		*/
		void checkMovings(vtkIntArray* index,
			vtkPolyData* mesh, ImageHandler::ImageType::Pointer image, int part, double &distance, int &nOfRight);
	
	};
}
#endif //REGISTATIONVALIDATIONVOLUMEMESH_H