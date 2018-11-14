/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTRATIONVALIDATIONVOLUMEMESH_H
#define REGISTRATIONVALIDATIONVOLUMEMESH_H

#include <ImageHandler.h>

#include <vtkIntArray.h>
#include <vtkPolyData.h>

#include <QString>
#include <QStringList>

using namespace image;

namespace mesh{

	/**
	* It computes the validation of the Ansys Morpher registration for the complete volume.
	*
	* For the reference bone, each node is label as part of the cortical bone, trabecular bone or marrow using the segmentation mask (classified image).
	* For each moving mesh, the nodes are compared with the reference one: if they belong to the same part of the bone,
	* they are labelled as "right", otherwise as "wrong". Average and standard deviation in % are given for the right positioned nodes
	* Average, standard deviation and standard error are given for the wrong positioned nodes.
	*/

	class RegistrationValidationVolumeMesh{

	public:

		// constructor/destructor
		RegistrationValidationVolumeMesh();
		~RegistrationValidationVolumeMesh();

		// accessors
		/**
		* Sets the reference mesh file name
		*/
		void SetReferenceMeshFileName (QString referenceMeshFileName) {_referenceMeshFileName = referenceMeshFileName;}
		/**
		* Sets the reference mask file name
		*/
		void SetReferenceMaskFileName (QString referenceMaskFileName) {_referenceMaskFileName = referenceMaskFileName;}
		/**
		* Sets the moving mesh file names
		*/
		void SetMovingMeshFileNames (QStringList movingMeshFileNames) {_movingMeshFileNames = movingMeshFileNames;}
		/**
		* Sets the moving mask file names
		*/
		void SetMovingMaskFileNames (QStringList movingMaskFileNames) {_movingMaskFileNames = movingMaskFileNames;}


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
	protected:

		// data members from accessors
		QString _referenceMeshFileName;
		QString _referenceMaskFileName;
		QStringList _movingMeshFileNames;
		QStringList _movingMaskFileNames;
	};
}
#endif //REGISTRATIONVALIDATIONVOLUMEMESH_H