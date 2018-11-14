/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHMORPHERVALIDATIONSURFACE_H
#define MESHMORPHERVALIDATIONSURFACE_H

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <QStringList>


namespace mesh{


	/**
	* It computes the validation of the Ansys Morpher registration for the outer surface.
	* 
	* It compares the original .stl file with the surface extracted from the volume .cdb file
	* (to do that it recalls the class MeshExtractOuterSurface). The comparison is done in two directions:
	* .cdb nodes -> .stl surface and .stl nodes -> .cdb surface. The max distance is considered in both cases
	* and the max bextween the two max is then considered as measure of the quality of the registration. (Hausdorff distance)
	* The average distance is also taken into account in the two cases.
	* 
	* ! The first five letters of the files must be correspondent for the two dataset. 
	* Ex. file name of the .stl mesh: 00077-ca-f-024-162-057_femur-l-Segmented-SLICES_affine.stl
	*     file name of the volume mesh: 00077VolumeResult.cdb
	*/


	class MeshMorpherValidationSurface{
		
	public:
		
		// constructor/destructor
		MeshMorpherValidationSurface();
		~MeshMorpherValidationSurface();

		// accessors
		/**
		* Sets the file names of the .stl meshes
		*/
		void SetStlFileNames (QStringList stlFileNames) {_stlFileNames = stlFileNames;}
		/**
		* Sets the file names of the .cdb meshes
		*/
		void SetVolumeMeshFileNames (QStringList volumeMeshFileNames) {_volumeMeshFileNames = volumeMeshFileNames;}
		/**
		* Gets the vector of all the Hausdorff distances
		*/
		vtkDoubleArray* GetDistanceVector () {return _HausdorffDistances;}
		/**
		* Gets the average of all the Hausdorff distances
		*/
		double GetDistanceAverage () {return _average;}
		/**
		* Gets the standard deviation of the Hausdorff distances
		*/
		double GetDistanceStandardDeviation () {return _standardDeviation;}
		/**
		* Gets the standard error of the Hausdorff distances
		*/
		double GetDistanceStandardError () {return _standardError;}
		
		// member functions
		/**
		* Computes the surface morphing validation 
		*/
		void Update();


	
	protected:

		// functions
		void findPointToSurfaceDistances(vtkPoints* points, vtkPolyData* surface, double & maxDistance);
	
		// data members from accessors
		QStringList _stlFileNames;
		QStringList _volumeMeshFileNames;
		vtkDoubleArray* _HausdorffDistances;
		double _average;
		double _standardDeviation; 
		double _standardError;

	};
}
#endif //MESHMORPHERVALIDATIONSURFACE_H