/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTRATIONVALIDATIONSURFACE_H
#define REGISTRATIONVALIDATIONSURFACE_H

#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <QStringList>


namespace registration{


	/**
	*
	* Abstract class to evaluate the surface registration
	* 
	* It loads the originals and the registereds and compares their surfaces.
	* ! The first five letters of the files must be correspondent for the two dataset. 
	* Ex. file name of the .stl mesh: 00077-ca-f-024-162-057_femur-l-Segmented-SLICES_affine.stl
	*     file name of the volume mesh: 00077VolumeResult.cdb
	*/


	class RegistrationValidationSurface{
		
	public:
		
		// constructor/destructor
		RegistrationValidationSurface();
		~RegistrationValidationSurface();

		// accessors
		/**
		* Sets the file names of the originals
		*/
		void SetOriginalFileNames (QStringList originalFileNames) {_originalFileNames = originalFileNames;}
		/**
		* Sets the file names of the registereds
		*/
		void SetRegisteredFileNames (QStringList registeredFileNames) {_registeredFileNames = registeredFileNames;}
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
		* Computes the registration surface validation 
		*/
		virtual void Update() = 0;


	
	protected:

		// functions
		void findPointToSurfaceDistances(vtkPoints* points, vtkPolyData* surface, double & maxDistance);
	
		// data members from accessors
		QStringList _originalFileNames;
		QStringList _registeredFileNames;
		vtkDoubleArray* _HausdorffDistances;
		double _average;
		double _standardDeviation; 
		double _standardError;

	};
}
#endif //REGISTRATIONVALIDATIONSURFACE_H