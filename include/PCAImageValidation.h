/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCAIMAGEVALIDATION_H
#define PCAIMAGEVALIDATION_H

#include <QStringList>
#include <vtkDoubleArray.h>
#include <itkImage.h>


namespace pca{

	/**
	* Abstract class for the validation of the PCA on images
	*/

	class PCAImageValidation{

	public:
		
		// constructor/destructor
		PCAImageValidation();
		~PCAImageValidation();
		
		typedef signed short VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;

		// accessors
		/**
		* Sets the file names of the training dataset
		*/
		void SetDataSetFileNames (QStringList fileNames) {_datasetFileNames = fileNames;}
		/**
		* Get the vector of distances between the training dataset and the evaluation set
		*/
		vtkDoubleArray* GetDistanceVector () {return _distanceVector;}
		/**
		* Get the average distance between the training dataset and the evaluation set
		*/
		double GetDistanceAverage () {return _distanceAverage;}
		/**
		* Get the standard deviation of the distances between the training dataset and the evaluation set
		*/
		double GetDistanceStandardDeviation () {return _distanceStandardDeviation;}
		/**
		* Get the standard error of the distances between the training dataset and the evaluation set
		*/
		double GetDistanceStandardError () {return _distanceStandardError;}

		// pure virtual function
		/**
		* Computes the calculation
		*/
		virtual void Update() = 0;

	protected:

		// data members from accessors
		QStringList _datasetFileNames;
		vtkDoubleArray* _distanceVector;
		double _distanceAverage;
		double _distanceStandardDeviation;
		double _distanceStandardError;

	};

}
#endif // PCAIMAGEVALIDATION_H