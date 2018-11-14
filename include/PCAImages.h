/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCAIMAGES_H
#define PCAIMAGES_H

#include <PCA.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <QStringList>


namespace pca{

	/**
	* It calculates the Statistical models on images, validates and create new instances.
	*
	* The Statistical Shape Model is computed on velocity Fields using itkFieldPCAShapeModelEstimator.h. 
	* Call SetImageFileNames() before ShapePCA(). (velocity fields are considered as float)
	*
	* The Statistical Intensity Model is computed on images using itkImagePCAShapeModelEstimator.h.
	* Call SetImageFileNames() before IntensityPCA(). (images are considered as float)
	* 
	* The Statistical Combined Model is computed on matrices using PCA.h.
	* Call SetVelocityFieldFileNames(), SetShapeAverageFileName(), SetShapeEVectorsFileNames(),
	* SetImageFileNames(), SetIntensityAverageFileName(), SetIntensityEVectorsFileNames(),
	* SetWFileName() before CombinedPCA().
	* 
	* The validation is performed recreating instances which are either part of the origina dataset (representation)
	* or not (Generalisation).
	* Call SetInstanceShapeFileNames(), SetShapeAverageFileName(), SetShapeEVectorsFileNames(),
	* SetInstanceIntensityFileNames(), SetIntensityAverageFileName(), SetIntensityEVectorsFileNames(),
	* SetCombinedEVectorsFileName(), SetCombinedEValuesFileName(),
	* SetNumberOfModes(), SetWFileName(), SetOutputFolder() before InstanceRecreation().
	*
	* The creation of new instaces is performed giving the file of weight for the modes.
	* Call SetShapeAverageFileName(), SetShapeEVectorsFileNames(),
	* SetIntensityAverageFileName(), SetIntensityEVectorsFileNames(),
	* SetCombinedEVectorsFileName(), SetCombinedEValuesFileName(), SetCombinedWeightsFileName()
	* SetNumberOfInstances(), SetNumberOfModes(),SetWFileName(), SetOutputFolder() before InstanceCreation().
	*/

	class PCAImages: public PCA{
	

	public:

		// image
		typedef float VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;

		// vector image
		typedef float FieldVoxelType;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 

		// readers
		typedef itk::ImageFileReader< ImageType > ImageReaderType;
		typedef itk::ImageFileReader< FieldType > FieldReaderType;
		
		// writers
		typedef itk::ImageFileWriter< ImageType > ImageWriterType;
		typedef itk::ImageFileWriter< FieldType > FieldWriterType;
		
		// iterators
		typedef itk::ImageRegionIterator< ImageType > ImageRegionIteratorType;
		typedef itk::ImageRegionIterator< FieldType > FieldRegionIteratorType;

		
		// constructor/destructor
		PCAImages();
		~PCAImages();

		// members
		/**
		* Sets the velocity field file names 
		*/
		void SetVelocityFieldFileNames(QStringList velocityFieldFileNames) { _velocityFieldFileNames = velocityFieldFileNames; }
		/**
		* Sets the shape eigenvectors file names 
		*/
		void SetShapeEVectorsFileNames(QStringList shapeEVectorsFileNames) { _shapeEVectorsFileNames = shapeEVectorsFileNames; }
		/**
		* Sets the intensity eigenvectors file names 
		*/
		void SetIntensityEVectorsFileNames(QStringList intensityEVectorsFileNames) { _intensityEVectorsFileNames = intensityEVectorsFileNames; }
		
		// members - new instance creation
		/**
		* Sets combined eigenvectors file names 
		*/
		void SetCombinedEVectorsFileNames(QStringList combinedEVectorsFileNames) { _combinedEVectorsFileNames = combinedEVectorsFileNames; }
		/**
		* Sets combined eigenvalues file names 
		*/
		void SetCombinedEValuesFileName(QString combinedEValuesFileName) { _combinedEValuesFileName = combinedEValuesFileName; }
		/**
		* Sets output folder
		*/
		void SetOutputFolder(QString outputFolder) { _outputFolder = outputFolder; }
			
		
		// overwritten virtual function
		void ShapePCA();
		void IntensityPCA();
		void CombinedPCA();
		void InstanceCreation();
		void InstanceRecreation();


	//private:
		void RecreateShapePca();
		

	protected:
		
		QStringList _velocityFieldFileNames;
		QStringList _shapeEVectorsFileNames;
		QStringList _intensityEVectorsFileNames;
		QString _combinedEValuesFileName;
		QStringList _combinedEVectorsFileNames;
		QString _outputFolder;
		
	};
}

#endif //PCAIMAGES_H 