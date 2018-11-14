/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCAMESH_H
#define PCAMESH_H

#include <PCA.h>

#include <vtkPCAAnalysisFilter.h>
#include <vtkPolyData.h>



namespace pca{

	
		/**
		* It computes the Principal Component Analysis for meshes.
		* 
		* SHAPE PCA
		* Meshes have to be aligned, not necessary in the center (0,0,0) (use RegistrationProcrustesAlignment before). 
		* Input file types can be ansys quadratic volume mesh (.cdb), abaqus quadratic volume mesh (.inp) and quadratic surface mesh (.stl)
		* Eigenvalues, eigenvectors and mean are saved as .txt files (the mean is saved also in the same format as the input mesh).
		* The implementation is done using vtkPCAAnalysisFilter
		* Use SetMeshFileNames() to set the input. To get the outputs refer to PCA.h
		*
		* INTENSITY PCA
		* Meshes and images must be in the same coordinate system. 
		* Image extrusion is performed (see Extrusion() of ImageHandler.h) and grey levels linearly interpolated (see GreyLevelsExtractor() from FemAssignerNode)
		* Eigenvalues, eigenvectors and mean are saved as .txt files
		* Use SetMeshFileNames() and SetImageFileNames() to set the input. To get the outputs refer to PCA.h
		*
		* COMBINED PCA
		* The output from the previous PCAs are used
		* Use ShapeAveragefileName() to load the average in order to know the number of points and coordinates.
		* Use SetShapeDatasetCoordinatesFileName(), SetShapeEValuesFileName() and SetShapeEVectorsFileName() to set the input for the shape part.
		* Use SetIntensityDatasetGreyLevelFileName(), SetIntensityEValuesFileName(), SetIntensityEVectorsFileName() to set the input for the intensity part.
		* Use SetNOfInstances() to set the number of instances that were used in the shape and combined PCAs.
		* To get the outputs refer to PCA.h
		*/

	
	class PCAMesh: public PCA{
	

	public:

		// constructor/destructor
		PCAMesh();
		~PCAMesh();
		
		// for the combined PCA
		// Set shape PCA file names
		/**
		* Sets the shape dataset coordinates file name
		*/
		void SetShapeDatasetCoordinatesFileName (QString fileName) { _shapeDatasetCoordinatesFileName = fileName; }
		/**
		* Sets the shape eigenvalues file name
		*/
		void SetShapeEVectorsFileName (QString fileName) { _shapeEVectorsFileName = fileName; }
		
		// Set intensity PCA file names (for the combined PCA)
		/**
		* Sets the intensity dataset greylevel file name
		*/
		void SetIntensityDatasetGreyLevelFileName (QString fileName) { _intensityDatasetCoordinatesFileName = fileName; }
		/**
		* Sets the intensity eigenvalues file name
		*/
		void SetIntensityEVectorsFileName (QString fileName) { _intensityEVectorsFileName = fileName; }
		
		// Set combined PCA file names (for the instances creation)
		/**
		* Sets the output folder
		*/
		void SetOutputFolder (QString fileName) { _outputFolder = fileName; }
		
			
		// overwritten virtual functions
		void ShapePCA();
		void IntensityPCA();
		void CombinedPCA();
		void InstanceCreation();
		void InstanceRecreation();
		void ShapeInstanceRecreation();
		
		// function to test vnl matrix pca
		void TestMatrixPCA();

	private:
		void LoadShapeAverage();
		void LoadShapeEValues();
		void LoadShapeEVectors();
		void LoadIntensityAverage();
		void LoadIntensityEValues();
		void LoadIntensityEVectors();
		void LoadCombinedEValues();
		void LoadCombinedEVectors();
		void LoadWandNofMeshes();
		void LoadModeNumbers(); 
		void LoadShapeToRecreate(QString fileName);
		void LoadIntensityToRecreate(QString fileName);
		void LoadWeightsForInstanceCreation();
		void WriteShapeRecreated(QString fileName);
		void WriteIntensityRecreated(QString fileName, vnl_vector<double> vector);
		
		
		vtkPolyData* _shapeAverageMesh;
		vnl_vector<double> _shapeAverage;
		vnl_matrix<double> _shapeEigenvectors;
		vnl_vector<double> _intensityAverage;
		vnl_matrix<double> _intensityEigenvectors;
		vnl_vector<double> _combinedEValues;
		vnl_matrix<double> _combinedEigenVectors;
		double _w;
		vnl_diag_matrix<double> _W;
		vnl_diag_matrix<double> _Winv;
		double _nOfOriginalMeshes;
		vtkPolyData* _meshToRecreate;
		vnl_vector<double> _shapeToRecreate;
		vnl_vector<double> _intensityToRecreate;
		vnl_matrix<double> _weights;
				

	protected:
		// data members of the public members
		QString _shapeDatasetCoordinatesFileName;
		QString _shapeEVectorsFileName;
		
		QString _intensityDatasetCoordinatesFileName;
		QString	_intensityEVectorsFileName;
		
		QString _outputFolder;

	};
}
#endif //PCAMESH_H 