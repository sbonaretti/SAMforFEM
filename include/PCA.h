/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCA_H
#define PCA_H

#include <QString>
#include <QStringList>

#include <vnl/vnl_vector.h> 
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_svd.h>

namespace pca{
	
	/**
	* Abstract class for the computation of the Principal Component Analysis.
	*
	* It is built using the same static functions that are in vtkPCAAnalysisFilter.
	* The input and output data are in vnl_vector and vnl_matrix formats for compatibility with other classes.
	*/
	
	class PCA{

	public:
		
		

		// constructor/destructor
		PCA();
		~PCA();

		
		// accessors
		
		/**
		* Sets training set images file names
		*/
		void SetImageFileNames (QStringList fileNames) { _imageFileNames = fileNames; }
		/**
		* Sets training set mesh file names
		*/
		void SetMeshFileNames (QStringList fileNames) { _meshFileNames = fileNames; }
		
		/**
		* Sets the shape average file name
		*/
		void SetShapeAverageFileName (QString fileName) { _shapeAverageFileName = fileName; }
		/**
		* Sets the shape eigenvalues file name
		*/
		void SetShapeEValuesFileName (QString fileName) { _shapeEValuesFileName = fileName; }
		
		/**
		* Sets the intensity average file name
		*/
		void SetIntensityAverageFileName (QString fileName) { _intensityAverageFileName = fileName; }
		/**
		* Sets the intensity eigenvalues file name
		*/
		void SetIntensityEValuesFileName (QString fileName) { _intensityEValuesFileName = fileName; }
		
		/**
		* Sets the b file name
		*/
		void SetBFileName (QString fileName) { _bFileName = fileName; }
		/**
		* Sets the combined eigenvalues file name
		*/
		void SetCombinedEValuesFileName (QString fileName) { _combinedEValuesFileName = fileName; }
		/**
		* Sets the combined eigenvalues file name
		*/
		void SetCombinedEVectorsFileName (QString fileName) { _combinedEVectorsFileName = fileName; }
		/**
		* Sets the W file name
		*/
		void SetWFileName (QString fileName) { _wFileName = fileName; }
		/**
		* Sets the nOfModes file name
		*/
		void SetModeNumbersFileName (QString fileName) { _modeNumbersFileName = fileName; }
		
		
		// for the protected members
		/**
		* Set the matrix on which the PCA is calculated. Columns are the observations 
		*/
		void SetMatrix(vnl_matrix<double> matrix) { _matrix = matrix; }
		/**
		* Gets the mean
		*/
		vnl_vector<double> GetMean() const { return _mean; }
		/**
		* Gets the eigenvalues
		*/
		vnl_vector<double> GetEValues() const { return _eValues; }
		/**
		* Gets the normalized eigenvalues
		*/
		vnl_vector<double> GetNormalizedEValues() const { return _normalizedEValues; }
		/**
		* Gets the eigenvectors
		*/
		vnl_matrix<double> GetEVectors() const { return _eVectors; }

		// for the creation of new instances
		/**
		* Sets the number of modes
		*/
		void SetNumberOfModes (int value) { _nOfModes = value; }
		/**
		* Sets the number of instances
		*/
		void SetNumberOfInstances (int value) { _nOfInstances = value; }
		/**
		* Sets the weights file name
		*/
		void SetCombinedWeightsFileName (QString fileName) { _combinedWeightsFileName = fileName; }
	

		// for the instance recreation (validation)
		/**
		* Sets the instance shape file name
		*/
		void SetInstanceShapeFileNames (QStringList fileNames) { _instanceShapeFileNames = fileNames; }
		/**
		* Sets the instance intensity file name
		*/
		void SetInstanceIntensityFileNames (QStringList fileNames) { _instanceIntensityFileNames = fileNames; }

	

		// pure virtual functions
		/**
		* Executes the shape PCA
		*/
		virtual void ShapePCA() = 0;
		/**
		* Executes the intensity PCA
		*/
		virtual void IntensityPCA() = 0;
		/**
		* Executes the combined PCA
		*/
		virtual void CombinedPCA() = 0;
		/**
		* Calculates the new instances
		*/
		virtual void InstanceCreation() = 0;
		/**
		* Calculates the parameters of given instances in order to recreate them
		*/
		virtual void InstanceRecreation() = 0;

		
	protected:
		
		// data members of the public functions
		QStringList _imageFileNames;
		QStringList _meshFileNames;
		
		QString _shapeAverageFileName;
		QString _shapeEValuesFileName;
		
		QString _intensityAverageFileName;
		QString _intensityEValuesFileName;
		
		QString _bFileName;
		QString _combinedEValuesFileName;
		QString	_combinedEVectorsFileName;
		QString _wFileName;
		QString _modeNumbersFileName;
		
		// creation of new instances and creation of existing ones (validation)
		int _nOfInstances;
		int _nOfModes;
		QString _combinedWeightsFileName;		
		QStringList _instanceShapeFileNames;
		QStringList _instanceIntensityFileNames;
		vnl_vector<double> _modeNumbers;

		
		// data members of the protected members
		vnl_matrix<double> _matrix;
		vnl_vector<double> _mean;
		vnl_vector<double> _eValues;
		vnl_vector<double> _normalizedEValues;
		vnl_matrix<double> _eVectors;

		/**
		* Calculates principle component analysis on a vertical matrix (simplified formula)
		*/
		void PrincipalComponentAnalysis();
		
	};
}
#endif // PCA_H
		