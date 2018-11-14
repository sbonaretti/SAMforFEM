/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <iostream>
#include <PCA.h>

#include <vtkMath.h>
#include <vtkInformation.h>

namespace pca{

	// static functions from vtkPCAAnalysisFilter
	static inline double** NewMatrix(int rows, int cols) {
		
		double *matrix = new double[rows*cols];
		double **m = new double *[rows];
		for(int i = 0; i < rows; i++) {
			m[i] = &matrix[i*cols];
		}
		return m;
	}
	
	static inline double* NewVector(int length){
		double *vec = new double[length];
		return vec;
	}


	static inline void SubtractMeanColumn(double **m, double *mean, int rows, int cols){
		int r,c;
		double csum;
		for (r = 0; r < rows; r++) {
			csum = 0.0F;
		for (c = 0; c < cols; c++) {
			csum += m[r][c];
		}
		// calculate average value of row
		csum /= cols;

		// Mean shape vector is updated
		mean[r] = csum;

		// average value is subtracted from all elements in the row
		for (c = 0; c < cols; c++) {
		  m[r][c] -= csum;
		}
		}
	}

	static inline void SmallCovarianceMatrix(double **a, double **c,
										 int arows, int acols) { 
		const int s = acols;

		// c must have size acols*acols (we assume this)
		for(int i = 0; i < acols; i++) {
		for(int j = 0; j < acols; j++) {
		  // Use symmetry
		  if (i <= j) {
			  c[i][j] = 0.0;
			  for(int k = 0; k < arows; k++) {
				  c[i][j] += a[k][i]*a[k][j];
			  }
			  c[i][j] /= (s-1);
			  c[j][i] = c[i][j];
		  }
		}
		}
	} 

	static inline void MatrixMultiply(double **a, double **b, double **c,
                                  int arows, int acols, 
                                  int brows, int bcols) {
	if(acols != brows)  {
		return; // acols must equal br otherwise we can't proceed
	}
  
	// c must have size arows*bcols (we assume this)
	for(int i = 0; i < arows; i++) {
		for(int j = 0; j < bcols; j++) {
			c[i][j] = 0.0;
			for(int k = 0; k < acols; k++) {
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	}


	
	static inline void NormaliseColumns(double **m, int rows, int cols){
		// Normalise all eigenvectors columns to have length 1
		for (int c = 0; c < cols; c++){
		double cl = 0;
		for (int r = 0; r < rows; r++){
			cl += m[r][c] * m[r][c];
		}
		cl = sqrt(cl);
	    
		// If cl == 0 something is rotten, dont do anything now
		if (cl != 0) {
			for (int r = 0; r < rows; r++) {
				m[r][c] /= cl;
			}
		}
		}
	}



	
	// constructor
	PCA::PCA(){
	}

	// destructor
	PCA::~PCA(){
	}

	// protected functions
	void PCA::PrincipalComponentAnalysis (){

		int r = _matrix.rows();
		int c = _matrix.cols();

		
		// if matrix is horizontal, it is turned to vertical
		if (r < c){
			_matrix = _matrix.transpose();
			r = _matrix.rows();
			c = _matrix.cols();
		}

		// matrix put in the form wanted by the static functions
		double **vtkD = NewMatrix(r, c);
		for (int i = 0; i < r; i++){
			for (int j = 0; j < c; j++){
				vtkD[i][j] = _matrix(i,j);
			}
		}

		// calculating and subtracting the mean
		double *meanshape = NewVector(r);
		SubtractMeanColumn(vtkD, meanshape, r, c);
		_mean.set_size(r);
		for (int i=0; i<r; i++){
			_mean(i) = meanshape[i];
		}
	
		// covariance matrix calculation
		double **vtkT = NewMatrix(c, c);
		SmallCovarianceMatrix(vtkD, vtkT, r, c);

		// matrix decomposition, eigenvalues vector and eigenvalues matrix allocation
		double *ev = NewVector(c);
		double **evecMat = NewMatrix(c, c);
		vtkMath::JacobiN(vtkT, c, ev, evecMat);

		// eigenvalues
		_eValues.set_size(c);
		for (int i=0; i<c; i++){
			_eValues(i) = ev[i];
		}

		// normalize eigenvalues
		double sum = _eValues.sum();
		_normalizedEValues.set_size(c);
		for (int i=0; i<c; i++)
			_normalizedEValues(i) = _eValues(i) / sum * 100;

		// compute eigenvecs of DD' instead of T which is D'D
		double **evecMat2 = NewMatrix(r, c);
		MatrixMultiply(vtkD, evecMat, evecMat2, r, c, c, c);

		// normalise eigenvectors
		NormaliseColumns(evecMat2, r, c);
		_eVectors.set_size(r,c);
		for (int i = 0; i < r; i++){
			for (int j = 0; j < c; j++){
				_eVectors(i,j) = evecMat2[i][j];
			}
		}


		/*
		////// vnl code for comparison //////
		std::cout << "vnl code" << std::endl;
		
		// mean column
		_mean.set_size(r);
		vnl_vector<double> temp;
		temp.set_size(c);
		for (int i=0; i<r; i++){
			temp = _matrix.get_row(i);
			_mean(i) = temp.mean();
		}
		for (int i=0; i<9; i++)
			std::cout << "mean " << _mean(i) << std::endl;
	
		// create matrix whose columns are the mean
		vnl_matrix<double> meanMatrix;
		meanMatrix.set_size(r,c);
		for (int i=0; i<c; i++)
			meanMatrix.set_column(i, _mean);

		// create matrix D 
		vnl_matrix<double> D;
		D.set_size(r,c);
		D = _matrix - meanMatrix;
		//for (int i=0; i<9; i++) std::cout << "D " << D(i,0) << " " << D(i,1) << " " << D(i,2) << std::endl;
		
		// covariance T
		vnl_matrix<double> T;
		T.set_size(c,c);
		vnl_matrix<double> Dtransp;
		Dtransp.set_size(c,r);
		Dtransp  = D.transpose();
		T = (Dtransp * D)/ (c-1);
		// for (int i=0; i<c; i++) std::cout << "T " << T(i,0) << " " << T(i,1) << " " << T(i,2) << std::endl;
		
		// svd 
		vnl_svd<double> svd (T);
		
		// eigenvalues
		vnl_diag_matrix<double> diag = svd.W();
		_eValues.set_size(c);
		_eValues = diag.diagonal();
		for (int i=0; i<c; i++)
			std::cout << "_eValues " << _eValues(i) << std::endl;

		
		// eigenvectors and their normalization missing
		
		std::cout << "end of vnl code" << std::endl;
		*/
	}

}