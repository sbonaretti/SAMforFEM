/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLWRITERMATRIX_H
#define VNLWRITERMATRIX_H

#include <VnlWriter.h>

#include <vnl/vnl_matrix.h>

namespace vnl{

	/**
	* Writes the vnl matrix into a text file. 
	*
	* The first two numbers represent the number of row and the number of columns
	*/

	class VnlWriterMatrix: public VnlWriter{

	public:
		
		// constructor/destructor
		VnlWriterMatrix();
		~VnlWriterMatrix();

		// accessor
		/**
		* Sets the vnl vector
		*/
		void SetVnlMatrix(vnl_matrix<double> matrix) {_matrix = matrix;}


		// overwritten virtual function
		/**
		* Writes the vnl matrix.
		* The matrix is written like an array (rows are written one after the other)
		*/
		void Update();
		
		// member
		/**
		* Writes the vnl matrix.
		* It keeps the matrix shape
		*/
		void MatrixShapeUpdate();

	protected:
		vnl_matrix<double> _matrix;

					
	};
}
#endif //VNLWRITERMATRIX_H