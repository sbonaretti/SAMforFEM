/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLREADERMATRIX_H
#define VNLREADERMATRIX_H

#include <VnlReader.h>

#include <vnl/vnl_matrix.h>

namespace vnl{

	/**
	* Reads the vnl matrix from a text file. 
	*
	* The first two numbers of the file must represent the number of row and the number of columns
	*/

	class VnlReaderMatrix: public VnlReader{

	public:
		
		// constructor/destructor
		VnlReaderMatrix();
		~VnlReaderMatrix();

		// accessor
		/**
		* Sets the vnl vector
		*/
		vnl_matrix<double> GetVnlMatrix() {return _matrix;}


		// overwritten virtual function
		void Update();

	protected:
		vnl_matrix<double> _matrix;

					
	};
}
#endif //VNLREADERMATRIX_H